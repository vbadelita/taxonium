import argparse
import orjson
import json
import pandas as pd
import datetime
import gc
import gzip
import sys
import os
import logging
import treeswift
from alive_progress import alive_it, alive_bar

from Bio import SeqIO

logging.getLogger('treetime').setLevel(logging.ERROR)

if __package__ is None or __package__ == "":
    sys.path.insert(
        0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from taxoniumtools import core_mutations
    from taxoniumtools import utils
    try:
        from taxoniumtools import _version
        version = _version.version
    except ImportError:
        version = "dev"
else:
    from . import core_mutations
    from . import utils
    try:
        from . import _version
        version = _version.version
    except ImportError:
        version = "dev"


def get_parser():
    parser = argparse.ArgumentParser(
        description=
        'Convert a Newick tree, fasta alignment, and metadata to Taxonium jsonl format'
    )
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        help='File path to input Newick tree')
    parser.add_argument(
        '-a',
        '--aln',
        type=str,
        help='File path(s) to input fasta alignment(s) (optional)',
        nargs='+')
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='File path for output Taxonium jsonl file')
    parser.add_argument('-m',
                        '--metadata',
                        type=str,
                        help='File path for input metadata file (CSV/TSV)')
    parser.add_argument(
        '-g',
        '--genbank',
        type=str,
        help='File path(s) for GenBank file(s) containing reference genome(s)',
        nargs='+')
    parser.add_argument(
        '-c',
        '--columns',
        type=str,
        help='Column names to include in the metadata, separated by commas',
        default="")
    parser.add_argument('-t',
                        '--title',
                        type=str,
                        help='A title for the tree.',
                        default=None)
    parser.add_argument(
        '--key_column',
        type=str,
        help=
        'The column in the metadata file which is the same as the names in the tree',
        default=None)
    parser.add_argument('-j',
                        '--config',
                        type=str,
                        help='Config file in JSON format')
    parser.add_argument(
        '--ancestral_reconstruction',
        type=str,
        choices=['treetime', 'parsimony'],
        default=None,
        help='Ancestral sequence reconstruction method for all segments: '
        '"treetime" (accurate, slow, memory-heavy) or "parsimony" (naive but fast). '
        'In config JSON, set "ancestral_reconstruction": "treetime"|"parsimony" per segment.'
    )
    return parser


def make_root_aa_mutations(ref_seq, loader, segment_name):
    """Generate X->ref_AA AAMutations for every codon in every gene.

    These are placed on the root so Taxonium knows the reference amino acid
    state for all positions in this segment.
    """
    aa_muts = []
    for codon_positions_by_index in loader.nuc_to_codon.values():
        pass  # we iterate below via genes
    seen_codons = set()
    for genome_pos, codons in loader.nuc_to_codon.items():
        for codon in codons:
            if codon in seen_codons:
                continue
            seen_codons.add(codon)
            codon_seq = [ref_seq[codon.positions[x]] for x in range(3)]
            codon_str = "".join(codon_seq)
            if codon.strand == -1:
                codon_str = core_mutations.complement(codon_str)
            ref_aa = core_mutations.codon_table.get(codon_str, "X")
            aa_muts.append(
                core_mutations.AAMutation(
                    gene=codon.gene,
                    one_indexed_codon=codon.codon_number + 1,
                    initial_aa="X",
                    final_aa=ref_aa,
                    nuc_for_codon=codon.positions[1]))
    return aa_muts


def make_missing_aa_mutations(loader):
    """Generate X->X AAMutations for every codon in every gene.

    These are placed on tip nodes that are absent from a segment's alignment,
    so Taxonium shows X (unknown) rather than inheriting the parent's inferred
    state.
    """
    aa_muts = []
    seen_codons = set()
    for genome_pos, codons in loader.nuc_to_codon.items():
        for codon in codons:
            if codon in seen_codons:
                continue
            seen_codons.add(codon)
            aa_muts.append(
                core_mutations.AAMutation(
                    gene=codon.gene,
                    one_indexed_codon=codon.codon_number + 1,
                    initial_aa="X",
                    final_aa="X",
                    nuc_for_codon=codon.positions[1]))
    return aa_muts


def parsimony_ancestral_reconstruction(tree, tip_label_to_node, aln_f, ref_seq,
                                       loader, segment_name):
    """Naive parsimony reconstruction.

    For each tip: compute absolute mutations vs reference (set of NucMutation).
    Post-order: internal node's absolute set = intersection of children that
    have data; children with no data are ignored.
    Pre-order: branch mutations = absolute set minus parent's absolute set,
    plus reversions (in parent but not in node). Then annotate AA mutations.
    """
    print(f"  Reading alignment and annotating tips for {segment_name}...")
    matched = 0
    for record in SeqIO.parse(aln_f, "fasta"):
        if record.id not in tip_label_to_node:
            continue
        matched += 1
        node = tip_label_to_node[record.id]
        tip_seq = str(record.seq).upper()
        abs_muts = {}
        for i, (r, t) in enumerate(zip(ref_seq, tip_seq)):
            if r != t and r != '-' and r != 'N' and t != '-' and t != 'N':
                abs_muts[i + 1] = (r, t)
        node._parsimony_abs = abs_muts

    print(f"  Matched {matched} tips for {segment_name}.")

    # Post-order: propagate up
    print(
        f"  Propagating mutations up tree (post-order) for {segment_name}...")
    for node in tree.traverse_postorder():
        if node.is_leaf():
            if not hasattr(node, '_parsimony_abs'):
                node._parsimony_abs = None  # missing from alignment
        else:
            child_sets = [
                c._parsimony_abs for c in node.children
                if c._parsimony_abs is not None
            ]
            if not child_sets:
                node._parsimony_abs = None
            elif len(child_sets) == 1:
                node._parsimony_abs = dict(child_sets[0])
            else:
                # Intersection: positions where ALL children agree on same mut
                common_positions = set(child_sets[0].keys())
                for cs in child_sets[1:]:
                    common_positions &= cs.keys()
                node._parsimony_abs = {
                    pos: child_sets[0][pos]
                    for pos in common_positions if all(
                        cs[pos] == child_sets[0][pos] for cs in child_sets[1:])
                }

    # Pre-order: convert absolute -> branch mutations, annotate AA.
    # At the boundary where a node has no data but its parent does (or it is the
    # root with no data), emit X->X AA mutations once. All descendants inherit X
    # by the normal parent-walk in Taxonium, so no per-tip mutations are needed.
    print(f"  Computing branch mutations (pre-order) for {segment_name}...")
    missing_aa_muts = make_missing_aa_mutations(loader) if loader else []
    for node in tree.traverse_preorder():
        if node._parsimony_abs is None:
            parent_has_data = (node.parent is not None
                               and node.parent._parsimony_abs is not None)
            is_root_with_no_data = node.parent is None
            if parent_has_data or is_root_with_no_data:
                node.aa_muts.extend(missing_aa_muts)
            continue

        my_abs = node._parsimony_abs
        par_abs = (node.parent._parsimony_abs or {}) if node.parent else {}

        branch_muts = []
        # Mutations gained on this branch
        for pos, (ref_nuc, mut_nuc) in my_abs.items():
            if pos not in par_abs:
                branch_muts.append(
                    core_mutations.NucMutation(one_indexed_position=pos,
                                               par_nuc=ref_nuc,
                                               mut_nuc=mut_nuc,
                                               chromosome=segment_name))
        # Reversions: in parent but not in this node
        for pos, (ref_nuc, mut_nuc) in par_abs.items():
            if pos not in my_abs:
                branch_muts.append(
                    core_mutations.NucMutation(one_indexed_position=pos,
                                               par_nuc=mut_nuc,
                                               mut_nuc=ref_nuc,
                                               chromosome=segment_name))

        node.nuc_mutations.extend(branch_muts)
        if loader and branch_muts:
            node.aa_muts.extend(
                core_mutations.get_mutations({},
                                             branch_muts,
                                             ref_seq,
                                             loader.nuc_to_codon,
                                             chromosome=segment_name))

    # Free the temporary absolute mutation sets
    for node in tree.traverse_preorder():
        if hasattr(node, '_parsimony_abs'):
            del node._parsimony_abs


def do_processing(input_tree,
                  output_file,
                  aln_files=None,
                  metadata_file=None,
                  genbank_files=None,
                  columns="",
                  title=None,
                  key_column="strain",
                  ancestral_reconstruction_flags=None):
    metadata_dict, metadata_cols = utils.read_metadata(metadata_file, columns,
                                                       key_column)

    config = {}
    if title is not None:
        config['title'] = title

    print("Loading tree...")
    try:
        import dendropy
        dt = dendropy.Tree.get(path=input_tree,
                               schema='newick',
                               preserve_underscores=True,
                               case_sensitive_taxon_labels=True)
        clean_newick = dt.as_string(schema='newick',
                                    suppress_internal_node_labels=True)
        tree = treeswift.read_tree_newick(clean_newick)

        # Save clean tree to temporary file for TreeAnc
        import tempfile
        fd, temp_tree_path = tempfile.mkstemp(suffix=".nwk")
        with os.fdopen(fd, 'w') as f:
            f.write(clean_newick)
        treetime_input = temp_tree_path
    except ImportError:
        tree = treeswift.read_tree_newick(open(input_tree).read())
        treetime_input = input_tree
        temp_tree_path = None

    # Optional sequence inference
    if aln_files:
        if not isinstance(aln_files, list):
            aln_files = [aln_files]
        if genbank_files and not isinstance(genbank_files, list):
            genbank_files = [genbank_files]

        if genbank_files and len(genbank_files) == 1 and len(aln_files) > 1:
            genbank_files = genbank_files * len(aln_files)
        elif not genbank_files:
            genbank_files = [None] * len(aln_files)

        if ancestral_reconstruction_flags is None:
            ancestral_reconstruction_flags = [None] * len(aln_files)
        elif not isinstance(ancestral_reconstruction_flags, list):
            ancestral_reconstruction_flags = [ancestral_reconstruction_flags
                                              ] * len(aln_files)

        # Initialize nodes
        for node in tree.traverse_preorder():
            node.nuc_mutations = []
            node.aa_muts = []

        print("Naming internal nodes...")
        for i, node in enumerate(
                core_mutations.preorder_traversal_internal(tree.root)):
            if not node.label:
                node.label = f"NODE_{i+1:07d}"
        if not tree.root.label:
            tree.root.label = "NODE_0000000"

        # Build label->node lookup for tips
        tip_label_to_node = {}
        for node in tree.traverse_leaves():
            if node.label:
                tip_label_to_node[node.label] = node

        all_gene_details = []

        for aln_idx, (aln_f, gb_f) in enumerate(zip(aln_files, genbank_files)):
            ar_mode = ancestral_reconstruction_flags[
                aln_idx]  # None, "treetime", "parsimony"
            segment_name = os.path.splitext(os.path.basename(aln_f))[0]
            mode_label = f" ({ar_mode})" if ar_mode else " (tips only)"
            print(
                f"Processing segment {aln_idx+1}/{len(aln_files)}: {segment_name}{mode_label}"
            )

            if ar_mode is None:
                # Tips-only mode: read FASTA directly, diff against reference,
                # annotate only leaf nodes. No TreeTime.
                if not gb_f:
                    print(
                        f"  Warning: tips_only without GenBank for {segment_name}, skipping AA annotation."
                    )

                ref_seq = None
                loader = None
                if gb_f:
                    loader = core_mutations.GenbankLoader(gb_f)
                    ref_seq = str(loader.genbank.seq)
                    all_gene_details.extend(list(loader.genes.keys()))

                # Stream alignment sequences one at a time to avoid loading all into memory
                print(f"  Annotating tips for {segment_name}...")
                matched = 0
                matched_tips = set()
                for record in SeqIO.parse(aln_f, "fasta"):
                    if record.id not in tip_label_to_node:
                        continue
                    matched += 1
                    matched_tips.add(record.id)
                    node = tip_label_to_node[record.id]
                    tip_seq = str(record.seq).upper()

                    if ref_seq is None:
                        continue

                    # Compute nuc mutations vs reference
                    tip_nuc_muts = []
                    for i, (r, t) in enumerate(zip(ref_seq, tip_seq)):
                        if r != t and r != '-' and r != 'N' and t != '-' and t != 'N':
                            mut = core_mutations.NucMutation(
                                one_indexed_position=i + 1,
                                par_nuc=r,
                                mut_nuc=t,
                                chromosome=segment_name)
                            tip_nuc_muts.append(mut)
                    node.nuc_mutations.extend(tip_nuc_muts)

                    # Compute AA mutations
                    if loader:
                        aa_muts = core_mutations.get_mutations(
                            {},
                            tip_nuc_muts,
                            ref_seq,
                            loader.nuc_to_codon,
                            chromosome=segment_name)
                        node.aa_muts.extend(aa_muts)

                print(
                    f"  Matched {matched}/{len(tip_label_to_node)} tips for {segment_name}."
                )

                # Mark missing tips with X->X so Taxonium shows unknown.
                # Place X->X only at the boundary node (highest ancestor whose
                # entire subtree is absent from this segment) rather than on
                # every missing tip, so descendants inherit X without needing
                # their own mutation entries.
                if loader:
                    missing_aa_muts = make_missing_aa_mutations(loader)
                    # Count matched tips in each node's subtree (postorder)
                    for node in tree.traverse_postorder():
                        if node.is_leaf():
                            node._seg_matched = 1 if node.label in matched_tips else 0
                        else:
                            node._seg_matched = sum(c._seg_matched
                                                    for c in node.children)
                    # Preorder: emit X->X at the transition boundary
                    for node in tree.traverse_preorder():
                        par_matched = (node.parent._seg_matched
                                       if node.parent else -1)
                        if node._seg_matched == 0 and par_matched != 0:
                            node.aa_muts.extend(missing_aa_muts)
                    for node in tree.traverse_preorder():
                        if hasattr(node, '_seg_matched'):
                            del node._seg_matched

                # Add root AA reference mutations so Taxonium knows the reference
                # amino acid state for all positions.
                if ref_seq and loader:
                    segment_root_aa_muts = make_root_aa_mutations(
                        ref_seq, loader, segment_name)
                    tree.root.aa_muts.extend(segment_root_aa_muts)

                del loader, ref_seq
                gc.collect()

            elif ar_mode == "parsimony":
                loader = None
                ref_seq = None
                if gb_f:
                    loader = core_mutations.GenbankLoader(gb_f)
                    ref_seq = str(loader.genbank.seq)
                    all_gene_details.extend(list(loader.genes.keys()))

                if ref_seq:
                    parsimony_ancestral_reconstruction(tree, tip_label_to_node,
                                                       aln_f, ref_seq, loader,
                                                       segment_name)

                    # Add root AA reference mutations
                    segment_root_aa_muts = make_root_aa_mutations(
                        ref_seq, loader, segment_name)
                    tree.root.aa_muts.extend(segment_root_aa_muts)
                else:
                    print(
                        f"  Warning: parsimony without GenBank for {segment_name}, skipping."
                    )

                del loader, ref_seq
                gc.collect()

            else:
                # Full TreeTime ancestral reconstruction mode
                from treetime import TreeAnc

                print(
                    f"  Running Ancestral Sequence Reconstruction with TreeTime for {segment_name}..."
                )
                ta = TreeAnc(tree=treetime_input,
                             aln=aln_f,
                             gtr='JC69',
                             verbose=0)
                ta.infer_ancestral_sequences()

                dummy = None
                if gb_f:
                    print(
                        f"  Loading GenBank annotations for {segment_name}...")
                    dummy = core_mutations.GenbankLoader(gb_f)
                    all_gene_details.extend(list(dummy.genes.keys()))

                ta_node_dict = {
                    node.name: node
                    for node in ta.tree.find_clades() if node.name is not None
                }

                print(
                    f"  Mapped sequences for {segment_name}. Extracting nucleotide mutations..."
                )
                for node in tree.traverse_preorder():
                    if node.parent:
                        my_seq = ta_node_dict[
                            node.
                            label].sequence if node.label in ta_node_dict else None
                        par_seq = ta_node_dict[
                            node.parent.
                            label].sequence if node.parent.label in ta_node_dict else None

                        if my_seq is not None and par_seq is not None:
                            for i, (p, m) in enumerate(zip(par_seq, my_seq)):
                                if p != m and p != '-' and p != 'N' and m != '-' and m != 'N':
                                    mut = core_mutations.NucMutation(
                                        one_indexed_position=i + 1,
                                        par_nuc=p,
                                        mut_nuc=m,
                                        chromosome=segment_name)
                                    node.nuc_mutations.append(mut)

                if gb_f:
                    print(
                        f"  Performing amino acid translation analysis for {segment_name}..."
                    )
                    root_seq = getattr(ta_node_dict[tree.root.label],
                                       'sequence', None)
                    if root_seq is None:
                        print(
                            f"  Root sequence missing from TreeTime for {segment_name}, defaulting to GenBank sequence."
                        )
                        root_seq = str(dummy.genbank.seq)
                    else:
                        root_seq = "".join(root_seq)

                    with alive_bar(tree.num_nodes(),
                                   title=f"  Annotating AA for {segment_name}"
                                   ) as pbar:
                        core_mutations.recursive_mutation_analysis(
                            tree.root, {},
                            root_seq,
                            dummy.cdses,
                            pbar,
                            dummy.nuc_to_codon,
                            chromosome=segment_name)

                    # Add root AA reference mutations for this segment
                    segment_root_aa_muts = make_root_aa_mutations(
                        root_seq, dummy, segment_name)
                    tree.root.aa_muts.extend(segment_root_aa_muts)

                    # Mark nodes absent from this segment with X->X at the
                    # boundary (highest ancestor whose entire subtree is absent)
                    # so descendants inherit X without per-tip entries.
                    missing_aa_muts = make_missing_aa_mutations(dummy)
                    for node in tree.traverse_postorder():
                        if node.is_leaf():
                            node._seg_matched = 1 if node.label in ta_node_dict else 0
                        else:
                            node._seg_matched = sum(c._seg_matched
                                                    for c in node.children)
                    for node in tree.traverse_preorder():
                        par_matched = (node.parent._seg_matched
                                       if node.parent else -1)
                        if node._seg_matched == 0 and par_matched != 0:
                            node.aa_muts.extend(missing_aa_muts)
                    for node in tree.traverse_preorder():
                        if hasattr(node, '_seg_matched'):
                            del node._seg_matched

                del ta, ta_node_dict
                gc.collect()

        config['gene_details'] = sorted(list(set(all_gene_details)))

        if temp_tree_path:
            os.remove(temp_tree_path)

    else:
        # No alignment, just basic tree structure
        for i, node in enumerate(
                core_mutations.preorder_traversal_internal(tree.root)):
            if not node.label:
                node.label = f"node_{i+1}"
        for node in tree.traverse_preorder():
            node.nuc_mutations = []
            node.aa_muts = []

    print("Assigning tips...")
    for node in tree.traverse_postorder():
        if node.is_leaf(): node.num_tips = 1
        else: node.num_tips = sum(child.num_tips for child in node.children)

    total_tips = tree.root.num_tips

    print("Ladderizing tree...")
    tree.ladderize(ascending=False)

    # We must ensure edge_length exists so set_x_coords works
    for node in tree.traverse_preorder():
        if node.edge_length is None:
            node.edge_length = 0.0001

    print("Setting coordinates...")
    utils.set_x_coords(tree.root, chronumental_enabled=False)
    utils.set_terminal_y_coords(tree.root)
    utils.set_internal_y_coords(tree.root)

    nodes_sorted_by_y = utils.sort_on_y(tree)

    all_aa_muts_objects = utils.get_all_aa_muts(tree.root)
    all_nuc_muts = utils.get_all_nuc_muts(tree.root)

    all_mut_inputs = all_aa_muts_objects + all_nuc_muts
    all_mut_objects = [
        utils.make_aa_object(i, thing)
        if thing.type == "aa" else utils.make_nuc_object(i, thing)
        for i, thing in enumerate(all_mut_inputs)
    ]

    input_to_index = {thing: i for i, thing in enumerate(all_mut_inputs)}

    config['num_tips'] = total_tips
    config['date_created'] = datetime.datetime.now().strftime("%Y-%m-%d")

    first_json = {
        "version": version,
        "mutations": all_mut_objects,
        "total_nodes": len(nodes_sorted_by_y),
        "config": config
    }

    node_to_index = {node: i for i, node in enumerate(nodes_sorted_by_y)}

    if "gz" in output_file:
        outfile = gzip.open(output_file, 'wb')
    else:
        outfile = open(output_file, 'wb')

    outfile.write(orjson.dumps(first_json) + b"\n")
    for node in alive_it(nodes_sorted_by_y, title="Converting JSON"):
        node_object = utils.get_node_object(node,
                                            node_to_index,
                                            metadata_dict,
                                            input_to_index,
                                            metadata_cols,
                                            chronumental_enabled=False)
        outfile.write(orjson.dumps(node_object) + b"\n")
    outfile.close()

    print(
        f"Done. Output written to {output_file}, with {len(nodes_sorted_by_y)} nodes."
    )


def main():
    parser = get_parser()
    args = parser.parse_args()

    config_data = {}
    if args.config:
        with open(args.config, 'r') as f:
            config_data = json.load(f)

    # CLI overrides config
    input_tree = args.input or config_data.get('input')
    output_file = args.output or config_data.get('output')
    metadata_file = args.metadata or config_data.get('metadata')
    columns = args.columns or config_data.get('columns', "")
    title = args.title or config_data.get('title')
    key_column = args.key_column or config_data.get(
        'metadata_id_column') or config_data.get('key_column', "strain")

    # For alignments and genbanks, we can take from segments or lists
    aln_files = args.aln
    if not aln_files:
        if 'segments' in config_data:
            aln_files = [s.get('alignment') for s in config_data['segments']]
        else:
            aln_files = config_data.get('aln')

    genbank_files = args.genbank
    if not genbank_files:
        if 'segments' in config_data:
            genbank_files = [s.get('genbank') for s in config_data['segments']]
        else:
            genbank_files = config_data.get('genbank')

    # ancestral_reconstruction: per-segment from config, or global from CLI flag
    # Values: None (tips only), "treetime", "parsimony"
    # Config may use true (legacy, treated as "treetime") or a string
    def _parse_ar(val):
        if val is True:
            return "treetime"
        if val in (False, None):
            return None
        return val  # already a string

    ancestral_reconstruction_flags = None
    if 'segments' in config_data:
        ancestral_reconstruction_flags = [
            _parse_ar(s.get('ancestral_reconstruction'))
            for s in config_data['segments']
        ]
    if args.ancestral_reconstruction:
        # CLI --ancestral_reconstruction overrides: set all segments
        n = len(aln_files) if aln_files else 0
        ancestral_reconstruction_flags = [args.ancestral_reconstruction] * n

    if not input_tree or not output_file:
        print(
            "Error: input and output are required (either via CLI or config file)"
        )
        sys.exit(1)

    do_processing(
        input_tree=input_tree,
        output_file=output_file,
        aln_files=aln_files,
        metadata_file=metadata_file,
        genbank_files=genbank_files,
        columns=columns,
        title=title,
        key_column=key_column,
        ancestral_reconstruction_flags=ancestral_reconstruction_flags)


if __name__ == "__main__":
    main()
