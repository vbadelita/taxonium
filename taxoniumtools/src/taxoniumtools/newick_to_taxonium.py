import argparse
import orjson
import json
import pandas as pd
import datetime
import gzip
import treeswift
from alive_progress import alive_it, alive_bar

import sys
import os
import logging

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
    return parser


def do_processing(input_tree,
                  output_file,
                  aln_files=None,
                  metadata_file=None,
                  genbank_files=None,
                  columns="",
                  title=None,
                  key_column="strain"):
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
        from treetime import TreeAnc
        if not isinstance(aln_files, list):
            aln_files = [aln_files]
        if genbank_files and not isinstance(genbank_files, list):
            genbank_files = [genbank_files]

        if genbank_files and len(genbank_files) == 1 and len(aln_files) > 1:
            genbank_files = genbank_files * len(aln_files)
        elif not genbank_files:
            genbank_files = [None] * len(aln_files)

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

        all_gene_details = []

        for aln_idx, (aln_f, gb_f) in enumerate(zip(aln_files, genbank_files)):
            segment_name = os.path.splitext(os.path.basename(aln_f))[0]
            print(
                f"Processing segment {aln_idx+1}/{len(aln_files)}: {segment_name}"
            )

            print(
                f"  Running Ancestral Sequence Reconstruction with TreeTime for {segment_name}..."
            )
            ta = TreeAnc(tree=treetime_input, aln=aln_f, gtr='JC69', verbose=0)
            ta.infer_ancestral_sequences()

            dummy = None
            if gb_f:
                print(f"  Loading GenBank annotations for {segment_name}...")
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
                root_seq = getattr(ta_node_dict[tree.root.label], 'sequence',
                                   None)
                if root_seq is None:
                    print(
                        f"  Root sequence missing from TreeTime for {segment_name}, defaulting to GenBank sequence."
                    )
                    root_seq = str(dummy.genbank.seq)
                else:
                    root_seq = "".join(root_seq)

                with alive_bar(
                        tree.num_nodes(),
                        title=f"  Annotating AA for {segment_name}") as pbar:
                    core_mutations.recursive_mutation_analysis(
                        tree.root, {},
                        root_seq,
                        dummy.cdses,
                        pbar,
                        dummy.nuc_to_codon,
                        chromosome=segment_name)

                # Add root mutations for this segment
                root_muts = []
                for i, character in enumerate(root_seq):
                    root_muts.append(
                        core_mutations.NucMutation(one_indexed_position=i + 1,
                                                   mut_nuc=character,
                                                   par_nuc="X",
                                                   chromosome=segment_name))

                segment_root_aa_muts = core_mutations.get_mutations(
                    {},
                    root_muts,
                    root_seq,
                    dummy.nuc_to_codon,
                    disable_check_for_differences=True,
                    chromosome=segment_name)
                tree.root.aa_muts.extend(segment_root_aa_muts)
                tree.root.nuc_mutations.extend(root_muts)

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

    if not input_tree or not output_file:
        print(
            "Error: input and output are required (either via CLI or config file)"
        )
        sys.exit(1)

    do_processing(input_tree=input_tree,
                  output_file=output_file,
                  aln_files=aln_files,
                  metadata_file=metadata_file,
                  genbank_files=genbank_files,
                  columns=columns,
                  title=title,
                  key_column=key_column)


if __name__ == "__main__":
    main()
