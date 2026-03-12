from collections import defaultdict
from dataclasses import dataclass
from Bio import SeqIO


def reverse_complement(input_string):
    return input_string.translate(str.maketrans("ATCGNagctn",
                                                "TAGCNtagcn"))[::-1]


def complement(input_string):
    return input_string.translate(str.maketrans("ATCGNagctn", "TAGCNtagcn"))


@dataclass(eq=True, frozen=True)
class AAMutation:
    gene: str
    one_indexed_codon: int
    initial_aa: str
    final_aa: str
    nuc_for_codon: int
    type: str = "aa"


@dataclass(eq=True, frozen=True)
class NucMutation:  #hashable
    one_indexed_position: int
    par_nuc: str
    mut_nuc: str
    chromosome: str = "chrom"
    type: str = "nt"


@dataclass(eq=True, frozen=True)
class Gene:
    name: str
    strand: int
    start: int  # zero-indexed
    end: int  # 0-indexed


@dataclass(eq=False)
class Codon:
    gene: Gene
    codon_number: int  # zero-indexed
    positions: dict  # zero-indexed positions e.g. {0:123,1:124,2:125}
    strand: int

    def __eq__(self, other):
        if isinstance(other, Codon):
            return self.gene == other.gene and self.codon_number == other.codon_number
        return False

    def __hash__(self):
        return hash((self.gene, self.codon_number))


def get_codon_table():
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    return dict(zip(codons, amino_acids))


codon_table = get_codon_table()


def get_gene_name(cds):
    """Returns gene if available, otherwise locus tag"""
    if "gene" in cds.qualifiers:
        return cds.qualifiers["gene"][0]
    elif "locus_tag" in cds.qualifiers:
        return cds.qualifiers["locus_tag"][0]
    else:
        raise ValueError(f"No gene name or locus tag for {cds}")


def get_genes_dict(cdses):
    genes = {}
    for cds in cdses:
        strand = cds.location.strand if cds.location.strand is not None else 1
        genes[get_gene_name(cds)] = Gene(get_gene_name(cds), strand,
                                         cds.location.start, cds.location.end)
    return genes


def get_mutations(past_nuc_muts_dict,
                  new_nuc_mutations_here,
                  seq,
                  nuc_to_codon,
                  disable_check_for_differences=False,
                  chromosome=None):
    by_codon = defaultdict(list)
    for mutation in new_nuc_mutations_here:
        if chromosome is not None and mutation.chromosome != chromosome:
            continue
        zero_indexed_pos = mutation.one_indexed_position - 1
        if zero_indexed_pos in nuc_to_codon:
            for codon in nuc_to_codon[zero_indexed_pos]:
                by_codon[codon].append(mutation)
    mutations_here = []
    for gene_codon, mutations in by_codon.items():
        initial_codon = [seq[gene_codon.positions[x]] for x in range(3)]
        relevant_past_muts = [(x, past_nuc_muts_dict[x])
                              for x in gene_codon.positions.values()
                              if x in past_nuc_muts_dict]
        flipped_dict = {
            position: offset
            for offset, position in gene_codon.positions.items()
        }
        for position, value in relevant_past_muts:
            initial_codon[flipped_dict[position]] = value
        final_codon = initial_codon.copy()
        for mutation in mutations:
            pos_in_codon = flipped_dict[mutation.one_indexed_position - 1]
            final_codon[pos_in_codon] = mutation.mut_nuc
        initial_codon = "".join(initial_codon)
        final_codon = "".join(final_codon)
        if gene_codon.strand == -1:
            initial_codon = complement(initial_codon)
            final_codon = complement(final_codon)
        initial_codon_trans = codon_table.get(initial_codon, "X")
        final_codon_trans = codon_table.get(final_codon, "X")
        if initial_codon_trans != final_codon_trans or disable_check_for_differences:
            mutations_here.append(
                AAMutation(gene=gene_codon.gene,
                           one_indexed_codon=gene_codon.codon_number + 1,
                           initial_aa=initial_codon_trans,
                           final_aa=final_codon_trans,
                           nuc_for_codon=gene_codon.positions[1]))
    for mutation in new_nuc_mutations_here:
        past_nuc_muts_dict[mutation.one_indexed_position -
                           1] = mutation.mut_nuc
    return mutations_here


def recursive_mutation_analysis(node,
                                past_nuc_muts_dict,
                                seq,
                                cdses,
                                pbar,
                                nuc_to_codon,
                                chromosome=None):
    pbar()
    new_nuc_mutations_here = node.nuc_mutations
    new_past_nuc_muts_dict = past_nuc_muts_dict.copy()

    current_aa_muts = get_mutations(new_past_nuc_muts_dict,
                                    new_nuc_mutations_here,
                                    seq,
                                    nuc_to_codon,
                                    chromosome=chromosome)

    if not hasattr(node, "aa_muts"):
        node.aa_muts = []
    node.aa_muts.extend(current_aa_muts)

    for child in node.children:
        recursive_mutation_analysis(child,
                                    new_past_nuc_muts_dict,
                                    seq,
                                    cdses,
                                    pbar,
                                    nuc_to_codon,
                                    chromosome=chromosome)


NUC_ENUM = "ACGT"


def preorder_traversal(node):
    yield node
    for clade in node.children:
        yield from preorder_traversal(clade)


def preorder_traversal_internal(node):
    yield node
    for clade in node.children:
        for x in preorder_traversal_internal(clade):
            if not x.is_leaf():
                yield x


def preorder_traversal_iter(node):
    return iter(preorder_traversal(node))


def find_cds(position, cdses):
    for cds in cdses:
        if cds.location.start <= position <= cds.location.end:
            return cds
    return None


def find_codon(position, cds):
    strand = cds.location.strand if cds.location.strand is not None else 1
    if strand == 1:
        codon_number = (position - cds.location.start) // 3
        codon_start = cds.location.start + codon_number * 3
        codon_end = codon_start + 3
    else:
        codon_number = (cds.location.end - position - 1) // 3
        codon_end = cds.location.end - codon_number * 3
        codon_start = codon_end - 3
    return codon_number, codon_start, codon_end


class GenbankLoader:

    def __init__(self, genbank_file):
        self.genbank = SeqIO.read(genbank_file, "genbank")
        self.cdses = [x for x in self.genbank.features if x.type == "CDS"]
        self.genes = get_genes_dict(self.cdses)
        by_everything = defaultdict(lambda: defaultdict(dict))
        total_lengths = {}
        for feature in self.cdses:
            gene_name = get_gene_name(feature)
            nucleotide_counter = 0
            for part in feature.location.parts:
                ranger = range(part.start,
                               part.end) if part.strand == 1 else range(
                                   part.end - 1, part.start - 1, -1)
                for genome_position in ranger:
                    cur_codon_number = nucleotide_counter // 3
                    cur_pos_in_codon = nucleotide_counter % 3
                    by_everything[gene_name][cur_codon_number][
                        cur_pos_in_codon] = genome_position
                    nucleotide_counter += 1
            total_lengths[gene_name] = nucleotide_counter
        nuc_to_codon = defaultdict(list)
        for feat_name, codons in by_everything.items():
            for codon_index, codon_dict in codons.items():
                codon_obj = Codon(feat_name, codon_index, codon_dict,
                                  self.genes[feat_name].strand)
                assert len(codon_dict) % 3 == 0
                for k, v in codon_dict.items():
                    nuc_to_codon[v].append(codon_obj)
        self.nuc_to_codon = nuc_to_codon
