# Changelog

Human-readable log of changes to the taxoniumtools codebase.

---

## 2026-03-16

### `newick_to_taxonium.py` — parsimony ancestral reconstruction fixes

**Files changed:** `newick_to_taxonium.py`

Three related issues in the parsimony mode meant that coloring by an amino acid residue in a segment that is missing from some tips (e.g. PB2 in an influenza tree where only ~48% of tips have PB2 sequences) showed `X` (unknown) for nearly all nodes rather than the correct reference residue.

**Root reference state was wrong.** The root node receives `X→ref_AA` AA mutations so Taxonium can resolve residues at positions where a node has no downstream branch mutation. Previously, these were generated via `get_mutations(..., disable_check_for_differences=True)`, which builds the initial codon from the reference sequence and therefore always produced same-AA mutations (e.g. `K→K`). Taxonium cannot use `K→K` to establish `K` as the reference state. Fixed by adding `make_root_aa_mutations()`, which directly constructs `AAMutation(initial_aa='X', final_aa=ref_AA)` for every codon in every gene of the segment.

**Nodes missing from a segment's alignment generated spurious reversion mutations.** In the pre-order pass, nodes with `_parsimony_abs = None` (absent from the alignment) were coerced to `{}` via `or {}`, making them look like "reference sequence". When their parent had inferred mutations, the node would generate reversion branch mutations back to reference, resulting in `K→X`-style AA mutations (because the codon table returned `X` for invalid codons — see next point). Fixed by skipping the branch mutation computation entirely for `None` nodes; they correctly inherit their parent's state.

**Soft-masked (lowercase) bases produced invalid codons.** Alignments may contain lowercase nucleotides (e.g. `a`, `g`) in soft-masked regions. These passed the `!= 'N'` and `!= '-'` filters and were stored in `_parsimony_abs`. When assembled into codons, strings like `'gaA'` are not in the codon table, returning `'X'`. Fixed by calling `.upper()` on tip sequences before comparison in both parsimony and tips-only modes.
