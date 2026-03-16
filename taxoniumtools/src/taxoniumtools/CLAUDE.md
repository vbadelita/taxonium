# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**taxoniumtools** is a Python package that converts phylogenetic trees into Taxonium's JSONL format for web visualization. It is part of the larger Taxonium monorepo (`/home/vlad/code/taxonium/`). Two main CLI entry points:

- `usher_to_taxonium` — converts UShER protobuf (.pb/.pb.gz) trees (legacy path, depends on UShER)
- `newick_to_taxonium` — converts Newick trees with FASTA alignments and GenBank annotations. **This is the actively developed path**, designed to avoid the UShER dependency. Key advantage: supports multi-segment genomes (e.g. influenza with 8 segments) by passing multiple alignment/GenBank pairs.

## Environment & Installation

Uses [pixi](https://pixi.sh) for environment management (conda-forge + bioconda channels). The `pixi.toml` is at `../../pixi.toml` (the taxoniumtools root).

```bash
# From taxoniumtools/ directory:
pixi install            # Set up conda environment
pixi run python -m taxoniumtools  # Run the package

# Or install editably with pip:
pip install -e .
```

Entry points are defined in `setup.cfg`. Version is managed by `setuptools_scm` (generates `_version.py` from git tags).

## Running newick_to_taxonium

The primary workflow uses a JSON config file (see `testtree/config.json` for an influenza example with 8 segments):

```bash
newick_to_taxonium -j testtree/config.json
```

The config specifies `input` (Newick tree), `output` (JSONL), `segments` (array of `{alignment, genbank}` pairs), `metadata`, `metadata_id_column`, `columns`, and `title`. CLI args override config values:

```bash
newick_to_taxonium -i tree.nwk -a aln1.fasta aln2.fasta -g ref1.gb ref2.gb -m metadata.csv -o out.jsonl.gz
```

## Architecture

The data pipeline converts tree → JSONL:

1. **Load tree**: Parse UShER protobuf (`ushertools.py` → `UsherMutationAnnotatedTree`) or Newick (`newick_to_taxonium.py` using treeswift + optionally dendropy for clean parsing)
2. **Annotate mutations**: Nucleotide mutations from the tree; amino acid mutations computed via `core_mutations.py` using GenBank CDS annotations
3. **Ancestral reconstruction** (newick path only): TreeTime `TreeAnc` for inferring ancestral sequences
4. **Layout**: `utils.py` computes x-coordinates (branch length or time) and y-coordinates (leaf ordering), normalizes x to 95th percentile
5. **Serialize**: Sort nodes by y, build mutation index, write gzipped JSONL (first line = header with config/mutations, subsequent lines = one node each)

### Key files

- **`usher_to_taxonium.py`** — CLI + `do_processing()` for UShER protobuf input. Supports Chronumental time trees, shearing, clade annotations.
- **`newick_to_taxonium.py`** — CLI + `do_processing()` for Newick input. The actively developed converter. Supports multi-segment genomes via a JSON config with a `segments` array (each entry has `alignment` + `genbank`). Each segment runs ancestral reconstruction independently and mutations are accumulated per-node across segments. See `testtree/config.json` for an example (influenza with 8 segments: HA, NA, NP, MP, NS, PA, PB1, PB2). Also accepts CLI args directly; CLI overrides config values. Contains `make_root_aa_mutations()` which generates `X→ref_AA` AA mutations for every codon in a segment's genes; these are placed on the root so Taxonium can resolve residues at all positions even for nodes that carry no branch mutations for that segment.
- **`ushertools.py`** — `UsherMutationAnnotatedTree` class: parses protobuf, expands condensed nodes, annotates mutations/clades, reconstructs root sequence, performs AA analysis.
- **`core_mutations.py`** — Frozen dataclasses (`AAMutation`, `NucMutation`, `Gene`, `Codon`), codon translation, `GenbankLoader` for parsing GenBank CDS features into codon lookup maps, `recursive_mutation_analysis()` for tree-wide AA annotation.
- **`utils.py`** — Metadata loading (TSV/CSV via pandas), coordinate computation (`set_x_coords`, `set_terminal_y_coords`, `set_internal_y_coords`), node JSON serialization (`get_node_object`), Chronumental integration.
- **`parsimony_pb2.py`** — Generated protobuf module for UShER's mutation-annotated tree format.

### JSONL output format

Line 1 (header): `{"version", "mutations": [...], "total_nodes", "config": {...}}`
Lines 2+: `{"name", "x_dist", "y", "mutations": [indices], "is_tip", "meta_*", "parent_id", "node_id", "num_tips"}`

Mutations are stored as a global array in the header; nodes reference them by index.

## Key patterns

- Tree traversal uses treeswift's `traverse_preorder()` / `traverse_postorder()` plus custom generators in `core_mutations.py`
- Mutations are frozen dataclasses (hashable) stored in sets for deduplication, then indexed for compact JSONL output
- `orjson` is used instead of `json` for fast binary serialization
- `alive_progress` bars are used throughout for progress reporting
- Metadata columns are prefixed with `meta_` in the output JSON
- `newick_to_taxonium.py` uses dendropy (if available) to clean the Newick before passing to treeswift, stripping internal node labels that can cause parse issues. It writes a temp file for TreeTime input.
- Multi-segment support: each segment runs ancestral reconstruction independently (treetime, parsimony, or tips-only); nuc/AA mutations are appended to nodes across segments. The `chromosome` field on `NucMutation` tracks which segment a mutation belongs to. Not all tips need to be present in every segment's alignment — nodes missing from a segment's alignment receive no branch mutations for that segment and inherit their ancestor's state.
- Parsimony ancestral reconstruction (`parsimony_ancestral_reconstruction()`): post-order pass infers absolute mutation sets per node by intersecting children that have data (children with no alignment match are skipped); pre-order pass converts absolute sets to branch mutations. Nodes with no alignment data (`_parsimony_abs = None`) are skipped in the pre-order pass entirely and inherit their parent's state. Alignment sequences are uppercased before comparison to handle soft-masked (lowercase) bases in alignments, which would otherwise produce invalid codons.
- Root reference state: each segment places `X→ref_AA` AA mutations on the root node via `make_root_aa_mutations()`. This allows Taxonium to display the correct residue for any node that has no downstream AA mutation for a given position, including tips that are absent from a segment's alignment.

## Monorepo context

The parent repo also contains JS/TS packages (use `npm`, never `yarn`):

- `taxonium_component` — React visualization component (Deck.gl). Run `npm run check-types` before committing.
- `taxonium_website` — Vite React app consuming the component
- `taxonium_data_handling` — shared JS library for parsing JSONL, filtering, exporting
- `taxonium_backend` — Node.js server for large trees

Formatting: `yapf` for Python, `prettier` for JS/JSON (configured via pre-commit hooks).

## Changelog

A human-readable log of codebase changes is maintained in `CHANGELOG.md`. **Append an entry every time you make a meaningful change** — include the date, files changed, and a plain-English summary of what was changed and why. Do not describe changes as bug fixes vs features; just describe what the code does now and the reasoning behind it.
