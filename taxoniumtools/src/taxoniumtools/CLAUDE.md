# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**taxoniumtools** is a Python package that converts phylogenetic trees into Taxonium's JSONL format for web visualization. It is part of the larger Taxonium monorepo (`/home/vlad/code/taxonium/`). Two main CLI entry points:

- `usher_to_taxonium` — converts UShER protobuf (.pb/.pb.gz) trees
- `newick_to_taxonium` — converts Newick trees (optionally with FASTA alignments and GenBank annotations)

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

## Architecture

The data pipeline converts tree → JSONL:

1. **Load tree**: Parse UShER protobuf (`ushertools.py` → `UsherMutationAnnotatedTree`) or Newick (`newick_to_taxonium.py` using treeswift + optionally dendropy for clean parsing)
2. **Annotate mutations**: Nucleotide mutations from the tree; amino acid mutations computed via `core_mutations.py` using GenBank CDS annotations
3. **Ancestral reconstruction** (newick path only): TreeTime `TreeAnc` for inferring ancestral sequences
4. **Layout**: `utils.py` computes x-coordinates (branch length or time) and y-coordinates (leaf ordering), normalizes x to 95th percentile
5. **Serialize**: Sort nodes by y, build mutation index, write gzipped JSONL (first line = header with config/mutations, subsequent lines = one node each)

### Key files

- **`usher_to_taxonium.py`** — CLI + `do_processing()` for UShER protobuf input. Supports Chronumental time trees, shearing, clade annotations.
- **`newick_to_taxonium.py`** — CLI + `do_processing()` for Newick input. Supports multi-segment genomes (multiple alignment + GenBank files), JSON config files with `segments` array.
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
- The `newick_to_taxonium.py` `main()` supports both CLI args and JSON config files (CLI takes precedence)

## Monorepo context

The parent repo also contains JS/TS packages (use `npm`, never `yarn`):

- `taxonium_component` — React visualization component (Deck.gl). Run `npm run check-types` before committing.
- `taxonium_website` — Vite React app consuming the component
- `taxonium_data_handling` — shared JS library for parsing JSONL, filtering, exporting
- `taxonium_backend` — Node.js server for large trees

Formatting: `yapf` for Python, `prettier` for JS/JSON (configured via pre-commit hooks).
