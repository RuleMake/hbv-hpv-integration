# Structural-variant catalogue

These sources establish the paired calls, breakends and donor signatures counted by the spatial analyses. Event eligibility is explicit, so record duplication and biological counting units remain separate.

## Source files

| File | Purpose |
|---|---|
| `build_and_audit_sv_call_catalogue.py` | Builds paired-call, breakend and donor-signature tables and reconciles earlier representations. |
| `sv_event_definitions.py` | Defines source-call identities, chromosome eligibility, interval geometry and class labels. |
| `test_sv_event_definitions_and_parser.py` | Exercises event definitions and source parsing with constructed examples. |

## Required inputs and layout

Per-sample PCAWG rearrangement files; sample-to-donor and aliquot mappings; integration-status metadata; archived variant representations; input manifests.

The expected catalogue directory contains `data/`, `provenance/`, `outputs/` and `runtime_deps/`. Source-call construction precedes representation comparisons and downstream spatial analyses.

See the release input inventory for data sources and `../source_map.csv` for source locations. Configure the required files, dependencies and output paths in a separate work directory before execution. Release verification checks syntax and file hashes without running these upstream analyses.
