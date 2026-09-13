# Changelog

Notable user-visible changes to liftAssess are recorded here.

The project is still alpha scientific software. Pre-release compatibility may change when the change is deliberate, documented, and scientifically justified.

## Unreleased

### Changed

- Replaced the legacy aggregate result model with orthogonal factual states, deterministic factual headlines, bounded interpretation, and explicit scope boundaries.
- Introduced schema-v2 single-locus JSON. The redesigned schema intentionally removes obsolete aggregate-result and mapping-selection fields rather than silently redefining them.
- Made human reporting progressively disclose material partial coverage, fragmentation, target discontinuity, multiple mappings, comparative relationships, and genomic context without converting them into a composite confidence score.

### Added

- Exact-resource chain indexing for scalable single-locus, reverse, flanking-interval, and batch lookup while retaining the original UCSC chain bytes as scientific provenance.
- Actual reverse-direction chain mapping as a result dimension distinct from UCSC reciprocal-best membership.
- Automatic centered 101-bp flanking-interval checks for 1-bp point queries, with explicit larger odd-width windows available through `--context-bases`.
- Paired standard liftOver versus all-chain inventory comparison and deterministic provenance-aware comparative relationships.
- Indexed BED3+ and 1-based interval-table batch assessment, including exact target collision/overlap relationships and shared COMPARATIVE resource traversal.
- Authoritative UCSC source-sequence validation for sequence names, aliases, and bounds before mapping is attempted.
- Version-bound UCSC/NCBI target sequence-role reporting when exact assembly metadata is available.
- UCSC Segmental Duplications source/target overlap context with exact provenance and an explicit descriptive-only interpretation boundary.

### Validation

- Completed the Milestone 23 held-out language/usability gate across all five pre-registered cases.
- Incorporated post-M23 outside-user feedback that unusual findings should explain why they matter and what a biologist can do next. Non-reciprocal reverse liftOver, multiple mappings, comparative evidence, and Segmental Duplications context now carry bounded interpretation and actionable follow-up without asserting a biological cause or preferred locus.
- Added direct UCSC Genome Browser navigation for unusual single-locus results and documented the scientific literature supporting the interpretation of non-reciprocal genome-build mappings.

## 0.1.0a1 - 2026-08-17

- First public alpha release.
- Shipped the initial UCSC chain/net/reciprocal-best evidence path, cache/acquisition workflow, provenance model, CLI, and schema-v1 aggregate-result reporting.
