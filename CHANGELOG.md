# Changelog

Notable user-visible changes to liftAssess are recorded here.

The project is still alpha scientific software. Pre-release compatibility may change when the change is deliberate, documented, and scientifically justified.

## Unreleased

### Changed

- Replaced the legacy aggregate `WELL_SUPPORTED` / `CONTESTED` / `INDETERMINATE` result model with orthogonal factual states, deterministic factual headlines, bounded interpretation, and explicit scope boundaries.
- Introduced schema-v2 single-locus JSON. The redesigned schema intentionally removes legacy aggregate `verdict`, verdict-derived `decision_reason`, and preferred-candidate semantics rather than silently redefining them.
- Made human reporting progressively disclose material partial coverage, fragmentation, target discontinuity, multiple projections, comparative relationships, and contextual observations without converting them into a composite confidence score.

### Added

- Exact-resource chain indexing for scalable single-locus, reverse, point-context, and batch lookup while retaining the original UCSC chain bytes as scientific provenance.
- Actual reverse-direction chain mapping as a result dimension distinct from UCSC reciprocal-best membership.
- Automatic centered 101-bp forward-chain context for 1-bp point queries, with explicit larger odd-width windows available through `--context-bases`.
- Paired ordinary-filtered versus all-chain inventory comparison and deterministic provenance-aware comparative relationships.
- Indexed BED3+ and 1-based interval-table batch assessment, including exact target collision/overlap relationships and shared COMPARATIVE resource traversal.
- Authoritative UCSC source-sequence preflight for sequence names, aliases, and bounds before mapping is attempted.
- Version-bound UCSC/NCBI target sequence role/context reporting when exact assembly metadata is available.
- Typed UCSC `genomicSuperDups` source/target overlap context with exact provenance and an explicit descriptive-only interpretation boundary.

### Validation

- Added held-out Milestone 23 real-case execution. One case exposed a presentation gap in the shared comparative interpretation; the presentation was corrected and the affected case rerun unchanged. Outside-user/domain review remains pending before the redesigned result language is considered release-ready.

## 0.1.0a1 - 2026-08-17

- First public alpha release.
- Shipped the initial UCSC chain/net/reciprocal-best evidence path, cache/acquisition workflow, provenance model, CLI, and schema-v1 aggregate-verdict reporting.
