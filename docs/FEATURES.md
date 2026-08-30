# liftAssess Feature Catalog

This document catalogs the capabilities implemented in the current liftAssess
codebase. It is a user-facing inventory, not the scientific specification or
development roadmap.

- [`DESIGN.md`](DESIGN.md) is authoritative for scientific semantics, invariants,
  result semantics, coordinate rules, and current scope.
- [`ROADMAP.md`](ROADMAP.md) is authoritative for implementation status and planned
  work.
- [`GETTING_STARTED.md`](GETTING_STARTED.md) explains how to use the current CLI.

The catalog deliberately separates **implemented behavior** from features that are
represented in the design or data model but are not yet implemented. This prevents a
planned evidence kind or future milestone from being mistaken for a current
capability.

## End-to-end assessment

liftAssess can currently run a complete single-locus UCSC assessment from:

1. a source UCSC database identifier;
2. a target UCSC database identifier; and
3. one source genomic interval.

The `assess-liftover` CLI composes resource discovery or verified cache reuse,
resource integrity checks, candidate generation, evidence extraction, deterministic
result-profile derivation, and report rendering. The companion `prepare-liftassess-index`
command explicitly builds a reusable local chain index from an already verified cache bundle;
it never contacts UCSC.

The common single-locus form is:

```text
assess-liftover SOURCE_DB TARGET_DB CHR:START-END
```

CLI coordinates are UCSC-style **1-based, inclusive**. Comma-grouped browser-style
coordinates such as `chr16:12,345-12,400` are accepted. liftAssess converts input at
the boundary to its canonical **0-based, half-open** internal representation.

Prepared indexes also support BED3-or-later and simple interval-table batch assessment:

```text
assess-liftover SOURCE_DB TARGET_DB --bed loci.bed
assess-liftover SOURCE_DB TARGET_DB --interval-table loci.tsv
```

BED coordinates remain native 0-based, half-open. Interval tables require a tab-delimited `sequence`, `start`, `end` header with optional `label` and use 1-based, inclusive coordinates like the single-locus CLI. Both normalize to the same canonical batch model. Batch execution is cache-only and index-only; it never falls back to a whole-chain scan.

## Candidate generation and mapping structure

The built-in UCSC engine implements:

- streaming UCSC chain parsing;
- forward- and reverse-orientation chain geometry;
- source-to-target interval projection;
- split candidate mappings represented as exact aligned segments;
- one target bounding interval for each candidate, kept distinct from the exact
  aligned segments so split mappings are not presented as continuous alignment;
- source-locus mapping coverage as `FULL` or `PARTIAL`;
- exact uncovered source intervals for partial mappings;
- chain-gap geometry through the requested locus;
- raw UCSC chain score as contextual evidence; and
- stable candidate IDs tied to the source record that produced each mapping;
- optional exact-resource chain indexing using 65,536-bp source-coordinate memberships and
  single-copy encounter-order records in independently compressed blocks; and
- transparent indexed lookup when a matching validated index is present; CLI assessment uses
  full-traversal fallback when it is absent or unusable, while lower-level library callers surface
  index-corruption errors for caller-directed recovery.

Candidate encounter order is preserved for reproducibility, but it is **not** a
candidate rank.

## Indexed batch assessment

liftAssess can assess BED3-or-later or simple interval-table record sets through one prepared exact-resource chain index. The current batch layer:

- rejects zero-width/empty BED intervals at input validation and validates 1-based inclusive table coordinates before normalization;
- preserves deterministic row IDs plus optional BED names or interval-table labels;
- preserves rows with zero candidate projections;
- derives exact target collisions separately from overlapping-but-offset projections across distinct input records;
- compares exact mapped target segments after adjacent coverage is canonicalized, never target bounding spans, using a target-local candidate sweep rather than an all-record-pair cross product;
- records the selected chain publication class, exact SHA-256 resource identity, and chain provenance;
- applies automatic 101-bp point context to one-base rows from either batch input form using the same prepared index, with `--context-bases` for a different odd-width point window and no widening of ordinary interval rows;
- keeps input-row and point-context relationships as separate scales, including explicit neighborhood-level target collisions and overlapping-but-offset context projections;
- exposes compact human and schema-v2 JSON batch reports with exact tested context intervals and per-record context run/not-run state;
- for a complete cached `COMPARATIVE` bundle, attaches ordinary-net and reciprocal-best-chain observations to every submitted-row candidate with one shared pass over each resource, without rescanning the indexed all-chain; and
- requires a usable prepared chain index, with no provider access, automatic index build, refresh, or whole-chain fallback.

`LIFTOVER-ONLY` batches remain chain-only. The current COMPARATIVE batch scope does not run the paired filtered-vs-all-chain inventory comparison or categorical comparative relationship classifier used by single-locus results; those dimensions are explicitly reported as not assessed. COMPARATIVE point-context candidates also remain chain-only in this milestone; net/reciprocal-best evidence is attached only to submitted rows rather than silently being inferred at the neighborhood scale. Batch input now receives authoritative cached source-name/bounds/alias preflight before chain assessment; a zero-candidate row therefore means that a valid source interval had no candidate in the selected chain index. One cached source metadata catalog is shared across the batch, and authoritative source sequence length bounds point-context clipping. Target sequence role/context is also cache-only in batch mode: when matching version-bound UCSC/NCBI role metadata is cached it is shared across records, and when it is absent the role dimension is explicitly `UNAVAILABLE` with no inference from sequence naming. Actual reverse mapping is still not re-run per batch row.

## Authoritative assembly metadata and typed context

Before single-locus mapping, liftAssess validates the submitted source sequence and interval against authoritative UCSC assembly metadata. `chromInfo` supplies canonical source names and lengths; exact `chromAlias` correspondences may support a suggestion but are never silently rewritten. Unrecognized sequences and out-of-bounds intervals fail preflight before mapping rather than becoming biological-looking no-projection results.

For target sequence role/context, liftAssess requires an exact versioned assembly binding from the UCSC assembly description before attaching an NCBI Datasets sequence report. Provider-native assembly unit, role, chromosome name, GenBank accession, RefSeq accession, and exact provenance are reported when available. If that binding cannot be established, the role dimension remains `UNAVAILABLE`; names such as `_alt`, `_random`, or `chrUn` are not used as role heuristics.

Single-locus assessment can also attach typed UCSC `genomicSuperDups` context. It checks exact source-query overlap and overlap against exact mapped target segments while preserving the paired interval, strand, provider UID, aligned-base count, fraction matching, resource identity, and provenance. Segmental-duplication overlap is descriptive context only: it does not penalize a projection, establish a mechanism, or prove biological correctness. Missing or unusable optional context degrades to `UNAVAILABLE` without destroying an already-valid primary mapping assessment.

Target-role metadata is cache-only in batch mode. Typed segmental-duplication context is currently a single-locus capability and is not silently implied for batch rows.

## Actual reverse-mapping context

liftAssess can attach actual reverse-direction mapping facts to each forward candidate.
This capability is distinct from UCSC reciprocal-best membership.

- every exact forward target segment is reversed independently; a fragmented candidate's
  bounding span is never queried across an unaligned gap;
- reverse execution consumes only the reverse chain, not net or reciprocal-best
  artifacts;
- the lower-level API can share one verified chain traversal across all segment queries
  or use a validated chain index;
- the automatic CLI is cache-only and index-only for reverse execution and uses the
  same chain publication class (`COMPARATIVE` all-chain or filtered `LIFTOVER-ONLY`) as
  the forward assessment;
- no matching cached reverse chain is `UNAVAILABLE`; a matching chain without a usable
  prepared index is `NOT_RUN`; and
- completed runs preserve original-source coverage, original-versus-elsewhere return
  relationships, exact-geometry reconstruction, resource identity, and provenance.

Normal assessment never downloads reverse resources, builds a reverse index, or starts
an exhaustive reverse-chain fallback implicitly.

## Point-query local context

For a 1-bp source query, the CLI automatically requests a centered 101-bp local-context
assessment when the exact forward chain has a prepared validated index. The context check:

- uses the same forward chain publication class as the point assessment;
- reports the exact tested source window and its actual width;
- clips at indexed source-sequence bounds rather than shifting the point away from center;
- evaluates chain projection count, source coverage, fragmentation, and target discontinuity;
- distinguishes mapped agreement from no projection at either tested scale, and reports newly
  revealed partial coverage, fragmentation, and target discontinuity as distinct facts together with
  whether the result changes materially with query scale; and
- never silently widens again to 1 kb, 10 kb, or another scale.

Automatic context is **forward chain only**. The point/context relationship is derived from
candidate identity, coverage, fragmentation, and target discontinuity; raw chain score is not used
as a rank or threshold. A `COMPARATIVE` point assessment may consume net and reciprocal-best
evidence at the point itself, but those resources are not re-run for the
101-bp neighborhood. The context therefore does not imply neighborhood-scale comparative support.
If no usable matching forward index is available, or the index cannot provide a safe source bound,
the context result is `NOT_RUN`; no extra whole-chain fallback is started.

The 101-bp width is a product default, not a confidence threshold or biological universal. For a
1-bp query, `--context-bases N` requests a different odd-width window explicitly. Ordinary interval
queries are not widened automatically.

## Comparative evidence

When a full comparative resource bundle is available, liftAssess can attach:

- net aligned-base (`ali`) observations;
- net duplicated-query-base (`qDup`) observations;
- net classification such as `top`, `syn`, or `nonSyn` when present in the resource;
- net hierarchy/depth context; and
- locus-specific reciprocal-best membership.

Reciprocal-best membership is categorical:

- `FULL` — all aligned source bases for the candidate are retained;
- `PARTIAL` — some, but not all, aligned source bases are retained; or
- `NONE` — none of the candidate's aligned source bases are retained.

`PARTIAL` and `NONE` are emitted only when liftAssess has been told that the checked
reciprocal-best material is complete for the relevant scope. An arbitrary partial scan
is not treated as evidence of absence.

Net and reciprocal-best observations retain provenance linking them to their exact upstream
resources. In automatic CLI runs, UCSC resources for one source/target database direction are
conservatively grouped as dependent evidence so they are not presented as independent
confirmation. That pair-level grouping does not verify one exact UCSC processing run.

### Evidence-kind implementation matrix

The public `EvidenceKind` enum includes current and planned vocabulary. The exact implementation status is:

| `EvidenceKind` | Current status | Current result use |
| --- | --- | --- |
| `MAPPING_COVERAGE` | Implemented | Factual profile/headline |
| `CHAIN_GAPS` | Implemented | Factual geometry/context |
| `CANDIDATE_RANK` | Not yet emitted | None yet |
| `TARGET_PLACEMENT` | Not yet emitted | None yet |
| `CHAIN_SCORE` | Implemented | Reported context |
| `ALIGNED_BASES` | Implemented | Reported comparative context |
| `DUPLICATED_QUERY_BASES` | Implemented | Reported comparative context |
| `NET_CLASSIFICATION` | Implemented | Reported comparative context |
| `NET_HIERARCHY` | Implemented | Reported comparative context |
| `RECIPROCAL_BEST_MEMBERSHIP` | Implemented for `COMPARATIVE` | Reported comparative relationship |
| `FLANKING_GENE_SYNTENY` | Not yet emitted | None yet |

`TARGET_PLACEMENT` refers to future interpretation of target-sequence role/context. It is distinct from the already implemented target coordinates and target bounding interval recorded for every candidate.

No evidence kind is converted into a numeric score or hidden weighted vote.

## Evidence-availability tiers

liftAssess reports evidence availability separately from the factual mapping result.

### `COMPARATIVE`

A complete v1 comparative bundle was discovered or loaded from cache. The bundle
contains:

- all-chain;
- ordinary classified net;
- syntenic net;
- reciprocal-best chain; and
- reciprocal-best net.

The current assessment engine consumes the all-chain, ordinary net, and reciprocal-best
chain when candidates exist. The syntenic net and reciprocal-best net are retained in
the bundle and report as retrieval context but are not currently parsed as assessment
evidence.

If chain projection produces no candidates, comparative resources are not needlessly
parsed; the report records which cached resources were actually consumed.

### `LIFTOVER-ONLY`

Only a UCSC liftOver chain is available. liftAssess can still generate candidates and
chain-derived evidence, but comparative net and reciprocal-best evidence is not
available.

These tiers describe **what could be checked**, not how confident liftAssess is.

## Deterministic factual result profile

The active result path does **not** assign `WELL_SUPPORTED`, `CONTESTED`, or `INDETERMINATE`, and it does not replace them with another aggregate verdict.

A dedicated derived `ResultProfile` sits between scientific candidate/evidence data and both renderers. It deterministically records currently available dimensions including:

- projection count;
- source coverage with exact numerator/denominator;
- mapped-segment count and uncovered source intervals;
- target bounding span and target gaps;
- orientation;
- maximum candidate source coverage when multiple projections exist;
- point-query local-context state, exact tested window, and factual point/context relationships;
- evidence tier and consumed resource roles; and
- explicit not-yet-assessed/not-run boundaries for later result dimensions.

The profile also derives a factual headline such as `NO CHAIN PROJECTION`, `ONE COMPLETE CHAIN PROJECTION`, `PARTIAL SOURCE COVERAGE`, `COMPLETE BUT DISCONTINUOUS PROJECTION`, or `MULTIPLE CHAIN PROJECTIONS`, plus a bounded deterministic interpretation.

Raw chain score, net `ali`, net `qDup`, net classification, net hierarchy, and reciprocal-best membership remain exact reported observations. They are not combined through arbitrary weights or thresholds, and shared UCSC provenance is preserved.

## Human-readable reporting

The default CLI output is a facts-first progressive summary containing:

- the deterministic factual headline;
- the source interval and coordinate convention;
- exact source coverage and target geometry needed to understand the result;
- evidence tier and consumed-resource context;
- bounded deterministic interpretation;
- relevant scope/identity boundaries; and
- the unconditional biological-correctness caveat.

Uncomplicated one-complete-projection cases stay compact. The current first slice expands automatically for partial source coverage, fragmented or target-discontinuous geometry, and multiple projections.

For `COMPARATIVE` results, the summary states that UCSC-derived comparative observations are conservatively treated as dependent and that exact shared processing-run provenance is not verified.

For single-locus assessment, `--details` emits the complete currently available factual dossier, including:

- every result-profile field and explicit scope boundary;
- candidate IDs and UCSC chain IDs where applicable;
- exact mapped segments, uncovered source intervals, target gaps, and orientation;
- every evidence observation;
- net hierarchy and reciprocal-best context;
- resource URLs, cache paths, retrieval times, sizes, checksums, and terms context;
- consumed-versus-unconsumed resource status; and
- the complete provenance dependency graph.

## Machine-readable JSON reporting

`--json` renders schema version 2. Single-locus JSON comes from the same completed report and derived `ResultProfile` used by the human-readable reports; batch JSON uses the separate `liftassess.ucsc_batch_result` report type from the indexed batch result. Neither is a second candidate-generation path.

Schema v2 includes:

- the UCSC database pair;
- source interval and explicit coordinate-system metadata;
- factual headline, bounded interpretation, projection count, and exact source-coverage summary;
- complete candidate result profiles;
- ordered candidate records with exact mapping segments and target bounding intervals;
- typed evidence values without verdict-derived supporting/contradicting roles;
- evidence tier and exact resource-consumption metadata;
- flattened provenance sources and dependency edges;
- explicit scope states for result dimensions not yet assessed; and
- the biological-correctness caveat.

Schema v2 intentionally does **not** preserve the legacy aggregate `verdict`, verdict-derived `decision_reason`, or preferred-candidate field.

All JSON genomic intervals use canonical **0-based, half-open** coordinates. The human-facing tier name `LIFTOVER-ONLY` is serialized as the enum token `LIFTOVER_ONLY` in JSON and is exposed as `EvidenceAvailabilityTier.LIFTOVER_ONLY` in Python. Status and progress remain on stderr, so stdout can be redirected directly to a JSON file.

`--json` and `--details` are mutually exclusive. `--details` is not yet implemented for batch mode.

## UCSC resource discovery

Automatic resource discovery:

- checks UCSC's published directory listings instead of assuming that a constructed URL
  exists;
- prefers a complete comparative bundle when all required resources are published;
- falls back to a verified liftOver-only chain when a complete comparative set is not
  available;
- checks both observed UCSC pair-directory layouts for directional reciprocal-best
  resources, including asymmetric publication layouts; and
- distinguishes a genuine absence from a provider/network failure so a transient
  transport error cannot silently downgrade evidence availability.

Automatic UCSC discovery is a convenience layer, not a requirement of the candidate/
evidence engine. Expert callers can supply local resources and provenance directly.

## Resource terms and transfer planning

Before automatic UCSC acquisition, liftAssess:

- classifies the provider terms applicable to each planned resource;
- distinguishes dedicated `liftOver/*.over.chain.gz` resources from comparative
  resources because their terms can differ;
- requires explicit acknowledgement of the displayed UCSC terms before provider
  metadata inspection or acquisition;
- performs body-free HTTP HEAD inspection after terms acknowledgement;
- preserves provider-advertised transfer metadata such as `Content-Length`,
  `Accept-Ranges`, `Last-Modified`, `ETag`, and `Content-Encoding` when present;
- displays the complete bundle transfer plan and cache destination; and
- requires a separate explicit transfer-plan acknowledgement before acquisition.

For non-interactive workflows, `--acknowledge-ucsc-terms` and
`--accept-transfer-plan` provide the same explicit acknowledgements without prompts.

## Acquisition, integrity, and cache

The cache/acquisition layer implements:

- caller-selected external cache directories;
- platform-appropriate default user-cache locations;
- SHA-256 content-addressed artifact storage;
- URL-to-artifact cache indexes;
- exact cached-byte re-verification before direct provider-artifact consumption;
- validated derived-chain identity can replace a redundant reread of the original chain during indexed assessment; normal queries authenticate a compact lookup catalog plus queried bin metadata and selected record blocks, while index build/rebuild still verifies the source bytes and records a full database checksum for deep verification;
- provider-published MD5 verification when an exact checksum entry is available;
- transfer-length validation when an exact identity-encoded length is available;
- atomic publication only after completed resources pass required verification;
- complete-or-error bundle acquisition rather than returning a partially acquired
  comparative bundle;
- convergence of identical bytes retrieved from different URLs onto one content
  artifact; and
- retrieval metadata retained separately from exact byte identity.

The default CLI is cache-first. A complete cached bundle is reused without contacting UCSC.
`--evidence-tier COMPARATIVE` or `--evidence-tier LIFTOVER-ONLY` requests one exact
publication class and disables automatic tier fallback for cache selection and discovery.
Cached provider artifacts other than the indexed source chain retain normal SHA-256 verification;
when an exact-resource chain index is already validated, that derived artifact carries the
source-chain identity for indexed lookup. Query-relevant bin membership/record-locator rows and
selected compressed blocks are verified without requiring full reads of either the unused original
chain or the large SQLite lookup database.

`--offline` guarantees zero provider access and fails if no complete verified bundle is
available locally.

`--refresh` deliberately bypasses cache-first reuse and checks/reacquires current
provider resources.

## Restart-safe resumable transfers

When UCSC provides the required transport guarantees, liftAssess can resume interrupted
HTTPS downloads. Persistent resume is enabled only when the resource has:

- an exact provider checksum;
- an exact identity-encoded size;
- advertised byte-range support; and
- a strong ETag.

Retained partials are bound to the exact URL, total size, and validator. Resume uses
`Range` plus `If-Range`. A changed validator, an unexpected full `200` response, or an
invalid `Content-Range` causes a clean restart instead of appending potentially
incompatible bytes.

Completed resumable data is copied to a private snapshot and independently verified
before atomic publication into the immutable content-addressed cache.

## Measured progress reporting

Interactive terminal runs provide measured byte progress for three different kinds of
work:

- resource transfer;
- cached-bundle SHA-256 verification; and
- assessment-time resource reading.

Progress is based on actual bytes, not an estimated biological/algorithmic completion
percentage. When the provider does not supply a trustworthy total size, liftAssess
shows byte counts without inventing a percentage.

Progress displays are terminal-only and are suppressed by `--quiet`. Redirected report
output therefore remains clean.

## Local-resource and programmatic APIs

The public Python package exposes expert-level boundaries in addition to the CLI.
Current capabilities include:

- `iter_chain_file()` and `iter_net_file()` for plain-text or gzip-compressed local
  UCSC resources; gzip is detected from the file bytes rather than the filename suffix,
  so content-addressed/renamed cached artifacts still parse correctly;
- `build_ucsc_candidates()` for parsed UCSC records;
- `build_ucsc_candidates_from_files()` for explicitly supplied local resources and
  provenance;
- `build_ucsc_candidates_from_cached_bundle()` for verified cached bundles;
- `build_result_profile()` for already-normalized candidates plus their evidence;
- `assess_ucsc_cached_bundle()` for end-to-end assessment of an acquired bundle;
- `discover_ucsc_resources()` for provider discovery;
- transfer planning and metadata inspection APIs;
- individual-resource and complete-bundle acquisition APIs;
- cached-bundle loading and verification APIs; and
- checksum, SHA-256 identity, and file-provenance helpers.

The public model includes typed assembly, interval, mapping-segment, candidate, evidence,
provenance, result-profile, resource, and report objects.

### Public package version

`liftassess.__version__` exposes the installed package version from distribution
metadata.

### Public Python functions

The package currently exports these callable boundaries:

- result/orchestration: `build_result_profile()`, `assess_ucsc_cached_bundle()`;
- candidate/file handling: `build_ucsc_candidates()`,
  `build_ucsc_candidates_from_files()`, `build_ucsc_candidates_from_cached_bundle()`,
  `iter_chain_file()`, `iter_net_file()`;
- resource discovery/planning/cache: `discover_ucsc_resources()`,
  `plan_ucsc_bundle_acquisition()`, `inspect_ucsc_resource()`,
  `inspect_ucsc_bundle_transfer_plan()`, `acquire_ucsc_resource()`,
  `acquire_ucsc_resource_bundle()`, `load_cached_ucsc_resource_bundle()`, and
  `ucsc_resource_terms()`; and
- identity/provenance: `compute_resource_checksum()`, `verify_resource_checksum()`,
  `sha256_identifier_for_file()`, and `provenance_source_for_file()`.

### Public Python types

The package currently exports these model/resource types and errors:

- core models: `AssemblyIdentifier`, `GenomicInterval`, `MappingSegment`,
  `MappingOrientation`, `NormalizedCandidate`, and `EvidenceAvailabilityTier`;
- result profile: `ResultProfile`, `CandidateResultProfile`,
  `CandidateReverseMappingProfile`, `ResultScopeProfile`, `FactualHeadline`,
  `ProjectionCountState`, `SourceCoverageState`, and explicit state enums for later
  dimensions;
- actual reverse mapping: `CandidateReverseMappingResult`, `ReverseSegmentResult`,
  `ReverseCheckState`, `ReverseRelationshipState`,
  `ReverseOriginalSourceCoverageState`,
  `build_reverse_mapping_results_from_cached_bundle()`,
  `build_reverse_mapping_results_from_cached_chain()`,
  `attach_reverse_mapping_context()`, and `attach_reverse_mapping_results()`; the CLI
  runs this automatically only when a matching reverse-direction chain and prepared
  reverse chain index are already cached;
- evidence: `EvidenceKind`, `EvidenceObservation`,
  `MappingCoverageStatus`, `MappingCoverageSummary`, `ChainGap`, `ChainGapSummary`,
  `NetHierarchySummary`, `ReciprocalBestMembershipStatus`,
  `ReciprocalBestMembershipSummary`, and `ReciprocalBestResourceCompleteness`;
- provenance/identity: `ProvenanceIdentifier`, `ProvenanceIdentifierKind`,
  `ProvenanceSource`, `ResourceChecksumAlgorithm`, `ResourceChecksumMismatchError`,
  and `ResourceIdentityMismatchError`;
- acquired resources: `CachedResource`, `CachedUCSCChainResource`,
  `CachedUCSCResourceBundle`,
  `ProviderChecksum`, `UCSCResourceBundle`, `UCSCResourceClass`,
  `UCSCBundleResourceRole`, `UCSCResourceTerms`, `UCSCRemoteResourceMetadata`,
  `UCSCBundleAcquisitionItem`, `UCSCBundleAcquisitionPlan`,
  `UCSCBundleTransferInspectionItem`, `UCSCBundleTransferInspection`,
  `UCSCResourceDiscoveryError`, `UCSCResourceAcquisitionError`,
  `UCSCResourceTermsAcknowledgementRequired`, and
  `UCSCBundleAcquisitionPlanAcknowledgementRequired`; and
- reporting/orchestration: `UCSCAssessmentResource`, `UCSCAssessmentReport`, and
  `ResourceReadProgressCallback`.

There is intentionally only one concrete candidate-generation engine in v1. The clean
normalized-candidate boundary is **not** a plugin registry or engine auto-discovery
system.

## Validation and auditability

The current repository includes automated coverage for, among other behaviors:

- forward and reverse mappings;
- split mappings and gaps;
- repeated net chain IDs;
- provenance dependency diamonds;
- reciprocal-best subsetting and completeness;
- resource discovery and asymmetric reciprocal-best layouts;
- cache integrity and corrupt-cache recovery;
- resumable-transfer safety;
- cache-first, offline, refresh, and acknowledgement CLI paths;
- summary, detailed, and JSON reporting semantics; and
- deterministic factual result-profile, headline, reverse-context, and scope semantics.

A real `canFam3` → `canFam4` comparative mechanical fixture is also maintained outside
the repository's normal test suite because the provider resources are multi-gigabyte.
It validates evidence extraction and deterministic software behavior, not biological
ground truth.

## Important current limitations and non-features

The following are **not currently implemented**, even when related concepts appear in
the design or model vocabulary:

- candidate-rank evidence with defined locus-scoped semantics;
- flanking-gene orthology/synteny evidence;
- freshly computed sequence identity from raw bases;
- a new alignment run such as minimap2 or lastz;
- a numeric composite confidence score;
- machine learning;
- automatic claims of orthology or biological truth;
- reverse evidence across many batch loci;
- implicit reverse-direction acquisition or automatic reverse-index construction during
  ordinary assessment;
- reproducible case manifests or portable resource packets;
- a completed truth-bearing historical-resolution locus;
- a second candidate-generation engine or plugin-management framework;
- hosted infrastructure; or
- general assembly alias/canonicalization resolution beyond the explicit UCSC database
  names/aliases needed at current boundaries.

Large chain resources can be prepared once with `prepare-liftassess-index`; matching later
assessments reuse the validated region-addressable chain index automatically. Without one,
assessment retains the original full-chain traversal. Net and reciprocal-best access remain separate
comparative-resource costs where profiling shows they are material. See
[`PERFORMANCE.md`](PERFORMANCE.md) for measured performance and current optimization priorities.
