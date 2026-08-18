# liftAssess Feature Catalog

This document catalogs the capabilities implemented in the current liftAssess
codebase. It is a user-facing inventory, not the scientific specification or
development roadmap.

- [`DESIGN.md`](DESIGN.md) is authoritative for scientific semantics, invariants,
  verdict definitions, coordinate rules, and v1 scope.
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
assessment, and report rendering.

The common form is:

```text
assess-liftover SOURCE_DB TARGET_DB CHR:START-END
```

CLI coordinates are UCSC-style **1-based, inclusive**. Comma-grouped browser-style
coordinates such as `chr16:12,345-12,400` are accepted. liftAssess converts input at
the boundary to its canonical **0-based, half-open** internal representation.

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
- stable candidate IDs tied to the source record that produced each mapping.

Candidate encounter order is preserved for reproducibility, but it is **not** a
candidate rank.

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

Net and reciprocal-best observations retain provenance linking them to their upstream
resources and shared alignment lineage. Multiple observations derived from a shared
upstream alignment are therefore not presented as independent confirmation.

### Evidence-kind implementation matrix

The public `EvidenceKind` enum includes current and planned v1 vocabulary. The exact
implementation status is:

| `EvidenceKind` | Current status | Current assessment role |
| --- | --- | --- |
| `MAPPING_COVERAGE` | Implemented | Verdict-driving |
| `CHAIN_GAPS` | Implemented | Context |
| `CANDIDATE_RANK` | Not yet emitted | None yet |
| `TARGET_PLACEMENT` | Not yet emitted | None yet |
| `CHAIN_SCORE` | Implemented | Context |
| `ALIGNED_BASES` | Implemented | Context |
| `DUPLICATED_QUERY_BASES` | Implemented | Context |
| `NET_CLASSIFICATION` | Implemented | Context |
| `NET_HIERARCHY` | Implemented | Context |
| `RECIPROCAL_BEST_MEMBERSHIP` | Implemented for `COMPARATIVE` | Verdict-driving |
| `FLANKING_GENE_SYNTENY` | Not yet emitted | None yet |

`TARGET_PLACEMENT` refers to future interpretation of target-sequence role/context. It
is distinct from the already implemented target coordinates and target bounding
interval recorded for every candidate.

## Evidence-availability tiers

liftAssess reports evidence availability separately from the assessment verdict.

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

## Deterministic assessment

v1 uses exactly three verdicts:

- `WELL_SUPPORTED`
- `CONTESTED`
- `INDETERMINATE`

The assessor is deterministic and deliberately has **no numeric confidence score**.
It currently drives verdicts from categorical source-locus mapping coverage and, in the
`COMPARATIVE` tier, reciprocal-best membership.

Raw chain score, net `ali`, net `qDup`, net classification, and net hierarchy are
reported as scientifically useful context, but v1 does not convert them into arbitrary
weights or thresholds.

Every assessment records one assessor-owned terminal `decision_reason`. The current
reason vocabulary is:

- `NO_CANDIDATES`
- `LIFTOVER_MULTIPLE_CANDIDATES`
- `LIFTOVER_SINGLE_FULL_MAPPING`
- `LIFTOVER_SINGLE_PARTIAL_MAPPING`
- `COMPARATIVE_MULTIPLE_MATERIAL_CANDIDATES`
- `COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_FULL`
- `COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_NONE`
- `COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_PARTIAL`
- `COMPARATIVE_SOLE_MATERIAL_PARTIAL`
- `COMPARATIVE_NO_MATERIAL_CANDIDATE`

A preferred candidate is reported only when the deterministic rules support one. A
`WELL_SUPPORTED` result still does not establish biological correctness.

## Human-readable reporting

The default CLI output is a concise summary containing:

- the source locus with its display coordinate convention;
- evidence availability;
- the verdict;
- the preferred candidate when one exists, otherwise candidate count/context;
- a plain-language explanation derived from the recorded `decision_reason`; and
- the unconditional biological-correctness caveat.

For `COMPARATIVE` assessments, the summary also states that comparative observations
are not assumed to be independent and points to the detailed outputs for dependency
provenance.

`--details` emits the full human-readable dossier, including:

- exact verdict and `decision_reason`;
- candidate IDs and UCSC chain IDs where applicable;
- exact mapped segments and orientation;
- every evidence observation;
- categorical assessment roles for observations;
- net hierarchy context;
- resource URLs, cache paths, retrieval times, sizes, checksums, and terms context;
- consumed-versus-unconsumed resource status; and
- the complete provenance dependency graph.

## Machine-readable JSON reporting

`--json` renders schema version 1 from the same completed assessment/report model used
by the human-readable reports. It is not a second assessment path.

Schema v1 includes:

- the UCSC database pair;
- source interval and explicit coordinate-system metadata;
- evidence tier, verdict, and required `decision_reason`;
- preferred-candidate reference when one exists;
- supporting and contradicting evidence references;
- ordered candidate records;
- exact mapping segments and target bounding intervals;
- typed evidence values and categorical assessment roles;
- resource consumption, retrieval, checksum, and provider-terms metadata;
- flattened provenance sources and dependency edges; and
- the biological-correctness caveat.

All JSON genomic intervals use canonical **0-based, half-open** coordinates. The
human-facing tier name `LIFTOVER-ONLY` is serialized as the enum token `LIFTOVER_ONLY`
in schema-v1 JSON and is exposed as `EvidenceAvailabilityTier.LIFTOVER_ONLY` in Python.
Status and progress remain on stderr, so stdout can be redirected directly to a JSON
file.

`--json` and `--details` are mutually exclusive.

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

Automatic UCSC discovery is a convenience layer, not a requirement of the assessor
core. Expert callers can supply local resources and provenance directly.

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
- exact cached-byte re-verification before reuse;
- provider-published MD5 verification when an exact checksum entry is available;
- transfer-length validation when an exact identity-encoded length is available;
- atomic publication only after completed resources pass required verification;
- complete-or-error bundle acquisition rather than returning a partially acquired
  comparative bundle;
- convergence of identical bytes retrieved from different URLs onto one content
  artifact; and
- retrieval metadata retained separately from exact byte identity.

The default CLI is cache-first. A complete verified cached bundle is reused without
contacting UCSC.

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
- `assess_candidates()` for already-normalized candidates;
- `assess_ucsc_cached_bundle()` for end-to-end assessment of an acquired bundle;
- `discover_ucsc_resources()` for provider discovery;
- transfer planning and metadata inspection APIs;
- individual-resource and complete-bundle acquisition APIs;
- cached-bundle loading and verification APIs; and
- checksum, SHA-256 identity, and file-provenance helpers.

The public model includes typed assembly, interval, mapping-segment, candidate,
evidence, provenance, assessment, resource, and report objects.

### Public package version

`liftassess.__version__` exposes the installed package version from distribution
metadata.

### Public Python functions

The package currently exports these callable boundaries:

- assessment/orchestration: `assess_candidates()`, `assess_ucsc_cached_bundle()`;
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
  `MappingOrientation`, `NormalizedCandidate`, `Assessment`, `Verdict`,
  `AssessmentDecisionReason`, and `EvidenceAvailabilityTier`;
- evidence: `EvidenceKind`, `EvidenceObservation`, `EvidenceReference`,
  `MappingCoverageStatus`, `MappingCoverageSummary`, `ChainGap`, `ChainGapSummary`,
  `NetHierarchySummary`, `ReciprocalBestMembershipStatus`,
  `ReciprocalBestMembershipSummary`, and `ReciprocalBestResourceCompleteness`;
- provenance/identity: `ProvenanceIdentifier`, `ProvenanceIdentifierKind`,
  `ProvenanceSource`, `ResourceChecksumAlgorithm`, `ResourceChecksumMismatchError`,
  and `ResourceIdentityMismatchError`;
- acquired resources: `CachedResource`, `CachedUCSCResourceBundle`,
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
- the complete deterministic decision-reason vocabulary.

A real `canFam3` → `canFam4` comparative mechanical fixture is also maintained outside
the repository's normal test suite because the provider resources are multi-gigabyte.
It validates evidence extraction and deterministic software behavior, not biological
ground truth.

## Important current limitations and non-features

The following are **not currently implemented**, even when related concepts appear in
the design or model vocabulary:

- candidate-rank evidence with defined locus-scoped semantics;
- target-sequence role/placement interpretation such as primary versus alternative or
  unplaced sequence;
- flanking-gene orthology/synteny evidence;
- freshly computed sequence identity from raw bases;
- a new alignment run such as minimap2 or lastz;
- a numeric composite confidence score;
- machine learning;
- automatic claims of orthology or biological truth;
- scalable batch assessment that reuses resource work across many loci;
- reproducible case manifests or portable resource packets;
- a completed truth-bearing historical-resolution locus;
- a second candidate-generation engine or plugin-management framework;
- hosted infrastructure; or
- general assembly alias/canonicalization resolution beyond the explicit UCSC database
  names/aliases needed at current boundaries.

Current single-locus assessment streams large comparative resources rather than using a
prebuilt genomic index, so large assembly pairs can be slow even for small loci. See
[`PERFORMANCE.md`](PERFORMANCE.md) for measured performance and current optimization
priorities.
