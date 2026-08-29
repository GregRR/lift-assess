# liftAssess

**liftAssess** is a Python command-line tool and library for assessing ambiguous genomic coordinate liftOver mappings between genome assemblies. It reports factual projection geometry, evidence availability, provenance, and bounded interpretation without collapsing those facts into a single confidence verdict.

> **Status:** Early public alpha development. Core candidate generation, evidence extraction, cached-bundle orchestration, a derived factual result profile, progressive human reporting, and schema-v2 JSON reporting are implemented and tested. The current redesign intentionally breaks the original alpha result schema; this remains early scientific software.

**New to liftAssess?** Start with [`GETTING_STARTED.md`](docs/GETTING_STARTED.md). See [`FEATURES.md`](docs/FEATURES.md) for the complete catalog of implemented capabilities, expert APIs, and current limitations.

## Why liftAssess exists

Coordinate liftOver tools answer an important question: *where can this interval map in another assembly?* They do not, by themselves, explain exact source coverage, split or discontinuous geometry, competing projections, which comparative resources were examined, or how those observations depend on one another.

Use liftAssess when genome assembly or genome build coordinate conversion produces one-to-many mappings, split mappings, surprising target spans, or other cases where a coordinate alone does not tell the whole story.

liftAssess is intended to sit downstream of or alongside coordinate-conversion tools and answer a different question:

> **What physically happened to this interval in the consumed mapping resources, what evidence was examined, and what does that evidence not establish?**

liftAssess is an **assessor, not a resolver**. It does not claim that a locus is biologically “correct,” and it does not produce a composite numeric confidence score.

## Result model

The current result model is facts-first rather than verdict-first. It derives orthogonal factual states from the scientific report, including projection count, exact source coverage, mapped-segment geometry, target discontinuity, orientation, evidence availability, resource consumption, and provenance.

A deterministic factual headline summarizes the dominant mapping event, for example `ONE COMPLETE CHAIN PROJECTION`, `PARTIAL SOURCE COVERAGE`, or `MULTIPLE CHAIN PROJECTIONS`. A bounded interpretation explains that event without turning it into a claim of biological truth or candidate correctness.

The result profile also carries actual reverse-mapping context when a matching reverse-direction chain and its prepared index are already available in the local cache. Reverse mapping is chain-only, is reported separately from reciprocal-best membership, and preserves exact original-source return geometry. The automatic CLI uses the same chain publication class as the forward assessment. If no matching reverse chain is cached, the check is `UNAVAILABLE`; if the chain exists but scalable indexed access is not prepared or is unusable, the check is `NOT_RUN`. An explicit forward `--refresh` also leaves reverse mapping `NOT_RUN` rather than mixing freshly reacquired forward resources with an unrefreshed reverse chain. Normal assessment does not silently contact UCSC, refresh reverse resources, build a reverse index, or start another exhaustive reverse-chain scan.

For 1-bp point queries, the profile additionally carries automatic local-context chain evidence when the prepared forward chain index is available. The default window is centered at 101 bp (boundary-clipped when needed), its exact tested interval is reported, and factual states distinguish mapped agreement, no projection at either tested scale, newly revealed partial coverage, fragmentation, target discontinuity, or other query-scale change. `--context-bases N` requests another odd-width point window explicitly. This context reuses the same forward chain publication class and does not re-run net/reciprocal-best evidence or fall back to another full chain scan. BED3+ batch input is also available through prepared chain indexes; COMPARATIVE batches attach ordinary-net and reciprocal-best-chain observations to submitted rows with one shared pass over each resource, while point-context remains chain-only. Target sequence role and typed external context remain later dimensions.

## Scientific principles

A few rules are central to the project:

- **Evidence provenance is part of the result.** Every observation records where it came from.
- **Dependent observations are not treated as independent confirmation.** For example, UCSC chain, net, and reciprocal-best evidence derived from the same upstream alignment remain linked through provenance.
- **Evidence availability is separate from factual result state.** `COMPARATIVE` and `LIFTOVER-ONLY` say what resource families were available, not whether a mapping is good, bad, correct, or confident.
- **Interpretation stays close to the evidence.** liftAssess reports what the data support without automatically promoting an evidence pattern to a specific biological mechanism such as paralogy or pseudogene status.
- **Internal coordinates are 0-based, half-open.** Input and output boundaries must state or explicitly convert coordinate conventions rather than guessing them.
- **Same-species assembly comparisons are the v1 operational envelope.** Differences between assemblies from different individuals may represent real structural variation and are not automatically errors.

## Current implementation

The current development code includes:

- typed assembly, interval, candidate, evidence, result-profile, and provenance models;
- a minimal streaming UCSC chain reader;
- forward- and reverse-strand chain geometry with 0-based half-open internal coordinates;
- chain-backed candidate projection, including split mappings;
- actual chain-only reverse-mapping context from an exact-class cached reverse chain, with explicit `RUN`, `UNAVAILABLE`, and `NOT_RUN` states and no implicit reverse acquisition or index build;
- automatic indexed 101-bp local chain context for 1-bp queries, with exact tested geometry, explicit query-scale findings, and `--context-bases` for a different odd-width window;
- mapping-coverage and chain-gap evidence;
- a minimal streaming UCSC net reader with hierarchy preservation;
- net evidence for aligned bases (`ali`), duplicated query bases (`qDup`), net classification, and hierarchy;
- dependency-aware provenance linking chain- and net-derived observations;
- reciprocal-best membership evidence with `FULL`, `PARTIAL`, and `NONE` states;
- explicit completeness requirements for reciprocal-best absence/partial evidence;
- UCSC resource discovery and evidence-availability tier detection, including verified reciprocal-best lookup across asymmetric pair-directory publication layouts;
- single-pass UCSC engine orchestration that generates candidates and attaches available net and
  reciprocal-best evidence without rescanning whole comparative resources per candidate;
- streaming local-file adapters for plain-text or gzip-compressed chain/net resources, including
  direct parser-to-engine orchestration for user-supplied or cached files;
- content-addressed SHA-256 provenance for exact local resource-file bytes, with separate optional
  verification of provider-published MD5 or SHA-256 checksums;
- UCSC acquisition into a caller-selected external cache, with explicit provider-terms acknowledgement
  before network access, exact-filename MD5 and HTTP-length verification when available, atomic
  content-addressed storage, retrieval metadata, offline cache reuse, and refresh;
- cache-first CLI reuse of complete verified bundles with zero provider access, plus explicit `--offline` and `--refresh` modes so offline analysis and freshness checks are never conflated;
- interactive assessment progress based on exact compressed bytes consumed during SHA-256-verified parsing, with Chain/Net/Reciprocal-best progress bars and numeric byte percentages;
- interactive cache-verification progress based on exact artifact bytes hashed across the required cached bundle, with one aggregate row that reaches 100% only after all required SHA-256 checks pass;
- interactive resource-transfer progress for fresh and resumable UCSC acquisition, using measured exact bytes, retained-prefix-aware resume state, explicit cache-hit labeling, and byte-only display rather than invented percentages when the provider size is unknown;
- explicit bundle transfer planning and complete-or-error bundle acquisition for discovered
  `COMPARATIVE` and `LIFTOVER-ONLY` resource sets, with a separate transfer-plan acknowledgement before
  any planned resource acquisition begins;
- terms-gated, body-free remote metadata inspection using HTTP HEAD with identity encoding requested,
  preserving provider-advertised `Content-Length`, `Accept-Ranges`, `Last-Modified`, `ETag`, and
  `Content-Encoding` without transferring resource bodies;
- restart-safe resumable HTTPS acquisition when UCSC publishes an exact checksum and advertises an identity-encoded size, byte-range
  support, and a strong ETag; retained partials are bound to that exact URL/size/validator and resumed with
  `Range` + `If-Range`, while contract mismatches restart fresh rather than splicing representations;
- regression coverage for forward/reverse mappings, split mappings, gaps, repeated net chain IDs, provenance diamonds, reciprocal-best subsetting, and resource-discovery failure modes.
- a real `canFam3`→`canFam4` comparative mechanical fixture, replayed through the production cached-bundle path from exact externally cached UCSC resources without committing provider data to the repository.
- a dedicated derived result-profile layer that validates durable candidate/evidence invariants and deterministically summarizes projection count, coverage, geometry, orientation, evidence availability, and scope boundaries without an aggregate verdict.
- progressive human rendering that stays compact for uncomplicated results and expands for currently detectable partial, fragmented/discontinuous, or multiple-projection states.
- schema-v2 JSON reporting from the same result profile, retaining exact candidates, evidence, resources, and provenance while removing the legacy `verdict`, verdict-derived `decision_reason`, and preferred-candidate fields.
- cached-bundle orchestration that connects acquired UCSC resources to candidate/evidence generation and result-profile derivation while preserving evidence tier, resource-consumption metadata, retrieval context, and shared provenance.

## Not implemented yet

The project now implements the common-case CLI, concise summary, human-readable detail, schema-versioned JSON reporting, and measured cache-verification/assessment/transfer progress paths. The post-hardening real comparative CLI and independent mechanical-fixture verifier have both completed successfully against the established external `canFam3`→`canFam4` cache. Planned work beyond the first public alpha includes:

- a future truth-bearing historical-resolution locus for the planned `canFam3.1`→`canFam6` sanity-check pedigree;
- optional flanking-gene orthology/synteny evidence;
- defensible candidate-rank and target-sequence context backed by explicit evidence/metadata sources rather than naming heuristics;
- reverse context across batch loci;
- reproducible case manifests and, where redistribution terms permit, portable resource packets.

`0.1.0a1` is the first public alpha. It should be treated as early scientific software under active development rather than as a mature or stable analysis platform.

## Evidence-availability tiers

liftAssess distinguishes **which evidence resources are available** from the factual mapping result.

- `COMPARATIVE` — a full comparative resource set is available, including chain/net context and reciprocal-best resources needed by the v1 UCSC evidence path.
- `LIFTOVER-ONLY` — only a liftOver chain resource is available, so candidate generation and chain-derived evidence can still proceed but comparative evidence is unavailable.

These are evidence-availability tiers, **not confidence tiers**.

## Architecture

```text
Candidate/evidence engine
        |
        v
NormalizedCandidate[] + evidence/provenance
        |
        v
UCSC scientific report
        |
        v
Derived ResultProfile
(factual states + headline + bounded interpretation)
        |                 |
        v                 v
Human renderer        Schema-v2 JSON
```

Candidate generation and evidence extraction remain separate from result-language synthesis. The derived result profile is the single boundary used by both human and machine renderers, so renderers do not independently rediscover result semantics.

The candidate/evidence engine is designed around normalized candidates plus explicit provenance. The current implementation intentionally has **one** concrete candidate-generation engine: an internal minimal UCSC chain/net reader. liftAssess does not include a plugin registry, engine auto-discovery, or other speculative plugin infrastructure.

Chains and nets have different responsibilities:

- **chains generate candidate mappings**;
- **nets annotate those candidates with comparative context**.

Net availability is therefore not required for basic chain-backed candidate generation.

## Installation

liftAssess requires Python 3.11 or newer. For the public alpha release, install the pre-release package from PyPI with:

```bash
python -m pip install --pre liftassess
```

Verify the installed command without contacting any external provider:

```bash
assess-liftover --help
```

For development from a source checkout, use `uv sync` as described under [Development setup](#development-setup).

## CLI workflow

After installation, the common-case command is:

```text
assess-liftover canFam3 canFam4 chrUn_JH373233:1845736-1845835
```

CLI loci use the familiar UCSC-style 1-based, inclusive display convention and are converted immediately to liftAssess's canonical 0-based, half-open internal representation.

Before resource acquisition, the command displays the applicable UCSC terms and requires explicit acknowledgement. It then performs body-free HEAD inspection of the exact transfer plan, displays provider-advertised resource sizes and the cache destination, and requires a separate transfer-plan acknowledgement before downloading or verifying cached resources. The displayed size is the provider resource-set size, not a promise that all bytes will be transferred: verified cache hits can avoid body transfer unless `--refresh` is used. `--acknowledge-ucsc-terms` and `--accept-transfer-plan` provide explicit non-interactive acknowledgements; `--cache-dir`, `--refresh`, `--offline`, `--details`, `--json`, `--quiet`, `--context-bases`, `--bed PATH`, and `--interval-table PATH` control cache placement, provider access, report format, progress output, local-context width for one-base queries, and batch input.

The default cache is the platform user cache (`~/Library/Caches/liftassess` on macOS, `%LOCALAPPDATA%\liftassess\Cache` on Windows, and `$XDG_CACHE_HOME/liftassess` or `~/.cache/liftassess` on other platforms). Before a single-locus scientific assessment begins, liftAssess now validates the submitted source sequence and interval against the UCSC database's authoritative `chromInfo` metadata. When available, exact `chromAlias` correspondences are used only to suggest the canonical UCSC name; submitted aliases are not silently rewritten. These small database-table artifacts use the same content-addressed cache and exact SHA-256 provenance as other external resources but are not mapping evidence. A normal online run discovers/acquires them when missing; `--offline` requires a verified cached `chromInfo` artifact before assessment can proceed.

The documented `canFam3`→`canFam4` example is the real comparative mechanical fixture and requires a complete approximately 2.50 GiB UCSC comparative bundle; initial acquisition and an unindexed assessment can therefore take substantial time depending on storage and CPU performance. Once that bundle is cached, `--offline` reuses it without contacting UCSC. Unindexed runs re-verify the source chain directly; when a validated exact-resource chain index exists, liftAssess verifies its compact lookup-integrity catalog, the queried bin metadata, and selected compressed record blocks instead of rereading either the unused multi-gigabyte source chain or the complete SQLite lookup database, while the other cached bundle artifacts retain their normal direct integrity checks.

For repeated work on a large assembly pair, an optional one-time preparation command builds a reusable chain index from the already verified local cache and **never contacts UCSC**:

```text
prepare-liftassess-index canFam3 canFam4
```

When both UCSC publication classes exist,
`assess-liftover --evidence-tier LIFTOVER-ONLY` can explicitly acquire/use the ordinary
filtered liftOver chain instead of the default COMPARATIVE-preferred selection. This is
useful when preparing the filtered-chain side of a paired comparison; provider terms and
transfer-plan acknowledgement remain unchanged.

Index construction parses the complete chain once and can take many minutes and several GiB of additional local cache space for very large resources. Later `assess-liftover` CLI runs automatically reuse the exact-resource index when present. The CLI falls back to the original verified full traversal when the index is absent or unusable; lower-level library calls surface `ChainIndexCorruptionError` so callers can choose their own recovery policy. The index is a derived acceleration artifact; scientific provenance continues to identify the original UCSC chain bytes. See [`docs/PERFORMANCE.md`](docs/PERFORMANCE.md) for measured examples and scope.

Batch assessment uses that same prepared exact-resource index and is intentionally stricter than the single-locus path. BED3-or-later input keeps native 0-based, half-open coordinates:

```text
assess-liftover SOURCE_DB TARGET_DB --bed loci.bed --json > batch.json
```

A simple tab-delimited interval table is also accepted with a required `sequence\tstart\tend` header and optional fourth `label` column. Table coordinates use the same 1-based, inclusive convention as the single-locus CLI and are normalized immediately to canonical 0-based, half-open intervals:

```text
sequence	start	end	label
chr1	101	200	region-a
chr2	500	500	point-b
```

```text
assess-liftover SOURCE_DB TARGET_DB --interval-table loci.tsv --json > batch.json
```

Both batch input forms support `-` for stdin. Batch output uses canonical 0-based, half-open intervals. Batch mode is cache-only and index-required in this milestone: it does not contact UCSC, refresh resources, build an index, or fall back to a whole-chain traversal. If the preferred exact publication class lacks a prepared index, prepare it explicitly with `prepare-liftassess-index`. `LIFTOVER-ONLY` batches remain chain-only. When a complete cached `COMPARATIVE` bundle is available, submitted-row candidates use the prepared all-chain index and then receive ordinary-net and reciprocal-best-chain observations from one shared full pass over each of those two comparatively small resources; the all-chain itself is not rescanned. The current batch COMPARATIVE scope does not run the paired filtered-vs-all-chain inventory comparison or categorical comparative relationship classifier used by single-locus results; those dimensions are reported as not assessed. Exact target collisions remain separate from overlapping-but-offset projections using exact mapped target segments. One-base batch rows automatically receive the same 101-bp point-context check as single-locus points through the prepared chain index; those context candidates remain chain-only, exact context-scale collisions are reported separately as neighborhood-level target collisions, and `--context-bases` selects another odd-width point window. Ordinary interval rows are not widened. Authoritative assembly-sequence name/alias preflight and reverse batch evidence are not yet assessed. Use `--evidence-tier` when the filtered `LIFTOVER-ONLY` publication class is required.

The default command emits a concise facts-first summary headed by a deterministic mapping headline. Uncomplicated results stay compact; currently detectable partial coverage, fragmented/discontinuous geometry, and multiple projections expand automatically with the material geometry. `COMPARATIVE` summaries also state that UCSC-derived observations are conservatively treated as dependent, not independent votes, and exact shared processing-run provenance is not verified. For single-locus runs, `--details` emits the complete currently available factual profile, exact mapped segments and gaps, evidence observations, resource URLs/checksums and consumed-vs-unconsumed status, scope boundaries, and the provenance dependency graph. `--json` emits schema version 2 using canonical 0-based, half-open interval objects; batch JSON uses the separate `liftassess.ucsc_batch_result` report type. Schema v2 intentionally omits the legacy aggregate `verdict`, verdict-derived `decision_reason`, and preferred-candidate fields. `--details` is not yet available with batch input; `--details` and `--json` remain mutually exclusive. All completed report modes retain the explicit caveat that coordinate/evidence observations do not establish biological correctness.

For machine use, stdout remains the JSON document while status/progress stays on stderr, so normal shell redirection is safe:

```text
assess-liftover canFam3 canFam4 chrUn_JH373233:1845736-1845835 --json > assessment.json
```

## Validation strategy

Validation uses two complementary tracks.

### Mechanical evidence fixture

`canFam3` ↔ `canFam4` is now the real mechanical fixture for verifying extraction of chain/net/reciprocal-best evidence. The selected source interval is `chrUn_JH373233:1845735-1845835` in 0-based, half-open coordinates (CLI display locus `chrUn_JH373233:1845736-1845835`). The production cached-bundle path reproduces 170 candidate mappings across 114 target sequences, including a reverse, split, net-annotated, reciprocal-best-supported candidate and contrasting partial/reciprocal-best-absent alternatives. A historical pre-redesign run on 2026-08-17 reported the legacy `CONTESTED` result and 170 candidates; that legacy label is retained only as benchmark history. The durable fixture facts are the candidate/evidence/provenance observations themselves. These assemblies come from different dogs, so the fixture establishes **mechanical correctness of evidence extraction and deterministic result derivation**, not biological ground truth.

The exact multi-gigabyte UCSC resources remain outside the repository. Developers with the acquired cache can replay the frozen verification with `scripts/verify_canFam3_canFam4_mechanical_fixture.py`.

### Historical-resolution fixture

A `canFam3.1` → `canFam6` pedigree is planned for a future truth-bearing sanity check because those references derive from the same individual. A specific independently established locus still needs to be identified before this becomes a real fixture.

The historical fixture will be a sanity check, not a calibration set. v1 has no numeric score or fitted threshold to calibrate.

## Development setup

The project currently uses `uv` for environment and dependency management.

```bash
uv sync
uv run pytest
uv run ruff check src tests
uv run ruff format --check src tests
uv run mypy --strict src tests
git diff --check
```

The package currently targets Python 3.11 or newer.

## Scientific transparency

Non-obvious genomic, coordinate, provenance, and evidence decisions should be documented directly in code comments and docstrings. When implementation behavior is materially based on a scientific paper, standard, or primary-source implementation/documentation, the relevant code should cite that source near the logic it supports.

The goal is for researchers to be able to inspect not only **what** liftAssess concluded, but also **what evidence was examined, where it came from, which observations share upstream sources, and what assumptions the implementation made**.

[`docs/DESIGN.md`](docs/DESIGN.md) is the authoritative design document and contains the detailed scientific rationale, result semantics, scope, invariants, validation plan, and open questions.

## External resources

liftAssess does not bundle UCSC chain/net resources or depend on the UCSC liftOver executable for its core logic. UCSC and other external resources remain subject to their providers' own licensing and usage terms; liftAssess's GPL-3.0-only license does not relicense those external files.

UCSC resource terms are not uniform simply because multiple resources use chain format. In particular, UCSC's dedicated `liftOver/*.over.chain.gz` files are subject to UCSC's liftOver chain-file terms, including non-commercial-use restrictions unless an applicable commercial license has been obtained. Comparative `vsTarget/` chain/net resources follow the terms published for their own download directory. The established `canFam3/vsCanFam4/` mechanical-fixture directory states that its files are freely available for public use.

The resolver and acquisition layers remain separate. The acquisition API can retrieve one explicitly requested UCSC resource or execute an explicit plan for a complete discovered resource bundle into a caller-selected cache outside the source tree. Planning is no-network and enumerates every required URL plus its provider-terms classification; bundle execution additionally requires explicit acknowledgement of that transfer plan before any resource acquisition begins. This is deliberately separate from terms acknowledgement. Dedicated `liftOver/*.over.chain.gz` files are identified separately because UCSC applies additional liftOver-chain restrictions. Provider `md5sum.txt` entries are verified when an exact filename entry exists, while liftAssess SHA-256 remains the canonical artifact identity. Verified cached URL→artifact reuse is intentionally available offline and does not claim remote freshness; callers request an explicit refresh to contact UCSC and reacquire current bytes. A separate body-free metadata-inspection step can query provider HTTP headers after explicit terms acknowledgement and before transfer-plan acknowledgement. It does not create cache artifacts or transfer resource bodies. A live canFam3→canFam4 check on 2026-08-14 verified exact `Content-Length` values for all five comparative resources and `Accept-Ranges: bytes`; a separate small-range probe verified `206 Partial Content`, exact `Content-Range`, stable `ETag`, and `If-Range` behavior without downloading the full chain.

The acquisition layer now uses that verified transport contract opportunistically. If UCSC publishes an exact checksum and HEAD supplies an exact identity-encoded size, explicit byte-range support, and a strong ETag, an interrupted download retains a partial whose cache path is bound to the exact URL, total size, and validator. A later attempt resumes with `Range` + `If-Range`; a `200` response, malformed/mismatched `Content-Range`, or changed validator is never appended and instead triggers a fresh transfer. Shared resumable partials are never promoted directly into the content-addressed store: completion is copied into a unique private snapshot and independently MD5/SHA-256 verified before atomic publication, so a concurrent stale writer cannot mutate a published artifact. If the required resume metadata is unavailable, acquisition keeps the existing non-resumable streaming behavior.

A live interrupted acquisition check on 2026-08-14 exercised the full resume path against the 5,403,921-byte `canFam3.canFam4.rbest.chain.gz`: the first transfer was stopped after 262,144 bytes, retry revalidated with HEAD and requested `Range: bytes=262144-` plus `If-Range`, and the completed file matched both UCSC's published MD5 and the previously measured liftAssess SHA-256 identity.

A fully cached bundle can now feed the existing file-backed candidate engine directly while preserving content-addressed SHA-256 provenance and explicit shared upstream alignment provenance. `LIFTOVER-ONLY` consumes its chain. A complete `COMPARATIVE` bundle retains all five acquired resources, but the current v1 engine consumes only the all-chain, ordinary classified net, and reciprocal-best chain. UCSC's current automation produces the optional syntenic net by filtering the ordinary net for synteny, so liftAssess does not substitute it for the ordinary net evidence stream; the syntenic net and reciprocal-best net remain available on the cached bundle as retrieval/provenance context. The bridge validates the bundle's UCSC source/target database names against explicitly recorded assembly names or aliases rather than performing general alias resolution.

Automatic UCSC discovery is intended as a convenience, not a permanent hard dependency. User-supplied resources are part of the v1 design and remain subject to their own provider terms.

## Project documentation

- [`GETTING_STARTED.md`](docs/GETTING_STARTED.md) — beginner-oriented CLI guide, output interpretation, cache/network modes, and common mistakes.
- [`FEATURES.md`](docs/FEATURES.md) — complete catalog of implemented capabilities, expert APIs, and current non-features.
- [`DESIGN.md`](docs/DESIGN.md) — authoritative scientific and architectural specification.
- [`ROADMAP.md`](docs/ROADMAP.md) — implementation history, current review state, and planned milestones.
- [`PERFORMANCE.md`](docs/PERFORMANCE.md) — measured runtime characteristics, profiling results, and current optimization priorities.
- [`REFERENCES.md`](docs/REFERENCES.md) — literature, provider/format documentation, and technical evidence used by the project.

## Citation

If you use liftAssess in research, please cite the software using the version and release-date metadata in [`CITATION.cff`](CITATION.cff).

## License

liftAssess is licensed under the **GNU General Public License v3.0 (GPL-3.0-only)**.

## Project scope

v1 is deliberately narrow. It does not attempt to perform new sequence alignment by default, build large species-support databases, assign biological orthology automatically, use machine learning, or produce numeric confidence scores.


The objective is simpler: make ambiguous liftOver results **inspectable, evidence-based, provenance-aware, and explicit about uncertainty**.
