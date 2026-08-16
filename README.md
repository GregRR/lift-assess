# liftAssess

**liftAssess** evaluates ambiguous genomic coordinate liftOver mappings using transparent, provenance- and dependency-aware evidence, reporting whether mappings are **well supported**, **contested**, or **indeterminate**.

> **Status:** Active development. Core candidate-generation and evidence-extraction components are implemented and tested, but the end-to-end assessment pipeline, CLI, and final reporting layer are still under construction.

## Why liftAssess exists

Coordinate liftOver tools answer an important question: *where can this interval map in another assembly?* They do not, by themselves, explain how much support an ambiguous result deserves when there are multiple candidate targets, split mappings, duplicated sequence, alternative or unplaced scaffolds, or disagreement between methods.

liftAssess is intended to sit downstream of or alongside coordinate-conversion tools and answer a different question:

> **What does the available evidence say about the competing mappings, and how dependent are those lines of evidence on one another?**

liftAssess is an **assessor, not a resolver**. It does not claim that a locus is biologically “correct.” A `WELL_SUPPORTED` verdict means that the available informative evidence favors one candidate without material contradiction; it does **not** establish biological truth.

## Verdicts

v1 uses exactly three qualitative verdicts:

- `WELL_SUPPORTED` — available informative evidence favors one candidate, with no material evidence contradicting it.
- `CONTESTED` — multiple candidates retain meaningful support, or informative evidence materially disagrees.
- `INDETERMINATE` — available evidence is insufficient, non-discriminating, or too mutually dependent to distinguish candidates.

liftAssess deliberately does **not** produce a composite numeric confidence score in v1.

## Scientific principles

A few rules are central to the project:

- **Evidence provenance is part of the result.** Every observation records where it came from.
- **Dependent observations are not treated as independent confirmation.** For example, UCSC chain, net, and reciprocal-best evidence derived from the same upstream alignment remain linked through provenance.
- **Evidence availability is separate from evidentiary support.** A richly resourced assembly pair can still yield an `INDETERMINATE` result, while a sparse evidence tier can sometimes support a clear conclusion.
- **Interpretation stays close to the evidence.** liftAssess reports what the data support without automatically promoting an evidence pattern to a specific biological mechanism such as paralogy or pseudogene status.
- **Internal coordinates are 0-based, half-open.** Input and output boundaries must state or explicitly convert coordinate conventions rather than guessing them.
- **Same-species assembly comparisons are the v1 operational envelope.** Differences between assemblies from different individuals may represent real structural variation and are not automatically errors.

## Current implementation

The current development code includes:

- typed assembly, interval, candidate, evidence, verdict, and provenance models;
- a minimal streaming UCSC chain reader;
- forward- and reverse-strand chain geometry with 0-based half-open internal coordinates;
- chain-backed candidate projection, including split mappings;
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
- explicit bundle transfer planning and complete-or-error bundle acquisition for discovered
  `COMPARATIVE` and `LIFTOVER_ONLY` resource sets, with a separate transfer-plan acknowledgement before
  any planned resource acquisition begins;
- terms-gated, body-free remote metadata inspection using HTTP HEAD with identity encoding requested,
  preserving provider-advertised `Content-Length`, `Accept-Ranges`, `Last-Modified`, `ETag`, and
  `Content-Encoding` without transferring resource bodies;
- restart-safe resumable HTTPS acquisition when UCSC publishes an exact checksum and advertises an identity-encoded size, byte-range
  support, and a strong ETag; retained partials are bound to that exact URL/size/validator and resumed with
  `Range` + `If-Range`, while contract mismatches restart fresh rather than splicing representations;
- regression coverage for forward/reverse mappings, split mappings, gaps, repeated net chain IDs, provenance diamonds, reciprocal-best subsetting, and resource-discovery failure modes.
- a real `canFam3`→`canFam4` comparative mechanical fixture, replayed through the production cached-bundle path from exact externally cached UCSC resources without committing provider data to the repository.

## Not implemented yet

The project is not yet an end-to-end user tool. Major v1 work still includes:

- the assessor logic that deterministically converts evidence into `WELL_SUPPORTED`, `CONTESTED`, or `INDETERMINATE`;
- a future CLI/user-cache default (the library currently requires an explicit caller-supplied cache root);
- orchestration from the resulting candidate evidence into an assessment verdict and report;
- the command-line interface;
- human-readable summary and detailed/JSON reports;
- a future truth-bearing historical-resolution locus for the planned `canFam3.1`→`canFam6` sanity-check pedigree;
- optional flanking-gene orthology/synteny evidence;
- a defensible definition of candidate-rank and target-placement evidence where those concepts can be supported without inventing unsupported heuristics.

Until those pieces are complete, the repository should be treated as a developing scientific software project rather than a released analysis tool.

## Evidence-availability tiers

liftAssess distinguishes **how much evidence can be checked** from the eventual verdict.

- `COMPARATIVE` — a full comparative resource set is available, including chain/net context and reciprocal-best resources needed by the v1 UCSC evidence path.
- `LIFTOVER_ONLY` — only a liftOver chain resource is available, so candidate generation and chain-derived evidence can still proceed but comparative evidence is unavailable.

These are evidence-availability tiers, **not confidence tiers**.

## Architecture

```text
Candidate-generation engine
        |
        v
NormalizedCandidate[] + provenance
        |
        v
Assessor core
(evidence extraction, dependency/provenance labeling, verdict)
        |
        v
Assessment report
(summary + detailed dossier)
```

The assessor core is designed to consume normalized candidates plus provenance without knowing how those candidates were generated.

v1 intentionally has **one** concrete candidate-generation engine: an internal minimal UCSC chain/net reader. The interface boundary is kept clean so another engine can be added later if a real need appears, but liftAssess does not currently include a plugin registry, engine auto-discovery, or other speculative plugin infrastructure.

Chains and nets have different responsibilities:

- **chains generate candidate mappings**;
- **nets annotate and evaluate those candidates**.

Net availability is therefore not required for basic chain-backed candidate generation.

## Planned user workflow

The intended common-case interface is approximately:

```text
assess-liftover canFam3 canFam4 chr16:12345-12400
```

This is a **planned interface**, not a currently released command.

CLI-typed loci are intended to use the familiar UCSC-style 1-based, inclusive display convention and be converted immediately to liftAssess's canonical 0-based, half-open internal representation.

The final report is intended to provide two levels of output:

1. a concise summary with the verdict, candidates, and major evidence;
2. a detailed dossier/JSON representation with provenance, dependencies, chain/net details, checksums, and coordinate conventions.

Every report should make clear that evidentiary support is not proof of biological correctness.

## Validation strategy

Two complementary validation tracks are planned.

### Mechanical evidence fixture

`canFam3` ↔ `canFam4` is now the real mechanical fixture for verifying extraction of chain/net/reciprocal-best evidence. The selected source interval is `chrUn_JH373233:1845735-1845835` in 0-based, half-open coordinates. The production cached-bundle path reproduces 170 candidate mappings across 114 target sequences, including a reverse, split, net-annotated, reciprocal-best-supported candidate and contrasting partial/reciprocal-best-absent alternatives. These assemblies come from different dogs, so the fixture establishes **mechanical correctness of evidence extraction**, not biological ground truth.

The exact multi-gigabyte UCSC resources remain outside the repository. Developers with the acquired cache can replay the frozen verification with `scripts/verify_canFam3_canFam4_mechanical_fixture.py`.

### Historical-resolution fixture

A `canFam3.1` → `canFam6` pedigree is planned for a future truth-bearing sanity check because those references derive from the same individual. A specific independently established locus still needs to be identified before this becomes a real fixture.

The historical fixture will be a sanity check, not a calibration set. v1 has no numeric score or fitted threshold to calibrate.

## Development setup

The project currently uses `uv` for environment and dependency management.

```bash
uv sync --extra test
uv run pytest
uvx ruff check src tests
uv run --extra test --with mypy mypy --strict src tests
git diff --check
```

The package currently targets Python 3.11 or newer.

## Scientific transparency

Non-obvious genomic, coordinate, provenance, and evidence decisions should be documented directly in code comments and docstrings. When implementation behavior is materially based on a scientific paper, standard, or primary-source implementation/documentation, the relevant code should cite that source near the logic it supports.

The goal is for researchers to be able to inspect not only **what** liftAssess concluded, but also **what evidence was examined, where it came from, which observations share upstream sources, and what assumptions the implementation made**.

`DESIGN.md` is the authoritative v1 design document and contains the detailed scientific rationale, scope, invariants, validation plan, and open questions.

## External resources

liftAssess does not bundle UCSC chain/net resources or depend on the UCSC liftOver executable for its core logic. UCSC and other external resources remain subject to their providers' own licensing and usage terms; liftAssess's GPL-3.0 license does not relicense those external files.

UCSC resource terms are not uniform simply because multiple resources use chain format. In particular, UCSC's dedicated `liftOver/*.over.chain.gz` files are subject to UCSC's liftOver chain-file terms, including non-commercial-use restrictions unless an applicable commercial license has been obtained. Comparative `vsTarget/` chain/net resources follow the terms published for their own download directory. The planned `canFam3/vsCanFam4/` mechanical-fixture directory currently states that its files are freely available for public use.

The resolver and acquisition layers remain separate. The acquisition API can retrieve one explicitly requested UCSC resource or execute an explicit plan for a complete discovered resource bundle into a caller-selected cache outside the source tree. Planning is no-network and enumerates every required URL plus its provider-terms classification; bundle execution additionally requires explicit acknowledgement of that transfer plan before any resource acquisition begins. This is deliberately separate from terms acknowledgement. Dedicated `liftOver/*.over.chain.gz` files are identified separately because UCSC applies additional liftOver-chain restrictions. Provider `md5sum.txt` entries are verified when an exact filename entry exists, while liftAssess SHA-256 remains the canonical artifact identity. Verified cached URL→artifact reuse is intentionally available offline and does not claim remote freshness; callers request an explicit refresh to contact UCSC and reacquire current bytes. A separate body-free metadata-inspection step can query provider HTTP headers after explicit terms acknowledgement and before transfer-plan acknowledgement. It does not create cache artifacts or transfer resource bodies. A live canFam3→canFam4 check on 2026-08-14 verified exact `Content-Length` values for all five comparative resources and `Accept-Ranges: bytes`; a separate small-range probe verified `206 Partial Content`, exact `Content-Range`, stable `ETag`, and `If-Range` behavior without downloading the full chain.

The acquisition layer now uses that verified transport contract opportunistically. If UCSC publishes an exact checksum and HEAD supplies an exact identity-encoded size, explicit byte-range support, and a strong ETag, an interrupted download retains a partial whose cache path is bound to the exact URL, total size, and validator. A later attempt resumes with `Range` + `If-Range`; a `200` response, malformed/mismatched `Content-Range`, or changed validator is never appended and instead triggers a fresh transfer. Shared resumable partials are never promoted directly into the content-addressed store: completion is copied into a unique private snapshot and independently MD5/SHA-256 verified before atomic publication, so a concurrent stale writer cannot mutate a published artifact. If the required resume metadata is unavailable, acquisition keeps the existing non-resumable streaming behavior.

A live interrupted acquisition check on 2026-08-14 exercised the full resume path against the 5,403,921-byte `canFam3.canFam4.rbest.chain.gz`: the first transfer was stopped after 262,144 bytes, retry revalidated with HEAD and requested `Range: bytes=262144-` plus `If-Range`, and the completed file matched both UCSC's published MD5 and the previously measured liftAssess SHA-256 identity.

A fully cached bundle can now feed the existing file-backed candidate engine directly while preserving content-addressed SHA-256 provenance and explicit shared upstream alignment provenance. `LIFTOVER_ONLY` consumes its chain. A complete `COMPARATIVE` bundle retains all five acquired resources, but the current v1 engine consumes only the all-chain, ordinary classified net, and reciprocal-best chain. UCSC's current automation produces the optional syntenic net by filtering the ordinary net for synteny, so liftAssess does not substitute it for the ordinary net evidence stream; the syntenic net and reciprocal-best net remain available on the cached bundle as retrieval/provenance context. The bridge validates the bundle's UCSC source/target database names against explicitly recorded assembly names or aliases rather than performing general alias resolution.

Automatic UCSC discovery is intended as a convenience, not a permanent hard dependency. User-supplied resources are part of the v1 design and remain subject to their own provider terms.

## Project documentation

- [`DESIGN.md`](DESIGN.md) — authoritative scientific and architectural specification.
- [`ROADMAP.md`](ROADMAP.md) — implementation history, current review state, and planned milestones.

## Citation

If you use liftAssess in research, please cite the software using the metadata in [`CITATION.cff`](CITATION.cff). Version-specific citation metadata will be added with the first formal release.

## License

liftAssess is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.

## Project scope

v1 is deliberately narrow. It does not attempt to perform new sequence alignment by default, build large species-support databases, assign biological orthology automatically, use machine learning, or produce numeric confidence scores.

The objective is simpler: make ambiguous liftOver results **inspectable, evidence-based, provenance-aware, and explicit about uncertainty**.