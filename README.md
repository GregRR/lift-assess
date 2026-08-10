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
- UCSC resource discovery and evidence-availability tier detection;
- regression coverage for forward/reverse mappings, split mappings, gaps, repeated net chain IDs, provenance diamonds, reciprocal-best subsetting, and resource-discovery failure modes.

## Not implemented yet

The project is not yet an end-to-end user tool. Major v1 work still includes:

- the assessor logic that deterministically converts evidence into `WELL_SUPPORTED`, `CONTESTED`, or `INDETERMINATE`;
- orchestration from assembly pair + locus through candidate generation, evidence extraction, and assessment;
- a practical strategy for obtaining and caching external resources without silently downloading very large comparative datasets;
- the command-line interface;
- human-readable summary and detailed/JSON reports;
- end-to-end validation against real assembly/locus fixtures;
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

`canFam3` ↔ `canFam4` is the planned mechanical fixture for verifying extraction of chain/net/reciprocal-best evidence. These assemblies come from different dogs, so this fixture is for **mechanical correctness of evidence extraction**, not biological ground truth.

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

liftAssess does not bundle UCSC chain/net resources or depend on the UCSC liftOver executable for its core logic. UCSC and other external resources remain subject to their providers' own licensing and usage terms.

Automatic UCSC discovery is intended as a convenience, not a permanent hard dependency. User-supplied resources are part of the v1 design.

## License

liftAssess is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.

## Project scope

v1 is deliberately narrow. It does not attempt to perform new sequence alignment by default, build large species-support databases, assign biological orthology automatically, use machine learning, or produce numeric confidence scores.

The objective is simpler: make ambiguous liftOver results **inspectable, evidence-based, provenance-aware, and explicit about uncertainty**.