# liftAssess Roadmap

This roadmap tracks implementation status and sequencing for **liftAssess**. It is a development-status document, not the scientific specification.

`DESIGN.md` remains authoritative for the project's problem definition, scientific invariants, coordinate semantics, verdict meanings, v1 scope, architecture, licensing constraints, and validation requirements. This file answers a different set of questions: **what has been built, what is being reviewed now, what comes next, and what is deliberately deferred?**

## Current status — 2026-08-12

liftAssess is an active, private development repository. It is not yet a released end-to-end analysis tool.

The project has moved past isolated parsers and models into an integrated UCSC evidence pipeline. The current code can:

- parse UCSC chain and net resources as streams;
- project source intervals into normalized candidate mappings, including reverse-strand and split mappings;
- extract mapping coverage, chain-gap, chain-score, net, and reciprocal-best evidence;
- preserve provenance relationships so observations derived from the same upstream alignment are not presented as independent confirmation;
- discover UCSC comparative and liftOver-only resource availability without silently treating transport failures as resource absence;
- stream plain or gzip-compressed local chain/net resources into the engine without materializing whole files;
- run a real UCSC liftOver chain through the local-file → parser → engine path;
- and, in the current provenance batch, identify exact local resource bytes with content-addressed SHA-256 provenance, verify those bytes again as they are streamed into parsing, and optionally verify provider-published checksums.

The current reciprocal-best discovery tree has **132 tests** after adding an orchestration-level regression for sibling-directory transport failures. The preceding 131-test tree passed pytest, Ruff, and strict mypy locally and a direct live UCSC resolver check; the final added test is synthetic and exercises error propagation only. Resource-provenance review separately reproduced the original hash-before-parse race, verified the corrected failure behavior, stress-tested gzip/raw-byte hashing, checked handle cleanup and streaming memory behavior, and found no remaining blocking provenance issues.

What liftAssess still cannot do is the most important user-facing part: it does **not yet assign `WELL_SUPPORTED`, `CONTESTED`, or `INDETERMINATE` verdicts**, and it does not yet provide the planned CLI/report workflow.

## Implementation history

### 1. Foundation and scientific data model — complete

Implemented:

- assembly and genomic-interval models;
- canonical 0-based, half-open internal coordinates;
- normalized candidate representation;
- mapping segments and orientation;
- evidence observations and evidence kinds;
- provenance sources, typed identifiers, and provenance graph relationships;
- evidence-availability tiers (`COMPARATIVE`, `LIFTOVER_ONLY`);
- exact v1 verdict vocabulary (`WELL_SUPPORTED`, `CONTESTED`, `INDETERMINATE`).

Important boundary established: an `Assessment` model may hold a verdict, but no numeric confidence score exists and no code should claim biological correctness.

### 2. UCSC chain parsing and candidate projection — complete

Implemented and reviewed:

- streaming chain parsing;
- forward- and reverse-strand coordinate handling;
- source-locus projection into candidate target geometry;
- split mappings represented as multiple aligned segments under one candidate;
- mapping coverage and uncovered-source intervals;
- source/destination chain-gap characterization;
- chain score as evidence.

Key design decision: the old/source assembly corresponds to UCSC chain target coordinates and the new/destination assembly to query coordinates. Reverse query coordinates are normalized before entering the core model.

### 3. UCSC net evidence — complete

Implemented and reviewed:

- streaming net parsing;
- hierarchy depth preservation;
- repeated chain IDs preserved rather than treated as unique net keys;
- aligned-bases (`ali`) evidence;
- duplicated-query-bases (`qDup`) evidence;
- net classification and hierarchy evidence;
- candidate matching against actual aligned segments rather than only a candidate's bounding target span.

Key design decision: chains generate candidates; nets annotate/evaluate them. Net availability is never required for chain-backed candidate generation.

### 4. Reciprocal-best membership — complete

Implemented and reviewed:

- exact geometry-based reciprocal-best matching rather than chain-ID equality;
- `FULL`, `PARTIAL`, and `NONE` membership states;
- aligned candidate bases as the denominator;
- explicit completeness requirements before absence/partial membership can be interpreted;
- internal-gap handling;
- provenance linking reciprocal-best evidence to its upstream alignment when appropriate.

The engine's candidate-relevance filter is deliberately lossless with respect to the downstream reciprocal-best matcher, so the caller's declared completeness scope is preserved rather than rewritten by orchestration.

### 5. UCSC resource discovery — complete

Implemented and reviewed:

- verification of actual UCSC directory contents instead of assuming constructed URLs exist;
- detection of full comparative resource sets versus liftOver-only resources;
- no partial-comparative state represented as `COMPARATIVE`;
- 404/resource absence distinguished from timeout, DNS, and server/transport failures;
- relative and absolute directory links normalized before resource comparison.

The resolver discovers and classifies resources only; it does **not** download them. A 2026-08-12 live-directory review exposed an asymmetric UCSC publication layout in which reciprocal-best files for a source→target pair were hosted under the sibling/reverse `target/vsSource/reciprocalBest/` directory while `source/vsTarget/` held the ordinary comparative resources. Milestone 10 implemented and live-verified that fallback without changing coordinate semantics.

### 6. Single-pass UCSC engine orchestration — complete

Implemented and reviewed:

- one public orchestration path from chain records to normalized candidates plus available evidence;
- one-shot chain/net/reciprocal-best iterators consumed safely;
- no whole-resource rescanning per candidate;
- repeated net fills preserved;
- reverse-strand and split mappings regression-tested through the orchestration boundary;
- no ranking, preferred-candidate selection, or verdict assignment in the engine.

This is the current boundary between candidate/evidence generation and the still-unimplemented assessor logic.

### 7. Local resource-file integration — complete

Implemented and reviewed:

- streaming local chain and net files into the existing parsers and engine;
- gzip detection by file magic bytes rather than filename suffix;
- support for both plain-text and gzip-compressed resources;
- file adapters that do not infer licensing, download resources, manufacture provenance, or assume reciprocal-best completeness.

This created the first real file → parser → candidate/evidence path.

### 8. Real LIFTOVER-ONLY smoke run — complete mechanical checkpoint

A real UCSC `canFam3ToCanFam4.over.chain.gz` resource was exercised on 2026-08-12 after its provider-published MD5 was verified.

Measured result:

- exact downloaded-byte SHA-256: `c79c9e7c2a3d546f7a9d7efe27cc8815da611d79adb0da4e4ff1556810f28f48`;
- test source interval, 0-based half-open: `chr1:12514-12534`;
- one chain-derived candidate: `chr1:660-680`;
- same orientation, one aligned segment;
- full 20/20 source-base coverage;
- no chain gaps;
- no assessment verdict computed.

This establishes real-file mechanical plumbing under `LIFTOVER_ONLY`. It is **not** the planned full comparative fixture and does not establish biological support or correctness.

### 9. Resource identity and integrity — complete

Implemented and reviewed:

- computes SHA-256 over the exact local resource-file bytes;
- uses that SHA-256 as the canonical file provenance identity;
- creates content-addressed file `ProvenanceSource` identities so identical bytes at different paths/names converge on the same artifact identity;
- keeps provider MD5/SHA-256 verification separate from scientific provenance identity;
- validates canonical `sha256:<64 lowercase hexadecimal characters>` provenance identifiers;
- requires callers creating file provenance to make the `derived_from` relationship explicit, including explicitly choosing `()` when no upstream source is claimed;
- verifies the SHA-256 again over the exact raw bytes streamed into the parser and raises on mismatch, preventing a file changed after provenance construction from silently producing evidence with stale provenance;
- verifies file-backed chain and net provenance can retain a shared upstream alignment relationship through the comparative evidence path.

Review conclusions retained as design guidance:

- exact on-disk bytes, including gzip/compression bytes, are the v1 artifact identity;
- provider MD5 is integrity metadata only, never provenance identity or evidentiary strength;
- canonical lowercase SHA-256 validation is intentional;
- content-addressed file IDs are appropriate internally, but human-facing reports should use labels/metadata rather than raw hash-heavy candidate identifiers.

The confirmed hash→parse mismatch was fixed by verifying the raw byte stream actually consumed by parsing against the provenance SHA-256 before a file-backed engine call can return. This preserves the current engine API; producing provenance only after a one-pass parse would require a broader orchestration/model redesign because evidence needs finalized provenance while it is being attached.

## Next milestones

### 10. Correct reciprocal-best discovery across asymmetric UCSC pair directories — complete

Measured against the live UCSC directories on 2026-08-12: `canFam3/vsCanFam4/` contains the ordinary comparative chain/net resources but no listed `reciprocalBest/` directory, while `canFam4/vsCanFam3/reciprocalBest/` contains both directional reciprocal-best chain/net files, including `canFam3.canFam4.rbest.chain.gz` and `canFam3.canFam4.rbest.net.gz`. A second checked pair (`rheMac10`/`hg38`) showed the same publication pattern. This is evidence for a real fallback location, not a claim that UCSC guarantees the layout universally.

Implemented in the resolver:

- checks the forward comparative directory's `reciprocalBest/` location first;
- if the exact directional files are absent, checks the sibling/reverse comparative directory's `reciprocalBest/`;
- accepts the fallback only when the exact `source.target.rbest.{chain,net}.gz` files are actually observed;
- keeps unit regressions synthetic, including the measured asymmetric layout and an opposite-direction-only negative case;
- leaves projection and reciprocal-best evidence geometry unchanged because the directional filename and file headers, not the hosting directory, define coordinate semantics.

Live verification completed on 2026-08-13: `discover_ucsc_resources("canFam3", "canFam4")` returned `COMPARATIVE`, with chain/net/syn-net URLs under `canFam3/vsCanFam4/` and the exact directional reciprocal-best URLs under `canFam4/vsCanFam3/reciprocalBest/`. The unit suite also includes an orchestration-level regression confirming that a transport failure during the sibling lookup propagates as `UCSCResourceDiscoveryError` rather than being misread as absence and silently downgrading the evidence tier.

### 11. Resource acquisition and cache

Goal: connect already-implemented UCSC resource discovery to exact local files without mixing network behavior, licensing acceptance, and scientific interpretation.

Planned responsibilities:

- retrieve discovered resources only after any required provider/licensing acknowledgement;
- retain source URL and retrieval metadata;
- verify provider checksums when available;
- compute liftAssess's canonical SHA-256 file identity;
- store downloads outside the source tree/release artifacts;
- avoid redundant transfers when the correct content-addressed artifact is already cached;
- support user-supplied local resources as an equal path, not a second-class fallback;
- use a practical bulk-transfer strategy for very large UCSC comparative resources.

**Repository note:** when the project introduces its first runtime/cache directory, `.gitignore` must be updated in the same batch.

### 12. Full `canFam3` ↔ `canFam4` comparative mechanical fixture

Goal: prove real extraction of the comparative evidence families already implemented.

The fixture should exercise, from real UCSC comparative resources after milestone 10's asymmetric reciprocal-best discovery correction:

- multiple candidate mappings where available;
- chain score;
- coverage and gaps;
- net `ali` and `qDup`;
- net classification/hierarchy;
- reciprocal-best membership;
- dependency-aware provenance showing which observations share one upstream alignment.

This fixture is mechanical validation only because canFam3 and canFam4 represent different dogs. It must not be described as biological ground truth.

A practical acquisition strategy is required before treating the multi-gigabyte full comparative chain as routine test data. Large provider resources should not be committed to the repository.

### 13. Assessor core and deterministic verdict logic

This is the largest remaining scientific implementation milestone.

Goal: transform normalized candidates plus provenance-aware evidence into exactly one of:

- `WELL_SUPPORTED`;
- `CONTESTED`;
- `INDETERMINATE`.

Constraints:

- no numeric composite confidence score;
- no automatic claim that a candidate is biologically correct;
- dependent observations cannot masquerade as independent confirmation;
- evidence availability and evidentiary support remain separate;
- one candidate can still be `CONTESTED` if materially contradictory evidence exists;
- qualitative rules must be deterministic enough that the same evidence set produces the same verdict.

This logic should be built from explicit evidence patterns and adversarial tests, not by ranking raw chain scores or inventing thresholds because they are convenient.

### 14. Assessment/report orchestration

Once verdict logic exists:

- build the full source + target + locus → resources → candidates → evidence → verdict path;
- create the final `Assessment` object from real engine output;
- ensure evidence tier is always displayed independently from verdict;
- expose provenance/dependency detail sufficient for scientific audit.

### 15. CLI and user-facing reports

Implement the planned common-case workflow:

```text
assess-liftover canFam3 canFam4 chr16:12345-12400
```

Required behavior:

- CLI locus strings use UCSC-style 1-based, inclusive display coordinates;
- conversion to 0-based, half-open occurs explicitly at the boundary;
- plain-language evidence-availability tier is shown before interpretation;
- default output is concise;
- `--details` exposes the evidence dossier;
- JSON output preserves coordinates, provenance, dependencies, resource hashes, and verdict semantics.

The README should gain real installation/usage examples only when these commands actually work.

### 16. First public alpha milestone

The repository is currently private. A good public-transition milestone is reached when:

- a researcher can run one documented end-to-end real UCSC assessment from the CLI;
- the comparative mechanical fixture works;
- provenance/resource identity is visible in output;
- evidence availability and support are clearly separated;
- README usage commands are real and reproducible;
- the project still clearly labels itself pre-release/alpha and states that well-supported does not mean biologically correct.

A formal versioned release can follow after the CLI/report workflow is coherent enough for an outside researcher to try without project-specific guidance.

## Remaining v1 work after the first public alpha

### Historical-resolution sanity fixture

Identify one concrete, independently established `canFam3.1` → `canFam6` locus from the same-individual Tasha pedigree. This should test sensible behavior under sparse evidence and later assembly resolution without being used as a calibration set.

### Assembly identity/canonicalization boundary

`AssemblyIdentifier` currently uses frozen-dataclass structural equality, including optional accession and aliases. That conservative behavior is safe today because `resources.py` still uses plain UCSC `source_db`/`target_db` strings and each engine call threads one `AssemblyIdentifier` instance through its candidate geometry.

The issue becomes live when the resolver/cache/CLI layer first converts those database strings into assembly objects that may then be compared with independently constructed objects carrying additional accession/alias metadata. Before that bridge lands, define explicit identity/canonicalization semantics rather than either relying indefinitely on structural equality or weakening `__eq__` speculatively.

### Candidate-rank and target-placement evidence

These remain intentionally deferred until a defensible locus-scoped definition exists. Raw chain-score ordering is not an acceptable substitute for a scientifically justified candidate-rank concept.

### Optional flanking-gene synteny context

Select and verify a real source, likely an external orthology/synteny provider, and model its provenance/dependencies explicitly. Do not assume orthology calls are independent merely because they come from a different API or file.

### Detailed/JSON schema stabilization

Finalize stable machine-readable report semantics once the assessor core and at least one real end-to-end assessment exist. Avoid freezing a schema around synthetic-only assumptions.

## Deliberately deferred beyond v1

The following should not distract from making the first version scientifically useful:

- plugin registry, entry-point discovery, or general engine configuration framework;
- automatic support for many candidate-generation engines before a real second engine exists;
- fresh minimap2/lastz alignment by default;
- machine-learning confidence models;
- composite numeric confidence scores;
- automatic biological orthology claims;
- large bundled species databases;
- hosted service infrastructure.

## Development and review discipline

Each substantial batch should continue to follow the current pattern:

1. implement one narrow capability;
2. cover its failure modes with tests;
3. pass pytest, Ruff, strict mypy, and `git diff --check`;
4. perform a targeted review before commit when the batch changes scientific, coordinate, provenance, completeness, or interpretation semantics;
5. verify reviewer findings against the actual code/design before accepting proposed fixes;
6. commit only after confirmed findings are resolved and gates pass again.

The purpose of this discipline is not to maximize test count. It is to make each scientific claim in the software reconstructable from explicit code, provenance, tests, and primary-source checks.
