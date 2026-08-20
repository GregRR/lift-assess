# liftAssess Roadmap

This roadmap tracks implementation status and sequencing for **liftAssess**. It is a development-status document, not the scientific specification.

[`DESIGN.md`](DESIGN.md) remains authoritative for the project's problem definition, scientific invariants, coordinate semantics, result-model semantics, current scope, architecture, licensing constraints, and validation requirements. This file answers a different set of questions: **what has been built, what is being reviewed now, what comes next, and what is deliberately deferred?**

## Current status — 2026-08-19

liftAssess `v0.1.0a1` was released on 2026-08-17 as the project's first public alpha. The project remains active scientific software under development and should not be treated as a mature or stable analysis platform.

The post-alpha 50-case real-world validation / UX program and subsequent design review are now complete enough to set the next implementation direction. The approved redesign removes the legacy aggregate verdict taxonomy from the target model, adopts a facts-first orthogonal result profile with progressive disclosure, deliberately permits a new pre-release schema version, and starts scalable/indexed resource-access work immediately after the first renderer/profile slice. The released alpha code remains the historical baseline until those milestones are implemented.

The project now has an integrated UCSC evidence pipeline plus a completed real comparative mechanical fixture. The current code can:

- parse UCSC chain and net resources as streams;
- project source intervals into normalized candidate mappings, including reverse-strand and split mappings;
- extract mapping coverage, chain-gap, chain-score, net, and reciprocal-best evidence;
- preserve provenance relationships so observations derived from the same caller-declared upstream alignment are not presented as independent confirmation;
- discover, inspect, acquire, cache, resume, and integrity-check UCSC comparative and liftOver-only resources;
- bridge a complete cached resource bundle directly into the file-backed candidate/evidence engine;
- reproduce a measured `canFam3`→`canFam4` comparative fixture through that production cached-bundle path while keeping all multi-gigabyte provider resources outside the repository;
- run a reviewed deterministic assessor over normalized coverage and reciprocal-best evidence without a numeric score;
- compose a complete cached UCSC bundle through candidate/evidence generation and the assessor into an auditable assessment report;
- and run the `assess-liftover` CLI cache-first or fully offline with measured cache-verification, transfer, and assessment-read progress.

The routine automated suite contains **316 tests**. The real comparative fixture is intentionally an external-cache integration verification rather than a routine pytest case because its five UCSC resources total 2,686,242,854 compressed bytes.

The deterministic assessor and assessment/report orchestration milestones are complete and reviewed. The common-case CLI, cache-first/offline execution, human-readable and JSON detailed reporting, measured cache-verification progress, measured transfer progress, and measured assessment-read progress are implemented. Milestone 15 and the pre-alpha semantic/output-hardening slice are complete and reviewed. Post-hardening real comparative CLI and independent-verifier runs completed successfully on 2026-08-17. A real-provider transfer-progress smoke check remains desirable but non-gating. Milestone 16 is complete. The first public alpha, `v0.1.0a1`, was tagged and published to PyPI on 2026-08-17.

## Implementation history

### 1. Foundation and scientific data model — complete

Implemented:

- assembly and genomic-interval models;
- canonical 0-based, half-open internal coordinates;
- normalized candidate representation;
- mapping segments and orientation;
- evidence observations and evidence kinds;
- provenance sources, typed identifiers, and provenance graph relationships;
- evidence-availability tiers (`COMPARATIVE`, `LIFTOVER-ONLY`);
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

This establishes real-file mechanical plumbing under `LIFTOVER-ONLY`. It is **not** the planned full comparative fixture and does not establish biological support or correctness.

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

### 11. Resource acquisition and cache — complete

Goal: connect already-implemented UCSC resource discovery to exact local files without mixing network behavior, licensing acknowledgement, and scientific interpretation.

Implemented in the first acquisition slice:

- acquires one explicitly requested resolver-produced UCSC URL at a time;
- requires a caller-selected cache root outside the source tree rather than inventing a repository-local cache location;
- surfaces UCSC's general licensing page plus the relevant comparison/liftOver directory terms reference and requires explicit acknowledgement **before any network access**;
- distinguishes dedicated `liftOver/*.over.chain.gz` resources from comparative `vsTarget/` resources so the restricted liftOver-chain case is not generalized to every chain-format file;
- reads the resource directory's `md5sum.txt` and verifies MD5 when an exact filename entry is published; an existing `md5sum.txt` without that filename is treated as checksum unavailable, not as corruption;
- streams downloads into a temporary file while computing liftAssess's canonical SHA-256 and any provider MD5, validates HTTP `Content-Length` when supplied, then publishes the artifact atomically under a content-addressed `artifacts/sha256/...` path;
- writes an atomic URL index retaining source URL, retrieval timestamp, SHA-256, provider checksum metadata when available, and applicable terms references;
- reuses a cached artifact only after re-verifying its local SHA-256; retains the provider checksum recorded at acquisition time without requiring a network freshness check;
- supports `refresh=True` to contact the provider and reacquire the current representation, which is especially relevant when the provider publishes no checksum for that exact file;
- cleans ordinary temporary downloads on failure; identical bytes acquired from different URLs converge on one artifact.

Implemented in the second acquisition slice:

- converts a verified `UCSCResourceBundle` into an inspectable, no-network transfer plan containing the exact resource roles/URLs required by its evidence tier and the terms classification for each URL;
- requires a separate `transfer_plan_acknowledged=True` before bundle execution can call the single-resource acquisition layer, so discovery alone cannot silently trigger a five-resource comparative transfer;
- acquires the plan sequentially through the already-reviewed single-resource cache path and returns a `CachedUCSCResourceBundle` only after every required item succeeds or is verified in cache;
- preserves the no-partial-comparative invariant at the returned-object boundary while deliberately retaining any immutable content-addressed artifacts successfully published before a later item fails, so a retry can reuse them;
- keeps `LIFTOVER-ONLY` as a one-chain local bundle and `COMPARATIVE` as the exact five-role set: chain, net, syntenic net, reciprocal-best chain, reciprocal-best net;
- binds each plan item's surfaced terms metadata back to its exact URL classification so an inspectable plan cannot display one terms class while executing another URL.

Implemented in the third acquisition slice:

- after explicit provider-terms acknowledgement, issues body-free HTTP HEAD requests for the exact
  URLs already present in an acquisition plan;
- requests identity encoding and preserves provider-advertised `Content-Length`, `Accept-Ranges`,
  `Last-Modified`, `ETag`, and `Content-Encoding` values rather than guessing from directory display
  text; a contrary non-identity content encoding is retained as metadata but excluded from transfer-size
  totals;
- reports a complete bundle byte total only when every planned resource advertises `Content-Length`, while still exposing the sum of individually known lengths;
- preserves the existing role/URL/directional-pair validation in the inspection result;
- does not treat `Accept-Ranges` or any other header as proof of resumable HTTP support and does not begin resource-body acquisition.

Implemented in the fourth acquisition slice:

- treats HEAD metadata as an opportunistic resume capability rather than a new acquisition prerequisite: if HEAD inspection fails or lacks the required fields, the existing fresh streaming path remains available;
- enables resumable HTTPS only when UCSC also publishes an exact checksum and the representation has an exact identity-encoded `Content-Length`, explicit `Accept-Ranges: bytes`, and a strong ETag; weak ETags and `Last-Modified` alone are not used as resume validators;
- retains interrupted partial bytes under a cache path derived from the source URL plus the exact total length and strong-ETag hash, so a changed provider representation cannot be appended to an old prefix by filename accident;
- resumes from the existing prefix with `Range: bytes=<offset>-` and `If-Range: <ETag>`, requiring `206 Partial Content` and a matching `Content-Range` before any resumed bytes are written;
- rejects a changed/missing resume contract by restarting through the fresh-transfer path rather than splicing incompatible representations;
- can publish a fully received retained partial on retry without another resource-body GET, while still re-running provider checksum lookup/HEAD and final SHA-256/provider-MD5 verification;
- never promotes the shared deterministic partial inode directly into the content-addressed store: a completing process snapshots its open partial into a unique private temporary file, recomputes provider MD5 and SHA-256 over that snapshot, and publishes only the private file, preventing a concurrent stale writer from mutating an already-published artifact;
- requests identity encoding for ordinary resource GETs as well, so cached bytes remain the provider resource representation even when no provider checksum is published.

Measured provider detail checked 2026-08-13: UCSC's `canFam3/vsCanFam4/md5sum.txt` publishes an MD5 for `canFam3.canFam4.all.chain.gz` and `canFam3.canFam4.syn.net.gz` but not `canFam3.canFam4.net.gz`; `canFam4/vsCanFam3/reciprocalBest/md5sum.txt` publishes MD5 values for both directional reciprocal-best chain/net files. Therefore the exact-filename checksum-optional behavior is required by the real fixture resources rather than being hypothetical.

Measured size context checked from the live UCSC directory listings on 2026-08-14: the planned canFam3→canFam4 comparative set includes a roughly 2.5 GB `all.chain.gz`, alongside a roughly 10 MB net, 9.1 MB syntenic net, 5.2 MB directional reciprocal-best chain, and 7.8 MB directional reciprocal-best net. This is why bundle planning and explicit acknowledgement were added before any user-facing automatic comparative transfer. The subsequent HEAD verification below supplied exact machine-readable lengths for all five resources.

Live provider verification completed 2026-08-14 for the five-resource canFam3→canFam4 comparative plan. HEAD returned exact `Content-Length` values for every resource, totaling 2,686,242,854 bytes; the 2,652,632,416-byte forward all-chain accounts for nearly the entire transfer. Every resource advertised `Accept-Ranges: bytes`, and the large chain supplied both `ETag` and `Last-Modified`. A separate small-range probe against that chain returned `206 Partial Content`, exact `Content-Range` and `Content-Length`, a stable ETag across adjacent requests, successful `If-Range`, and byte-identical reconstruction of adjacent ranges versus one combined range. These checks transferred only a few KiB of the chain and did not execute the full comparative acquisition.

End-to-end resume behavior was then measured on 2026-08-14 against the 5,403,921-byte directional reciprocal-best chain. The first GET was intentionally interrupted after 262,144 bytes; retry performed a fresh HEAD and requested `Range: bytes=262144-` with the same strong ETag in `If-Range`. The completed resource matched UCSC MD5 `03bc68aa1c8ce4582cff71a8813f0b6b` and liftAssess SHA-256 `34f4061fd29e7720c7eb2adc1ea8299e86f21f08e18f97f9ba468cf8b466690c`, and the retained partial was removed after publication.

Focused review then reproduced a concurrent-writer corruption bug in the original resumable finalization path: directly renaming the shared partial into the artifact store allowed a second process holding the same inode open to mutate the published artifact after `os.replace`. The corrected path snapshots the shared partial into a private temporary file through the completing process's open descriptor, recomputes provider MD5 and SHA-256 over that private snapshot, and publishes only the snapshot. A deterministic regression reproduces the corruption against the pre-fix implementation and verifies that the published artifact remains immutable after the fix.

Implemented in the fifth acquisition slice:

- bridges a complete `CachedUCSCResourceBundle` directly into the existing file-backed candidate engine while preserving the lower-level user-supplied-file API unchanged;
- builds content-addressed file provenance from the SHA-256 identities already recorded by acquisition, avoiding an extra pre-parse full-file hash while retaining the parser's exact-byte SHA-256 verification on every consumed resource;
- requires caller-supplied upstream alignment provenance, so acquisition metadata is never mistaken for evidence independence or alignment-process provenance;
- validates cached bundle `source_db`/`target_db` strings against only the source/target assembly's explicit name or aliases, avoiding speculative general alias resolution;
- maps `LIFTOVER-ONLY` to its chain-only engine path and `COMPARATIVE` to the all-chain, ordinary classified net, and reciprocal-best chain with `COMPLETE_RESOURCE` semantics;
- deliberately retains the syntenic net and reciprocal-best net on the five-resource comparative bundle without parsing them as current v1 engine inputs. UCSC's current automation confirms that `*.syn.net.gz` is a `netFilter -syn` derivative of the ordinary net, so substituting it would discard non-syntenic placements rather than add evidence.

This closes the acquisition/cache milestone itself. The future CLI's default cache location, progress/refresh controls, and large-transfer confirmation belong to milestone 15. Presentation of retrieval metadata in the final assessment belongs to milestone 14.

**Repository note:** this implementation deliberately does **not** create a runtime/cache directory in the repository, so no `.gitignore` change is required. If a future development workflow introduces any repo-local runtime/cache path, `.gitignore` must change in the same batch.

### 12. Full `canFam3` ↔ `canFam4` comparative mechanical fixture — complete

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

The acquisition path now has measured size preflight plus restart-safe resumable HTTPS, so the multi-gigabyte all-chain can be acquired deliberately for this fixture without treating it as routine unit-test data. Large provider resources must remain outside the repository and must not be committed as fixtures.

Measured fixture setup on 2026-08-14:

- acquired the complete five-resource `canFam3`→`canFam4` comparative bundle into an external cache: 2,686,242,854 compressed bytes total;
- verified the 2,652,632,416-byte all-chain as SHA-256 `f10a6b48b5461bb8378ffaff311fb7355b1910511131ce0f5df5402c4db67519` and provider MD5 `8dd10ad24f866e8eb88b5442b1e26742`;
- verified the ordinary net as SHA-256 `c889134e95ff82741c0092b1673b3e5fe0125aa82f611ee97d91e468b56d51ac`;
- verified the syntenic net as SHA-256 `39ba8ca12f935755ced5eaa555b9b476460c615cc7d6c0f122ac43194a9fabce`;
- reused the previously measured reciprocal-best chain identity SHA-256 `34f4061fd29e7720c7eb2adc1ea8299e86f21f08e18f97f9ba468cf8b466690c`;
- verified the reciprocal-best net as SHA-256 `cd4891e80eaa8b8625620162306c50fc291eab6912227d78c566d9dae7fe716e`.

The first direct parse of that real all-chain exposed a format-compatibility gap: its decompressed prefix contains UCSC `#`/`##` lastz/axtChain metadata before the first `chain` record. Direct inspection of the cached bytes confirmed gzip magic `1f8b` and the leading metadata prefix. The chain parser now ignores metadata/comment lines only while looking for the next record header and still rejects them inside a chain record, preserving strict block parsing while accepting that measured UCSC prefix. The same outside-record rule intentionally permits comments between complete records and after the final record; those placements are parser policy covered by synthetic tests, not a claim that they were observed in this real resource.

After that compatibility fix, the exploratory fixture scanner streamed the complete 2,652,632,416-byte cached all-chain without another parse error, then completed its ordinary-net and reciprocal-best passes. This is end-to-end mechanical validation of reading that exact real chain resource, not biological ground truth and not evidence that every synthetically supported outside-record comment placement occurs in the provider file.

Full comparative fixture verification completed 2026-08-16 through the public `build_ucsc_candidates_from_cached_bundle()` path, using the exact cached resource identities listed above and no network access. The selected source locus is `chrUn_JH373233:1845735-1845835` in canonical 0-based, half-open coordinates. Measured production-path results:

- 170 chain-derived candidates across 114 distinct target sequences;
- chain 573 maps in reverse orientation to `chr35:925644-925938`, with two aligned segments, full 100/100 source-base coverage, one target-side chain gap, chain score 16,617,372, net `ali=3603`, `qDup=4098`, `nonSyn` classification at hierarchy depth 7, and `FULL` reciprocal-best membership covering 100/100 candidate source bases;
- chain 5170 maps to `chrUn_MU018764v1:171661-171760` with partial 99/100 source-base coverage, one gap, and `NONE` reciprocal-best membership;
- chain 2692 maps to `chrUn_JAAHUQ010000602v1:62326-62622` with full 100/100 source-base coverage, one gap, and `NONE` reciprocal-best membership;
- an independent preflight over the small reciprocal-best chain counted 3, 1, and 2 relevant source/target/orientation chain records for chains 573, 5170, and 2692 respectively, matching each production `chains_examined` result;
- primary chain 573 segment geometry was measured as source `1845735-1845808` → target `925865-925938` and source `1845808-1845835` → target `925644-925671`, with a target-side gap `chr35:925671-925865` at source boundary 1845808;
- chain, net, and reciprocal-best file provenance preserved one caller-declared shared upstream alignment ancestor. This verifies provenance wiring/dependency handling; it does not independently infer that common ancestry from the resource bytes;
- reciprocal-best completeness remains the caller's `COMPLETE_RESOURCE` claim from the prior complete acquisition, while the production file path independently SHA-256-verifies every consumed raw stream before returning candidates;
- no assessment verdict was computed and no biological ground-truth claim was made.

The reproducible verifier is `scripts/verify_canFam3_canFam4_mechanical_fixture.py`. It requires the already-acquired external cache and deliberately does not download or commit UCSC bulk resources. This completes the mechanical comparative fixture milestone; the next implementation milestone is the assessor core and deterministic verdict logic.

### 13. Assessor core and deterministic verdict logic — complete

This was the largest scientific implementation milestone remaining at that stage.

The deterministic assessor core is implemented and reviewed. It transforms normalized
candidates plus provenance-aware evidence into exactly one of `WELL_SUPPORTED`, `CONTESTED`, or
`INDETERMINATE` without a numeric score or biological-correctness claim.

The alpha verdict-driving rules are intentionally limited to categorical, locus-specific evidence
whose direction is already explicit: source-locus mapping coverage and, for `COMPARATIVE`,
reciprocal-best membership. Raw chain score, `ali`, `qDup`, net classification, and net hierarchy
remain report context and do not silently become weights or thresholds.

The implementation also:

- keeps evidence availability separate from verdict strength, including allowing a single full
  `LIFTOVER-ONLY` mapping to be `WELL_SUPPORTED`;
- treats multiple sparse-tier candidates as `CONTESTED` instead of ranking them by chain score;
- treats a single full comparative candidate with reciprocal-best `NONE` as `CONTESTED`, while
  `PARTIAL` remains `INDETERMINATE` because v1 has no threshold for deciding when partial
  self-consistency disagreement becomes material;
- never counts multiple observations as votes, and rejects duplicate verdict-driving observations
  rather than allowing duplication to strengthen a result;
- validates source geometry, coverage denominators, and reciprocal-best denominators at the
  assessor boundary so a malformed normalized candidate cannot receive a plausible-looking verdict;
- sets a preferred candidate only for `WELL_SUPPORTED`.

Adversarial tests cover sparse and comparative ambiguity, partial mappings, internal contradictory
evidence, irrelevant high-valued score/qDup context, missing/duplicate verdict evidence, multiple
fully retained candidates, 3+-candidate cases, both exhaustive reciprocal-best completeness bases,
and inconsistent normalized geometry. Focused review corrected the sparse-tier evidence-role
classification and narrowed sole-candidate `PARTIAL` reciprocal-best semantics to `INDETERMINATE`.
It also confirmed that `COMPLETE_CANDIDATE_SUBSET` is exhaustive for the generated candidates, not a
weaker arbitrary partial scan.

The final candidate-equivalence review item is also resolved. A real audit of the frozen
`canFam3`→`canFam4` fixture found zero exact-geometry and zero same-bounding duplicate groups among
its 170 candidates, but synthetic valid-chain checks demonstrated that distinct chain records can
still project one assessed locus to identical local coordinate geometry. v1 therefore treats
candidate multiplicity as hypothesis-level rather than record-level: the assessor rejects distinct
IDs with equivalent canonical local mapping geometry instead of silently merging their
provenance/evidence or manufacturing a `CONTESTED` verdict. Adjacent collinear segment partitions
are canonicalized for this check; equal target bounds with genuinely different internal mapping
geometry remain distinct candidates.

### 14. Assessment/report orchestration — complete

Implemented and reviewed:

- build the full source + target + locus → resources → candidates → evidence → verdict path;
- create the final `Assessment` object from real engine output;
- ensure evidence tier is always displayed independently from verdict;
- expose provenance/dependency detail sufficient for scientific audit;
- surface cached retrieval metadata (source URLs, retrieval timestamps, provider checksum metadata, and terms references) alongside file/evidence provenance without claiming that unconsumed bundle resources were assessed.

### 15. CLI and user-facing reports — complete

Implemented the planned common-case workflow:

```text
assess-liftover canFam3 canFam4 chrUn_JH373233:1845736-1845835
```

Implemented so far:

- UCSC database identifiers and CLI loci are parsed at an explicit boundary, with 1-based inclusive display coordinates converted immediately to canonical 0-based half-open intervals;
- concise human-readable summaries show evidence availability before verdict interpretation and always retain the biological-correctness caveat;
- the `assess-liftover` console entry point now composes discovery, terms review, HEAD transfer inspection, separate transfer-plan acknowledgement, cached acquisition, assessment orchestration, and summary rendering;
- the CLI chooses a platform user-cache default while preserving `--cache-dir`; complete verified bundles are reused cache-first without provider access, `--offline` guarantees zero network access, and `--refresh` explicitly forces a fresh provider check/acquisition;
- interactive assessment progress reports measured compressed bytes consumed for Chain, Net, and Reciprocal-best inputs with a visual bar, numeric percentage, and byte counts; `--quiet` suppresses it;
- interactive cache verification reports one measured aggregate SHA-256 row across the required cached bundle, reusing the same progress rendering primitives and withholding 100% until every required artifact passes integrity; `--quiet` suppresses it;
- interactive UCSC acquisition reports measured per-resource transfer progress, starts resumable resources at their retained validator-bound prefix, labels verified cache hits as cache reuse rather than downloaded bytes, and avoids invented percentages when exact transfer size is unavailable; `--quiet` and non-TTY stderr suppress the display;
- interactive acknowledgements are the default, with explicit `--acknowledge-ucsc-terms` and `--accept-transfer-plan` flags for non-interactive use;
- automatic UCSC runs use a conservative pair-level shared lineage node only for provenance dependency grouping, while exact consumed-file identity remains SHA-256-addressed beneath it;
- the first real CLI smoke run completed 2026-08-16 against the established external `canFam3`→`canFam4` comparative fixture, reporting `COMPARATIVE`, 170 candidates, and `CONTESTED` for display locus `chrUn_JH373233:1845736-1845835`;
- an independent fixture cross-check derives `CONTESTED` directly from the extracted public evidence without using production verdict logic, identifying 138 material candidates and agreeing with the production assessor;
- `--details` emits the full human-readable evidence dossier: exact mapped segments, categorical verdict-evidence roles, every observation, resource retrieval/checksum context, consumed-vs-unconsumed status, and the complete provenance dependency graph without implying candidate rank or independent confirmation;
- `--json` emits schema version 1 from the same assessment/report model, preserving canonical 0-based half-open coordinates, structured evidence roles/values, candidate order without ranking semantics, resource consumption/checksums/terms, provenance dependency edges, and the biological-correctness caveat.

Milestone 15 closure review completed:
- transfer-progress implementation and terminal semantics were reviewed against the code/tests with no confirmed defects; automated coverage now includes a fresh transfer whose real response omits `Content-Length`, proving the acquisition callback preserves an unknown total end to end rather than inventing a percentage;
- a real-provider transfer-progress smoke check remains desirable but non-gating and should not force a multi-gigabyte reacquisition solely for UI validation.

### 15.5. Pre-alpha semantic and output hardening — complete

This focused correctness/compatibility slice was discovered during the final project-level review of Milestone 15. It did not reopen transfer-progress work.

Implemented and reviewed:

- fixed the confirmed concise-summary bug in which an `INDETERMINATE` comparative assessment with multiple raw candidates but only one material partial candidate can be described as though the evidence failed to distinguish the candidates; the remaining uncertainty is incomplete source-locus coverage, not unresolved material multiplicity;
- added one required assessor-owned categorical `decision_reason` to every `Assessment`, using the ten exhaustive/mutually-exclusive terminal conditions specified in the pre-alpha design baseline rather than mixing biological findings with evidence-rule names;
- split the final comparative fallback into `COMPARATIVE_SOLE_MATERIAL_PARTIAL` versus `COMPARATIVE_NO_MATERIAL_CANDIDATE` using the already-defined material-candidate predicate, and regression-tested the exact boundary (`PARTIAL` coverage with reciprocal-best `FULL`/`PARTIAL` is material; reciprocal-best `NONE` is not);
- made reporting consume the recorded decision reason instead of reconstructing assessor semantics from candidate count, verdict, or evidence values;
- required every `Assessment` construction to supply a decision reason, handle every declared decision-reason enum member explicitly in reason/verdict/reporting mappings without wildcard fallback, and test that assessor coverage reaches the complete declared reason vocabulary; branch-specific regression tests must still verify that each semantic boundary selects the correct reason;
- included the required decision reason in human detail and JSON output and finalized schema v1 before alpha; private pre-alpha schema v1 remains mutable, while the first public alpha freezes it as an external compatibility surface and later incompatible structural/semantic changes require a new schema version;
- added a concise dependence qualification only for `COMPARATIVE` summaries: comparative observations are not assumed to be independent, with dependency provenance available in `--details` / `--json`. `LIFTOVER-ONLY` output should not receive that context-free qualification.

At that historical alpha-hardening point, the slice preserved the legacy three-verdict model, two evidence tiers, no-score policy, and assessor-not-resolver boundary. The routine suite now contains 316 passing tests; Ruff, Ruff formatting, strict mypy, package build, and `git diff --check` all passed. A post-hardening offline CLI run on 2026-08-17 reported `COMPARATIVE`, 170 candidates, and `CONTESTED` with the new dependence qualification and biological-correctness caveat. The independently derived mechanical-fixture verifier then passed against the same cached bundle, deriving `CONTESTED` from 138 material candidates and matching the production assessor without calling its verdict logic.

### 16. First public alpha milestone

The public-transition milestone is defined by:

- a researcher can run one documented end-to-end real UCSC assessment from the CLI;
- the comparative mechanical fixture works;
- provenance/resource identity is visible in output;
- evidence availability and support are clearly separated;
- README usage commands are real and reproducible;
- the project clearly labels itself alpha and states that well-supported does not mean biologically correct.

As of 2026-08-17, Milestone 16 is complete. v0.1.0a1 was tagged and published to PyPI as the first public alpha, and the documented external PyPI install path was verified successfully. The public-alpha compatibility period began with this release; the post-50-case owner decision below deliberately allows the upcoming redesign to break that early-alpha schema rather than preserve obsolete result semantics.

## Post-alpha redesign after the 50-case program

The first public alpha proved that liftAssess can discover/acquire UCSC resources, build candidates,
extract comparative evidence, preserve provenance, and produce end-to-end CLI/JSON output. The
completed 50-case real-world validation / UX program then established that the legacy aggregate
verdict interface is not an adequate primary explanation of what happened to a locus.

The corpus was deliberately enriched for difficult/support-question cases and controls. It is not a
prevalence, sensitivity, specificity, or general accuracy benchmark. The `38 YES / 10 MOSTLY / 2 NO`
language replay is same-corpus design evidence, not held-out validation. Likewise, counts of legacy
`WELL_SUPPORTED` misses demonstrate that the old label can conceal known failure modes in this
corpus; they are not estimates of how often arbitrary real-world mappings have hidden problems.

### Approved post-alpha policy changes

The owner review on 2026-08-19 closes the main design-policy questions:

- remove `WELL_SUPPORTED`, `CONTESTED`, and `INDETERMINATE` from the target result model; do not
  introduce a replacement one-word aggregate verdict;
- use orthogonal factual states, evidence/provenance, deterministic factual headlines, and bounded
  interpretation;
- deliberately make a pre-release machine-schema compatibility break rather than carrying the
  legacy verdict schema through the redesign;
- use progressive disclosure: keep the full result profile structured, keep uncomplicated terminal
  output compact, and expand materially unusual results;
- begin indexing/shared-traversal work immediately after the first factual renderer slice while
  assembly metadata/preflight and contextual-evidence work may proceed in parallel;
- use categorical, provenance-aware comparative relationships with no hidden numeric weighting;
  human output must explain *how* mixed/conflicting evidence differs rather than print an opaque
  `MIXED` label alone;
- initially add automatic centered 101-bp local context to 1-bp point queries once the scalable
  resource-access path makes it practical;
- make BED/simple interval-table input first-class with batch support, with batch relationships in a
  separate result layer;
- treat nonzero exit status as usage/input/operational failure, not as a scientific ambiguity/no-
  projection classifier;
- allow constrained BED12/custom-track export and Genome Browser links as visualization/navigation,
  not evidence; and
- begin the difficult-region pilot now, modeling each source as typed, provenance-bearing context
  rather than one generic warning. Start with UCSC segmental-duplication context where source,
  terms, and assembly coverage are verified; evaluate GIAB/excluderanges categories separately.

### 17. Factual result profile, new schema, and progressive renderer

**Goal:** replace the target aggregate-verdict interface without rewriting candidate generation.

Implement a dedicated derived result-profile/view-model layer over the existing scientific report
and new composite-analysis results. The profile should represent input validity, projection count,
source coverage, continuity/geometry, target role, orientation, reverse result, query-scale context,
comparative relationships, batch relationships, typed external context, evidence tier/resource
consumption, and provenance/dependence.

First-slice work:

- define the new machine schema version and explicitly document the alpha-v1 compatibility break;
- remove legacy verdict/`decision_reason`/preferred-candidate semantics from the target model rather
  than carrying them forward solely for compatibility;
- derive literal factual headlines from evidence already computed;
- add coverage/fragmentation and large-region summaries, including maximum candidate source coverage,
  uncovered bases/spans, segment count, target gaps, and alternatives;
- implement progressive-disclosure default output plus complete detail/JSON output;
- preserve the six common-use lenses as scope/model concepts without forcing six invariant terminal
  lines on every clean result;
- reserve explicit input/preflight states in the result-profile/schema so Milestone 18 can populate
  them from authoritative metadata;
- retain the no-score, assessor-not-resolver, provenance-dependence, coordinate, and target-bounding-
  span invariants;
- keep Genome Browser/locus links as optional navigation aids;
- keep a compact profile-vector string deferred unless real workflow demand appears.

The 50-case work does **not** justify a candidate-generation rewrite. Existing candidate semantics
should change only in response to a separately measured defect.

### 18. Start scalable resource access; assembly metadata and difficult-region context in parallel

**Trigger:** begin this work immediately after the first Milestone-17 renderer/profile slice is
working. Do not wait until reverse/neighborhood/batch features have already multiplied scans.

Measured performance background (2026-08-17): on the tested implementation, a 2.47-GiB
canFam3→canFam4 all-chain single-locus run took roughly 10.5 minutes on an Apple M4 regardless of
whether the final legacy result was complex or simple, while profiling showed Python-level chain
parsing/object construction dominating. The smaller-resource probes scaled roughly with chain size,
and mapped comparative processing added material additional cost. These measurements establish the
problem, not the final architecture.

Required work:

- prototype region-addressable/indexed/shared-traversal approaches and benchmark them against the
  frozen performance probes;
- preserve exact coordinate, candidate, evidence-completeness, resource-identity, and provenance
  semantics;
- add authoritative assembly-sequence metadata for source-name validation, bounds, aliases, and
  target role; chain-file names alone are not sufficient;
- reject invalid source names and out-of-range coordinates before scientific mapping; expose
  reusable preflight metadata for later BED/batch intake;
- begin the typed difficult-region pilot, first against the duplication/paralogy cases using UCSC
  segmental-duplication context when its source/terms/assembly coverage are verified;
- evaluate GIAB stratifications and relevant `excluderanges` categories separately rather than
  treating them as interchangeable.

No particular index design is frozen. Local interval index, header index, compact parsed
representation, shared traversal, or another design remains eligible until measured prototype
results select one.

### 19. Actual reverse-mapping context

Add reverse assessment as its own structured result dimension.

- distinguish `returns to original source`, `returns elsewhere`, `unavailable`, and `not run`;
- never relabel current UCSC reciprocal-best membership as actual reverse mapping;
- treat non-reciprocity as context, not automatic proof that the forward projection is wrong;
- reuse the scalable resource-access architecture rather than launching another exhaustive scan per
  candidate/query.

### 20. Point neighborhood / multi-scale context

Add automatic local context for 1-bp point queries:

- centered 101-bp window (±50 bases when source bounds permit);
- exact tested window always reported;
- 101 bp described as a product/context default, not a confidence threshold or universal biological
  scale;
- ordinary interval queries are not automatically widened;
- explicit larger-context controls remain available;
- a point/context disagreement is reported directly and does not silently trigger recursive 1-kb/
  10-kb widening.

Corpus B provides six matched clean 101-bp controls showing that widening can remain neutral when
local geometry is uncomplicated; difficult cases show that context can also expose fragmentation or
cross-record relationships. Revisit the default after held-out/outside-user testing.

### 21. Filtered/all-chain comparison and comparative relationships

Make ordinary filtered liftOver versus all-chain candidate inventory explicit comparative context.

Initial categorical interpretation must support at least:

- comparative evidence favors one placement;
- comparative evidence does not separate placements; and
- comparative evidence is mixed/conflicting.

The first accepted `favors one placement` pattern is the B14-style relationship: multiple
full all-chain placements, exactly one retained by the ordinary filtered chain, that same
placement top-net + full reciprocal-best, and no competing full placement with equivalent
categorical top-net + full reciprocal-best support.

Do not create hidden weights from chain score, `ali`, `qDup`, net hierarchy, or reciprocal-best
membership. `ali`/`qDup` remain descriptive until a separately justified deterministic rule exists.
When evidence is mixed/conflicting, human output must list the material relationship—what the
filtered chain retained, what net/rbest supports, and what conflicts—instead of stopping at a vague
label. Shared UCSC alignment lineage remains explicit; these observations are not independent votes.

### 22. BED/table batch input and cross-record relationships

Add first-class BED/simple interval-table batch assessment only on top of shared traversal/indexed
resource access.

- validate BED semantics at intake, including zero-width/empty intervals;
- retain the existing single-locus CLI syntax;
- represent exact target collisions separately from overlapping-but-offset projections;
- keep batch relationships in a batch result layer rather than candidate-level evidence;
- support neighborhood-level collision relationships when point-context analysis is available;
- preserve per-record exact resources/provenance and deterministic results;
- do not implement batch as a naive outer loop that reparses multi-gigabyte resources per row.

Optional BED12/custom-track export may visualize one candidate whose blocks share a target sequence;
it must not collapse multiple candidates/target sequences or replace source-coverage reporting.

### 23. Held-out result-language and outside-user gate

Before describing the redesigned language as validated or release-ready:

- run a small set of real cases not used to derive the language;
- include uncomplicated controls and difficult/ambiguous examples;
- exercise the 101-bp point-context default and typed contextual observations;
- obtain outside-user/domain feedback on whether the factual headline, expanded unusual-case output,
  comparative explanations, and scope boundaries answer the practical question without implying
  biological certainty.

This gate tests communication/usefulness, not a numeric confidence model. The existing 50 cases
remain same-corpus design evidence.

### Parallel / non-blocking research items

These remain useful but should not disrupt the sequence above:

- identify one concrete independently established `canFam3.1` → `canFam6` same-individual locus to
  turn the historical-resolution pedigree into a truth-bearing sanity fixture;
- define portable case manifests and any later byte-containing packets under provenance and
  redistribution-term constraints;
- select and verify a real optional flanking-gene synteny/orthology source;
- build a curated worked-examples gallery from the 50 cases for onboarding, clearly labeled as
  examples rather than independent validation;
- revisit candidate-rank evidence only when a defensible locus-scoped semantics exists.

## Deliberately deferred beyond the current redesign

The following should not distract from making the first version scientifically useful:

- plugin registry, entry-point discovery, or general engine configuration framework;
- automatic support for many candidate-generation engines before a real second engine exists;
- a second mapping-evidence source/provider until a concrete source is verified to come from a genuinely independent upstream mapping/alignment process and the cross-engine hypothesis-equivalence semantics are defined; when that exists, independent agreement/disagreement becomes an evidence capability rather than a reason to add speculative plugin infrastructure;
- fresh minimap2/lastz alignment by default;
- machine-learning confidence models;
- composite numeric confidence scores;
- automatic biological orthology claims;
- large bundled species databases;
- hosted service infrastructure;
- compact profile-vector strings unless a real workflow demonstrates need for a second compact
  compatibility surface;
- VCF-specific vocabulary/API commitments until a VCF data model exists.

## Development and review discipline

Each substantial batch should continue to follow the current pattern:

1. implement one narrow capability;
2. cover its failure modes with tests;
3. pass pytest, Ruff, strict mypy, and `git diff --check`;
4. perform a targeted review before commit when the batch changes scientific, coordinate, provenance, completeness, or interpretation semantics;
5. verify reviewer findings against the actual code/design before accepting proposed fixes;
6. commit only after confirmed findings are resolved and gates pass again.

The purpose of this discipline is not to maximize test count. It is to make each scientific claim in the software reconstructable from explicit code, provenance, tests, and primary-source checks.
