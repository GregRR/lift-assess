# LiftOver Ambiguity Assessor — liftAssess - Design Baseline

Status: post-alpha design baseline; redesign approved after the 50-case real-world validation / UX
program. This document remains the authoritative scientific and architectural specification, not a
development-status log. See [`ROADMAP.md`](ROADMAP.md) for implementation history, current review
state, and planned milestones. The released v0.1.0a1 implementation still contains the legacy
three-verdict model; historical sections may describe that behavior, but the target design defined
here no longer uses an aggregate verdict taxonomy.

## 1. Problem

Coordinate liftover (UCSC liftOver, CrossMap, rtracklayer, NCBI Remap, BCFtools/liftover,
Liftoff, and others) converts coordinates between assemblies but does not explain confidence
when the result is ambiguous: multiple candidate targets, a placement on an unplaced/alt
scaffold instead of a named chromosome, a duplicated-region mapping, a split interval, or two
tools disagreeing on the same input. Existing tools compete on conversion completeness (fewer
dropped or mis-lifted records). We have not identified an existing tool that ships a transparent,
evidence-based explanation for the records that remain genuinely ambiguous after conversion.
That gap is what this tool fills — it sits downstream of / alongside a converter, not in
competition with one.

### Evidence the problem is real (Gate 1 — validated)
- NAR Genomics and Bioinformatics (2020): benchmark of 6 liftover tools (UCSC liftOver,
  rtracklayer, CrossMap, NCBI Remap, flo, segment_liftover) on WGBS/ChIP-seq data. Failures
  concentrated in mapping to alternative chromosomes.
- Bioinformatics (2024), BCFtools/liftover paper: benchmark of 6 variant-liftover tools against
  45,595,458 real 1000 Genomes SNVs. Dropped/discordant record counts: Transanno 31,436
  (mostly multi-region mapping), CrossMap/VCF 46,607, Picard/LiftoverVcf 21,503, Genozip/DVCF
  21,068, GenomeWarp 22,251, BCFtools/liftover 12,582 (best of the six). Documents tool-vs-tool
  discordance on the same input.
- Liftoff (Bioinformatics, 2021): built specifically because standard liftover strategies split
  or mis-map gene/exon structure when an interval isn't contiguous in the target assembly.

### Evidence people hit this and ask the exact question (Gate 2 — validated)
- Biostars: canFam4 (dog) liftOver disagreement between an unplaced scaffold and chr16 for a
  miRNA-region SNP set — "which locus is the true ortholog."
- Biostars: human liftOver landing on one paralogous FAM72/SRGAP2 neighborhood while the
  historical investigation favored another; responders used flanking genes and reciprocal mapping
  to interpret the discrepancy — exactly the kind of context this tool should expose without
  claiming biological truth.
- Bioconductor support: rtracklayer liftOver silently mapped a SNP to an unrelated "random"
  position; maintainer had to explain liftOver has no knowledge of what's being mapped.
- Bioconductor support: rtracklayer and UCSC's own web liftOver gave different results (one
  clean region vs. 1,072 fragments) for the identical input and chain file, due to differing
  internal gap-handling — a case where even picking "the standard tool" doesn't resolve which
  answer to trust.
- GitHub (GIAB remap repo): CrossMap-based liftover introduced reference/alternate allele
  inversions caught only during downstream validation.

### Important scope caution (do not overstate)
The ~12,582 SNVs BCFtools/liftover still drops are evidence the underlying problem is real, not
12,582 cases this assessor could rescue. Many failures (chain gaps, representation changes) have
no competing candidate to evaluate — nothing for an assessor to do. The assessor's actual target
population is narrower: multi-mapped loci, alt/chrUn placements, duplicated-region mappings,
split mappings, and tool-vs-tool discordance — cases where there is evidence to examine.

## 2. Philosophy

- **Assessor, not resolver.** Never outputs "correct locus = X." It reports what physically happened
  to the requested interval, what evidence supports those observations, how that evidence depends on
  upstream sources, and the strongest bounded interpretation those facts justify. A clean coordinate
  projection is not a claim of biological truth or identity.
- **No aggregate verdict taxonomy.** The current design does not use `WELL_SUPPORTED`, `CONTESTED`,
  `INDETERMINATE`, or a replacement one-word verdict. Those labels are legacy v0.1.0a1 behavior.
  The redesigned result consists of orthogonal factual states, evidence/provenance, a deterministic
  factual headline, and bounded interpretation.
- **No composite numeric score.** A number like "0.91" claims a precision that has not been validated
  against ground truth. Show categorical facts and relationships instead. Do not hide weighted
  combinations of chain score, `ali`, `qDup`, net hierarchy, reciprocal-best membership, or other
  observations behind prose labels.
- **Evidence carries provenance, not just a result.** Every observation records its upstream source
  (e.g. "UCSC canFam3→canFam6 alignment," "Ensembl Compara release X," "independent minimap2
  alignment run locally"). This lets the report detect when multiple observations trace back to the
  same underlying alignment and must not be presented as independent confirmation. Chain score, net
  type, qDup, and reciprocal-best membership, when derived from the same UCSC alignment, are related
  observations, not independent votes.
  - Nets are derived from chains (UCSC's own pipeline: chains → chainNet → net). Not independent.
  - UCSC's standard reciprocal-best chains are typically produced by swapping/filtering the same
    alignment lineage, not by an independent second alignment run. Reciprocal-best membership is a
    self-consistency observation. Actual reverse mapping is a separate result dimension (§4).
  - Typed digest model for provenance identifiers (do not misapply refget to non-sequence files):
    - reference sequence → GA4GH refget sequence digest (SQ....)
    - sequence collection → GA4GH SeqCol identity
    - chain / net / GTF / any other file → plain `sha256:<bytes>`
  - Independence is established from provenance case by case, never assigned to an evidence type in
    advance.
- **Interpretation stays one rung below evidence.** E.g. `qDup` + `nonSyn` may support the literal
  description *"duplication-associated nonsyntenic placement"*; they do not by themselves prove a
  specific mechanism such as paralogy, processed pseudogene, or assembly error. Promote a named
  mechanism only when additional evidence supports it.
- **Neutral structural facts stay neutral.** Reverse orientation, interchromosomal projection,
  alternate targets, and non-reciprocity can be important context but are not automatic error or
  quality labels.
- **Same species remains the operational envelope for the present UCSC workflow.** Different
  individuals' assemblies may legitimately differ biologically, so disagreement is not automatically
  error. Same-individual pairs remain especially valuable for truth-bearing validation fixtures.

## 3. Scientific invariants

A short checklist to test any future implementation decision against:

1. Never count observations that share an upstream alignment or source as independent confirmation.
2. Absence of evidence for a competing candidate is not proof that no competing candidate exists.
3. Evidence availability is not confidence. `COMPARATIVE` and `LIFTOVER-ONLY` describe what could be
   examined, not how "good" a result is.
4. Different-individual assemblies may legitimately differ biologically; a disagreement between them
   is not automatically an error.
5. A factual headline describes the mapping/evidence state, not biological truth.
6. Biological interpretation sits one rung below evidence and is promoted to a named mechanism only
   when independently supported.
7. Input validity precedes scientific assessment. Invalid sequence names, out-of-range coordinates,
   and empty future BED intervals must not degrade into biological-looking no-projection results.
8. A target bounding span is a summary only; it must never imply continuity when the aligned segments
   are fragmented.
9. Candidate encounter order is reproducibility-only, not rank. Raw chain score is not a substitute
   for a scientifically defined candidate-rank concept.
10. Result dimensions are not votes and are not collapsed into a hidden score or single verdict.

## 4. Result model and deterministic interpretation

### 4.1 No aggregate verdict

The redesigned result model intentionally has **no single aggregate verdict**. The legacy alpha
labels `WELL_SUPPORTED`, `CONTESTED`, and `INDETERMINATE` are removed from the target model rather
than renamed or preserved as a parallel policy layer.

A completed result is a structured profile of factual dimensions plus evidence/provenance and a
bounded deterministic interpretation. Human rendering may summarize that profile, but it must not
recreate a hidden aggregate verdict under another name.

### 4.2 Orthogonal result dimensions

The profile must be able to represent, as applicable:

1. **Input validity** — valid input; unrecognized source sequence; out-of-range coordinate; empty
   input; or another explicit preflight failure.
2. **Projection count** — no chain projection; one chain projection; multiple chain projections.
3. **Source coverage** — complete or partial coverage with exact covered/source-base counts and
   uncovered source intervals.
4. **Continuity / geometry** — contiguous, fragmented, target-discontinuous, exact mapped segments,
   target bounding span, source gaps, and target gaps.
5. **Target role** — primary, alternate, unplaced, or another authoritative assembly role when a
   defined metadata source supports the label.
6. **Orientation** — same or reverse; mixed across a candidate set when relevant.
7. **Reverse mapping context** — actual reverse assessment, distinct from reciprocal-best
   membership. Record check state (`NOT_RUN`, `UNAVAILABLE`, or `RUN`) separately from the
   observed return relationship. A completed run records whether reverse projections are absent,
   touch only the original aligned source geometry, land only elsewhere, or do both; exact
   returned/source-base coverage and exact reconstruction of the original aligned geometry remain
   separate factual fields. For fragmented forward candidates, reverse each exact mapped target
   segment rather than the target bounding span across unaligned gaps.
8. **Query-scale context** — point result compared with its automatic 101-bp context and any explicit
   larger context requested by the user.
9. **Comparative relationship** — filtered/all-chain agreement, additional all-chain placements, and
   categorical relationships among filtered-chain, net, reciprocal-best, and other comparative
   observations.
10. **Batch relationship** — exact target collision, overlapping target projection, neighborhood-level
    collision, or no detected relationship within the assessed batch.
11. **Typed external context** — resource-specific segmental-duplication/difficult-region,
    assembly-history, variant identity, gene/transcript identity, or other optional observations when
    explicitly implemented.
12. **Evidence tier / resource consumption** — `LIFTOVER-ONLY` or `COMPARATIVE`, exactly which
    resources were consumed, and what was not assessed.
13. **Provenance / dependence** — exact resource identities and dependency relationships.

These dimensions can coexist. No single field is expected to encode the whole scientific story.

### 4.3 Deterministic factual headlines

Human output leads with a literal mapping/event headline derived mechanically from the profile. The
headline vocabulary should describe what happened, for example:

- `NO CHAIN PROJECTION`
- `ONE COMPLETE CHAIN PROJECTION`
- `PARTIAL SOURCE COVERAGE`
- `PARTIAL AND FRAGMENTED PROJECTION`
- `COMPLETE BUT DISCONTINUOUS PROJECTION`
- `MULTIPLE CHAIN PROJECTIONS`
- `SOURCE INTERVAL SPLITS ACROSS MULTIPLE PROJECTIONS`

Exact final token spelling is a schema/renderer implementation detail, but headline semantics must be
factual and testable. Large-region output gives source coverage/fragmentation precedence over raw
candidate count.

### 4.4 Comparative interpretation

Comparative interpretation is categorical and provenance-aware, not weighted. The initial model must
support at least:

- comparative evidence **favors one placement**;
- comparative evidence **does not separate placements**; and
- comparative evidence is **mixed/conflicting**.

The first accepted `favors one placement` pattern is the B14-style case: multiple full all-chain
placements exist; exactly one is retained by ordinary filtered liftOver; that same placement is
top-net and has full reciprocal-best membership; and no competing full placement has equivalent
categorical top-net + full reciprocal-best support.

`ali` and `qDup` remain reported observations initially, not hidden weights or thresholds.

A human-facing result must never stop at a vague statement such as `COMPARATIVE EVIDENCE IS MIXED`.
It must state *how* the evidence is mixed: which placement is retained by the filtered chain, which
has top-net support, which has reciprocal-best membership, and which observations conflict or fail
to separate candidates. Shared UCSC lineage remains explicit; these are not independent votes.

### 4.5 Point-query local context

Once the scalable resource-access path in §7 exists, a 1-bp point query automatically receives a
centered 101-bp local-context assessment (±50 bases when source bounds permit). The exact tested
window is always reported.

The 101-bp default is a product/context choice, not a confidence threshold, biological universal, or
calibrated scale. Ordinary interval queries are not automatically widened. Users may request larger
context explicitly. If the 101-bp result differs materially from the point result, report that
relationship; do not silently recurse through increasingly large windows.

### 4.6 Common-use evidence boundaries

The result model preserves six common-use lenses that help prevent overclaiming:

- sequence-coordinate projection;
- named variant / rsID identity;
- gene/transcript identity;
- interval coverage/structure;
- batch relationships; and
- file/downstream workflow.

These are model/scope concepts, not six mandatory terminal lines. Detailed/JSON output must preserve
clear `NOT TESTED`, `NOT CHECKED`, or `NOT ASSESSED` boundaries whenever a domain was not evaluated.
Default human output uses progressive disclosure (§9).

## 5. Coordinate semantics

Not addressed in earlier drafts of this design and worth getting right before any code exists.

- **Canonical internal representation: 0-based, half-open intervals** — matches both the BED
  convention and the chain file convention, which is what the tool primarily consumes.
- A genomic interval is an **unstranded span in forward-reference coordinates**. Mapping
  orientation is a separate relation on a candidate mapping: `SAME` when source and target are
  aligned in the same orientation, `REVERSE` when they are aligned in opposite orientations.
  Do not attach alignment orientation to the interval itself. For a reverse-strand UCSC chain
  query, normalize the query coordinates back to forward-reference coordinates before creating
  the target interval.
- A split mapping produced by one chain remains **one candidate**, not several independent
  candidates. Preserve its exact aligned structure as ordered source/target mapping segments.
  `NormalizedCandidate.target_interval` is the smallest forward-reference span containing all
  mapped target segments; it is a summary/bounding span only and must never imply that bases
  between separated segments are aligned.
- The core interval type can represent an empty half-open span, but v1 chain candidate generation
  does not assign semantics to zero-length source intervals. Reject them rather than implicitly
  borrowing UCSC liftOver's special zero-width handling; add such support only after the desired
  point/insertion semantics are specified explicitly.
- Every input parser converts explicitly into this representation at the boundary. Nothing
  internal ever guesses a convention.
- CLI-typed locus strings (e.g. `chr16:12345-12400`) follow the common UCSC-style display
  convention (1-based, inclusive) — matching what a researcher would type or paste from the
  browser — and are converted to internal 0-based half-open immediately on parse.
- BED input retains its native 0-based half-open convention as-is. VCF input retains its native
  1-based convention as-is, with explicit conversion at the same boundary.
- All internal computation, comparison, and evidence extraction happens in the canonical 0-based
  half-open representation.
- Any coordinate shown in output states which convention it's displayed in. No implicit guessing
  on the way out either.

## 6. Current target scope

**In scope — target mapping/evidence capabilities:**

- **Mapping structure:** source-locus coverage, exact mapped segments, uncovered source intervals,
  source/target chain gaps, target bounding spans, orientation, and candidate multiplicity.
  `FULL`/complete coverage means every requested source base is represented by an aligned mapping
  segment; partial coverage means one or more requested source bases are not aligned. Preserve exact
  geometry rather than reducing it to a quality-like number.
- **UCSC chain/net comparative evidence:** chain score, aligned bases (`ali`), duplicated query bases
  (`qDup`), net classification/hierarchy, and reciprocal-best membership. Reciprocal-best membership
  remains locus-specific to the candidate's aligned source bases and can be `FULL`, `PARTIAL`, or
  `NONE` only after a complete relevant resource scope has actually been checked. It is a
  self-consistency observation unless provenance demonstrates an independent upstream lineage.
- **Filtered-chain versus all-chain comparison:** ordinary directional liftOver output and the broader
  all-chain candidate inventory become explicit comparative context when both are available.
- **Actual reverse mapping context:** a reverse assessment is a separate capability from
  reciprocal-best membership. It preserves check availability, no-projection state, original-source
  overlap/coverage, elsewhere returns, and exact original-geometry reconstruction as distinct facts.
  Fragmented forward mappings are reversed through their exact mapped target segments; a target
  bounding span must never manufacture a reverse query across an unaligned gap.
- **Point local context:** 1-bp queries gain the approved automatic 101-bp neighborhood assessment
  once indexing/shared traversal makes the additional work practical.
- **Assembly-sequence metadata/preflight:** source sequence validity, source bounds, authoritative
  aliases, and target sequence role belong to one coherent metadata capability rather than ad-hoc
  naming heuristics. Chain-file names alone are not authoritative proof that an assembly sequence is
  valid.
- **Batch relationships:** exact collisions and overlapping target projections are relationships among
  records and belong in a separate batch result layer. BED/simple interval-table input becomes
  first-class with batch support.
- **Typed contextual evidence:** optional difficult-region/context resources are represented as
  resource-specific, provenance-bearing observations. The first active pilot is UCSC segmental-
  duplication context; GIAB stratifications and `excluderanges` categories must be evaluated
  separately rather than collapsed into one generic warning.
- **Optional genomic-context evidence:** flanking-gene orthology/synteny or other context may be added
  when a real source is selected and its provenance/dependencies are explicit.

**Explicitly out of the base mapping layer:**

- Freshly computed sequence identity from raw bases.
- A new alignment run (minimap2/lastz) by default.
- A numeric composite confidence score or hidden weighted evidence formula.
- Machine learning.
- Automatic claims of biological orthology or "correct locus" identity.
- Automatic dbSNP/variant identity unless a separate identity module is explicitly invoked.
- Automatic gene/transcript equivalence unless a separate annotation module is explicitly invoked.
- VCF allele/normalization semantics until a VCF-specific data model exists.
- Downstream-tool/file diagnosis unless a dedicated workflow module or user-supplied evidence
  establishes it.
- Large bundled species-support databases.
- Hosted infrastructure.

Candidate rank remains deferred until a defensible locus-scoped definition exists. Raw chain-score
order, UCSC output order, or encounter order must not be silently promoted into scientific rank.

## 7. Architecture

```
Candidate-generation / evidence engine
        │
        ▼
NormalizedCandidate[] + evidence + provenance
        │
        ▼
Scientific report / composite-analysis results
        │
        ▼
Derived result profile (orthogonal factual states)
        │
        ├──────────────► schema-versioned machine output
        │
        ▼
Progressive-disclosure human renderer
```

- The scientific core knows nothing about how candidates were produced beyond the normalized
  candidate/evidence/provenance contract. Source-specific mechanical observations that cannot be
  reconstructed losslessly after normalization must therefore be extracted while the source record
  is still available and carried forward with provenance. Chain-gap evidence is one such case.
- **The result profile is a dedicated derived boundary.** It converts scientific report/composite
  analysis data into the orthogonal factual states defined in §4. Human renderers consume that
  profile rather than rediscovering semantics independently, and user-facing prose does not belong
  in core candidate/evidence dataclasses.
- **The current implementation has exactly one candidate-generation/evidence engine**: an internal,
  minimal chain/net reader against the documented UCSC chain/net format. Not pyliftover
  (point-converter only, no block/gap detail) and not the UCSC liftOver binary (see §8, licensing).
- **Chains generate candidate mappings; nets annotate and evaluate those candidates.** These are
  distinct responsibilities. For UCSC-generated liftOver map chains specifically, the old/source
  assembly is the chain target (`t*`) side and the new/destination assembly is the query (`q*`)
  side; candidate generation follows that documented direction explicitly rather than inferring it
  from assembly names. The LIFTOVER-ONLY evidence tier has no net file and must still be
  able to generate and report candidates from chain data alone. The implementation may keep both in
  one reader package, but the two responsibilities stay conceptually — and ideally structurally —
  separate, so net availability is never accidentally required for candidate generation itself.
- **UCSC engine orchestration consumes each external record stream once.** Candidate generation
  streams chains once; net records are filtered during one pass to fills associated with generated
  candidate chain IDs; reciprocal-best chains are filtered during one exhaustive pass to the
  source/target sequence and orientation combinations relevant to those candidates. This avoids
  rescanning large comparative resources once per candidate or materializing whole-genome files
  merely to make iterators reusable. Per-candidate matching still applies the stricter geometry and
  provenance rules defined above; the orchestration layer does not rank candidates or manufacture
  an aggregate result verdict.

- **Region-addressable/shared resource access becomes enabling architecture immediately after the
  first result-profile/renderer slice.** Milestone-18 prototype/benchmark work selected a reusable
  chain index with 65,536-bp source-coordinate bin memberships and encounter-order chain records
  stored exactly once in independently compressed blocks. The derived index is keyed to the exact
  source-chain SHA-256 and remains an acceleration artifact, not replacement scientific evidence.
  Automatic reverse mapping, point-context comparison, filtered/all-chain comparison, expanded
  comparative analysis, and batch assessment must not ship by multiplying the current exhaustive
  whole-resource scan. Indexed/preprocessed forms must preserve exact coordinate,
  evidence-completeness, resource-identity, provenance, and reproducible candidate-order semantics.
  A validated exact-resource index may be used transparently for chain-record access when present;
  the CLI falls back to the existing verified full traversal when that derived acceleration
  artifact is absent or unusable. Lower-level library callers receive `ChainIndexCorruptionError`
  for query-time index corruption and choose their own recovery policy. Integrity work must follow
  the bytes actually used: index build/rebuild verifies
  the complete original chain and records a full database SHA-256 for explicit deep verification;
  normal indexed assessment verifies a compact lookup catalog, validates each queried 65,536-bp
  bin's membership/record-locator rows against catalogued SHA-256 digests, and verifies selected
  compressed record blocks before parsing. It therefore does not reread either the unused original
  chain or the complete SQLite lookup database on every query. The other cached bundle artifacts
  retain their normal direct SHA-256 checks. This distinction matters on slow storage: a measured
  cold cache-bundle verification on an older iMac HDD took 95.976 seconds, while the same bundle
  verification on the M4 Mac mini took about 1.2 seconds. The integrity contract therefore remains
  exact for the query-relevant derived data without forcing storage-sensitive whole-artifact reads.
  Index construction remains an explicit user action rather than an implicit first-query pause.
- **Local resource-file integration stays a thin streaming boundary.** Plain-text and gzip-compressed
  chain/net files are opened locally and fed directly through the existing parsers into the UCSC
  engine; they are not materialized as whole-file objects in memory or rescanned once per candidate.
  This layer does not download files or infer provider permissions. Exact local file artifacts are
  identified by SHA-256 over their raw on-disk bytes and represented through the existing
  `ProvenanceIdentifierKind.SHA256` model; the content-addressed file provenance node derives its
  structural source ID from that same digest rather than from a filename or path. Provider-published
  checksums (including MD5 where that is what the provider publishes) are separate transfer-integrity
  checks and are not promoted into liftAssess provenance identity. Upstream alignment/process
  provenance remains caller-supplied because file bytes cannot reveal how the artifact was produced.
  The acquisition/cache layer computes SHA-256 during transfer while recording source URL, retrieval
  metadata, applicable provider terms, and any provider checksum verification. Single-resource
  acquisition and complete bundle planning/execution are implemented against a caller-supplied cache
  root; no repository-local cache is created. After explicit provider-terms acknowledgement, a
  separate metadata-inspection step may issue body-free HTTP HEAD requests to expose
  provider-advertised `Content-Length` and transport headers before transfer-plan acknowledgement.
  Live provider checks on 2026-08-14 verified exact HEAD sizes plus byte-range/`If-Range` behavior for
  the canFam3→canFam4 comparative resources. The acquisition layer now uses those semantics
  opportunistically: when an exact provider checksum, identity-encoded `Content-Length`, explicit byte-range support, and a strong
  `ETag` are available, interrupted HTTPS transfers retain a URL/size/validator-bound partial and resume
  with `Range` + `If-Range`. A response that fails the expected 206/`Content-Range`/validator contract is
  never appended to that prefix and instead restarts through the fresh streaming path. The shared partial
  itself is never promoted into the content-addressed store: completion is copied through the process's
  open handle into a unique private snapshot, provider MD5 and SHA-256 are recomputed over that snapshot,
  and only the private file is atomically published. This prevents another process holding the shared
  partial inode open from mutating an already-published artifact. When the required metadata is absent,
  acquisition remains non-resumable rather than inventing a weaker validator.
  A fully acquired cached bundle can now bridge directly into the same file-backed engine without
  re-downloading resources or pre-hashing them a second time: cache-recorded SHA-256 values create
  the file provenance nodes, while the parser still hashes every consumed raw stream and rejects a
  post-acquisition mutation before candidates return. The bridge requires caller-supplied upstream
  alignment provenance and validates the bundle's UCSC database strings only against explicit
  assembly names/aliases; it does not create a general alias resolver. For `COMPARATIVE`, all five
  acquired files remain on the bundle, but the v1 engine consumes the all-chain, ordinary net, and
  reciprocal-best chain. UCSC's current `doBlastzChainNet.pl` produces the ordinary net through
  `chainNet`/`netSyntenic` and `netClass`, while the optional `*.syn.net.gz` is a subsequent
  `netFilter -syn` derivative; substituting the syntenic net would therefore discard exactly the
  non-syntenic placements that the assessor needs to preserve. The reciprocal-best net likewise
  remains retained provider/retrieval context while membership is computed from the published
  reciprocal-best chain geometry. Primary implementation references:
  https://raw.githubusercontent.com/ucscGenomeBrowser/kent/refs/heads/master/src/hg/utils/automation/doBlastzChainNet.pl
  and https://raw.githubusercontent.com/ucscGenomeBrowser/kent/refs/heads/master/src/hg/utils/automation/doRecipBest.pl.
- **Reciprocal-best completeness metadata survives orchestration filtering unchanged.** The
  engine's source-sequence, target-sequence, and orientation filter is lossless for candidate
  membership: it removes only reciprocal-best chains that the downstream geometry matcher would
  reject for that candidate anyway. Therefore a caller's `COMPLETE_RESOURCE` or
  `COMPLETE_CANDIDATE_SUBSET` claim continues to describe the scope that was exhaustively checked
  and is not rewritten merely because the engine materializes a candidate-relevant subset
  internally. `chains_examined` is relevance-filtered audit context, not evidence strength or
  proof that the caller's completeness claim is true.
- **Net fills are not one-to-one with chains.** UCSC net `fill` records may carry the associated
  chain ID, but the same chain ID can occur on more than one fill in one net (including at
  different hierarchy positions). Candidate annotation must therefore use the chain association
  together with target-side overlap against the candidate's actual aligned source segments and
  preserve the hierarchy context of every matching fill; it must not treat `chain_id` as a unique
  net-record key or choose a single "best" fill. UCSC defines a fill as a portion of a chain, and
  its own `netToAxt` implementation subsets that chain using the fill's target start/end. This is
  also directly visible in UCSC's published net-format example.
- **Net metrics stay fill-scoped.** `ali`, `qDup`, net classification, and hierarchy depth are
  preserved exactly for each relevant fill rather than aggregated across repeated fills or
  reinterpreted as locus-specific quantities. Hierarchy depth is context, not confidence.
  Observations from one fill share a fill-level provenance node derived from the net resource;
  the net resource must share upstream provenance with the chain/alignment source so these
  observations cannot be presented as independent confirmation.
- **Net query coordinates are preserved but not yet used for genomic comparisons.** The current
  parser stores each net record's reported query start/span and relative orientation, but v1
  candidate matching intentionally uses only the target-side fill span plus the candidate's exact
  aligned source segments. Do not treat the stored net query coordinates as forward-reference
  genomic coordinates until their strand/origin semantics have been verified against the kent
  chainNet implementation. This is a deferred verification item, not an inferred convention.
- Reconstruct candidates independently from chain/net resources rather than only grading a
  converter's already-filtered output (a converter may have discarded a relevant alternative
  before the assessor ever sees the locus). A converter's flag ("this looks suspicious") is a
  trigger to investigate, not the candidate set itself.
- Pluggable is a real requirement for the *interface boundary* only. It does **not** mean a
  plugin manager, registration API, entry points, config system, or dynamic discovery in v1. No
  piece of "extensibility" belongs in v1 unless it makes the one real engine or the assessor
  itself simpler, clearer, or more testable right now. Add adapters (UCSC liftOver, pyliftover,
  CrossMap, BCFtools/liftover) later, against concrete requirements, when a genuine second engine
  exists.
- If multiple engines are ever run against the same locus, report their agreement/disagreement
  explicitly, and if they share an underlying data source, label that agreement as methodological
  consistency, not independent evidence — same provenance discipline as everywhere else. A future
  second mapping-evidence source counts as genuinely independent only after its upstream
  mapping/alignment lineage has been verified as independent of the first source and the project has
  defined how equivalent local mapping hypotheses from different engines are recognized; provider
  plurality alone is not evidence independence.
- **Assembly identifiers**: v1's CLI and resource resolver accept UCSC database identifiers for
  automatic discovery (e.g. `canFam3`, `canFam4`). The internal assembly representation should
  still be structured to record provider, accession, and known aliases (e.g. the UCSC db name
  `canFam3` versus the biological/NCBI lineage name `CanFam3.1`) even though v1 does not build a
  general alias-resolution system. Leaving this unrecorded would let exactly the kind of hidden
  naming ambiguity this project exists to catch creep into the project's own data model.
- **Resource resolver**: given source + target assembly names, discover and verify which evidence
  tier is actually available. UCSC's Golden Path layout is a reasonable heuristic for where to
  look first, but availability must be confirmed by checking, not inferred purely from
  constructing a plausible URL — the directory conventions are not a guaranteed, permanent
  contract:
  1. Look for `source/vsTarget/` comparative chain, net, and syntenic-net resources. For
     reciprocal-best, check that directory's `reciprocalBest/` first; if the exact directional
     `source.target.rbest.{chain,net}.gz` files are absent, also check the sibling/reverse
     `target/vsSource/reciprocalBest/` publication location and accept it only when those exact
     filenames are actually observed. A complete set → **COMPARATIVE** evidence tier.
  2. Else look for `source/liftOver/sourceToTarget.over.chain.gz` → **LIFTOVER-ONLY** tier.
  3. Else → unavailable; require user-supplied resources.
  4. Always tell the user which tier is in play, in plain language, before showing results.
  5. Always accept user-supplied chain/net resources directly — UCSC is a convenient default
     provider, not a hard dependency.
- **Resource acquisition/cache**: discovery and retrieval remain separate operations. v1 acquisition
  requires explicit acknowledgement of the applicable UCSC/general and directory-specific terms before
  network access, distinguishes restricted `liftOver/*.over.chain.gz` resources from comparative
  `vsTarget/` resources, verifies provider checksum metadata when an exact filename entry is available,
  validates transport length metadata when supplied, and stores exact downloaded bytes outside the source
  tree by liftAssess SHA-256. A provider checksum is an integrity check, not provenance identity. Cache
  reuse is based on the retained URL index plus exact content identity and must not be described as
  proof that the remote URL is unchanged; cached provider artifacts retain local SHA-256
  verification except that a validated exact-resource derived index may carry the source-chain
  identity without rereading unused original chain bytes. An explicit refresh path is required. Whole-bundle
  retrieval is a second explicit boundary: construct a no-network plan that enumerates every resource
  required by the discovered evidence tier, surface each resource's terms classification, and require a
  separate acknowledgement of that transfer plan before acquiring any item. A returned cached bundle is
  complete for its evidence tier; if a later item fails, no partial bundle object is returned, while any
  already-published content-addressed artifacts remain valid cache entries for retry. A separate body-free
  metadata-inspection step can record provider-advertised size and transport headers for the exact plan
  URLs after explicit provider-terms acknowledgement. Live provider verification established the required
  HEAD and byte-range behavior for the current comparative fixture. The acquisition layer now retains and
  resumes partial HTTPS transfers only when an exact provider checksum plus the required strong
  validator/size/range metadata are available; otherwise it uses the fresh streaming path.
- **Call these "evidence-availability tiers," not "confidence tiers."** How much evidence exists
  and what factual/interpretive conclusions it supports are separate questions (see invariant 3, §3).
- **Target-sequence role/context must come from a defined metadata source.** Prefer authoritative
  per-sequence assembly metadata when available. For assemblies represented by NCBI, the genome
  sequence report exposes sequence role, assembly unit, chromosome name, GenBank/RefSeq accessions,
  and UCSC-style sequence name, providing an explicit bridge from assembly metadata to UCSC names.
  Provider-specific naming patterns may be used only as an explicitly labeled fallback when such
  metadata is unavailable. Sequence role/context is descriptive evidence and must not silently
  become an error/quality rule or a claim that ambiguity was biologically caused by duplication or
  an alternate sequence.
- **Batch assessment must reuse resource work across loci.** A batch interface must not be a naive
  outer loop that reparses multi-gigabyte comparative resources once per locus. It should either
  evaluate many intervals during shared resource traversal or use an indexed/preprocessed local
  representation while preserving the same single-locus mapping/evidence semantics and per-locus
  provenance.
- **Portable case packets are constrained by redistribution terms.** A reproducible case manifest
  may record the schema-versioned result, exact resource SHA-256 identities, source URLs,
  retrieval/checksum/terms metadata, and provenance graph. An archive may embed byte-identical
  cached resource files only when their applicable redistribution terms permit it; local cache
  possession alone does not authorize rebundling provider data.

## 8. Licensing constraints

- **Do not treat every UCSC file in chain format as having the same license.** UCSC's general
  Genome Browser licensing page identifies the dedicated liftOver chain files as the data-file
  exception to its otherwise broadly reusable Genome Browser data. The README in a UCSC
  `liftOver/` directory states that its `*.over.chain.gz` files are available free for
  non-commercial use by independent researchers and nonprofit organizations; other use requires
  a commercial license. Downloading or using those files constitutes acceptance of UCSC's EULA,
  and redistribution must include the applicable UCSC README/license terms.
- Comparative resources under `source/vsTarget/` are a distinct publication class and must follow
  the terms published for that resource directory rather than automatically inheriting the
  `liftOver/*.over.chain.gz` restriction merely because some files use chain format. For the
  established `canFam3` ↔ `canFam4` mechanical fixture, UCSC's `canFam3/vsCanFam4/` README explicitly
  states that all files in that directory are freely available for public use. Future comparative
  assembly pairs must still retain and respect their own provider README/terms rather than
  generalizing from this one directory.
- UCSC's liftOver **program/source** is separately licensed from most kent command-line utilities;
  UCSC lists the `src/hg/liftOver` source directory under its non-commercial UC license and offers
  commercial licensing for liftOver. liftAssess therefore does not depend on or redistribute the
  UCSC liftOver implementation for its core logic.
- **Provider network access is explicitly gated.** The resource resolver may identify either a
  comparative URL or a restricted `liftOver/*.over.chain.gz` URL. Provider-contacting inspection and
  acquisition therefore refuse network access until the caller explicitly acknowledges
  review of UCSC's general and relevant directory-specific terms. Restricted liftOver chains are
  surfaced distinctly because UCSC states that downloading or using them indicates EULA acceptance and
  limits free use to the described non-commercial/nonprofit cases unless an applicable commercial
  license exists. No-network discovery/planning itself does not require terms acknowledgement, and HEAD
  inspection does not substitute for the separate transfer-plan acknowledgement required to download
  resource bodies.
- Do not bundle or mirror UCSC chain/net resources in the liftAssess source tree, package,
  release artifacts, or fixtures. Keep downloaded provider data in a separate local cache,
  preserve source URL, applicable terms/README provenance, retrieval metadata, and content
  checksum, and keep liftAssess's own GPL-licensed code independent of those external files.
  The acquisition layer now supports individual resources plus explicit complete-bundle transfer plans
  in a caller-selected external cache. Planning itself performs no network access, and bundle execution
  requires a separate transfer-plan acknowledgement before any resource acquisition starts. Terms-gated
  HEAD inspection can expose exact remote size before transfer. When an exact provider checksum, strong ETag, exact identity-encoded
  size, and byte-range support are available, interrupted HTTPS transfers can resume safely; otherwise the
  implementation falls back to a fresh streaming transfer rather than inventing a weaker resume contract.
- Treat UCSC as one external, terms-bound evidence provider. Always accept user-supplied resources;
  their licensing remains the user's/provider's responsibility and must not be represented as
  covered by liftAssess's GPL license.

Primary UCSC terms/checksum behavior checked through 2026-08-13:
- Genome Browser licensing: https://genome.ucsc.edu/license/
- canFam3 liftOver README/terms: https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/
- canFam3/vsCanFam4 comparative README: https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/
- canFam3/vsCanFam4 provider MD5 metadata: https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt
- canFam4/vsCanFam3 reciprocal-best MD5 metadata: https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/reciprocalBest/md5sum.txt
- kent source license: https://github.com/ucscGenomeBrowser/kent/blob/master/LICENSE

## 9. Output format

The redesigned output uses one structured result profile and multiple renderers. Human and machine
output are different representations of the same scientific/composite-analysis facts, not separate
assessment paths.

### Default human summary — progressive disclosure

Default terminal output is **facts-first**:

1. preflight/input failure, if any;
2. deterministic factual mapping headline;
3. the few measured facts needed to understand it;
4. evidence tier plus material consumed/context resources;
5. bounded deterministic interpretation;
6. relevant scope/identity boundaries; and
7. a pointer to complete details/machine output.

Uncomplicated results remain compact. The renderer expands when materially unusual states are
present, initially including partial coverage, fragmented/discontinuous geometry, multiple
projections, alternate/unplaced targets, reverse disagreement, point/101-bp disagreement,
filtered/all-chain disagreement, comparative conflict/non-separation, batch collision/overlap, or
typed difficult-region context.

The renderer does **not** print six invariant `NOT TESTED` lines on every clean result. The six
common-use lenses remain represented in the profile/details so that absent domains cannot be
mistaken for assessed evidence.

For large intervals, source coverage and fragmentation lead the story. Use a mechanically defined
maximum candidate source coverage unless a separate rank concept has been explicitly defined.

### Detailed human dossier (`--details`)

Detailed output exposes the complete result profile and supporting evidence, including:

- exact mapped segments and target bounding spans;
- source coverage and uncovered source intervals;
- source/target chain gaps;
- orientation and target-role context;
- every relevant evidence observation;
- filtered/all-chain and comparative relationships when available;
- reciprocal-best membership and actual reverse-result context as distinct concepts;
- point/neighborhood comparisons when run;
- typed external-context observations and exact provenance;
- batch relationships when assessing multiple records;
- resource retrieval/checksum/terms context and actual engine consumption; and
- the provenance dependency graph.

Candidate encounter order is preserved for reproducibility but is not presented as rank.

### Machine-readable output (`--json`)

The redesign deliberately introduces a **new schema version** rather than preserving the alpha-v1
verdict schema. Schema v1 remains historical v0.1.0a1 behavior and may be broken by this pre-release
redesign. Do not retain `verdict`, verdict-derived `decision_reason`, or preferred-candidate semantics
merely for compatibility with the obsolete aggregate model.

The new schema is built around the derived result profile plus exact candidates/evidence/resources/
provenance. Exact field names are implementation work, but the structure must preserve:

- schema version and report type;
- source/target assembly identifiers and exact source interval;
- orthogonal result-profile dimensions from §4;
- exact candidate geometry and evidence;
- resource consumption, retrieval/checksum/terms metadata, and SHA-256 identities;
- typed external-context observations;
- composite-analysis results (reverse, neighborhood, comparative, batch) when run;
- provenance/dependency edges; and
- explicit scope boundaries for untested identity/workflow domains.

Every genomic interval is emitted in canonical 0-based, half-open coordinates with an explicit
coordinate-system field. Target bounding spans must remain named/described as bounding spans so
fragmented mappings cannot be mistaken for continuous alignment. JSON field/array encounter order
is reproducibility context unless a field explicitly defines ordering semantics.

The compatibility break must be documented in release notes when the redesigned schema ships. A
migration layer is not current scope unless real users establish a need.

### Comparative wording

Comparative categorical states are summaries of structured relationships, not votes. If evidence is
mixed/conflicting, human output must identify the material conflict instead of merely printing a
label such as `COMPARATIVE EVIDENCE IS MIXED`. `ali` and `qDup` remain descriptive observations
until an explicit, validated interpretation rule says otherwise.

### Process status and automation

Invalid input, usage errors, and operational failures may return nonzero process status. A valid
completed scientific result—including no chain projection, multiple projections, or unresolved
comparative evidence—is not itself a process failure. Automation should branch on structured result
fields, not prose or exit status for scientific interpretation.

### Optional navigation/export

- Genome Browser/locus links may be emitted as navigation aids when the required assembly and
  coordinates are known. A link is not new scientific evidence or an identity check.
- BED12/custom-track export may represent one candidate whose mapped blocks are legally expressible
  on one target sequence. It must not collapse multiple target sequences/candidates into one BED12
  object and does not replace source-coverage reporting.
- A compact profile-vector string is deferred; structured JSON is the machine interface.

### CLI/resource behavior retained from the alpha implementation

The common single-locus CLI remains cache-first. A complete verified local bundle can be assessed
without provider access; `--offline` makes zero network access explicit, `--refresh` deliberately
checks/reacquires current provider resources, `--cache-dir` overrides the user cache, and `--quiet`
suppresses nonessential progress. Interactive transfer/cache-verification/assessment progress is
measured from actual bytes and is terminal-only; it must not fabricate biological/algorithmic
completion percentages or ETAs.

For automatic UCSC runs, consumed UCSC files remain conservatively grouped under their shared
upstream lineage for dependency purposes while each exact file is identified independently by
SHA-256. This prevents chain/net/reciprocal-best observations from reading as independent
confirmation without claiming more provider-process history than the downloaded bytes establish.

**Legacy implementation note:** v0.1.0a1 emitted schema v1 with
`WELL_SUPPORTED`/`CONTESTED`/`INDETERMINATE`, `decision_reason`, and the legacy default summary.
Those facts remain valid implementation history, not current target policy.

## 10. Validation strategy — fixtures, real-world corpus, and held-out UX gate

- **Mechanical evidence fixture** — proves correct extraction (score, coverage, `ali`, `qDup`,
  net type, hierarchy). Use `canFam3` ↔ `canFam4`: confirmed to have a full, published
  chain/net/reciprocalBest comparison at UCSC (documented lastz params: `M=254`; axtChain
  `minScore=3000`; processed through chainNet/netSyntenic/netClass). As measured from UCSC's
  live listings on 2026-08-12, the canFam3-referenced chain/net/syn-net files are under
  `canFam3/vsCanFam4/`, while the directional `canFam3.canFam4.rbest.{chain,net}.gz` files are
  published under the sibling `canFam4/vsCanFam3/reciprocalBest/` directory alongside the reverse
  direction. That hosting asymmetry is publication layout, not a coordinate-semantics change.
  Tasha (canFam3) and Mischka (canFam4, German Shepherd) are different individuals, so no
  ground-truth claim is needed or made — this fixture is purely about extraction correctness.
  - **Real-file smoke check completed 2026-08-12 (not the full comparative fixture).** The
    restricted `canFam3ToCanFam4.over.chain.gz` liftOver resource matched UCSC's published MD5
    `15123263dbe4f2c1eb670a98c9b0acf2`; the exact downloaded bytes had SHA-256
    `c79c9e7c2a3d546f7a9d7efe27cc8815da611d79adb0da4e4ff1556810f28f48`. The implemented
    file → parser → engine path mapped the 0-based half-open source interval
    `chr1:12514-12534` to one candidate at `chr1:660-680` with full 20/20 source-base coverage
    and no chain gaps. No assessment verdict was computed. This establishes real-file mechanical
    plumbing under the sparse LIFTOVER-ONLY path; it does not validate comparative evidence or
    biological support.
  - **Full comparative mechanical fixture completed 2026-08-16.** The selected canonical
    0-based half-open source interval is `chrUn_JH373233:1845735-1845835`. Running the exact
    externally cached five-resource bundle through `build_ucsc_candidates_from_cached_bundle()`
    produced 170 candidates across 114 target sequences. Chain 573 provides the primary
    evidence-rich mechanical case: reverse orientation to `chr35:925644-925938`, two aligned
    segments, full 100/100 source coverage, one target-side chain gap, chain score 16,617,372,
    `ali=3603`, `qDup=4098`, `nonSyn` net classification at depth 7, and `FULL` reciprocal-best
    membership. Chains 5170 and 2692 provide contrasting partial-coverage and reciprocal-best-
    absent cases. The production path SHA-256-verifies the consumed resource bytes and preserves
    one caller-declared shared upstream alignment ancestor across chain/net/reciprocal-best file
    provenance. The original completed fixture run verified extraction and provenance wiring only.
    The reproducible verifier also contains a legacy-alpha cross-check that independently derives
    the then-expected `CONTESTED` verdict from extracted candidate evidence without calling
    `assess_candidates()`; the post-hardening real-data rerun on 2026-08-17 derived that legacy
    result from 138 material candidates and matched the production alpha assessor. This remains
    historical regression coverage for v0.1.0a1 semantics, not a target-model requirement. When
    Milestone 17 replaces the aggregate verdict model, the same fixture should verify the new
    factual-profile/comparative-relationship semantics instead. No biological ground-truth claim is
    made. The verifier lives at `scripts/verify_canFam3_canFam4_mechanical_fixture.py`; the UCSC bulk
    files remain external.
- **Historical-resolution fixture pedigree (not yet a concrete fixture)** — identifies the right
  assembly pair for proving the report behaves sensibly against a known resolution; the specific
  truth-bearing locus within that pair has not yet been identified. `canFam3.1` → `canFam6`
  (Dog10K_Boxer_Tasha_1.0) is the right pedigree: **same individual** (Tasha) as the original
  reference and CanFam3.1, long-read assembly, independently documented to close >23,000
  CanFam3.1 gaps and improve contiguity >100-fold. Confirmed that UCSC only publishes the
  lightweight `canFam6ToCanFam3.over.chain.gz` for this pair — no full `vsCanFam6` comparative
  directory under canFam3 exists. Whatever fixture is eventually built from this pedigree will
  therefore also exercise graceful degradation under sparse evidence, which is realistic: most
  real target users (non-model-organism, thin-resourced pairs) will be in this situation more
  often than the fully-resourced case. See §13 for the open item to turn this pedigree into an
  actual fixture.
  - This fixture is a **sanity check**, not a calibration set. It must not influence any
    threshold or scoring logic (v1 has none to calibrate).
  - **Pedigree caution, confirmed and generalizable:** do not assume sequential assembly version
    numbers form one improving single-individual lineage. canFam4 = Mischka (German Shepherd),
    canFam5 = Zoey (Great Dane), canFam6 = Tasha (Boxer, same as canFam1/2/3). Three different
    individuals across four version numbers. Verify same-individual provenance explicitly before
    treating any later assembly's placement as ground truth for any species.

### Post-alpha 50-case real-world design corpus

The completed 50-case program is design/UX evidence, not a general accuracy benchmark. It was
intentionally enriched for historically troublesome liftover/support cases and controls. Its main
load-bearing conclusions for the current design are:

- a clean single-chain legacy result can omit material context such as non-reciprocity,
  duplication/paralogy, neighborhood fragmentation, or cross-record collision;
- coverage/fragmentation explains large-region threshold behavior more clearly than raw candidate
  count;
- point widening has a measured clean background in six matched 101-bp controls and can also expose
  local structure missed by a 1-bp query;
- filtered chains can conceal additional all-chain placements;
- comparative net/reciprocal-best observations can sometimes materially distinguish placements,
  but their shared provenance must remain explicit;
- batch collisions are a separate evidence scale; and
- input/preflight failures must not masquerade as biological no-projection outcomes.

The `38 YES / 10 MOSTLY / 2 NO` language replay is same-corpus design acceptance, not independent
validation. Likewise, counts of legacy `WELL_SUPPORTED` misses in the enriched corpus are not
population prevalence estimates.

### Held-out language/usability gate

Before the redesigned result language is described as validated or release-ready, test the actual
implemented output on a small set of real cases that were not used to derive the language. Include
both uncomplicated controls and difficult/ambiguous cases. Outside-user/domain feedback should
also exercise representative problem cases, including the automatic 101-bp point context and typed
context observations. This gate evaluates whether the result is understandable and appropriately
bounded; it is not a requirement to calibrate a numeric confidence score.

## 11. Three-gate status (as of this baseline)

| Gate | Question | Status |
|---|---|---|
| 1 | Is this needed? | **Yes.** Two independent peer-reviewed benchmarks + a dedicated tool (Liftoff) exist because this class of problem is real and unresolved by conversion accuracy alone. |
| 2 | Are people asking questions it answers? | **Yes.** Recurring across years, species (human, dog), and use cases (variant, annotation, epigenomic interval), in both Q&A forums and tool support threads. |
| 3 | Can it be extremely helpful and easy to use? | **Plausibly yes, with a measured redesign need.** The alpha works end to end and the 50-case corpus identified concrete interpretation/UX failures that the current redesign addresses. Held-out cases and outside-user/domain feedback remain required before calling the redesigned language validated. |

## 12. Sequencing

The first public alpha proved the end-to-end UCSC workflow, but the 50-case program showed that the
legacy aggregate-verdict interface is too easy to overread and omits important explanatory context.
The post-alpha redesign therefore proceeds in this order:

1. **Factual result profile + renderer.** Replace the target aggregate-verdict model with the
   orthogonal profile in §4, define the new schema version, add deterministic headlines,
   coverage/fragmentation summaries, evidence boundaries, and progressive disclosure using facts
   already computed. Do not begin by rewriting candidate generation.
2. **Start indexing/shared traversal immediately after the first renderer slice.** Prototype and
   benchmark a region-addressable/reusable resource path. In parallel, add authoritative assembly
   metadata/preflight and begin the typed difficult-region pilot (initially UCSC segmental
   duplications where source/terms/assembly coverage are verified).
3. **Actual reverse mapping context.** Keep it distinct from reciprocal-best membership.
4. **Point neighborhood / multi-scale context.** Add automatic centered 101-bp context for 1-bp
   queries plus explicit larger-window controls.
5. **Comparative expansion/asymmetry.** Pair ordinary filtered liftOver with all-chain inventory and
   expose categorical comparative relationships. Mixed/conflicting output must explain the actual
   evidence relationship rather than stop at an opaque label.
6. **Batch relationships.** Add BED/simple interval-table input, shared-traversal assessment, and a
   separate batch layer for exact collision/overlap relationships.
7. **Held-out language/usability gate.** Exercise the implemented redesign on unseen real cases and
   outside users/domain reviewers before describing it as validated/release-ready.

Reverse, neighborhood, filtered/all-chain expansion, and batch must not ship by multiplying the
current exhaustive whole-resource scan. The exact indexing/preprocessing implementation remains a
measured architecture decision rather than a design-policy commitment.

## 13. Open and deferred design items

The main product-policy questions from the post-alpha review are resolved. Remaining items are
implementation/evidence questions or deliberately deferred domains:

- Extend reusable/indexed access to additional resource families only where profiling demonstrates a
  material need. Chain-index preparation is now an explicit cache-only user action via
  `prepare-liftassess-index`; normal assessment never incurs an implicit first-query build and retains
  full traversal when no usable derived index is present.
- Finalize exact result-profile field/API names and the new JSON schema layout; the dedicated derived
  profile boundary and deliberate alpha compatibility break are already decided.
- Extend the accepted initial comparative classifier only with explicit deterministic semantics and
  tests. Do not introduce hidden weighting of `ali`, `qDup`, chain score, or related observations.
- Finalize BED/simple-table CLI syntax, batch file schema, browser-link UX, and optional BED12/export
  command surface within the constraints in §§4, 6, and 9.
- Pilot UCSC segmental-duplication context against the motivating real cases. Evaluate GIAB
  stratifications and relevant `excluderanges` categories separately before promoting them into
  implemented evidence families.
- Identify and document at least one concrete CanFam3.1 locus whose later placement in canFam6 is
  independently established, turning the historical-resolution pedigree in §10 into an actual
  truth-bearing fixture.
- Decide the source for optional flanking-gene synteny context and its fallback behavior when no
  ortholog table exists.
- Define reproducible case manifests and any later byte-containing portable packets under the
  provenance and redistribution-term constraints in §7 and §8.
- Add a genuinely independent mapping-evidence source only after verifying independent upstream
  mapping/alignment lineage and defining cross-engine local-hypothesis equivalence.
- Candidate rank remains deferred. Do not substitute raw chain score, UCSC output order, or encounter
  order for a defensible locus-scoped rank definition.
- Variant identity, gene/transcript identity, VCF normalization/allele semantics, assembly-history
  interpretation, and downstream workflow diagnosis remain separate optional evidence domains rather
  than implicit consequences of coordinate mapping.
- Plugin registries, automatic many-engine support, default fresh alignments, machine-learning
  confidence models, composite numeric scores, large bundled species databases, and hosted service
  infrastructure remain deliberately deferred until real requirements justify them.

