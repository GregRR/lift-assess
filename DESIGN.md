# LiftOver Ambiguity Assessor — liftAssess - Design Baseline (v1)

Status: v1 design baseline; implementation in progress. This document is a starting point for
implementation, not a constraint against revising it - building will expose mistakes here, and
it should change when it does.

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
- Biostars: human liftOver landing on a paralogous region (FAM72C/SRGAP2D) instead of the
  correct one (FAM72B/SRGAP2C); resolved manually by checking flanking genes and reciprocal
  mapping — exactly the evidence types this tool automates.
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

- **Assessor, not resolver.** Never outputs "correct locus = X." Always evidence + a verdict
  about support, never a claim of biological truth. Documentation should state, near-verbatim:
  *"Well supported does not mean correct."*
- **Three verdicts only:** `WELL_SUPPORTED`, `CONTESTED`, `INDETERMINATE`. Resist adding more, or
  a candidate-level fourth label — nuance belongs in the evidence detail, not the verdict label.
  See §4 for operational definitions.
- **No composite numeric score in v1.** A number like "0.91" claims a precision that hasn't been
  validated against ground truth. Show categorical evidence instead.
- **Evidence carries provenance, not just a result.** Every observation records its upstream
  source (e.g. "UCSC canFam3→canFam6 alignment," "Ensembl Compara release X," "independent
  minimap2 alignment run locally"). This lets the report detect when multiple observations trace
  back to the same underlying alignment and must not be presented as independent confirmation.
  Chain score, net type, qDup, and reciprocal-best membership, when derived from the same UCSC
  alignment, count as **one** line of evidence, not four.
  - Nets are derived from chains (UCSC's own pipeline: chains → chainNet → net). Not independent.
  - UCSC's standard "reciprocal-best" chains are typically produced by **swapping** the original
    alignment, not by an independent second alignment run. Reciprocal agreement using UCSC's own
    paired chains is closer to a self-consistency check than a second opinion. Label it as such.
    Genuine independence requires a freshly computed reverse alignment (e.g. minimap2), done only
    when explicitly requested.
  - Typed digest model for provenance identifiers (do not misapply refget to non-sequence files):
    - reference sequence → GA4GH refget sequence digest (SQ....)
    - sequence collection → GA4GH SeqCol identity
    - chain / net / GTF / any other file → plain `sha256:<bytes>`
  - Independence is a property established from provenance on a case-by-case basis, never a
    category assigned to an evidence type in advance (see §6, "genomic-context evidence").
- **Interpretation stays one rung below evidence.** E.g. `qDup` + `nonSyn` net type together is
  labeled *"duplication-associated nonsyntenic placement"* (what the evidence shows), not
  *"likely paralog/processed pseudogene"* (one specific cause among several: paralog, segmental
  duplication, transposed duplication, assembly artifact, processed pseudogene). Promote to a
  named interpretation only when additional supporting evidence exists (e.g. an annotated
  intronless copy of a multi-exon gene).
- **Same species is the v1 operational envelope; same individual is a validation criterion, not
  a usage restriction.** The assessor should work across different individuals' assemblies within
  a species — it should simply warn that divergent placements may reflect real structural
  variation, not error. Same-individual assembly pairs are what's needed to construct backtest
  fixtures with known ground truth (see §9). This also matches UCSC's own documented guidance:
  liftOver "was only designed to work between different assemblies of the same organism" and
  should not be used with `-multiple` for cross-species or fragmented/poor-quality assemblies.

## 3. Scientific invariants

A short checklist to test any future implementation decision against:

1. Never count observations that share an upstream alignment or source as independent
   confirmation of each other.
2. Absence of evidence for a competing candidate is not proof that no competing candidate exists.
3. Evidence availability (how much could be checked) is not confidence (how strong the answer
   is) — track them separately. A COMPARATIVE evidence tier can still be `INDETERMINATE`; a
   sparse tier can still be `WELL_SUPPORTED`.
4. Different-individual assemblies may legitimately differ biologically; a disagreement between
   them is not automatically an error.
5. A verdict describes evidentiary support, not biological truth. "Well supported" does not mean
   "correct."
6. Biological interpretation sits one rung below evidence and is only promoted to a named
   mechanism (paralog, pseudogene, etc.) when independently supported, not inferred from a single
   evidence pattern.

## 4. Verdict definitions (operational)

- **`WELL_SUPPORTED`** — the available informative evidence favors one candidate, and no material
  evidence contradicts it.
- **`CONTESTED`** — two or more candidates retain meaningful, non-negligible support, or
  informative evidence sources materially disagree with each other.
- **`INDETERMINATE`** — the available evidence is insufficient, non-discriminating between
  candidates, or too mutually dependent (shared provenance) to distinguish candidates.

These are qualitative and deterministic in intent, not numeric thresholds: the same evidence set
applied against this text should reliably produce the same label. They do not eliminate all
judgment — "material" and "meaningful" still require interpretation when implemented — but they
give every future decision a fixed text to be tested against, rather than leaving the boundary to
be reinvented case by case.

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

## 6. v1 Scope

**In scope — evidence types:**
- Mapping structure: coverage (full/partial locus), chain gaps through the locus, candidate
  rank, target chromosome/scaffold placement.
- Chain/net evidence: chain score, aligned bases (`ali`), duplicated query bases (`qDup`), net
  classification (`top`/`syn`/`nonSyn`/etc.), net hierarchy, reciprocal-best membership.
- Optional genomic-context evidence: flanking-gene orthology/synteny, when available, with its
  provenance and evidence dependence evaluated explicitly rather than assumed independent — an
  orthology call may itself incorporate alignment/synteny evidence from a related pipeline.

**Explicitly out of scope for v1:**
- Freshly computed sequence identity from raw bases.
- A new alignment run (minimap2/lastz) by default.
- A numeric composite confidence score.
- Machine learning.
- Automatic claims of orthology.
- Large species-support databases.
- Hosted infrastructure.

## 7. Architecture

```
Candidate-generation engine
        │
        ▼
NormalizedCandidate[] + provenance
        │
        ▼
Assessor core (evidence extraction, dependency/provenance labeling, verdict)
        │
        ▼
Assessment report (summary + detailed dossier)
```

- The assessor core knows nothing about how candidates were produced. It consumes only the
  normalized candidate representation and provenance.
- **v1 implements exactly one engine**: an internal, minimal chain/net reader against the
  documented UCSC chain/net format. Not pyliftover (point-converter only, no block/gap detail) —
  and not the UCSC liftOver binary (see §8, licensing).
- **Chains generate candidate mappings; nets annotate and evaluate those candidates.** These are
  distinct responsibilities. For UCSC-generated liftOver map chains specifically, the old/source
  assembly is the chain target (`t*`) side and the new/destination assembly is the query (`q*`)
  side; candidate generation follows that documented direction explicitly rather than inferring it
  from assembly names. The LIFTOVER-ONLY evidence tier has no net file and must still be
  able to generate and report candidates from chain data alone. v1 may implement both in one
  reader package, but the two responsibilities stay conceptually — and ideally structurally —
  separate, so net availability is never accidentally required for candidate generation itself.
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
  consistency, not independent evidence — same provenance discipline as everywhere else.
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
  1. Look for a full `source/vsTarget/` comparative directory (chain, net, synNet,
     reciprocalBest) → **COMPARATIVE** evidence tier.
  2. Else look for `source/liftOver/sourceToTarget.over.chain.gz` → **LIFTOVER-ONLY** tier.
  3. Else → unavailable; require user-supplied resources.
  4. Always tell the user which tier is in play, in plain language, before showing results.
  5. Always accept user-supplied chain/net resources directly — UCSC is a convenient default
     provider, not a hard dependency.
- **Call these "evidence-availability tiers," not "confidence tiers."** How much evidence exists
  and how strong the resulting verdict is are orthogonal (see invariant 3, §3).

## 8. Licensing constraints

- UCSC chain files are free to use/download/link/redistribute for **non-commercial use only**
  (per UCSC's own licensing page); redistribution must retain the original README/license terms.
- UCSC's liftOver **program itself** (not just the chain data) is separately listed among the
  restricted-license kent source directories — not covered by the general MIT license that
  covers most kent utilities. Commercial use requires a separate license from UCSC.
- Design response: do not bundle or mirror UCSC chain/net files; fetch live from UCSC into a
  local cache with the user accepting UCSC's terms directly (not laundered through this tool's
  license); do not depend on UCSC's liftOver binary/source for core logic — use the internal
  chain reader instead; keep this tool's own code fully open source and independent of UCSC's
  restricted components; treat UCSC as one external, terms-bound evidence provider.

## 9. Output format

Two layers, so rigor and glanceability aren't in tension:

**Summary (default, ~5 lines):**
```
Assessment: WELL SUPPORTED
Preferred candidate: chr16:...
Alternative: chrUn_...

Why:
  + stronger UCSC alignment placement
  + syntenic net context
  + no duplication-associated placement detected

This does not establish biological correctness.
```

**Detail (via `--details` or JSON):** full provenance, chain IDs, net hierarchy, evidence
dependency notes (what shares a source with what), checksums.

CLI target for the common case:
```
assess-liftover canFam3 canFam4 chr16:12345-12400
```
Everything else (resource discovery, evidence tier, candidate generation) happens automatically.
Expert users can override resources explicitly. Engine selection becomes a real option only once
a second engine exists — v1 has exactly one (§7).

## 10. Validation strategy — two fixtures, different jobs

- **Mechanical evidence fixture** — proves correct extraction (score, coverage, `ali`, `qDup`,
  net type, hierarchy). Use `canFam3` ↔ `canFam4`: confirmed to have a full, published
  chain/net/reciprocalBest comparison at UCSC (documented lastz params: `M=254`; axtChain
  `minScore=3000`; processed through chainNet/netSyntenic/netClass). Tasha (canFam3) and Mischka
  (canFam4, German Shepherd) are different individuals, so no ground-truth claim is needed or
  made — this fixture is purely about extraction correctness.
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

## 11. Three-gate status (as of this baseline)

| Gate | Question | Status |
|---|---|---|
| 1 | Is this needed? | **Yes.** Two independent peer-reviewed benchmarks + a dedicated tool (Liftoff) exist because this class of problem is real and unresolved by conversion accuracy alone. |
| 2 | Are people asking questions it answers? | **Yes.** Recurring across years, species (human, dog), and use cases (variant, annotation, epigenomic interval), in both Q&A forums and tool support threads. |
| 3 | Can it be extremely helpful and easy to use? | **Plausibly yes**, via automatic resource discovery, visible evidence tiers, and the two-layer summary/detail report. Untested until built. |

## 12. Sequencing


**LiftOver ambiguity assessor** — need and demand are now considered validated rather
than merely plausible. Remaining uncertainty is implementation, not justification.

v1 should not be shaped around any hypothetical downstream adopter or integration partner — the
pluggable engine boundary (§7) already preserves that flexibility for later without costing
anything now.

## 13. Open items for whenever this is picked back up

- Implement the internal chain/net reader against the documented format, with candidate
  generation and net annotation kept as separate responsibilities (§7).
- Build the `canFam3`/`vsCanFam4` mechanical fixture end to end.
- Identify and document at least one concrete CanFam3.1 locus whose later placement in canFam6
  is independently established (e.g. traceable via the Dog10K assembly paper's gap-closure or
  SNV-array mapping data), turning the historical-resolution pedigree in §10 into an actual
  fixture.
- Decide the exact `--details` / JSON schema.
- Decide the source for optional flanking-gene synteny context (e.g. Ensembl Compara) and its
  fallback behavior when no ortholog table exists for a species pair.
- Write the resource resolver against UCSC's Golden Path layout, with real existence checks
  rather than constructed-path assumptions; confirm behavior when neither tier is available and
  only user-supplied resources are given.
