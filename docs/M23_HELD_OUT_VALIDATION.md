# Milestone 23 held-out result-language validation

**Status:** INTERNAL CASE EXECUTION COMPLETE — H04 influenced presentation; outside review pending
**Selection date:** 2026-08-29
**Execution date:** 2026-08-29
**Purpose:** Milestone 23 held-out real-case and outside-user/domain gate

## Evidence boundary

This gate evaluates whether the implemented liftAssess result language is useful and scientifically
bounded on real cases that were not used to derive that language. It is not an accuracy benchmark,
prevalence study, sensitivity/specificity study, or calibration exercise.

The 50-case Corpus A/B program and the A03/A04 UCSC segmental-duplication pilot are excluded from this
held-out set. To reduce leakage further, the selected support/forum sources below are not sources used
in Corpus A/B.

The cases and acceptance criteria in this file were frozen before running liftAssess on any selected
coordinate. Do not replace a case because its output is inconvenient, and do not add follow-up queries
merely to move a result toward the historical expectation.

Historical reports are comparison evidence, not automatic biological ground truth. Where a source
mentions an rsID/database placement, that identifier-aware placement remains a separate evidence domain
from coordinate projection.

## Pre-registered cases

### H01 — rs200923174 coordinate-convention control

- **Role:** uncomplicated human point/control; coordinate-convention and scope-language check.
- **Source:** Biostars, “Tool: Converting Genome Coordinates From One Genome Version To Another”.
- **Source URL:** https://www.biostars.org/p/65558/
- **Query:** hg19 → hg38, `chr22:16287557-16287557` (1-based inclusive).
- **Historical comparison:** the discussion identifies a 0-based/1-based input error and reports that
  the correctly expressed 1-based hg19 point maps to hg38 `chr22:15690406`.
- **Why selected before execution:** supplies a clean point case from a source not used in Corpus A/B,
  with an explicit practical confusion that good coordinate wording should prevent.
- **What is being tested:** factual point headline, automatic 101-bp context, source/target conventions,
  reverse/context scope, and whether the output avoids treating rsID identity as established by chain
  geometry.

### H02 — mixed-build / out-of-bounds source point

- **Role:** invalid-input/preflight control.
- **Source:** Reddit r/bioinformatics, “How to perform liftover from 38 to 37 for GWAS summary
  statistics?”.
- **Source URL:** https://www.reddit.com/r/bioinformatics/comments/17y8f8b/
- **Query:** hg38 → hg19, `chr2:242193706-242193706` (1-based inclusive).
- **Historical comparison:** Hail rejected this as outside GRCh38 chr2; the user later reported that
  sampled members of the same rejected set appeared to be hg19 loci, raising a source-build mix-up.
- **Why selected before execution:** directly tests the redesigned requirement that invalid source
  coordinates stop before scientific mapping interpretation rather than becoming biological-looking
  “no projection” results.
- **What is being tested:** authoritative source bounds, user-facing preflight wording, nonzero failure
  behavior, and absence of a scientific mapping headline for invalid input.

### H03 — surprising telomeric interchromosomal point

- **Role:** difficult human point; local-context and typed-context exercise.
- **Source:** same Reddit discussion as H02.
- **Source URL:** https://www.reddit.com/r/bioinformatics/comments/17y8f8b/
- **Query:** hg38 → hg19, `chr10:10709-10709` (1-based inclusive).
- **Historical comparison:** the user reports that both Hail and UCSC liftOver map the GRCh38 point to
  GRCh37 chr18 near position 10905. Community replies suggest telomeric/repetitive complexity, but that
  mechanism is not treated as established evidence here.
- **Why selected before execution:** a real user found the chromosome change surprising; it tests whether
  liftAssess can describe a non-obvious projection without declaring it erroneous or biologically
  correct.
- **What is being tested:** automatic 101-bp context, actual reverse mapping, interchromosomal and
  orientation wording, typed UCSC segmental-duplication context when available, target-role context, and
  explicit separation between measured structure and unverified mechanism.

### H04 — CanFam3.1 LHX3 pituitary-dwarfism variant

- **Role:** genuinely new real canine COMPARATIVE case.
- **Primary database source:** OMIA variant `omia.variant:608`, LHX3.
- **Source URL:** https://omia.org/variant/omia.variant%3A608/
- **Publication:** Voorbij et al., 2011, PLoS ONE, PMID 22132174.
- **Query:** canFam3 → canFam4, `chr9:49252491-49252493` (1-based inclusive).
- **Source basis:** OMIA records CanFam3.1 `NC_006591.3:g.49252491_49252493dup` for the LHX3-related
  pituitary-dwarfism record. The publication establishes the LHX3 disease/variant context; OMIA notes
  that the genomic CanFam3.1 coordinates were supplied subsequently.
- **Why selected before execution:** it is outside the B12–B14 DoGA/OMIA loci used to shape comparative
  language and therefore exercises COMPARATIVE reporting on a new real canine locus.
- **What is being tested:** filtered/all-chain relationship, categorical net/reciprocal-best explanation,
  provenance/dependence language, target-role/context availability, and whether a comparative conclusion
  stays below biological variant identity/correctness.
- **No target coordinate is pre-declared as correct.** The case tests explanatory evidence, not agreement
  with an independently adjudicated CanFam4 locus.

### H05 — rs138257042 asymmetric chr22/chr14 mapping

- **Role:** difficult human point; variant-identity boundary and typed-context stress case.
- **Source:** same Biostars discussion as H01.
- **Source URL:** https://www.biostars.org/p/65558/
- **Query:** hg38 → hg19, `chr22:15528888-15528888` (1-based inclusive).
- **Historical comparison:** the discussion reports UCSC liftOver returning hg19 chr14 near 19378323,
  while an rsID-aware GRCh37 placement is reported on chr22 near 16449075. The database placement is
  treated as identifier-aware comparison evidence, not as proof that one coordinate projection is the
  biologically correct locus.
- **Why selected before execution:** this is an unseen locus/source relative to Corpus A/B and directly
  tests whether the redesigned output prevents a coordinate projection from being mistaken for variant
  identity. It is related to the chr22/chr14 duplication failure class represented by motivating case
  A04, so it must not be presented as independent evidence that the mechanism generalizes.
- **What is being tested:** automatic 101-bp context, reverse mapping, typed duplication context, scope
  boundaries, and whether the output explains the evidence conflict without resolving biological
  identity.

## Execution rule

Run the five baseline queries exactly as registered above. The automatic capabilities shipped by
liftAssess may run normally. Additional manual follow-ups are allowed only when an observed result leaves
two or more concrete hypotheses that a narrowly targeted query can distinguish. Record the reason before
running such a follow-up.

Do not substitute a case after execution. A provider/resource outage may leave an optional dimension
`UNAVAILABLE`; that is itself part of the observed gate result and should be recorded rather than worked
around silently.

## Case-level review questions

For each case, record:

1. What factual event does the headline say occurred?
2. Do the detailed geometry and structured fields support that headline exactly?
3. Does progressive disclosure expose the material unusual facts without making the clean cases noisy?
4. Are point-context, reverse, comparative, target-role, batch, and typed-context scope states truthful?
5. Does any sentence imply variant identity, gene identity, uniqueness, causal mechanism, or biological
   correctness beyond the evidence consumed?
6. Would the output have answered the practical confusion in the historical source, or at least made the
   remaining evidence gap explicit?
7. Did optional enrichment failure, if any, preserve the already-valid primary coordinate assessment?

## Internal execution record

All five pre-registered baseline cases were executed without substitution on 2026-08-29. No blocking
scientific-correctness or evidence-boundary failure remains from the internal pass. H04 did expose a
release-worthy presentation gap: the shared interpretation string reported multiple projections but did
not surface the already-computed `FAVORS_ONE_PLACEMENT` comparative relationship. That finding directly
influenced implementation, so this five-case set is **not** described as untouched held-out validation.
The interpretation was corrected without changing candidate generation or comparative classification,
and H04 was rerun unchanged after the full native test/lint/type-check gate passed.

### H01 result — coordinate-convention control

- **Observed:** one complete hg19→hg38 projection to `chr22:15690406`, matching the historical coordinate.
- **Local/reverse context:** the automatic 101-bp window mapped 101/101 bases contiguously through the
  same chain; reverse mapping returned only to the original hg19 source point.
- **Typed context:** source and target segmental-duplication overlaps were present, but remained
  descriptive and did not alter the coordinate interpretation.
- **Adjudication:** no blocker. The output distinguishes clean coordinate geometry from unassessed rsID
  identity even when typed context is present.

### H02 result — invalid source coordinate

- **Observed:** authoritative hg38 `chr2` length was 242,193,529 while the requested point was
  `chr2:242193706`.
- **Behavior:** source preflight stopped the run before mapping, printed that mapping was not attempted,
  and returned process exit status 1.
- **Adjudication:** no blocker. Invalid input did not become a biological-looking no-projection result.

### H03 result — telomeric interchromosomal projection

- **Observed:** hg38 `chr10:10709` projected completely to hg19 `chr18:10905`; the 101-bp context also
  mapped 101/101 bases contiguously and agreed with the point.
- **Reverse/context evidence:** actual reverse mapping was `ELSEWHERE_ONLY`. The source point overlapped
  three hg38 `genomicSuperDups` rows, including a chr10↔chr18 row with `fracMatch=0.979191`; target
  overlap was also assessed. These observations were reported as descriptive context rather than a
  causal explanation.
- **Target role:** `UNAVAILABLE` because the hg19 assembly description does not provide the exact
  versioned NCBI assembly binding required by the target-role model. The primary mapping continued.
- **Adjudication:** no blocker. The output reproduces the surprising chromosome change while separating
  clean local geometry, non-reciprocity, duplication context, and unverified mechanism.

### H04 result — new canine COMPARATIVE case

- **Observed:** the 3-bp CanFam3 locus produced seven complete all-chain projections. One canFam4 chr9
  placement was retained by the ordinary filtered chain, represented by a depth-1 top-net fill, and had
  full reciprocal-best membership; none of the six other complete placements had the same categorical
  support pattern. The comparative relationship was therefore `FAVORS_ONE_PLACEMENT`.
- **Evidence boundary:** filtered-chain, net, and reciprocal-best observations remained grouped as
  provenance-dependent UCSC-derived evidence rather than independent votes. Variant and gene identity
  remained unassessed, and the result did not establish a biological locus.
- **Optional dimensions:** reverse mapping was `UNAVAILABLE` because no cached reverse-direction chain
  with matching COMPARATIVE publication class was available and UCSC was not contacted. Typed
  segmental-duplication context was also `UNAVAILABLE`; neither absence changed the primary comparative
  assessment. The exact reverse-unavailability reason is currently present in run-status text while the
  durable dossier/JSON carries only the `UNAVAILABLE` state; this is a non-blocking release UX follow-up.
- **Implementation influence:** the first run's top-level interpretation said only that multiple chain
  projections existed even though the detailed comparative section already reported
  `FAVORS_ONE_PLACEMENT`. The interpretation was changed to state that available categorical comparative
  evidence favors one placement while preserving the no-biological-locus boundary. Regression coverage
  verifies the result profile, detailed dossier, and schema-v2 JSON.
- **Rerun:** the unchanged H04 query produced the same seven placements and same categorical evidence
  relationships after the presentation fix; the new top-level interpretation surfaced the comparative
  conclusion.
- **Adjudication:** the original presentation gap is resolved for the current candidate. Because H04
  caused the change, this case is implementation-influencing evidence rather than untouched validation.

### H05 result — rs138257042 asymmetric chr22/chr14 mapping

- **Observed:** hg38 `chr22:15528888` projected completely to hg19 `chr14:19378323`; the 101-bp context
  also mapped 101/101 bases contiguously and agreed with the point.
- **Reverse/context evidence:** reverse mapping was `ELSEWHERE_ONLY`. The source point overlapped six hg38
  segmental-duplication rows, including one with `fracMatch=0.996176`; the hg19 target point overlapped
  four rows, including one with `fracMatch=0.996027`.
- **Identity boundary:** target role was `UNAVAILABLE` under the strict hg19 assembly-binding rule, and
  named-variant identity remained explicitly unassessed.
- **Adjudication:** no blocker. The output exposes coordinate projection, local agreement, non-reciprocity,
  and duplication context without resolving the historical rsID-aware placement conflict.

### Internal gate disposition

The internal five-case pass therefore has no unresolved blocker. It did produce one implementation-
influencing H04 presentation correction and one non-blocking reverse-unavailability explanation follow-up.
Milestone 23 remains open until the outside-user/domain packet is reviewed against the current candidate.

### Pending UX observations from independent AI review

A separate AI review of the H01/H03/H04/H05 dossiers produced the following potential release UX improvements. These are **advisory observations, not adopted requirements or design decisions**. Keep the M23 outside-review candidate frozen and adjudicate the outside-user/domain feedback before deciding which, if any, to implement:

1. surface the actual reverse-mapping relationship more prominently near the top of detailed output, especially `ELSEWHERE_ONLY`, without folding it into the factual headline or turning it into a confidence verdict;
2. report same-sequence versus interchromosomal projection explicitly as a neutral geometric relationship;
3. preserve the structured reason for reverse-mapping `UNAVAILABLE` in the durable result/JSON rather than only in run-status text;
4. reconsider the human-facing word `categorical` in comparative interpretation because ordinary-language readers may hear it as “absolute” rather than “discrete relationship class”;
5. make scope-versus-result grammar consistent so availability states such as `ASSESSED`/`UNAVAILABLE` are not mixed ad hoc with conclusions such as `FAVORS_ONE_PLACEMENT`;
6. consider rendering point-context agreement as the point and tested 101-bp context mapping through the same forward chain, avoiding language that could sound like independent corroboration;
7. make the distinction between actual reverse mapping and precomputed UCSC reciprocal-best membership visually and terminologically unmistakable; and
8. consider a bounded reverse-orientation note for downstream strand-dependent sequence/allele use that explicitly says liftAssess did not transform or validate that downstream data.

The same review also proposed stronger interpretations that are **not** carried forward: chain-score-based confidence/ranking, claims that UCSC-derived observations are independent corroboration, causal duplication/mechanism conclusions from contextual overlap, or a replacement aggregate confidence/verdict layer. Those suggestions conflict with the current scientific model and previously reviewed invariants.

## Blocking failure criteria

Milestone 23 does not pass until any observed blocker is resolved and the affected held-out case is rerun.
A blocker includes:

- invalid input being rendered as a biological-looking mapping result;
- a factual headline contradicted by the detailed geometry/structured result;
- a material point/neighborhood disagreement being hidden or described incorrectly;
- typed contextual overlap being presented as a penalty, proof of error, or causal mechanism;
- comparative relationships being materially misstated, collapsed into hidden weighting, or presented
  as independent votes despite shared provenance;
- a scope state implying that evidence was assessed when it was not;
- an optional context/enrichment failure discarding an otherwise valid primary assessment;
- coordinate-convention, source/target, or strand wording that would reasonably lead a user to act on the
  wrong physical interval; or
- outside-user/domain feedback identifying a material scientific overclaim or a practical explanation
  failure that would make the result unsafe or substantially misleading.

Minor wording preferences do not fail the gate by themselves. Any correction made because of the
held-out set must be documented as such; the held-out set must then be described as having influenced the
implementation rather than as untouched validation evidence.

## Outside-user/domain packet

After the internal case review, send the outputs for H01, H03, H04, and H05 that represent the current
release candidate to at least one outside user or domain-informed reviewer who did not help derive the
result language. For H04, use the post-correction rerun and disclose that this case exposed and influenced
the comparative-interpretation wording. Include H02 as well if its preflight behavior or wording is
surprising. Do not present the packet as untouched held-out validation.

Ask the reviewer, without first supplying our preferred interpretation:

1. What do you think physically happened to the queried interval?
2. What evidence in the output makes you think that?
3. What, if anything, would you do next before using the projected coordinate?
4. Does any wording sound like a stronger claim than the evidence supports?
5. Is it clear which questions liftAssess did **not** assess (for example rsID/gene identity or biological
   correctness)?
6. For COMPARATIVE and typed-context cases, is it clear that related UCSC observations are contextual or
   provenance-dependent rather than independent votes?

Record the reviewer’s role/background only at the level needed to interpret the feedback; do not turn the
gate into an endorsement request.

## Gate completion rule

Milestone 23 passes when:

- all five pre-registered cases have been executed and adjudicated against the questions above;
- no unresolved blocking failure remains;
- the outside-user/domain packet has been reviewed;
- any blocking outside feedback has been resolved and the affected cases rerun; and
- the final record explicitly states whether the held-out cases remained untouched or influenced any
  implementation/language changes.
