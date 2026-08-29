# Milestone 23 held-out result-language validation

**Status:** PRE-REGISTERED — cases selected before liftAssess execution
**Selection date:** 2026-08-29
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

After the internal case review, send the **unchanged outputs** for H01, H03, H04, and H05 to at least one
outside user or domain-informed reviewer who did not help derive the result language. Include H02 as well
if its preflight behavior or wording is surprising.

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
