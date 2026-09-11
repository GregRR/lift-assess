# Milestone 23 held-out result-language validation

**Status:** COMPLETE — Milestone 23 gate passed
**Selection date:** 2026-08-29
**Initial execution date:** 2026-08-29
**Current-renderer verification:** 2026-09-06
**Outside complex-case review and final rerun:** 2026-09-06
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

All five pre-registered baseline cases were executed without substitution on 2026-08-29. The held-out
set subsequently influenced implementation in two presentation-only ways: outside-user feedback on H01
showed that the original default output was not understandable enough, and H04 exposed the need for a
compact multiple-mapping/comparative explanation. Candidate generation, comparative classification, and
scientific evidence semantics were not changed by those renderer corrections. The five-case set is
therefore **not** described as untouched held-out validation.

The current default renderer was verified across H01-H05 on 2026-09-06 using real reruns where material
output changed plus regression coverage for wording-only changes. The final reviewer-feedback renderer
patch passed 582 tests, Ruff lint, and strict mypy; after formatter-only normalization, Ruff lint/formatting
and `git diff --check` were clean. H03, H04, and H05 were then rerun from the real cached resources and
showed the intended final output. No unresolved scientific-correctness, evidence-boundary, or usability
blocker remains.

### H01 result — coordinate-convention control

- **Observed:** one hg19→hg38 liftOver mapping to `chr22:15690406`, matching the historical coordinate.
- **Reverse/flanking interval:** the centered 101-bp interval mapped 101/101 bases through the same chain,
  and reverse liftOver returned exactly to the original hg19 point.
- **Typed context:** both source and mapped coordinates overlapped the UCSC Segmental Duplications track.
  The overlap was presented descriptively and did not alter the mapping result.
- **Current human output:** `ONE LIFTOVER MAPPING`; ordinary primary-assembly metadata is not promoted into
  key findings, while the completed metadata check remains visible under `CHECKS PERFORMED`.
- **Adjudication:** no blocker. The output distinguishes coordinate conversion from unassessed variant,
  gene, transcript, uniqueness, or other biological identity.

### H02 result — invalid source coordinate

- **Observed:** UCSC hg38 chr2 has 242,193,529 bases; the requested point was 177 bp beyond the sequence
  end.
- **Behavior:** source preflight stopped before liftOver, reported `INVALID SOURCE COORDINATE`, stated that
  liftOver was not attempted, and returned process exit status 1.
- **Adjudication:** no blocker. Invalid input does not become a biological-looking no-mapping result.

### H03 result — telomeric interchromosomal mapping

- **Observed:** hg38 `chr10:10709` mapped completely to hg19 `chr18:10905`; the centered 101-bp interval
  also mapped 101/101 bases through the same forward chain.
- **Reverse liftOver:** hg19 `chr18:10905` mapped to hg38 `chr18:10905`, not to the original hg38 chr10
  coordinate.
- **Typed context:** both source and mapped coordinates overlapped the UCSC Segmental Duplications track.
  One overlapping hg38 record pairs the source chr10 region with an hg38 chr18 region. The renderer does
  not claim that this annotation explains the liftOver relationship.
- **Assembly metadata:** exact version-bound NCBI metadata for the hg19 target assembly was unavailable;
  mapping continued without inferring a target sequence role.
- **Current human output:** `INTERCHROMOSOMAL LIFTOVER MAPPING`, with reverse liftOver, flanking-interval
  behavior, and Segmental Duplications context stated separately. Because reverse liftOver returns
  elsewhere under `LIFTOVER_ONLY`, the footer also points to `--evidence-tier COMPARATIVE` as optional
  broader all-chain context, explicitly qualified as available only when the assembly pair supports it.
- **Adjudication:** no blocker. The chromosome change is explicit without being labeled erroneous or
  biologically correct.

### H04 result — new canine COMPARATIVE case

- **Observed:** the 3-bp CanFam3 locus has seven complete mappings in the UCSC all-chain alignments; every
  mapping covers 3/3 input bases and the mappings span four canFam4 chromosomes.
- **Comparative UCSC evidence:** the standard canFam3→canFam4 liftOver chain retains one of those seven
  mappings, `canFam4 chr9:49251380-49251382`. That mapping is also represented by a top-level net fill
  and all 3/3 input bases are present in the reciprocal-best chain. None of the other six complete
  mappings have the same combination.
- **Evidence dependence:** the standard liftOver chain, net, and reciprocal-best chain are explicitly
  described as related UCSC alignment evidence rather than independent confirmations.
- **Optional dimensions:** reverse liftOver was unavailable from the cached resources used for this run;
  UCSC Segmental Duplications context was also unavailable. Neither absence changed the comparative
  mapping facts. After outside review, reverse-liftOver unavailability was elevated into `KEY FINDINGS`
  so a skimming user does not have to infer that the round-trip check was missing.
- **Implementation influence:** H04 drove the compact multiple-mapping/comparative renderer. The default
  output now reports the seven complete mappings, shows the one mapping distinguished by the consumed
  UCSC comparative resources, and leaves the other six coordinates to `--details` rather than dumping
  them into ordinary terminal output.
- **Adjudication:** no blocker. The output says that the comparative evidence distinguishes one mapping
  within the assessed alignment resources; it does not declare that mapping biologically correct.

### H05 result — rs138257042 asymmetric chr22/chr14 mapping

- **Observed:** hg38 `chr22:15528888` mapped completely to hg19 `chr14:19378323`; the centered 101-bp
  interval mapped 101/101 bases through the same chain.
- **Reverse liftOver:** hg19 `chr14:19378323` mapped to hg38 `chr14:18601846`, not to the original hg38
  chr22 coordinate.
- **Typed context:** both source and mapped coordinates overlapped the UCSC Segmental Duplications track.
  One overlapping hg38 record pairs the source chr22 region with an hg38 chr14 region. The renderer does
  not infer a causal mechanism from that annotation.
- **Assembly metadata:** exact version-bound NCBI metadata for the hg19 target assembly was unavailable;
  mapping continued without inferring a target sequence role.
- **Identity boundary:** the limitations state that the coordinate result does not establish uniqueness
  or preservation of the same variant, gene, transcript, or other biological feature. Variant identity
  was not assessed.
- **Current follow-up guidance:** because reverse liftOver returns elsewhere under `LIFTOVER_ONLY`, the
  footer points to `--evidence-tier COMPARATIVE` for optional standard-chain/all-chain comparison when
  that evidence tier is available for the assembly pair.
- **Adjudication:** no blocker. H05 demonstrates the intended boundary between chain-based coordinate
  conversion and identifier-aware variant evidence.

### Internal gate disposition

The current renderer has no unresolved M23 blocker across H01-H05. H01 and H04 influenced presentation,
and the final outside review caused one H04 prominence fix plus conditional H03/H05 navigation guidance.
The set therefore remains implementation-influencing evidence rather than untouched validation. None of
these corrections changed candidate generation, comparative classification, or scientific evidence
semantics.

One non-blocking structured-output follow-up remains: when reverse liftOver is unavailable, the live
status text can explain why while the durable dossier/JSON may retain only the unavailable state. This
does not prevent the primary mapping assessment from completing and is not an M23 release blocker.

## Outside-user/domain feedback and final adjudication

Outside feedback directly influenced the current renderer. The reviewer is an experienced
bioinformatician and frequent liftOver user with UCSC Genome Browser/annotation experience; the record
does not treat that feedback as an endorsement.

- On the original H01-style output, the reviewer said the result was not understandable enough to count
  as human-readable. That feedback directly motivated the first default-output renderer slice.
- On the revised H01 output, the reviewer no longer needed terminology explained and immediately reasoned
  about the Segmental Duplications context. The reviewer independently suggested evaluating UCSC Self
  Chain and asked to see more complex cases. Self Chain remains deliberately deferred to Milestone 25 /
  `v0.3.0a1`; it was not pulled into the M23 renderer gate.
- The final outside packet contained the current H03, H04, and H05 default outputs without supplying a
  preferred interpretation first. The reviewer found H03 and H05 immediately understandable: the
  interchromosomal headline, prominent reverse-liftOver result, complete flanking-interval mapping, and
  descriptive-only Segmental Duplications context conveyed the structural facts without implying that the
  mapped coordinate preserved variant, gene, transcript, or other biological identity. The reviewer
  specifically found that H05 maintained the coordinate-versus-variant-identity boundary.
- For H04, the reviewer found the seven complete all-chain mappings and the one standard-chain mapping with
  top-level net and reciprocal-best membership clear and appropriately neutral. The explicit statement that
  chain, net, and reciprocal-best evidence are related rather than independent confirmations also satisfied
  the intended provenance boundary.
- The reviewer identified one material usability gap in H04: reverse liftOver was unavailable, but that fact
  appeared only in run-status text and `CHECKS PERFORMED`, not in the prominent result narrative. The
  renderer was corrected so unavailable reverse liftOver is now a `KEY FINDINGS` item for the multiple-
  mapping case.
- The reviewer also suggested that non-reciprocal H03/H05 `LIFTOVER_ONLY` results should point users toward
  broader all-chain comparison. The renderer now adds a bounded `--evidence-tier COMPARATIVE` suggestion
  only for relevant LIFTOVER_ONLY reverse results and explicitly says `when available`; ordinary clean
  reciprocal mappings do not receive the extra guidance.

After those two reviewer-driven changes, H03, H04, and H05 were rerun against the real cached resources.
H04 now states reverse-liftOver unavailability under `KEY FINDINGS`, while H03/H05 show the conditional
COMPARATIVE follow-up guidance. The reruns exposed no new overclaim, scope error, or practical explanation
failure.

**Outside-review adjudication:** the frozen outside-user/domain criterion is satisfied. The feedback found
no unresolved scientific overclaim and no remaining explanation failure that would make the current output
unsafe or substantially misleading.

## Post-completion H03 feedback received during Milestone 24

Additional H03 feedback arrived after the Milestone 23 gate had closed. It does not retroactively reopen
Milestone 23, but it is recorded because it informs release-preparation presentation and future evidence
priorities.

- **Default-output density:** the reviewer still found H03 too wordy and summarized the key practical
  takeaway as essentially "doesn't map back." Milestone 24 therefore uses H03 as a compactness case in the
  public-language audit. Shortening must preserve the source and mapped coordinates, complete forward
  coverage, reverse destination, exact flanking interval and coverage, Segmental Duplications context,
  evidence boundaries, and relevant follow-up guidance. `--details` and `--json` retain the complete
  structured result.
- **Why the mapping is non-reciprocal:** the reviewer asked whether the explanation could involve a short
  or boundary-adjacent chain alignment, a repeat and its repeat class, exonic/coding versus intronic or
  intergenic context, or assembly component/new-sequence history. The evidence collected for H03 does not
  answer those questions. In particular, Segmental Duplications overlap remains descriptive context and is
  not a causal explanation.
- **Feature disposition:** these questions are post-release evidence candidates, not Milestone 24 features.
  Point/gap-boundary context is already represented in the Milestone 25 plan. Broader alignment-extent/edge
  context, repeat annotation, gene/transcript feature context, and assembly component/history context are
  retained as priority-to-be-determined candidates. One outside-user request establishes genuine demand
  to retain the ideas, not an automatic priority or release commitment.

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

The final packet consisted of H03, H04, and H05 as one comprehension and scope-boundary check. It was
reviewed on 2026-09-06, produced the two presentation findings recorded above, and was followed by real
reruns of all three affected cases. H01 and H04 had already influenced renderer implementation before the
final packet, so this milestone is explicitly **not** described as untouched held-out validation.

## Gate completion rule

Milestone 23 passes when:

- all five pre-registered cases have been executed and adjudicated against the questions above;
- no unresolved blocking failure remains;
- the outside-user/domain packet has been reviewed;
- any blocking outside feedback has been resolved and the affected cases rerun; and
- the final record explicitly states whether the held-out cases remained untouched or influenced any
  implementation/language changes.

**Final disposition:** all completion conditions are met. Milestone 23 passes. The five pre-registered cases
were all executed and adjudicated; the outside complex-case packet was reviewed; the resulting H04 and
H03/H05 presentation findings were resolved and rerun; no blocking failure remains; and this record
explicitly preserves the implementation influence of the held-out set.
