# Getting Started with liftAssess

This guide is for someone who wants to run liftAssess without first learning the
internal architecture of UCSC chain/net files.

For the complete implemented capability list, see [`FEATURES.md`](FEATURES.md). For the
scientific result semantics and invariants, see [`DESIGN.md`](DESIGN.md).

## 1. What liftAssess does

A normal liftOver operation asks:

> Where can this source interval map in the target assembly?

liftAssess asks a different follow-up question:

> What happened to this interval in the mapping resources, why do the reported observations matter, what can I inspect next, and what does the evidence not establish?

The current result is a factual profile rather than a single quality label. It reports dimensions such as source-input validity, number of mappings, exact source coverage, split/discontinuous alignment geometry, orientation, reverse liftOver when prepared reverse resources are available, target sequence role, genomic context, evidence availability, resource consumption, and provenance.

The default human report begins with a deterministic factual headline such as `ONE LIFTOVER MAPPING`, `PARTIAL LIFTOVER MAPPING`, or `MULTIPLE LIFTOVER MAPPINGS`. Routine results stay compact. When an unusual observation is important, its explanation is kept next to that finding and the report gives a bounded next step.

liftAssess does **not** claim that a target is the biologically correct locus, and it does not produce a numeric confidence score.

## 2. What you give the command

The CLI needs three positional arguments:

```text
assess-liftover SOURCE_DB TARGET_DB LOCUS
```

### `SOURCE_DB`

The UCSC database identifier for the assembly your coordinates are currently on, for
example `canFam3`.

### `TARGET_DB`

The UCSC database identifier for the assembly you want to assess mappings onto, for
example `canFam4`.

Direction matters. `canFam3 canFam4` and `canFam4 canFam3` are different comparisons.

### `LOCUS`

One source interval in this form:

```text
chromosome:start-end
```

Example:

```text
chr1:10000001-10000100
```

Comma-grouped browser-style coordinates are also accepted:

```text
chr16:12,345-12,400
```

CLI coordinates are **1-based and inclusive**, matching familiar UCSC display
coordinates. The start and end you type both belong to the interval.

Internally, liftAssess converts the locus immediately to **0-based, half-open**
coordinates. You normally do not need to think about the internal convention unless
you inspect JSON or use the Python API.

The CLI retains this single-locus form and also accepts batch input through `--bed` or `--interval-table`. BED rows retain their native 0-based, half-open coordinates. The simple interval-table form requires a `sequence`, `start`, `end` header (plus optional `label`) and uses the same 1-based, inclusive coordinates as the single-locus CLI. Both forms normalize into the same indexed batch engine; batch mode does not silently start a whole-chain fallback.

## 3. Install liftAssess

liftAssess requires Python 3.11 or newer. For the public alpha, install the
pre-release package from PyPI:

```bash
python -m pip install --pre liftassess
assess-liftover --help
```

The rest of this guide assumes the installed `assess-liftover` command. If you are
working from a source checkout instead, run:

```bash
uv sync
uv run assess-liftover --help
```

and prefix the assessment commands below with `uv run`.

## 4. Your first assessment

Start with the ordinary directional liftOver chain rather than the much larger
comparative resource set:

```bash
assess-liftover \
  canFam3 \
  canFam4 \
  chr1:10000001-10000100 \
  --evidence-tier LIFTOVER-ONLY
```

`LIFTOVER-ONLY` is an evidence-availability class, **not a lower-confidence mode**.
During the 2026-08-30 release-readiness check, UCSC advertised a 5.4 MiB
`canFam3ToCanFam4.over.chain.gz` transfer for this exact example. Provider resources can
change, so review the displayed transfer plan rather than treating that size as a
permanent guarantee.

After learning the basic workflow, the real comparative mechanical fixture exercises
the complete comparative path:

```bash
assess-liftover \
  canFam3 \
  canFam4 \
  chrUn_JH373233:1845736-1845835
```

That comparative fixture currently requires an approximately 2.50 GiB compressed
resource set and can require substantial runtime. See [`PERFORMANCE.md`](PERFORMANCE.md)
before using it as a speed test.

### Choose a resource workflow

The CLI does not require the same network behavior for every run. Choose the workflow
that matches where your resources already live:

| What you want to do | Recommended path | UCSC/provider access |
| --- | --- | --- |
| Let liftAssess obtain the ordinary directional liftOver chain | `--evidence-tier LIFTOVER-ONLY` | Used only when the required resources are not already cached |
| Let liftAssess obtain comparative chain/net/reciprocal-best resources | default selection or `--evidence-tier COMPARATIVE` | Used when the comparative resources are not already cached |
| Repeat work from the verified cache | rerun the same command | Cache-first; missing metadata/context may still be acquired unless `--offline` is used |
| Guarantee zero provider access | add `--offline` | None; required resources must already be cached |
| Force a freshness check | add `--refresh` | Provider access required |
| Assess local BED/TSV records | `--bed` or `--interval-table` after `prepare-liftassess-index` | None during batch assessment |
| Use chain/net files you already have outside the liftAssess cache | public Python API | None automatically; you provide files and provenance explicitly |

`LIFTOVER-ONLY` and `COMPARATIVE` describe which UCSC alignment resources are used.
They are not low- and high-confidence modes. The explicit-local-file API is separate from
the CLI: there is currently no `--chain-file` or `--net-file` option. Full examples are in
[Usage scenarios](#15-usage-scenarios).

A first automatic single-locus run goes through these high-level stages:

1. **Validate the source interval.** Authoritative UCSC `chromInfo` metadata validates
   the submitted source sequence and bounds before mapping is attempted; exact
   `chromAlias` matches may be suggested but are never silently rewritten.
2. **Check the local mapping-resource cache.** If a complete verified mapping resource set
   already exists, liftAssess can use it without contacting UCSC for mapping resources.
3. **Discover UCSC mapping resources** if no reusable mapping resource set is available.
4. **Show the applicable UCSC mapping-resource terms.** You must explicitly acknowledge
   them before provider metadata inspection or acquisition for those resources.
5. **Inspect transfer metadata.** liftAssess uses body-free HTTP HEAD requests to show
   provider-advertised sizes and other transfer information when available.
6. **Show the transfer plan.** You separately accept the exact mapping-resource plan
   before acquisition begins.
7. **Acquire and verify mapping resources.** Verified artifacts are stored outside the
   source tree in the liftAssess cache.
8. **Prepare optional target sequence-role context.** When an exact versioned UCSC/NCBI assembly
   binding is available, target sequence role is attached; inability to obtain
   this optional context does not invalidate the mapping assessment.
9. **Assess the locus.** Mappings and the evidence available for the selected
   resource class are extracted and passed to deterministic result-profile derivation.
   `LIFTOVER-ONLY` uses chain evidence; `COMPARATIVE` can additionally attach supported
   net/reciprocal-best observations.
10. **Attach available post-assessment context.** Depending on prepared/cached resources
    and run mode, liftAssess can attach the standard liftOver-versus-all-chain comparison, flanking-interval
    context, reverse liftOver, and UCSC Segmental Duplications observations.
    These dimensions report clearly when a check was not run or its resources were unavailable and do not become
    confidence votes or biological-mechanism claims.
11. **Print the report.** The default is the concise human-readable summary.

Cancelling either acknowledgement stops before the corresponding network/resource
operation.

## 5. How to read the default result

A routine successful mapping now has this general shape. Exact coordinates and optional
context depend on the assembly pair, locus, selected resources, and prepared cache:

```text
* ONE LIFTOVER MAPPING *

Source:
    canFam3 chr1:10000001-10000100

Mapped interval:
    canFam4 chr1:10027740-10027839

= KEY FINDINGS =

Mapping:
    One liftOver chain maps 100/100 bp.

= LIMITATIONS =

This result does not establish that the mapped coordinate represents the same variant,
gene, transcript, or other biological feature.

Details:
    Use --details for alignment blocks, gaps, evidence, and provenance;
    use --json for machine-readable output.
```

Unusual findings expand in place. For example, if reverse liftOver does not return to
the source locus, the report explains immediately under that finding why the lack of a
one-to-one reciprocal correspondence matters for biological interpretation. If a mapping
overlaps UCSC Segmental Duplications, that finding likewise carries its own interpretation
boundary. A later `NEXT STEP` section can then combine those observations into concrete
follow-up actions such as Genome Browser review or target-assembly-specific feature
verification.

The default summary keeps each reported item label on its own line and indents the
corresponding value by four spaces so the result remains easy to scan.

### Factual headline

Names the observed mapping event. It is deterministic result language, not a confidence rating or biological judgment.

### `Source` and mapping coverage

Confirm exactly what interval was assessed and how many requested source bases are represented. For unusual results, the summary expands with uncovered source intervals, mapped-segment counts, target gaps, or multiple mappings as needed.

### `Target`

Shows the target coordinate for a single mapping. When a mapping contains multiple mapped segments, the report explicitly identifies the displayed target interval as a **bounding span** rather than implying continuous alignment.

### `Evidence`

Tells you what kind of evidence liftAssess examined for this assembly pair and which resource roles were consumed.

- `COMPARATIVE` means comparative UCSC resources were available for the current evidence path.
- `LIFTOVER-ONLY` means only directional chain mapping evidence was available.

These are evidence-availability concepts, **not confidence tiers**. For `COMPARATIVE`, the report also warns that UCSC-derived observations are conservatively treated as dependent, not independent votes, and exact shared processing-run provenance is not verified.

### `Reverse liftOver`

When a matching reverse-direction chain and its prepared index are already in the local cache, liftAssess maps each exact target alignment segment back toward the source assembly. This asks whether the forward mapping participates in a one-to-one reciprocal coordinate relationship under the available forward and reverse chain mappings. Fragmented forward mappings are checked segment-by-segment; the target bounding span is never used to manufacture a query across an unaligned gap. Net and reciprocal-best resources are not required for this reverse geometry.

Why this matters: if a mapped coordinate does **not** return to the original locus, the source and mapped loci do not have a one-to-one reciprocal correspondence under those available chain mappings. The forward mapping is not necessarily wrong, but coordinate conversion alone is then weaker evidence for treating the target coordinate as the same variant, gene, transcript, or other biological feature. This interpretation is grounded in the genome-build conversion literature summarized in [`REFERENCES.md`](REFERENCES.md).

Reverse liftOver is distinct from UCSC reciprocal-best membership. If no matching reverse chain is cached, the check is unavailable. If the matching chain is cached but its prepared index is absent or unusable, the check is not run; normal assessment does not fall back to a surprise exhaustive scan. During an explicit `--refresh`, reverse liftOver is also not run rather than combining freshly reacquired forward resources with an unrefreshed reverse chain. Prepare the matching reverse-direction index explicitly with `prepare-liftassess-index TARGET_DB SOURCE_DB --evidence-tier COMPARATIVE` or `--evidence-tier LIFTOVER-ONLY`. Normal assessment does not silently download reverse resources or build an index.

### Flanking-interval context for point queries

For a 1-bp query, liftAssess automatically requests a centered 101-bp flanking interval from the same
prepared forward chain index used for fast mapping lookup. The summary always states the exact
source interval actually tested. Near a sequence boundary, the tested window may be shorter than 101
bases because it is clipped rather than shifted.

This automatic flanking-interval check is **forward chain only**. Its point/flanking-interval relationship is derived
from mapping identity and structural geometry, not from chain-score ranking or a hidden
threshold. For a `COMPARATIVE` point, net and reciprocal-best evidence may be available for the point
itself, but those resources are not re-run
for the flanking interval. If the matching forward chain index is missing, unusable, or cannot provide a
safe source bound, the flanking-interval check reports that it was not run and does not start another whole-chain scan.

The 101-bp default is contextual evidence, not a confidence threshold. To request a different larger odd-width
window for a point, use for example:

```bash
assess-liftover hg19 hg38 chr1:120904787-120904787 --context-bases 1001
```

Ordinary interval queries are not widened automatically, and liftAssess never recursively expands a
point from 101 bp to 1 kb or 10 kb because the first context result looks unusual.

### Target sequence role and genomic context

When UCSC's assembly description provides an exact versioned NCBI assembly accession, liftAssess can attach provider-native target sequence role from the matching NCBI Datasets sequence report. If the exact binding cannot be established, the role remains unavailable; liftAssess does not guess from names such as `_alt`, `_random`, or `chrUn`.

Single-locus runs can also report exact source and mapped-target overlap with UCSC's assembly-specific `genomicSuperDups` table. Those observations are descriptive context only. They do not penalize a mapping, prove that duplication caused the mapping behavior, or establish which locus is biologically correct.

Both context sources preserve exact resource provenance. If optional context is missing or unusable, the already-valid primary coordinate assessment continues and the unavailable dimension is reported explicitly.

### `Interpretation`

A deterministic sentence that stays close to the measured geometry/evidence. It does not choose a biologically correct locus.

### `Scope`

States important evidence boundaries so an untested identity question is not implied to have been answered by coordinate mapping.

## 6. Progressive disclosure

Routine one-mapping results stay compact. The current progressive renderer surfaces material facts when present, including:

- partial source coverage;
- fragmented or target-discontinuous mapping geometry;
- multiple liftOver mappings;
- reverse liftOver relationships;
- point-versus-flanking-interval relationships;
- comparative mapping relationships;
- target sequence role; and
- genomic-context observations.

For large intervals and multiple mappings, source coverage leads the story. The report gives actual measured coverage; it does not apply a built-in 90% or other quality threshold.

Reverse liftOver context is reported when the matching prepared reverse index is available. For 1-bp point queries, liftAssess also requests automatic 101-bp flanking-interval context when the prepared forward chain index is available. Batch mode reports cross-record exact target collisions and overlapping-but-offset mappings from indexed chain records, and one-base rows from either supported batch format receive the same automatic flanking-interval context from that index. Context-scale exact collisions are reported separately as flanking-interval target collisions. COMPARATIVE batches now attach shared net/reciprocal-best evidence to submitted rows. Target sequence-role metadata is cache-only in batch mode; reverse batch evidence and genomic context are not currently assessed there and are not silently implied by the output.

## 7. Ask for the full human-readable report

Use `--details` when you need to understand or audit the evidence behind the summary:

```bash
assess-liftover \
  canFam3 canFam4 chr1:10000001-10000100 \
  --evidence-tier LIFTOVER-ONLY \
  --offline \
  --details
```

The detailed report includes:

- the complete currently available factual result profile;
- authoritative source-validation metadata and provenance;
- explicit states for result dimensions that were not run or are not yet assessed;
- every mapping and exact mapped segment;
- exact uncovered source intervals and target gaps;
- mapping orientation;
- target sequence role when available;
- genomic-context observations and provenance when available;
- every evidence observation;
- chain/net/reciprocal-best detail;
- resource URLs, cache paths, retrieval metadata, and checksums;
- which cached resources were actually consumed by the engine; and
- the provenance dependency graph showing which observations share upstream sources.

Use this mode when the compact progressive summary omits detail needed for scientific review.

## 8. Get JSON for scripts and pipelines

Use `--json` for the complete schema-v2 machine-readable report:

```bash
assess-liftover \
  canFam3 canFam4 chr1:10000001-10000100 \
  --evidence-tier LIFTOVER-ONLY \
  --offline \
  --json > assessment.json
```

The JSON document goes to **stdout**. Status and progress messages go to **stderr**, so
normal shell redirection does not mix progress text into the JSON file. Schema v2 carries
the same factual profile, exact mappings/evidence, resources, and provenance used by the
human renderer. It intentionally omits obsolete aggregate-result and mapping-selection fields from
the earlier alpha schema. Detailed and JSON reports include local
cache paths as run context, so inspect that metadata before publishing a report if local
filesystem paths are information you do not want to share.

Important coordinate difference:

- CLI and human-readable display intervals are 1-based, inclusive.
- JSON interval objects are 0-based, half-open and explicitly state their coordinate
  system.

`--json` and `--details` cannot be used together.

## 9. Understand the cache before downloading large resources

The default cache is outside the repository:

- macOS: `~/Library/Caches/liftassess`
- Windows: `%LOCALAPPDATA%\liftassess\Cache`
- Linux/other Unix-like systems: `$XDG_CACHE_HOME/liftassess`, falling back to
  `~/.cache/liftassess`

To use another location:

```bash
assess-liftover \
  canFam3 canFam4 chr1:10000001-10000100 \
  --evidence-tier LIFTOVER-ONLY \
  --cache-dir /path/to/liftassess-cache
```

Cache reuse is verified: liftAssess hashes the stored artifact bytes before accepting a
complete resource set for assessment.

The cache is content-addressed by SHA-256. The local path is run context; the digest is
the exact artifact identity recorded in provenance.

## 10. Run with guaranteed zero provider access

Once the needed resource set is cached, use `--offline`:

```bash
assess-liftover \
  canFam3 canFam4 chr1:10000001-10000100 \
  --evidence-tier LIFTOVER-ONLY \
  --offline
```

`--offline` is a guarantee, not merely a preference. liftAssess will fail instead of
contacting UCSC if it cannot find a complete verified local resource set for that direction.

A normal cache-first run also avoids UCSC when a complete verified mapping resource set is present,
but `--offline` makes that requirement explicit.

## 11. Deliberately check current provider resources

Use `--refresh` when you specifically want liftAssess to contact UCSC and reacquire the
current resource bytes instead of accepting cache-first reuse:

```bash
assess-liftover \
  canFam3 canFam4 chr1:10000001-10000100 \
  --evidence-tier LIFTOVER-ONLY \
  --refresh
```

`--refresh` and `--offline` are mutually exclusive because they request opposite
network behavior.

## 12. Use liftAssess non-interactively

Interactive terms and transfer-plan confirmations are the default. For a script or
other non-interactive run, explicit flags can supply those acknowledgements:

```bash
assess-liftover \
  canFam3 canFam4 chr1:10000001-10000100 \
  --evidence-tier LIFTOVER-ONLY \
  --acknowledge-ucsc-terms \
  --accept-transfer-plan \
  --json > assessment.json
```

Use these flags only when the applicable provider terms and planned transfer are
actually acceptable for the workflow. They skip the prompts; they do not change the
provider terms or relicense external resources.

## 13. Reduce terminal noise

Use `--quiet` to suppress nonessential status and measured progress displays:

```bash
assess-liftover \
  canFam3 canFam4 chr1:10000001-10000100 \
  --evidence-tier LIFTOVER-ONLY \
  --quiet
```

`--quiet` does **not** silently accept provider terms or a transfer plan. Required
acknowledgements still remain unless supplied explicitly with the acknowledgement
flags.

## 14. Complete CLI option reference

| Option | What it does | When to use it |
| --- | --- | --- |
| `-h`, `--help` | Shows the command syntax and option help | Quick command reference |
| `--cache-dir PATH` | Uses a specific resource cache | Shared storage, testing, or keeping large resources on another disk |
| `--evidence-tier {COMPARATIVE,LIFTOVER-ONLY}` | Requires one exact UCSC mapping-resource class instead of the default COMPARATIVE-preferred selection | Reproducible resource selection or preparing the standard-liftOver side of a comparison |
| `--offline` | Guarantees zero provider access and requires a complete verified cache | Reproducible offline analysis or restricted network environments |
| `--refresh` | Contacts UCSC and reacquires current resources instead of cache-first reuse | Explicit freshness checks |
| `--acknowledge-ucsc-terms` | Supplies the explicit terms acknowledgement without a prompt | Non-interactive workflows after terms review |
| `--accept-transfer-plan` | Supplies the separate transfer-plan acknowledgement without a prompt | Non-interactive workflows after reviewing the planned transfer |
| `--bed PATH` | Reads BED3-or-later batch input using native 0-based, half-open coordinates; `-` reads stdin | Indexed BED batch assessment |
| `--interval-table PATH` | Reads a tab-delimited `sequence/start/end[/label]` batch table using 1-based, inclusive coordinates; `-` reads stdin | Spreadsheet-style indexed batch assessment |
| `--details` | Prints the full single-locus human-readable evidence/resource/provenance report; not yet available with batch input | Scientific inspection and debugging |
| `--json` | Prints schema-v2 machine-readable output | Scripts, archives, downstream analysis |
| `--quiet` | Suppresses nonessential terminal progress/status | Logs, scripts, or less terminal output |
| `--context-bases N` | Uses another odd-width flanking-interval window for 1-bp point queries, including one-base batch rows | Explicitly testing a point at a scale other than the automatic 101-bp default |

`--offline` and `--refresh` are mutually exclusive. `--details` and `--json` are also
mutually exclusive.

## 15. Usage scenarios

The examples below separate **where the mapping resources come from** from **how you want
the result rendered**. You can add `--details`, `--json`, `--quiet`, `--context-bases`, or a
custom `--cache-dir` to the compatible CLI scenarios without changing their scientific
meaning.

### Scenario A: first run, let liftAssess obtain the standard UCSC liftOver chain

Use the ordinary directional liftOver chain when you want the smallest mapping-resource
workflow:

```bash
assess-liftover \
  SOURCE_DB TARGET_DB CHR:START-END \
  --evidence-tier LIFTOVER-ONLY
```

If the required files are not cached, liftAssess displays the applicable UCSC terms,
inspects provider transfer metadata, shows the exact transfer plan, and asks for the two
explicit acknowledgements before acquisition. Verified files are stored in the liftAssess
cache for reuse.

### Scenario B: first run, use comparative UCSC alignment resources

Request the comparative resource set explicitly:

```bash
assess-liftover \
  SOURCE_DB TARGET_DB CHR:START-END \
  --evidence-tier COMPARATIVE
```

Or omit `--evidence-tier` to use the default COMPARATIVE-preferred selection when a
complete comparative set is available. Comparative assessment can add ordinary-net and
reciprocal-best observations and, for single loci with the required prepared indexes, a
standard-liftOver-versus-all-chain comparison. These related UCSC resources are not treated
as independent votes.

### Scenario C: repeat an assessment from verified cache

Rerun the same command:

```bash
assess-liftover SOURCE_DB TARGET_DB CHR:START-END
```

The default CLI is cache-first. Mapping resources that are already complete and verified
are reused instead of being downloaded again. A normal cache-first run can still contact a
provider for missing metadata or optional context; use `--offline` when zero provider
access is a requirement.

### Scenario D: guarantee a completely local CLI run

```bash
assess-liftover \
  SOURCE_DB TARGET_DB CHR:START-END \
  --evidence-tier LIFTOVER-ONLY \
  --offline
```

`--offline` guarantees zero provider access. The command fails rather than silently
contacting UCSC or another provider if required cached resources are missing.

Use `--cache-dir PATH` in both the acquisition run and the later offline run when you keep
liftAssess resources on shared storage or another disk.

### Scenario E: deliberately refresh from UCSC

```bash
assess-liftover \
  SOURCE_DB TARGET_DB CHR:START-END \
  --evidence-tier LIFTOVER-ONLY \
  --refresh
```

`--refresh` deliberately contacts UCSC and reacquires current resource bytes instead of
accepting cache-first reuse. It is mutually exclusive with `--offline`.

### Scenario F: prepare an index for repeated local work

After the assembly-pair resources have been acquired and verified at least once:

```bash
prepare-liftassess-index SOURCE_DB TARGET_DB
```

This command is cache-only: it never contacts UCSC and does not require a locus. It parses
the selected chain once to build a reusable local index bound to that chain's exact
SHA-256 identity. Large chains can take many minutes and several GiB of additional cache
space to prepare, so index construction is explicit rather than an automatic first-query
pause.

Select the exact resource class when both are cached:

```bash
prepare-liftassess-index SOURCE_DB TARGET_DB --evidence-tier COMPARATIVE
prepare-liftassess-index SOURCE_DB TARGET_DB --evidence-tier LIFTOVER-ONLY
```

Automatic reverse liftOver requires a matching reverse-direction index. For a forward
`SOURCE_DB` → `TARGET_DB` assessment, prepare `TARGET_DB` → `SOURCE_DB` with the same
resource class if you want the reverse check to run.

### Scenario G: assess a local BED file as a batch

BED3-or-later input keeps native **0-based, half-open** coordinates:

```bash
assess-liftover \
  SOURCE_DB TARGET_DB \
  --bed loci.bed
```

Use `--bed -` for stdin. Batch mode is deliberately cache-only and index-only. It never
contacts UCSC, runs `--refresh`, builds an index automatically, or falls back to a full
chain traversal.

### Scenario H: assess a local interval table as a batch

For spreadsheet-style input, use a tab-delimited table with a required header:

```text
sequence	start	end	label
chr1	101	200	region-a
chr2	500	500	point-b
```

```bash
assess-liftover \
  SOURCE_DB TARGET_DB \
  --interval-table loci.tsv
```

The table uses **1-based, inclusive** coordinates like the single-locus CLI. `start == end`
is therefore a valid one-base point. Use `--interval-table -` for stdin.

For both batch formats, prepare the selected chain index first. Without an explicit
`--evidence-tier`, batch mode prefers a complete cached COMPARATIVE resource set with a
prepared all-chain index and otherwise uses an available LIFTOVER-ONLY chain. COMPARATIVE
batches attach ordinary-net and reciprocal-best observations to submitted rows but do not
run the single-locus standard liftOver-versus-all-chain comparison or reverse liftOver per row.

### Scenario I: use an explicit local UCSC chain file through Python

If you already have a UCSC chain file and do not want liftAssess to discover or acquire it,
use the public Python API. There is currently no `--chain-file` CLI option.

The local-file boundary accepts plain-text or gzip-compressed resources. It computes and
verifies SHA-256 provenance for the exact bytes consumed. This is a lower-level integration
path: it does not automatically perform the CLI's UCSC source-metadata validation, provider
resource discovery/acquisition, reverse-index preparation, target sequence-role lookup, or
Segmental Duplications acquisition. Callers can compose those public APIs separately when
needed. `GenomicInterval` uses canonical **0-based, half-open** coordinates:

```python
from pathlib import Path

from liftassess import (
    AssemblyIdentifier,
    EvidenceAvailabilityTier,
    GenomicInterval,
    build_result_profile,
    build_ucsc_candidates_from_files,
    provenance_source_for_file,
)

source = AssemblyIdentifier(name="hg38", provider="UCSC")
target = AssemblyIdentifier(name="hg19", provider="UCSC")
interval = GenomicInterval(source, "chr10", 10708, 10709)
chain_path = Path("/data/hg38ToHg19.over.chain.gz")

chain_provenance = provenance_source_for_file(
    chain_path,
    label="local hg38 to hg19 liftOver chain",
    derived_from=(),
)

mappings = build_ucsc_candidates_from_files(
    interval,
    chain_path,
    target_assembly=target,
    chain_provenance=chain_provenance,
)

profile = build_result_profile(
    interval,
    mappings,
    evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    consumed_resource_roles=("CHAIN",),
)

print(profile.interpretation)
```

The exact public API name `build_ucsc_candidates_from_files()` retains the package's
internal model vocabulary; the returned objects represent the liftOver mappings described
throughout the user documentation.

`derived_from=()` means that you are making **no claim about shared upstream alignment
provenance**. If you know that several local resources derive from one upstream alignment
process, create a shared `ProvenanceSource` and pass it as `derived_from=(alignment,)` for
each file instead of pretending the files are independent evidence.

### Scenario J: use local chain, net, and reciprocal-best files

The same Python boundary accepts comparative resources you already have:

```python
from pathlib import Path

from liftassess import (
    AssemblyIdentifier,
    GenomicInterval,
    ProvenanceSource,
    ReciprocalBestResourceCompleteness,
    build_ucsc_candidates_from_files,
    provenance_source_for_file,
)

source = AssemblyIdentifier(name="canFam3", provider="UCSC")
target = AssemblyIdentifier(name="canFam4", provider="UCSC")
interval = GenomicInterval(source, "chr9", 49_252_490, 49_252_493)

alignment = ProvenanceSource(
    source_id="local-canFam3-canFam4-alignment",
    label="shared upstream canFam3 to canFam4 alignment",
)

chain_path = Path("/data/canFam3.canFam4.all.chain.gz")
net_path = Path("/data/canFam3.canFam4.net.gz")
rbest_path = Path("/data/canFam3.canFam4.rbest.chain.gz")

mappings = build_ucsc_candidates_from_files(
    interval,
    chain_path,
    target_assembly=target,
    chain_provenance=provenance_source_for_file(
        chain_path,
        label="local all-chain alignment",
        derived_from=(alignment,),
    ),
    net_path=net_path,
    net_provenance=provenance_source_for_file(
        net_path,
        label="local net alignment",
        derived_from=(alignment,),
    ),
    reciprocal_best_chain_path=rbest_path,
    reciprocal_best_provenance=provenance_source_for_file(
        rbest_path,
        label="local reciprocal-best chain",
        derived_from=(alignment,),
    ),
    reciprocal_best_completeness=(ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE),
)
```

Supplying local files bypasses automatic provider discovery, terms prompting, download
planning, and cache acquisition. It does **not** waive the provider's license or usage
terms; the caller remains responsible for using those files under the applicable terms.
The API also cannot infer upstream alignment/process provenance from file bytes alone.

### Scenario K: choose the report form

For a complete machine-readable single-locus or batch report:

```bash
assess-liftover \
  SOURCE_DB TARGET_DB CHR:START-END \
  --json > assessment.json
```

For the complete human-readable single-locus evidence/provenance report:

```bash
assess-liftover \
  SOURCE_DB TARGET_DB CHR:START-END \
  --details
```

`--details` is not yet implemented for batch input. `--details` and `--json` are mutually
exclusive. Status and progress go to stderr while JSON goes to stdout, so ordinary shell
redirection is safe.

## 16. Common mistakes

### Swapping source and target

The first assembly is where the input coordinates currently live. The second is where
mappings are assessed. Reverse the arguments and you are asking a different
question.

### Using the wrong coordinate convention

The CLI expects 1-based, inclusive coordinates. If your interval is already 0-based,
half-open, convert it before putting it on the command line. For the same physical
interval, add one to the 0-based start and keep the half-open end as the inclusive
end.

### Treating `COMPARATIVE` as "high confidence"

It is only an evidence-availability tier. Read the factual mapping profile and exact
evidence separately.

### Treating one complete mapping as biological truth

A complete chain mapping means every requested source base is represented in that
chain relationship. It is not an orthology call, identity check, uniqueness claim, or
proof that the locus is biologically correct.

### Assuming mapping order is rank

Current mapping order is retained for reproducibility. liftAssess does not yet define
or emit mapping-rank evidence.

### Expecting `--offline` to work before the resources are cached

It will fail by design rather than contact UCSC.

### Expecting `--quiet` to bypass acknowledgements

It suppresses nonessential status/progress only. Provider terms and transfer-plan
acknowledgements remain explicit.

### Reading progress as percent of the scientific algorithm

Progress bars measure exact bytes transferred, hashed, or read. They are not an ETA or
an estimate of how close the scientific computation is to completion.

## 17. Current performance warning

Without a prepared chain index, single-locus assessment still streams and parses the complete
chain resource. Measured work on the `canFam3` → `canFam4` pair showed that this can dominate runtime
even for a small, ordinary locus. `prepare-liftassess-index` turns that repeated whole-chain parse
into a one-time explicit preparation step; later mapping lookup is region-addressable while the
original chain remains the scientific provenance source.

Net and reciprocal-best evidence are still read through their existing resource paths when required,
so indexed chain lookup does not imply that every comparative-evidence operation is already indexed.
See [`PERFORMANCE.md`](PERFORMANCE.md) for measured benchmark results and the current architecture.

## 18. Current scientific/use envelope

For v1, liftAssess is designed around **same-species assembly comparisons**. Assemblies
from different individuals can contain real structural differences, so a disagreement
between mappings is not automatically an error.

The current tool also does **not**:

- compute fresh sequence identity from the raw genomes;
- run a new aligner;
- infer orthology automatically;
- use machine learning;
- produce a numeric confidence score;
- provide flanking-gene synteny evidence yet;
- define mapping-rank evidence yet; or
- attach reverse evidence across batches yet.

See [`FEATURES.md`](FEATURES.md) for the complete implemented/non-implemented catalog.

## 19. Expert Python use

The CLI is the easiest entry point. The package also exposes lower-level Python APIs
for callers that already have normalized mappings, local chain/net resources, or a
verified cached resource set.

The main boundaries are:

- `build_result_profile()` — derive the factual profile from normalized mappings and evidence;
- `build_ucsc_candidates_from_files()` — build mappings from explicit local UCSC
  resources plus provenance;
- `assess_ucsc_cached_bundle()` — assess a verified cached resource set;
- `discover_ucsc_resources()` — discover provider resources; and
- the resource planning/acquisition, cache, checksum, and provenance helpers documented
  in [`FEATURES.md`](FEATURES.md).

These APIs are intended for expert integration. There is currently no plugin registry
or automatic engine-discovery framework.

## 20. Where to go next

- [`FEATURES.md`](FEATURES.md) — complete current feature catalog and non-features.
- [`DESIGN.md`](DESIGN.md) — authoritative scientific semantics and architecture.
- [`ROADMAP.md`](ROADMAP.md) — implementation history and upcoming work.
- [`PERFORMANCE.md`](PERFORMANCE.md) — measured performance and optimization priorities.
- [`REFERENCES.md`](REFERENCES.md) — scientific, provider, format, and implementation
  sources used by the project.
