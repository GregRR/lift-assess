# Getting Started with liftAssess

This guide is for someone who wants to run liftAssess without first learning the
internal architecture of UCSC chain/net files.

For the complete implemented capability list, see [`FEATURES.md`](FEATURES.md). For the
scientific result semantics and invariants, see [`DESIGN.md`](DESIGN.md).

## 1. What liftAssess does

A normal liftOver operation asks:

> Where can this source interval map in the target assembly?

liftAssess asks a different follow-up question:

> What physically happened to this interval in the consumed mapping resources, what evidence was examined, and what does that evidence not establish?

The current result is a factual profile rather than a single quality verdict. It reports dimensions such as projection count, exact source coverage, split/discontinuous geometry, orientation, actual reverse-mapping context when cached reverse resources are available, evidence availability, resource consumption, and provenance.

The default human report begins with a deterministic factual headline such as `ONE COMPLETE CHAIN PROJECTION`, `PARTIAL SOURCE COVERAGE`, or `MULTIPLE CHAIN PROJECTIONS`. It then gives the few facts needed to understand that headline and a bounded interpretation.

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

A known real mechanical-fixture command is:

```bash
assess-liftover \
  canFam3 \
  canFam4 \
  chrUn_JH373233:1845736-1845835
```

This example is useful for validating the complete comparative path, but it is **not a
small first download**. The current `canFam3` → `canFam4` comparative bundle is about
2.50 GiB compressed, and the current implementation streams large resources during
single-locus assessment. See [`PERFORMANCE.md`](PERFORMANCE.md) before using this pair
as a speed test.

A first automatic run goes through these stages:

1. **Check the local cache.** If a complete verified bundle already exists, liftAssess
   can use it without contacting UCSC.
2. **Discover UCSC resources** if no reusable cache bundle is available.
3. **Show the applicable UCSC terms.** You must explicitly acknowledge them before
   provider metadata inspection or acquisition.
4. **Inspect transfer metadata.** liftAssess uses body-free HTTP HEAD requests to show
   provider-advertised sizes and other transfer information when available.
5. **Show the transfer plan.** You separately accept the exact resource plan before
   acquisition begins.
6. **Acquire and verify resources.** Verified artifacts are stored outside the source
   tree in the liftAssess cache.
7. **Assess the locus.** Candidate mappings and evidence are extracted and passed to
   deterministic result-profile derivation.
8. **Check prepared reverse-direction resources.** If the query produced candidates and a reverse-direction chain of the same publication class plus its prepared index are already cached, liftAssess runs an actual chain-only reverse assessment. No matching reverse chain is `UNAVAILABLE`; a matching chain without a usable prepared index is `NOT_RUN`. Normal assessment does not trigger an implicit provider request, index build, or exhaustive reverse-chain scan.
9. **Print the report.** The default is the concise human-readable summary.

Cancelling either acknowledgement stops before the corresponding network/resource
operation.

## 5. How to read the default result

A concise uncomplicated report has this general shape:

```text
* ONE COMPLETE CHAIN PROJECTION *
Source:
    chr1:10000001-10000100 (1-based inclusive)
Source coverage:
    100/100 bases
Target:
    chr1:10027740-10027839 (1-based inclusive; same orientation)
Reverse mapping:
    exactly reconstructs the original aligned source geometry
Evidence:
    COMPARATIVE — ...
Interpretation:
    ...
Scope:
    coordinate projection/structure assessed; variant/gene identity not assessed.
Details:
    use --details for the full profile/evidence or --json for schema v2.
This does not establish biological correctness.
```

The default summary keeps each reported item label on its own line and indents the
corresponding value by four spaces so the result remains easy to scan.

### Factual headline

Names the observed mapping event. It is deterministic result language, not a confidence rating or biological verdict.

### `Source` and `Source coverage`

Confirm exactly what interval was assessed and how many requested source bases are represented. For unusual results, the summary expands with uncovered source intervals, mapped-segment counts, target gaps, or multiple projections as needed.

### `Target`

Shows the target coordinate for a single projection. When a projection contains multiple mapped segments, the report explicitly identifies the displayed target interval as a **bounding span** rather than implying continuous alignment.

### `Evidence`

Tells you what kind of evidence liftAssess examined for this assembly pair and which resource roles were consumed.

- `COMPARATIVE` means comparative UCSC resources were available for the current evidence path.
- `LIFTOVER-ONLY` means only directional chain mapping evidence was available.

These are evidence-availability concepts, **not confidence tiers**. For `COMPARATIVE`, the report also warns that UCSC-derived observations are conservatively treated as dependent, not independent votes, and exact shared processing-run provenance is not verified.

### `Reverse mapping`

When a reverse-direction chain of the same publication class as the forward assessment and its prepared chain index are already in the local cache, liftAssess reverses each exact mapped target segment through that chain and reports whether those projections reconstruct the original aligned source geometry, land elsewhere, do both, or produce no reverse projection. Fragmented forward mappings are reversed segment-by-segment; the target bounding span is never used to manufacture a query across an unaligned gap. Net and reciprocal-best artifacts are not required for this reverse geometry.

This is distinct from UCSC reciprocal-best membership. If no matching reverse chain is cached, the result reports `UNAVAILABLE`. If the matching chain is cached but its prepared index is absent or unusable, the result reports `NOT_RUN`; normal assessment does not fall back to a surprise exhaustive scan. During an explicit `--refresh`, reverse mapping is also `NOT_RUN` rather than combining freshly reacquired forward resources with an unrefreshed reverse chain; reverse-direction refresh is never implicit. Prepare the exact reverse chain class explicitly with `prepare-liftassess-index TARGET_DB SOURCE_DB --evidence-tier COMPARATIVE` or `--evidence-tier LIFTOVER-ONLY`, matching the forward result. Normal assessment does not silently download reverse resources or build an index.

### `Local context`

For a 1-bp query, liftAssess automatically requests a centered 101-bp neighborhood from the same
prepared forward chain index used for fast candidate lookup. The summary always states the exact
source window actually tested. Near a sequence boundary, the tested window may be shorter than 101
bases because it is clipped rather than shifted.

This automatic neighborhood is **forward chain only**. Its point/context relationship is derived
from projection identity and structural geometry, not from chain-score ranking or a hidden
threshold. For a `COMPARATIVE` point, net and reciprocal-best evidence may be available for the point
itself, but those resources are not re-run
for the neighborhood. If the matching forward chain index is missing, unusable, or cannot provide a
safe source bound, local context reports `NOT_RUN` and does not start another whole-chain scan.

The 101-bp default is context, not a confidence threshold. To request a different larger odd-width
window for a point, use for example:

```bash
assess-liftover hg19 hg38 chr1:120904787-120904787 --context-bases 1001
```

Ordinary interval queries are not widened automatically, and liftAssess never recursively expands a
point from 101 bp to 1 kb or 10 kb because the first context result looks unusual.

### `Interpretation`

A deterministic sentence that stays close to the measured geometry/evidence. It does not choose a biologically correct locus.

### `Scope`

States important evidence boundaries so an untested identity question is not implied to have been answered by coordinate projection.

## 6. Progressive disclosure

Routine one-complete-projection results stay compact. The current first result-profile slice expands automatically when it can already detect:

- partial source coverage;
- fragmented or target-discontinuous projection geometry; or
- multiple chain projections.

For large intervals and multiple projections, source coverage leads the story. The report gives actual measured coverage; it does not apply a built-in 90% or other quality threshold.

Actual reverse-mapping context is reported when the matching prepared reverse index is available. For 1-bp point queries, liftAssess also requests automatic 101-bp local chain context when the prepared forward chain index is available. Batch mode reports cross-record exact target collisions and overlapping-but-offset projections from indexed chain candidates, and one-base rows from either supported batch format receive the same automatic point context from that index. Context-scale exact collisions are reported separately as neighborhood-level target collisions. COMPARATIVE batches now attach shared net/reciprocal-best evidence to submitted rows; reverse batch evidence, target-role metadata, and typed contextual evidence remain later roadmap capabilities. Those later checks are not silently implied by the current output.

## 7. Ask for the full human-readable dossier

Use `--details` when you need to understand or audit the evidence behind the summary:

```bash
assess-liftover \
  canFam3 canFam4 chrUn_JH373233:1845736-1845835 \
  --details
```

The detailed report includes:

- the complete currently available factual result profile;
- explicit states for result dimensions that were not run or are not yet assessed;
- every candidate and exact mapped segment;
- exact uncovered source intervals and target gaps;
- mapping orientation;
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
  canFam3 canFam4 chrUn_JH373233:1845736-1845835 \
  --json > assessment.json
```

The JSON document goes to **stdout**. Status and progress messages go to **stderr**, so
normal shell redirection does not mix progress text into the JSON file. Schema v2 carries
the same factual profile, exact candidates/evidence, resources, and provenance used by the
human renderer. It intentionally omits the legacy aggregate `verdict`, verdict-derived
`decision_reason`, and preferred-candidate field. Detailed and JSON reports include local
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
  --cache-dir /path/to/liftassess-cache
```

Cache reuse is verified: liftAssess hashes the stored artifact bytes before accepting a
complete bundle for assessment.

The cache is content-addressed by SHA-256. The local path is run context; the digest is
the exact artifact identity recorded in provenance.

## 10. Run with guaranteed zero provider access

Once the needed resource bundle is cached, use `--offline`:

```bash
assess-liftover \
  canFam3 canFam4 chr1:10000001-10000100 \
  --offline
```

`--offline` is a guarantee, not merely a preference. liftAssess will fail instead of
contacting UCSC if it cannot find a complete verified local bundle for that direction.

A normal cache-first run also avoids UCSC when a complete verified bundle is present,
but `--offline` makes that requirement explicit.

## 11. Deliberately check current provider resources

Use `--refresh` when you specifically want liftAssess to contact UCSC and reacquire the
current resource bytes instead of accepting cache-first reuse:

```bash
assess-liftover \
  canFam3 canFam4 chr1:10000001-10000100 \
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
| `--offline` | Guarantees zero provider access and requires a complete verified cache | Reproducible offline analysis or restricted network environments |
| `--refresh` | Contacts UCSC and reacquires current resources instead of cache-first reuse | Explicit freshness checks |
| `--acknowledge-ucsc-terms` | Supplies the explicit terms acknowledgement without a prompt | Non-interactive workflows after terms review |
| `--accept-transfer-plan` | Supplies the separate transfer-plan acknowledgement without a prompt | Non-interactive workflows after reviewing the planned transfer |
| `--bed PATH` | Reads BED3-or-later batch input using native 0-based, half-open coordinates; `-` reads stdin | Indexed BED batch assessment |
| `--interval-table PATH` | Reads a tab-delimited `sequence/start/end[/label]` batch table using 1-based, inclusive coordinates; `-` reads stdin | Spreadsheet-style indexed batch assessment |
| `--details` | Prints the full single-locus human-readable evidence/resource/provenance dossier; not yet available with batch input | Scientific inspection and debugging |
| `--json` | Prints schema-v2 machine-readable output | Scripts, archives, downstream analysis |
| `--quiet` | Suppresses nonessential terminal progress/status | Logs, scripts, or less terminal output |

`--offline` and `--refresh` are mutually exclusive. `--details` and `--json` are also
mutually exclusive.

## 15. Common workflows

### Normal first run

```bash
assess-liftover SOURCE_DB TARGET_DB CHR:START-END
```

Discover resources if needed, review terms and the transfer plan, cache verified
artifacts, and print the summary.

### Repeat a previous analysis using the cache

```bash
assess-liftover SOURCE_DB TARGET_DB CHR:START-END
```

If a complete verified bundle is already cached, the default path uses it without
contacting UCSC.

### Prepare a reusable chain index for repeated work

After the assembly-pair resources have been acquired and verified at least once:

```bash
prepare-liftassess-index SOURCE_DB TARGET_DB
```

This command is cache-only: it never contacts UCSC and does not require a locus. It verifies the
existing cached bundle, then parses the complete source chain once to build a reusable local index
bound to that chain's exact SHA-256 identity. Large resources can take many minutes and several GiB
of additional cache space to prepare, so this is an explicit action rather than an automatic
first-query pause. Later `assess-liftover` runs use the matching index automatically when present. A
validated index establishes the exact source-chain identity for indexed lookup, so the original
chain file is not redundantly reread merely to rehash unused bytes. Normal indexed queries verify a
compact lookup catalog, authenticate the exact genomic-bin membership/record-locator rows they use,
and verify selected compressed record blocks rather than hashing the complete SQLite lookup database
on every run; the other cached bundle artifacts retain their normal direct integrity verification.
Without a usable index, the original
full-verification/full-traversal behavior remains unchanged.

Use the same custom cache location when applicable:

```bash
prepare-liftassess-index SOURCE_DB TARGET_DB --cache-dir PATH
```

When both UCSC chain publication classes are cached, select the exact one explicitly:

```bash
prepare-liftassess-index SOURCE_DB TARGET_DB --evidence-tier COMPARATIVE
prepare-liftassess-index SOURCE_DB TARGET_DB --evidence-tier LIFTOVER-ONLY
```

If the exact chain class is not cached yet, acquire it through the normal assessment workflow
first. `--evidence-tier` disables automatic tier fallback, so this can retrieve the filtered
chain even when a complete comparative bundle is also published:

```bash
assess-liftover SOURCE_DB TARGET_DB CHR:START-END --evidence-tier LIFTOVER-ONLY
```

The usual UCSC terms and transfer-plan acknowledgement still apply.

Automatic reverse mapping requires the reverse-direction index whose publication class
matches the forward assessment.

A valid existing index is reused. `--rebuild` explicitly discards and regenerates only the derived
index; the original verified UCSC resource remains untouched.

### Assess a batch with the prepared index

BED3-or-later input keeps native **0-based, half-open** coordinates:

```bash
assess-liftover SOURCE_DB TARGET_DB --bed loci.bed
```

The BED parser ignores blank/comment/`track`/`browser` lines, preserves the optional fourth-column name as a label, and rejects zero-width or reversed intervals before assessment. Use `--bed -` to read BED from stdin.

For spreadsheet-style input, use a simple tab-delimited interval table with a required header:

```text
sequence	start	end	label
chr1	101	200	region-a
chr2	500	500	point-b
```

```bash
assess-liftover SOURCE_DB TARGET_DB --interval-table loci.tsv
```

The table's `start` and `end` are **1-based and inclusive**, matching the single-locus CLI; `start == end` is therefore a valid one-base point. The optional `label` column must be declared in the header. Table rows are normalized immediately to canonical 0-based, half-open intervals, which is the coordinate form shown in batch JSON and the current batch human report. Blank lines and `#` comments are ignored. Use `--interval-table -` to read the table from stdin.

Batch mode is deliberately cache-only and index-only. It does not contact UCSC, honor `--refresh`, build an index automatically, or fall back to a full chain traversal if the prepared index is missing or unusable. Prepare the exact selected chain class first:

```bash
prepare-liftassess-index SOURCE_DB TARGET_DB --evidence-tier COMPARATIVE
prepare-liftassess-index SOURCE_DB TARGET_DB --evidence-tier LIFTOVER-ONLY
```

Without an explicit `--evidence-tier`, batch mode prefers a complete cached COMPARATIVE bundle with a prepared all-chain index and otherwise uses the available `LIFTOVER-ONLY` chain class. `LIFTOVER-ONLY` batches assess the chain only. COMPARATIVE batches generate submitted-row candidates from the prepared all-chain index, then scan the ordinary net once and reciprocal-best chain once across the complete submitted candidate collection; the indexed all-chain is not rescanned. This preserves the single-locus net/reciprocal-best evidence semantics without a per-row whole-resource loop. The current batch COMPARATIVE scope does not run the paired filtered-vs-all-chain inventory comparison or categorical comparative relationship classifier used by single-locus results; both are reported as not assessed. Authoritative assembly-sequence name/alias preflight is not yet available, and reverse mapping is not re-run per batch row. A zero-candidate row therefore means that the selected chain index produced no candidate for the submitted label/interval; it is not an authoritative claim that the submitted sequence name is a valid assembly sequence. Exact target collisions and positive but non-identical target overlaps are reported as separate cross-record relationships derived from exact mapped target segments, never from target bounding spans. Relationship discovery uses a target-local candidate sweep rather than enumerating every pair of input rows. One-base batch rows also receive an automatic centered 101-bp point-context query through the same prepared chain index. Use `--context-bases N` to request another odd-width point window. Ordinary interval rows are not widened. Submitted-row relationships and context-scale relationships stay separate; exact equality at the context scale is reported as `NEIGHBORHOOD_LEVEL_TARGET_COLLISION`, while offset overlaps remain `OVERLAPPING_TARGET_PROJECTIONS`. Point-context candidates remain forward-chain-only even in COMPARATIVE batches; net/reciprocal-best observations are not silently promoted to the neighborhood scale. If the index cannot provide the conservative source bound needed to define a point window, that row's context is `NOT_RUN` rather than guessed.

`--json` emits schema-v2 `liftassess.ucsc_batch_result` output suitable for downstream automation. `--details` is not yet implemented for batch mode and fails explicitly rather than implying that the full batch dossier exists.

### Require offline reproducibility

```bash
assess-liftover SOURCE_DB TARGET_DB CHR:START-END --offline
```

Fail rather than use the network.

### Save the complete machine-readable report

```bash
assess-liftover \
  SOURCE_DB TARGET_DB CHR:START-END \
  --json > assessment.json
```

### Inspect the complete factual profile and evidence

```bash
assess-liftover \
  SOURCE_DB TARGET_DB CHR:START-END \
  --details
```

## 16. Common mistakes

### Swapping source and target

The first assembly is where the input coordinates currently live. The second is where
candidate mappings are assessed. Reverse the arguments and you are asking a different
question.

### Using the wrong coordinate convention

The CLI expects 1-based, inclusive coordinates. If your interval is already 0-based,
half-open, convert it before putting it on the command line. For the same physical
interval, add one to the 0-based start and keep the half-open end as the inclusive
end.

### Treating `COMPARATIVE` as "high confidence"

It is only an evidence-availability tier. Read the factual mapping profile and exact
evidence separately.

### Treating one complete projection as biological truth

A complete chain projection means every requested source base is represented in that
chain relationship. It is not an orthology call, identity check, uniqueness claim, or
proof that the locus is biologically correct.

### Assuming candidate order is rank

Current candidate order is retained for reproducibility. liftAssess does not yet define
or emit candidate-rank evidence.

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
into a one-time explicit preparation step; later candidate lookup is region-addressable while the
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
- define candidate-rank evidence yet; or
- attach reverse evidence across batches yet.

See [`FEATURES.md`](FEATURES.md) for the complete implemented/non-implemented catalog.

## 19. Expert Python use

The CLI is the easiest entry point. The package also exposes lower-level Python APIs
for callers that already have normalized candidates, local chain/net resources, or a
verified cached resource bundle.

The main boundaries are:

- `build_result_profile()` — derive the factual profile from normalized candidates and evidence;
- `build_ucsc_candidates_from_files()` — build candidates from explicit local UCSC
  resources plus provenance;
- `assess_ucsc_cached_bundle()` — assess a verified cached bundle;
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
