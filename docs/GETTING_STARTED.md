# Getting Started with liftAssess

This guide is for someone who wants to run liftAssess without first learning the
internal architecture of UCSC chain/net files.

For the complete implemented capability list, see [`FEATURES.md`](FEATURES.md). For the
scientific rules behind the verdicts, see [`DESIGN.md`](DESIGN.md).

## 1. What liftAssess does

A normal liftOver operation asks:

> Where can this source interval map in the target assembly?

liftAssess asks a different follow-up question:

> What does the available evidence say about those candidate mappings?

It can report one of three results:

- `WELL_SUPPORTED` — the available informative evidence favors one candidate without
  material contradiction under the v1 rules;
- `CONTESTED` — more than one candidate remains materially supported, or the available
  evidence materially disagrees; or
- `INDETERMINATE` — the available evidence is not sufficient or discriminating enough
  for either of the other conclusions.

`WELL_SUPPORTED` does **not** mean "biologically proven correct." liftAssess assesses
support in the available evidence; it is not a truth oracle and it does not output a
numeric confidence score.

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

The current CLI assesses one locus at a time. It does not yet accept a BED file or a
batch of loci.

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
   the deterministic assessor.
8. **Print the report.** The default is the concise human-readable summary.

Cancelling either acknowledgement stops before the corresponding network/resource
operation.

## 5. How to read the default result

A concise report has this general shape:

```text
Source locus: chr1:10000001-10000100 (1-based inclusive)
Evidence availability: COMPARATIVE — mapping plus comparative evidence available
Assessment: WELL SUPPORTED
Preferred candidate: chr1:10027740-10027839 (1-based inclusive; same orientation)
Why: ...
...
This does not establish biological correctness.
```

The lines answer different questions.

### `Source locus`

Confirms exactly what interval was assessed and states the display coordinate
convention.

### `Evidence availability`

Tells you what kind of evidence liftAssess could examine for this assembly pair.

- `COMPARATIVE` means the full comparative resource bundle was available.
- `LIFTOVER-ONLY` means only chain mapping evidence was available.

This is **not a confidence rating**. A `COMPARATIVE` run can still be
`INDETERMINATE`, and a `LIFTOVER-ONLY` run can sometimes be `WELL_SUPPORTED`.

### `Assessment`

The final categorical verdict: `WELL_SUPPORTED`, `CONTESTED`, or `INDETERMINATE`.

### `Preferred candidate`

Appears only when the assessor's deterministic rules support one candidate. If the
rules do not support a preferred candidate, the report gives candidate context instead
of pretending that candidate order is a ranking.

### `Why`

A plain-language rendering of the exact assessor `decision_reason`. The detailed and
JSON reports expose the machine-readable reason code itself.

### Final caveat

Every report retains the reminder that evidence support does not establish biological
correctness.

## 6. The three verdicts in practical terms

### `WELL_SUPPORTED`

Read this as:

> Under the evidence that was available and the v1 categorical rules, one candidate is
> favored without material contradiction.

Do **not** read it as:

> liftAssess proved this is the biologically correct locus.

### `CONTESTED`

Read this as:

> There is a real reason not to collapse the result to one unqualified mapping.

Typical v1 causes include multiple material candidates or a full mapping that conflicts
with reciprocal-best evidence.

### `INDETERMINATE`

Read this as:

> The evidence does not justify either a well-supported or contested conclusion under
> the current rules.

Examples include no candidate, an incomplete single mapping in liftOver-only evidence,
or comparative evidence that is mixed but not sufficient for a stronger categorical
conclusion.

## 7. Ask for the full human-readable dossier

Use `--details` when you need to understand or audit the evidence behind the summary:

```bash
assess-liftover \
  canFam3 canFam4 chrUn_JH373233:1845736-1845835 \
  --details
```

The detailed report includes:

- the exact `decision_reason` code;
- every candidate and exact mapped segment;
- mapping orientation;
- every evidence observation;
- whether each observation is supporting, contradicting, both, or context;
- chain/net/reciprocal-best detail;
- resource URLs, cache paths, retrieval metadata, and checksums;
- which cached resources were actually consumed by the engine; and
- the provenance dependency graph showing which observations share upstream sources.

Use this mode when a concise verdict is not enough for scientific review.

## 8. Get JSON for scripts and pipelines

Use `--json` for the complete schema-versioned machine-readable report:

```bash
assess-liftover \
  canFam3 canFam4 chrUn_JH373233:1845736-1845835 \
  --json > assessment.json
```

The JSON document goes to **stdout**. Status and progress messages go to **stderr**, so
normal shell redirection does not mix progress text into the JSON file. Detailed and
JSON reports include local cache paths as run context, so inspect that metadata before
publishing a report if local filesystem paths are information you do not want to share.

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
| `--details` | Prints the full human-readable evidence/resource/provenance dossier | Scientific inspection and debugging |
| `--json` | Prints schema-v1 machine-readable output | Scripts, archives, downstream analysis |
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

### Inspect why a result was contested or indeterminate

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

It is only an evidence-availability tier. Always read the verdict separately.

### Treating a preferred candidate as biological truth

A preferred candidate means the v1 evidence rules support that mapping over the
alternatives that remained material. It is not an orthology call or proof that the
locus is biologically correct.

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
an estimate of how close the assessment logic is to a verdict.

## 17. Current performance warning

Single-locus assessment currently streams and parses large comparative resources rather
than using a prebuilt genomic index. Measured work on the `canFam3` → `canFam4` pair
showed that large-resource parsing can dominate runtime even for a small, ordinary
locus.

That means:

- a small genomic interval does not necessarily mean a fast run;
- cached data avoids repeat downloading but does not yet avoid resource parsing; and
- `--offline` changes provider access, not the core parsing cost.

See [`PERFORMANCE.md`](PERFORMANCE.md) for measured benchmark results and the planned
optimization direction.

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
- optimize multi-locus batch assessment yet.

See [`FEATURES.md`](FEATURES.md) for the complete implemented/non-implemented catalog.

## 19. Expert Python use

The CLI is the easiest entry point. The package also exposes lower-level Python APIs
for callers that already have normalized candidates, local chain/net resources, or a
verified cached resource bundle.

The main boundaries are:

- `assess_candidates()` — assess normalized candidates;
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
