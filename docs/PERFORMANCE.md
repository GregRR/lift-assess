# liftAssess Performance

This page records an initial performance characterization of the current streaming UCSC
implementation. It is intended to make present-day runtime expectations and known bottlenecks
visible; it is **not** a performance guarantee or a claim that these measurements generalize to
all assembly pairs, loci, operating systems, or Python versions.

## Benchmark scope

Measurements below were collected on 2026-08-17 from Git commit
`e1e0472c059bd4ec214793c508c97dd4636de99f` using:

- Apple Mac mini with Apple M4, 10 CPU cores (4 performance + 6 efficiency), and 16 GB RAM;
- macOS 26.5.2 (25F84), arm64;
- Python 3.14.7;
- uv 0.12.5.

The published environment intentionally omits machine-specific identifiers such as serial numbers,
hardware UUIDs, local usernames, and cache paths.

All timed runs used already-acquired, verified UCSC resources and `--offline`, so no provider
network access or resource download was included. The timed command **did** include liftAssess's
normal local cache verification before assessment. Wall/user/system time and maximum resident set
size were recorded with macOS `/usr/bin/time -l`:

```text
/usr/bin/time -l \
  uv run assess-liftover SOURCE_DB TARGET_DB LOCUS --offline \
  2>&1 | tee benchmark.txt
```

Because stderr was piped through `tee`, these timed runs were non-TTY runs and liftAssess's
interactive progress bars were suppressed as designed.

## Measured runs

CLI loci below use the documented 1-based, inclusive display convention. The
`canFam6` → `mm39` pair is intentionally used here only as a computational performance probe
with a much smaller comparative resource set. It is cross-species and therefore outside
liftAssess's same-species scientific operational envelope. These benchmarks predate the
Milestone-17 result redesign, so the table preserves the legacy alpha result labels only as
historical run metadata; they are not current result fields and are not used as scientific
validation.

| Assembly pair | Source locus | Historical alpha result | Candidate context | All-chain size | Real time | User time | Max RSS |
|---|---|---|---:|---:|---:|---:|---:|
| `canFam3` → `canFam4` | `chrUn_JH373233:1845736-1845835` | `CONTESTED` | 170 | 2.47 GiB | 641.17 s | 639.14 s | 47.8 MiB |
| `canFam3` → `canFam4` | `chrUn_JH373233:1845736-1845835` | `CONTESTED` | 170 | 2.47 GiB | 625.71 s | 624.12 s | 46.4 MiB |
| `canFam3` → `canFam4` | `chr1:10000001-10000100` | `WELL_SUPPORTED` | preferred candidate reported | 2.47 GiB | 640.74 s | 638.09 s | 51.6 MiB |
| `canFam6` → `mm39` | `chr1:10000001-10000100` | `INDETERMINATE` | 0 | 157.0 MiB | 42.50 s | 42.14 s | 170.1 MiB |
| `canFam6` → `mm39` | `chr6:43585033-43585132` | `WELL_SUPPORTED` | preferred candidate reported | 157.0 MiB | 103.12 s | 102.39 s | 178.4 MiB |

The two repeated `canFam3` → `canFam4` mechanical-fixture runs differed by 2.4%; their mean wall
time was 633.44 seconds (10:33.44).

The `canFam3` → `canFam4` all-chain resource is exactly 2,652,632,416 compressed bytes in the
frozen mechanical fixture. The `canFam6` → `mm39` transfer plan advertised a 157.0 MiB all-chain
resource. Comparing the simple `canFam3` → `canFam4` locus with the zero-candidate
`canFam6` → `mm39` locus, the chain was about 16 times smaller and the complete offline command
was about 15 times faster. This is consistent with chain-resource size being a major runtime
factor for the current implementation; five measurements are not enough to establish an exact
linear scaling law.

The two `canFam6` → `mm39` runs also show that resource size is not the whole story. Once a mapping
candidate existed and comparative evidence was consumed, wall time rose from 42.50 to 103.12
seconds on the same assembly pair. Candidate-dependent comparative processing can therefore be material even when the chain
resource is much smaller.

Maximum RSS did not scale monotonically with chain size in these runs. The cross-pair memory
difference was not investigated further and should not be treated as an established memory-scaling
model.

## CPU profile

A separate `cProfile` run used the real `canFam3` → `canFam4` mechanical-fixture locus. Profiling
instrumentation substantially increased runtime, so that run is used only to identify hotspots,
not as a wall-time benchmark.

The profile recorded more than 13.2 billion function calls. The dominant cumulative path was
UCSC chain iteration/parsing. Representative counts included:

- about 2.43 billion calls to the chain integer parser;
- about 752.6 million chain-block constructions;
- about 752.6 million chain-block validation/post-initialization calls;
- about 786.9 million string `split` calls;
- about 853.9 million string `strip` calls.

Direct gzip decompression time was small relative to the Python parsing/object-construction path in
that profile. Net-evidence attachment was also much smaller than chain parsing in the profiled
mechanical-fixture run. These observations identify Python-level whole-chain parsing as the primary
hotspot for that workload, not storage transfer or gzip decompression alone.

In the unprofiled timed runs, process CPU time was approximately equal to wall time. For example,
the first 641.17-second run used 639.14 seconds of user CPU plus 1.65 seconds of system CPU. That
is consistent with this workload being CPU-bound and effectively using about one CPU core during
the dominant work; it is not a claim about every liftAssess workload.

## Current interpretation

The measurements support three scoped conclusions about the implementation at this commit:

1. **Whole-chain traversal is a major single-locus cost.** A 170-candidate contested locus and a
   simple well-supported locus on the same large assembly pair took essentially the same time.
2. **Chain-resource size materially affects runtime.** Moving from the 2.47 GiB
   `canFam3` → `canFam4` chain to the 157.0 MiB `canFam6` → `mm39` chain reduced a
   zero-candidate/simple-scan workload by roughly an order of magnitude.
3. **Comparative evidence still adds candidate-dependent work.** On the smaller assembly pair,
   the mapped locus took materially longer than the zero-candidate locus.

These results strengthen an architectural requirement already present in
[`DESIGN.md`](DESIGN.md): scalable assessment should not multiply full comparative-resource
parsing by locus count. An indexed/preprocessed local representation or shared traversal across
many loci may improve both single-locus latency and future batch throughput, provided exact
resource identity, provenance, coordinate semantics, and assessment behavior remain unchanged.

Optimization should be measured in this order rather than assuming that parallelism is the first
answer:

1. avoid unnecessary whole-resource parsing per query where a scientifically equivalent indexed,
   preprocessed, or shared-traversal approach can be demonstrated;
2. profile and optimize the hot chain parser/object-construction path;
3. profile candidate-dependent comparative evidence lookup after the structural improvements;
4. evaluate safe multicore execution if meaningful CPU-bound work remains and deterministic
   results/provenance can be preserved.

See [`ROADMAP.md`](ROADMAP.md) for planned performance and batch work.

## Informal historical hardware comparison

An earlier development run of the same `canFam3` → `canFam4` mechanical-fixture locus on an older
iMac recorded 2,308.36 seconds (38:28.36) of wall time, versus a 633.44-second mean for the two M4
runs above, or about 3.6 times longer.

That older run was **not** captured under the same frozen benchmark protocol: its exact commit,
Python/tool environment, machine specification, and cached resource identities were not recorded
together with the timing. It is therefore an informal development observation, not a controlled
cross-hardware benchmark and not a performance claim about liftAssess releases.

## Milestone 18 chain-index prototype results

A second performance investigation on 2026-08-20 used the same M4 Mac mini and exact
2.47-GiB `canFam3` → `canFam4` all-chain resource to select the reusable chain-access
architecture for Milestone 18. These measurements are engineering benchmarks, not scientific
validation. Every indexed/selective result below was required to reproduce the exact candidate
tuple from the full parser for both the 170-candidate mechanical fixture and the ordinary
single-candidate control.

| Access path | 170-candidate locus | 1-candidate locus | Notes |
|---|---:|---:|---|
| Full traversal | 626.376 s | 632.612 s | Current parser baseline |
| Selective materialization after header scan | 348.127 s | 359.255 s | Still scans the full 2.47-GiB resource |
| Per-source-sequence shard | 10.339 s | 7.774 s | Reuses a one-time sequence partition |
| 65,536-bp bin index over sequence shards | 0.309 s | 0.101 s | 948 and 1 chain records selected, respectively |
| Single-copy blocked store over the bin index | 0.350 s | 0.104 s | Same lookup semantics with substantially smaller derived storage |
| Full-resource blocked/bin index | 1.353 s | 0.106 s | 948 and 1 chain records selected; exact candidate equivalence |

The full-resource prototype indexed **33,486,862 chain records** into **35,067,017 bin
memberships** and **8,346 independently compressed blocks**. The index itself took 1,234.746
seconds (20.6 minutes) of active build time and occupied 3.844 GiB (1.365 GiB SQLite metadata plus
2.478 GiB compressed record blocks).

The enclosing `/usr/bin/time` process initially reported 2,242.47 seconds (37.4 minutes), but that
wall-clock figure is **not** a valid compute-time measurement. A later `pmset -g log` check showed
that the M4 entered sleep twice during the benchmark, for 318 seconds and 692 seconds respectively
(1,010 seconds total). That almost exactly explains the 1,007.724-second difference between the
external wall clock and the script's build timer. The 37.4-minute value should therefore not be used
as an index-build performance result.

Follow-up M4 measurements also ruled out cache-integrity hashing as the source of that gap. SHA-256
verification of the 1.365-GiB index database took 1.27 seconds immediately after `sudo purge` and
0.60 seconds on an immediate warm rerun. Loading/verifying the complete cached `canFam3` →
`canFam4` COMPARATIVE bundle took 1.232 seconds after `sudo purge` and 1.192 seconds on a warm
rerun on that M4/SSD system.

A later cold-cache check on the older development iMac, whose cache resides on a comparatively slow
HDD, measured **95.976 seconds** inside `load_cached_ucsc_resource_bundle()` for the same complete
COMPARATIVE bundle after `sudo purge`. That measurement is hardware/storage-specific, not a general
liftAssess runtime, but it demonstrates that unconditional multi-gigabyte rehashing can dominate an
indexed query on slower storage even when cryptographic computation itself is cheap. The indexed
production path therefore keeps the integrity contract while avoiding a redundant reread of the
original chain: index build/rebuild verifies the full source chain; subsequent indexed assessment
verifies a compact lookup-integrity catalog, authenticates the membership/record-locator rows for
the exact queried genomic bins, and validates selected compressed blocks; the other cached bundle
artifacts continue to receive their normal direct verification. The full SQLite database SHA-256 is
retained for explicit deep verification but is not a prerequisite for every indexed query. This
prevents the integrity mechanism from merely moving the slow-storage bottleneck from the original
2.53-GiB chain to the roughly 1.365-GiB lookup database.

A targeted two-sequence blocked-store experiment reduced the corresponding derived representation
from about 732 MiB to about 311 MiB while retaining 0.10–0.35-second queries. A competing design
that duplicated complete record payloads into compressed genomic-bin frames had less predictable
storage behavior, especially for long chains crossing many bins, and was not selected.

A focused iMac follow-up measured the still-exhaustive comparative evidence parsers after chain
indexing had made candidate lookup cheap. Parsing the 10.03-MiB ordinary net produced 847,501
records in **8.174 seconds**; parsing the 5.15-MiB reciprocal-best chain produced 7,570 records in
**4.841 seconds**, for **13.015 seconds combined**. These are parser timings on the older iMac, not
M4-comparable end-to-end timings. They establish that comparative scans are material once chain
lookup is indexed, but they do not by themselves select another persistent index format. Future
reverse, point-context, comparative-expansion, and batch work must share or index this traversal
rather than multiplying the same whole-resource scans per derived query.

**Measured architecture decision:** use 65,536-bp source-coordinate bin memberships for spatial
selection and store each serialized chain record exactly once in encounter-order compressed blocks.
Bind the derived artifact to the canonical SHA-256 of the original chain, preserve original record
order for reproducibility, and keep the existing parser/projection logic authoritative for selected
records. The production path may reuse a validated index when present and fall back to the original
verified traversal otherwise. Initial index construction is deliberately explicit through
`prepare-liftassess-index SOURCE_DB TARGET_DB`: the command uses only an already verified local
bundle and never turns an ordinary first query into an implicit multi-minute preprocessing job.

For future long-running macOS benchmarks, run the timed command under `caffeinate` (unless sleep is
itself part of the experiment) so system idle sleep cannot silently inflate wall-clock timing.
