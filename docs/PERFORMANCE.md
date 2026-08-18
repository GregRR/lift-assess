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
liftAssess v1's same-species scientific operational envelope; its verdicts below are not used
as scientific validation.

| Assembly pair | Source locus | Result | Candidate context | All-chain size | Real time | User time | Max RSS |
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
