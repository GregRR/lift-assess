"""Offline clean-install smoke test for liftAssess release artifacts."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import cast

import liftassess
from liftassess.resource_cache import _write_url_index, ucsc_resource_terms

_SOURCE_DB = "canFam3"
_TARGET_DB = "canFam4"
_CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
    "canFam3ToCanFam4.over.chain.gz"
)
_CHROM_INFO_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/chromInfo.txt.gz"
)


def _gzip_bytes(text: str) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as handle:
        handle.write(text.encode("utf-8"))
    return buffer.getvalue()


def _publish_cached_resource(
    cache_root: Path,
    *,
    source_url: str,
    data: bytes,
) -> None:
    digest = hashlib.sha256(data).hexdigest()
    artifact = cache_root / "artifacts" / "sha256" / digest[:2] / digest
    artifact.parent.mkdir(parents=True, exist_ok=True)
    artifact.write_bytes(data)

    index_key = hashlib.sha256(source_url.encode("utf-8")).hexdigest()
    _write_url_index(
        cache_root / "by-url" / f"{index_key}.json",
        source_url=source_url,
        retrieved_at="2026-08-29T00:00:00Z",
        sha256=digest,
        size_bytes=len(data),
        provider_checksum=None,
        terms=ucsc_resource_terms(source_url),
    )


def _run_installed_cli(cache_root: Path) -> dict[str, object]:
    executable = Path(sys.executable).with_name("assess-liftover")
    if not executable.is_file():
        raise RuntimeError(f"installed CLI entry point is missing: {executable}")

    subprocess.run(
        [str(executable), "--help"],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    completed = subprocess.run(
        [
            str(executable),
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(cache_root),
            "--offline",
            "--evidence-tier",
            "LIFTOVER-ONLY",
            "--json",
            "--quiet",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    raw_payload: object = json.loads(completed.stdout)
    if not isinstance(raw_payload, dict):
        raise TypeError("installed CLI JSON output is not an object")
    return cast(dict[str, object], raw_payload)


def _assert_release_payload(payload: dict[str, object]) -> None:
    if payload.get("schema_version") != 2:
        raise RuntimeError("installed CLI did not emit schema-v2 JSON")
    if "verdict" in payload:
        raise RuntimeError("schema-v2 JSON unexpectedly contains legacy verdict")

    source = payload.get("source_assembly")
    target = payload.get("target_assembly")
    profile = payload.get("result_profile")
    if not isinstance(source, dict) or source.get("name") != _SOURCE_DB:
        raise RuntimeError("unexpected source assembly in release-smoke JSON")
    if not isinstance(target, dict) or target.get("name") != _TARGET_DB:
        raise RuntimeError("unexpected target assembly in release-smoke JSON")
    if not isinstance(profile, dict):
        raise TypeError("release-smoke JSON is missing result_profile")
    if profile.get("headline") != "ONE_COMPLETE_CHAIN_PROJECTION":
        raise RuntimeError("release-smoke JSON has unexpected factual headline")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected-version", required=True)
    args = parser.parse_args()
    expected_version = cast(str, args.expected_version)

    if liftassess.__version__ != expected_version:
        raise RuntimeError(
            "installed package version does not match release version: "
            f"{liftassess.__version__!r} != {expected_version!r}"
        )

    with tempfile.TemporaryDirectory(prefix="liftassess-release-smoke-") as temp_dir:
        cache_root = Path(temp_dir) / "cache"
        _publish_cached_resource(
            cache_root,
            source_url=_CHROM_INFO_URL,
            data=_gzip_bytes("chr1\t1000\t/gbdb/canFam3/canFam3.2bit\n"),
        )
        _publish_cached_resource(
            cache_root,
            source_url=_CHAIN_URL,
            data=_gzip_bytes(
                "chain 100 chr1 1000 + 100 120 chrA 2000 + 500 520 1\n20\n\n"
            ),
        )
        payload = _run_installed_cli(cache_root)
        _assert_release_payload(payload)

    print(
        "release smoke passed: "
        f"liftassess {liftassess.__version__}, CLI entry point, schema-v2 JSON"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
