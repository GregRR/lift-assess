from __future__ import annotations

import gzip
import hashlib
import json
from collections.abc import Callable, Mapping
from datetime import UTC, datetime
from email.message import Message
from io import BytesIO
from pathlib import Path
from types import TracebackType
from typing import Self
from urllib.error import HTTPError, URLError
from urllib.request import Request

import pytest

from liftassess import (
    ResourceChecksumAlgorithm,
    ResourceChecksumMismatchError,
    UCSCResourceAcquisitionError,
    UCSCResourceClass,
    UCSCResourceTermsAcknowledgementRequired,
    acquire_ucsc_resource,
    iter_chain_file,
    provenance_source_for_file,
    ucsc_resource_terms,
)
from liftassess.resource_cache import _acquire_ucsc_resource


class _Response(BytesIO):
    def __init__(self, data: bytes, *, content_length: int | None = None) -> None:
        super().__init__(data)
        self._content_length = len(data) if content_length is None else content_length

    def getheader(self, name: str, default: str | None = None) -> str | None:
        if name.casefold() == "content-length":
            return str(self._content_length)
        return default

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        self.close()


def _opener(
    payloads: Mapping[str, bytes],
    *,
    calls: list[str] | None = None,
) -> Callable[[Request], _Response]:
    def open_url(request: Request) -> _Response:
        url = request.full_url
        if calls is not None:
            calls.append(url)
        if url not in payloads:
            raise HTTPError(url, 404, "not found", hdrs=Message(), fp=None)
        return _Response(payloads[url])

    return open_url


def _fixed_now() -> datetime:
    return datetime(2026, 8, 13, 23, 30, tzinfo=UTC)


def _md5_line(data: bytes, filename: str) -> bytes:
    digest = hashlib.md5(data, usedforsecurity=False).hexdigest()
    return f"{digest}  {filename}\n".encode()


def test_liftover_terms_are_distinguished_from_comparative_terms() -> None:
    liftover = ucsc_resource_terms(
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
        "canFam3ToCanFam4.over.chain.gz"
    )
    comparative = ucsc_resource_terms(
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )

    assert liftover.resource_class is UCSCResourceClass.LIFTOVER_CHAIN
    assert liftover.restricted_liftover_chain is True
    assert liftover.directory_terms_url.endswith("/canFam3/liftOver/")
    assert comparative.resource_class is UCSCResourceClass.COMPARATIVE
    assert comparative.restricted_liftover_chain is False
    assert comparative.directory_terms_url.endswith("/canFam3/vsCanFam4/")

    reciprocal = ucsc_resource_terms(
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/"
        "reciprocalBest/canFam3.canFam4.rbest.chain.gz"
    )
    assert reciprocal.resource_class is UCSCResourceClass.COMPARATIVE
    assert reciprocal.directory_terms_url.endswith("/canFam4/vsCanFam3/")


def test_terms_acknowledgement_is_required_before_network_access(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
        "canFam3ToCanFam4.over.chain.gz"
    )

    def forbidden(_: Request) -> _Response:
        raise AssertionError("network must not be touched before acknowledgement")

    with pytest.raises(
        UCSCResourceTermsAcknowledgementRequired,
        match="EULA acceptance",
    ):
        _acquire_ucsc_resource(
            url,
            tmp_path,
            terms_acknowledged=False,
            open_url=forbidden,
            now=_fixed_now,
        )


def test_download_verifies_provider_md5_and_writes_content_addressed_cache(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"realistic compressed-resource bytes"
    expected_sha256 = hashlib.sha256(data).hexdigest()

    result = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=_opener(
            {
                checksum_url: _md5_line(data, "canFam3.canFam4.net.gz"),
                url: data,
            }
        ),
        now=_fixed_now,
    )

    assert result.cache_hit is False
    assert result.path == (
        tmp_path / "artifacts" / "sha256" / expected_sha256[:2] / expected_sha256
    )
    assert result.path.read_bytes() == data
    assert result.sha256 == f"sha256:{expected_sha256}"
    assert result.size_bytes == len(data)
    assert result.retrieved_at == "2026-08-13T23:30:00Z"
    assert result.provider_checksum is not None
    assert result.provider_checksum.algorithm is ResourceChecksumAlgorithm.MD5
    assert result.provider_checksum.source_url == checksum_url

    index_files = list((tmp_path / "by-url").glob("*.json"))
    assert len(index_files) == 1
    index = json.loads(index_files[0].read_text(encoding="utf-8"))
    assert index["source_url"] == url
    assert index["sha256"] == expected_sha256
    assert index["size_bytes"] == len(data)
    assert index["terms"]["directory_terms_url"].endswith("/canFam3/vsCanFam4/")


def test_verified_cache_hit_is_offline_and_retains_original_provider_md5(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"cached resource"
    payloads = {
        checksum_url: _md5_line(data, "canFam3.canFam4.net.gz"),
        url: data,
    }
    first_calls: list[str] = []
    first = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=_opener(payloads, calls=first_calls),
        now=_fixed_now,
    )
    assert first_calls == [checksum_url, url]

    def offline(_: Request) -> _Response:
        raise AssertionError("verified cache hit must not require network access")

    second = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=offline,
        now=_fixed_now,
    )

    assert second.cache_hit is True
    assert second.path == first.path
    assert second.retrieved_at == first.retrieved_at
    assert second.provider_checksum == first.provider_checksum


def test_same_bytes_from_two_urls_converge_on_one_artifact(tmp_path: Path) -> None:
    base = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    first_url = f"{base}first.net.gz"
    second_url = f"{base}second.net.gz"
    checksum_url = f"{base}md5sum.txt"
    data = b"same bytes from two provider URLs"
    md5 = hashlib.md5(data, usedforsecurity=False).hexdigest()
    md5sum = f"{md5}  first.net.gz\n{md5}  second.net.gz\n".encode()
    opener = _opener(
        {
            checksum_url: md5sum,
            first_url: data,
            second_url: data,
        }
    )

    first = _acquire_ucsc_resource(
        first_url,
        tmp_path,
        terms_acknowledged=True,
        open_url=opener,
        now=_fixed_now,
    )
    second = _acquire_ucsc_resource(
        second_url,
        tmp_path,
        terms_acknowledged=True,
        open_url=opener,
        now=_fixed_now,
    )

    assert first.path == second.path
    artifact_files = [
        path
        for path in (tmp_path / "artifacts" / "sha256").rglob("*")
        if path.is_file()
    ]
    assert len(artifact_files) == 1
    assert len(list((tmp_path / "by-url").glob("*.json"))) == 2



def test_cached_extensionless_gzip_artifact_integrates_with_file_provenance_and_parser(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.all.chain.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    chain_text = (
        b"chain 100 chr1 1000 + 10 20 chrA 1000 + 30 40 1\n"
        b"10\n\n"
    )
    data = gzip.compress(chain_text, mtime=0)

    acquired = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=_opener(
            {
                checksum_url: _md5_line(data, "canFam3.canFam4.all.chain.gz"),
                url: data,
            }
        ),
        now=_fixed_now,
    )
    provenance = provenance_source_for_file(
        acquired.path,
        label="cached comparative chain",
        derived_from=(),
    )
    records = tuple(iter_chain_file(acquired.path))

    assert acquired.path.suffix == ""
    assert provenance.identifiers[0].value == acquired.sha256
    assert len(records) == 1
    assert records[0].chain_id == 1


def test_provider_md5_mismatch_does_not_publish_partial_artifact(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )

    with pytest.raises(ResourceChecksumMismatchError, match="checksum mismatch"):
        _acquire_ucsc_resource(
            url,
            tmp_path,
            terms_acknowledged=True,
            open_url=_opener(
                {
                    checksum_url: b"0" * 32 + b"  canFam3.canFam4.net.gz\n",
                    url: b"different bytes",
                }
            ),
            now=_fixed_now,
        )

    artifacts = tmp_path / "artifacts"
    indexes = tmp_path / "by-url"
    assert not artifacts.exists() or not any(path.is_file() for path in artifacts.rglob("*"))
    assert not indexes.exists() or not list(indexes.glob("*.json"))
    assert not list((tmp_path / "tmp").glob("*.part"))


def test_missing_md5sum_is_allowed_but_transport_failure_is_not(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"resource without provider checksum"

    result = _acquire_ucsc_resource(
        url,
        tmp_path / "missing",
        terms_acknowledged=True,
        open_url=_opener({url: data}),
        now=_fixed_now,
    )
    assert result.provider_checksum is None

    def failing(request: Request) -> _Response:
        if request.full_url == checksum_url:
            raise URLError("simulated metadata transport failure")
        return _Response(data)

    with pytest.raises(UCSCResourceAcquisitionError, match="checksum metadata"):
        _acquire_ucsc_resource(
            url,
            tmp_path / "failure",
            terms_acknowledged=True,
            open_url=failing,
            now=_fixed_now,
        )


def test_non_object_cache_index_is_treated_as_cache_miss(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"resource bytes"
    calls: list[str] = []
    opener = _opener(
        {
            checksum_url: _md5_line(data, "canFam3.canFam4.net.gz"),
            url: data,
        },
        calls=calls,
    )

    first = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=opener,
        now=_fixed_now,
    )
    index_path = next((tmp_path / "by-url").glob("*.json"))
    index_path.write_text('["not", "an", "object"]', encoding="utf-8")
    calls.clear()

    second = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=opener,
        now=_fixed_now,
    )

    assert second.cache_hit is False
    assert second.path == first.path
    assert calls == [checksum_url, url]


def test_corrupt_cached_artifact_is_replaced_from_provider(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"correct resource bytes"
    opener = _opener(
        {
            checksum_url: _md5_line(data, "canFam3.canFam4.net.gz"),
            url: data,
        }
    )
    first = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=opener,
        now=_fixed_now,
    )
    first.path.write_bytes(b"corrupt cache bytes")

    second = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=opener,
        now=_fixed_now,
    )

    assert second.cache_hit is False
    assert second.path.read_bytes() == data
    assert second.sha256 == first.sha256



def test_existing_md5sum_without_exact_filename_entry_is_checksum_unavailable(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"net bytes whose provider directory omits an exact MD5 entry"

    result = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=_opener(
            {
                checksum_url: _md5_line(
                    b"different file", "canFam3.canFam4.all.chain.gz"
                ),
                url: data,
            }
        ),
        now=_fixed_now,
    )

    assert result.provider_checksum is None
    assert result.path.read_bytes() == data


def test_refresh_forces_new_transfer_when_provider_checksum_is_unavailable(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    first_data = b"first provider bytes"
    first = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=_opener({checksum_url: b"", url: first_data}),
        now=_fixed_now,
    )

    second_data = b"updated provider bytes"
    second = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        refresh=True,
        open_url=_opener({checksum_url: b"", url: second_data}),
        now=_fixed_now,
    )

    assert first.cache_hit is False
    assert second.cache_hit is False
    assert first.path != second.path
    assert second.path.read_bytes() == second_data



def test_truncated_download_fails_when_content_length_is_known(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"truncated bytes"

    def truncated(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(b"")
        return _Response(data, content_length=len(data) + 50)

    with pytest.raises(UCSCResourceAcquisitionError, match="incomplete UCSC resource"):
        _acquire_ucsc_resource(
            url,
            tmp_path,
            terms_acknowledged=True,
            open_url=truncated,
            now=_fixed_now,
        )

    assert not list((tmp_path / "tmp").glob("*.part"))
    assert not (tmp_path / "by-url").exists()


def test_download_transport_failure_cleans_partial_file(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )

    def failing(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(b"")
        raise URLError("simulated download failure")

    with pytest.raises(UCSCResourceAcquisitionError, match="failed to download"):
        _acquire_ucsc_resource(
            url,
            tmp_path,
            terms_acknowledged=True,
            open_url=failing,
            now=_fixed_now,
        )

    assert not list((tmp_path / "tmp").glob("*.part"))
    assert not (tmp_path / "by-url").exists()


def test_public_acquisition_signature_requires_explicit_acknowledgement(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )

    with pytest.raises(UCSCResourceTermsAcknowledgementRequired):
        acquire_ucsc_resource(url, tmp_path, terms_acknowledged=False)
