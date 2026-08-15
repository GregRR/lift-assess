from __future__ import annotations

import gzip
import hashlib
import json
from collections.abc import Callable, Mapping
from datetime import UTC, datetime
from email.message import Message
from io import BytesIO
from os import PathLike
from pathlib import Path
from threading import Event, Thread
from types import TracebackType
from typing import Self
from urllib.error import HTTPError, URLError
from urllib.request import Request

import pytest

from liftassess import (
    CachedResource,
    CachedUCSCResourceBundle,
    EvidenceAvailabilityTier,
    ResourceChecksumAlgorithm,
    ResourceChecksumMismatchError,
    UCSCBundleAcquisitionItem,
    UCSCBundleAcquisitionPlan,
    UCSCBundleAcquisitionPlanAcknowledgementRequired,
    UCSCBundleResourceRole,
    UCSCBundleTransferInspection,
    UCSCBundleTransferInspectionItem,
    UCSCRemoteResourceMetadata,
    UCSCResourceAcquisitionError,
    UCSCResourceBundle,
    UCSCResourceClass,
    UCSCResourceTermsAcknowledgementRequired,
    acquire_ucsc_resource,
    acquire_ucsc_resource_bundle,
    inspect_ucsc_bundle_transfer_plan,
    inspect_ucsc_resource,
    iter_chain_file,
    plan_ucsc_bundle_acquisition,
    provenance_source_for_file,
    ucsc_resource_terms,
)
from liftassess.resource_cache import (
    _acquire_ucsc_resource,
    _acquire_ucsc_resource_bundle,
    _inspect_ucsc_bundle_transfer_plan,
    _inspect_ucsc_resource,
)


class _Response(BytesIO):
    def __init__(
        self,
        data: bytes,
        *,
        content_length: int | None = None,
        headers: Mapping[str, str] | None = None,
        status: int = 200,
    ) -> None:
        super().__init__(data)
        self.status = status
        self._headers = {
            key.casefold(): value for key, value in (headers or {}).items()
        }
        if "content-length" not in self._headers:
            length = len(data) if content_length is None else content_length
            self._headers["content-length"] = str(length)

    def getheader(self, name: str, default: str | None = None) -> str | None:
        return self._headers.get(name.casefold(), default)

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
    assert first_calls == [checksum_url, url, url]

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
    assert calls == [checksum_url, url, url]


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


class _InterruptingResponse(_Response):
    def __init__(
        self,
        prefix: bytes,
        *,
        total_size: int,
        etag: str,
    ) -> None:
        super().__init__(
            prefix,
            headers={
                "Content-Length": str(total_size),
                "ETag": etag,
            },
        )
        self._first_read = True

    def read(self, size: int | None = -1) -> bytes:
        if self._first_read:
            self._first_read = False
            return super().read(size)
        raise URLError("simulated interrupted body")


def _resumable_head_headers(data: bytes, etag: str) -> dict[str, str]:
    return {
        "Content-Length": str(len(data)),
        "Accept-Ranges": "bytes",
        "ETag": etag,
        "Last-Modified": "Fri, 14 Aug 2026 12:00:00 GMT",
    }


def test_interrupted_resumable_download_retains_prefix_and_resumes_with_if_range(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"0123456789abcdef"
    prefix = data[:6]
    etag = '"stable-v1"'
    first_requests: list[Request] = []

    def first_open(request: Request) -> _Response:
        first_requests.append(request)
        if request.full_url == checksum_url:
            return _Response(_md5_line(data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(data, etag))
        return _InterruptingResponse(prefix, total_size=len(data), etag=etag)

    with pytest.raises(UCSCResourceAcquisitionError, match="partial state was retained"):
        _acquire_ucsc_resource(
            url,
            tmp_path,
            terms_acknowledged=True,
            open_url=first_open,
            now=_fixed_now,
        )

    partials = list((tmp_path / "partials").rglob("*.part"))
    assert len(partials) == 1
    assert partials[0].read_bytes() == prefix
    assert not (tmp_path / "by-url").exists()
    assert first_requests[-1].get_header("Range") is None
    assert first_requests[-1].get_header("Accept-encoding") == "identity"

    second_requests: list[Request] = []

    def second_open(request: Request) -> _Response:
        second_requests.append(request)
        if request.full_url == checksum_url:
            return _Response(_md5_line(data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(data, etag))
        assert request.get_header("Range") == f"bytes={len(prefix)}-"
        assert request.get_header("If-range") == etag
        tail = data[len(prefix) :]
        return _Response(
            tail,
            status=206,
            headers={
                "Content-Length": str(len(tail)),
                "Content-Range": f"bytes {len(prefix)}-{len(data) - 1}/{len(data)}",
                "ETag": etag,
            },
        )

    result = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=second_open,
        now=_fixed_now,
    )

    assert result.path.read_bytes() == data
    assert result.sha256 == f"sha256:{hashlib.sha256(data).hexdigest()}"
    assert not list((tmp_path / "partials").rglob("*.part"))
    assert second_requests[-1].get_header("Range") == f"bytes={len(prefix)}-"


def test_completed_partial_is_published_on_retry_without_another_resource_get(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"complete bytes before connection teardown"
    etag = '"stable-v1"'

    def interrupted_after_complete(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(_md5_line(data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(data, etag))
        return _InterruptingResponse(data, total_size=len(data), etag=etag)

    with pytest.raises(UCSCResourceAcquisitionError, match="partial state was retained"):
        _acquire_ucsc_resource(
            url,
            tmp_path,
            terms_acknowledged=True,
            open_url=interrupted_after_complete,
            now=_fixed_now,
        )

    partial = next((tmp_path / "partials").rglob("*.part"))
    assert partial.stat().st_size == len(data)
    retry_requests: list[Request] = []

    def retry(request: Request) -> _Response:
        retry_requests.append(request)
        if request.full_url == checksum_url:
            return _Response(_md5_line(data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(data, etag))
        raise AssertionError("a complete verified partial must not GET the resource again")

    result = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=retry,
        now=_fixed_now,
    )

    assert result.path.read_bytes() == data
    assert [request.get_method() for request in retry_requests] == ["GET", "HEAD"]


def test_changed_etag_starts_fresh_instead_of_splicing_stale_partial(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    old_data = b"old-representation"
    new_data = b"new-representation"
    assert len(old_data) == len(new_data)

    def old_open(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(_md5_line(old_data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(old_data, '"v1"'))
        return _InterruptingResponse(old_data[:5], total_size=len(old_data), etag='"v1"')

    with pytest.raises(UCSCResourceAcquisitionError):
        _acquire_ucsc_resource(
            url,
            tmp_path,
            terms_acknowledged=True,
            open_url=old_open,
            now=_fixed_now,
        )

    requests: list[Request] = []

    def new_open(request: Request) -> _Response:
        requests.append(request)
        if request.full_url == checksum_url:
            return _Response(_md5_line(new_data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(new_data, '"v2"'))
        assert request.get_header("Range") is None
        return _Response(
            new_data,
            headers={
                "Content-Length": str(len(new_data)),
                "ETag": '"v2"',
            },
        )

    result = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=new_open,
        now=_fixed_now,
    )

    assert result.path.read_bytes() == new_data
    resource_gets = [
        request for request in requests if request.full_url == url and request.get_method() == "GET"
    ]
    assert len(resource_gets) == 1
    assert resource_gets[0].get_header("Range") is None


def test_resume_returning_200_is_not_appended_and_restarts_fresh(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"0123456789abcdef"
    prefix = data[:5]
    etag = '"v1"'

    def seed(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(_md5_line(data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(data, etag))
        return _InterruptingResponse(prefix, total_size=len(data), etag=etag)

    with pytest.raises(UCSCResourceAcquisitionError):
        _acquire_ucsc_resource(
            url, tmp_path, terms_acknowledged=True, open_url=seed, now=_fixed_now
        )

    requests: list[Request] = []

    def restart(request: Request) -> _Response:
        requests.append(request)
        if request.full_url == checksum_url:
            return _Response(_md5_line(data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(data, etag))
        if request.get_header("Range") is not None:
            return _HeadResponse(
                {
                    "Content-Length": str(len(data)),
                    "ETag": etag,
                }
            )
        return _Response(data)

    result = _acquire_ucsc_resource(
        url, tmp_path, terms_acknowledged=True, open_url=restart, now=_fixed_now
    )

    assert result.path.read_bytes() == data
    resource_gets = [
        request for request in requests if request.full_url == url and request.get_method() == "GET"
    ]
    assert len(resource_gets) == 2
    assert resource_gets[0].get_header("Range") == f"bytes={len(prefix)}-"
    assert resource_gets[0].get_header("If-range") == etag
    assert resource_gets[1].get_header("Range") is None


def test_bad_resume_content_range_restarts_fresh_without_splicing(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"0123456789abcdef"
    prefix = data[:5]
    etag = '"v1"'

    def seed(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(_md5_line(data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(data, etag))
        return _InterruptingResponse(prefix, total_size=len(data), etag=etag)

    with pytest.raises(UCSCResourceAcquisitionError):
        _acquire_ucsc_resource(
            url, tmp_path, terms_acknowledged=True, open_url=seed, now=_fixed_now
        )

    requests: list[Request] = []

    def restart(request: Request) -> _Response:
        requests.append(request)
        if request.full_url == checksum_url:
            return _Response(_md5_line(data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(data, etag))
        if request.get_header("Range") is not None:
            tail = data[len(prefix) :]
            return _Response(
                tail,
                status=206,
                headers={
                    "Content-Length": str(len(tail)),
                    "Content-Range": f"bytes 0-{len(tail) - 1}/{len(data)}",
                    "ETag": etag,
                },
            )
        return _Response(data)

    result = _acquire_ucsc_resource(
        url, tmp_path, terms_acknowledged=True, open_url=restart, now=_fixed_now
    )

    assert result.path.read_bytes() == data
    resource_gets = [
        request for request in requests if request.full_url == url and request.get_method() == "GET"
    ]
    assert len(resource_gets) == 2
    assert resource_gets[0].get_header("Range") is not None
    assert resource_gets[1].get_header("Range") is None


def test_concurrent_resumable_writer_cannot_mutate_published_artifact(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    good_data = b"AAAAAAAAAAZZZZZZZZZZ"
    divergent_data = b"AAAAAAAAAAXXXXXXXXXX"
    prefix = good_data[:10]
    etag = '"stable-v1"'

    def seed_partial(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(_md5_line(good_data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(good_data, etag))
        return _InterruptingResponse(prefix, total_size=len(good_data), etag=etag)

    with pytest.raises(UCSCResourceAcquisitionError, match="partial state was retained"):
        _acquire_ucsc_resource(
            url,
            tmp_path,
            terms_acknowledged=True,
            open_url=seed_partial,
            now=_fixed_now,
        )

    partial = next((tmp_path / "partials").rglob("*.part"))
    assert partial.read_bytes() == prefix

    divergent_read_started = Event()
    allow_divergent_tail = Event()
    divergent_errors: list[ResourceChecksumMismatchError] = []

    class _PausedDivergentResponse(_Response):
        def __init__(self) -> None:
            tail = divergent_data[len(prefix) :]
            super().__init__(
                b"",
                status=206,
                headers={
                    "Content-Length": str(len(tail)),
                    "Content-Range": (
                        f"bytes {len(prefix)}-{len(divergent_data) - 1}/"
                        f"{len(divergent_data)}"
                    ),
                    "ETag": etag,
                },
            )
            self._read_number = 0

        def read(self, size: int | None = -1) -> bytes:
            self._read_number += 1
            if self._read_number == 1:
                divergent_read_started.set()
                if not allow_divergent_tail.wait(timeout=5):
                    raise TimeoutError("test did not release divergent writer")
                return divergent_data[len(prefix) :]
            return b""

    def divergent_open(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(_md5_line(good_data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(good_data, etag))
        assert request.get_header("Range") == f"bytes={len(prefix)}-"
        assert request.get_header("If-range") == etag
        return _PausedDivergentResponse()

    def run_divergent_writer() -> None:
        try:
            _acquire_ucsc_resource(
                url,
                tmp_path,
                terms_acknowledged=True,
                open_url=divergent_open,
                now=_fixed_now,
            )
        except ResourceChecksumMismatchError as exc:
            divergent_errors.append(exc)

    divergent_thread = Thread(target=run_divergent_writer)
    divergent_thread.start()
    assert divergent_read_started.wait(timeout=5)

    def good_open(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(_md5_line(good_data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            return _HeadResponse(_resumable_head_headers(good_data, etag))
        assert request.get_header("Range") == f"bytes={len(prefix)}-"
        assert request.get_header("If-range") == etag
        tail = good_data[len(prefix) :]
        return _Response(
            tail,
            status=206,
            headers={
                "Content-Length": str(len(tail)),
                "Content-Range": (
                    f"bytes {len(prefix)}-{len(good_data) - 1}/{len(good_data)}"
                ),
                "ETag": etag,
            },
        )

    result = _acquire_ucsc_resource(
        url,
        tmp_path,
        terms_acknowledged=True,
        open_url=good_open,
        now=_fixed_now,
    )
    assert result.path.read_bytes() == good_data

    allow_divergent_tail.set()
    divergent_thread.join(timeout=5)
    assert not divergent_thread.is_alive()
    assert len(divergent_errors) == 1
    assert isinstance(divergent_errors[0], ResourceChecksumMismatchError)

    expected_sha256 = hashlib.sha256(good_data).hexdigest()
    assert result.sha256 == f"sha256:{expected_sha256}"
    assert result.path.name == expected_sha256
    assert result.path.read_bytes() == good_data
    assert hashlib.sha256(result.path.read_bytes()).hexdigest() == result.path.name


def test_missing_provider_checksum_does_not_enable_persistent_resume_state(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    def open_url(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(b"")
        assert request.get_method() == "GET"
        raise URLError("simulated fresh transfer interruption")

    with pytest.raises(UCSCResourceAcquisitionError, match="failed to download"):
        _acquire_ucsc_resource(
            url, tmp_path, terms_acknowledged=True, open_url=open_url, now=_fixed_now
        )

    assert not (tmp_path / "partials").exists()
    assert not list((tmp_path / "tmp").glob("*.part"))


def test_weak_etag_does_not_enable_persistent_resume_state(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"weak validator bytes"

    def open_url(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(b"")
        if request.get_method() == "HEAD":
            return _HeadResponse(
                {
                    "Content-Length": str(len(data)),
                    "Accept-Ranges": "bytes",
                    "ETag": 'W/"weak-v1"',
                    "Last-Modified": "Fri, 14 Aug 2026 12:00:00 GMT",
                }
            )
        raise URLError("simulated fresh transfer interruption")

    with pytest.raises(UCSCResourceAcquisitionError, match="failed to download"):
        _acquire_ucsc_resource(
            url, tmp_path, terms_acknowledged=True, open_url=open_url, now=_fixed_now
        )

    assert not (tmp_path / "partials").exists()
    assert not list((tmp_path / "tmp").glob("*.part"))


def test_head_failure_falls_back_to_existing_fresh_streaming_download(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"fallback bytes"
    requests: list[Request] = []

    def open_url(request: Request) -> _Response:
        requests.append(request)
        if request.full_url == checksum_url:
            return _Response(_md5_line(data, "canFam3.canFam4.net.gz"))
        if request.get_method() == "HEAD":
            raise URLError("HEAD unavailable")
        return _Response(data)

    result = _acquire_ucsc_resource(
        url, tmp_path, terms_acknowledged=True, open_url=open_url, now=_fixed_now
    )

    assert result.path.read_bytes() == data
    assert [request.get_method() for request in requests] == ["GET", "HEAD", "GET"]
    assert not (tmp_path / "partials").exists()


def test_fresh_download_rejects_non_identity_content_encoding(tmp_path: Path) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    checksum_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt"
    )
    data = b"provider bytes"

    def open_url(request: Request) -> _Response:
        if request.full_url == checksum_url:
            return _Response(b"")
        if request.get_method() == "HEAD":
            return _HeadResponse({})
        assert request.get_header("Accept-encoding") == "identity"
        return _Response(data, headers={"Content-Encoding": "gzip"})

    with pytest.raises(UCSCResourceAcquisitionError, match="identity encoding"):
        _acquire_ucsc_resource(
            url, tmp_path, terms_acknowledged=True, open_url=open_url, now=_fixed_now
        )

    assert not (tmp_path / "by-url").exists()
    assert not list((tmp_path / "tmp").glob("*.part"))


def test_public_acquisition_signature_requires_explicit_acknowledgement(
    tmp_path: Path,
) -> None:
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )

    with pytest.raises(UCSCResourceTermsAcknowledgementRequired):
        acquire_ucsc_resource(url, tmp_path, terms_acknowledged=False)



def _comparative_bundle() -> UCSCResourceBundle:
    return UCSCResourceBundle(
        source_db="canFam3",
        target_db="canFam4",
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        chain_url=(
            "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
            "canFam3.canFam4.all.chain.gz"
        ),
        net_url=(
            "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
            "canFam3.canFam4.net.gz"
        ),
        syntenic_net_url=(
            "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
            "canFam3.canFam4.syn.net.gz"
        ),
        reciprocal_best_chain_url=(
            "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/"
            "reciprocalBest/canFam3.canFam4.rbest.chain.gz"
        ),
        reciprocal_best_net_url=(
            "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/"
            "reciprocalBest/canFam3.canFam4.rbest.net.gz"
        ),
    )


def _liftover_bundle() -> UCSCResourceBundle:
    return UCSCResourceBundle(
        source_db="canFam3",
        target_db="canFam4",
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain_url=(
            "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
            "canFam3ToCanFam4.over.chain.gz"
        ),
    )


def _cached_for_url(url: str, cache_root: Path, *, cache_hit: bool) -> CachedResource:
    terms = ucsc_resource_terms(url)
    digest = hashlib.sha256(url.encode("utf-8")).hexdigest()
    return CachedResource(
        path=cache_root / "artifacts" / "sha256" / digest[:2] / digest,
        source_url=url,
        retrieved_at="2026-08-14T07:30:00Z",
        sha256=f"sha256:{digest}",
        size_bytes=len(url),
        provider_checksum=None,
        terms=terms,
        cache_hit=cache_hit,
    )


def test_comparative_bundle_plan_is_complete_inspectable_and_no_network() -> None:
    bundle = _comparative_bundle()
    plan = plan_ucsc_bundle_acquisition(bundle)

    assert plan.source_db == "canFam3"
    assert plan.target_db == "canFam4"
    assert plan.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
    assert tuple(item.role for item in plan.items) == (
        UCSCBundleResourceRole.CHAIN,
        UCSCBundleResourceRole.NET,
        UCSCBundleResourceRole.SYNTENIC_NET,
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
        UCSCBundleResourceRole.RECIPROCAL_BEST_NET,
    )
    assert tuple(item.url for item in plan.items) == (
        bundle.chain_url,
        bundle.net_url,
        bundle.syntenic_net_url,
        bundle.reciprocal_best_chain_url,
        bundle.reciprocal_best_net_url,
    )
    assert all(
        item.terms.resource_class is UCSCResourceClass.COMPARATIVE
        for item in plan.items
    )
    assert plan.items[-1].terms.directory_terms_url.endswith("/canFam4/vsCanFam3/")


def test_liftover_bundle_plan_surfaces_restricted_chain_terms() -> None:
    plan = plan_ucsc_bundle_acquisition(_liftover_bundle())

    assert tuple(item.role for item in plan.items) == (UCSCBundleResourceRole.CHAIN,)
    assert plan.items[0].terms.resource_class is UCSCResourceClass.LIFTOVER_CHAIN
    assert plan.items[0].terms.restricted_liftover_chain is True


def test_bundle_transfer_plan_acknowledgement_precedes_any_resource_acquisition(
    tmp_path: Path,
) -> None:
    plan = plan_ucsc_bundle_acquisition(_comparative_bundle())

    def forbidden(
        url: str,
        cache_root: str | PathLike[str],
        *,
        terms_acknowledged: bool,
        refresh: bool = False,
    ) -> CachedResource:
        raise AssertionError(
            f"resource acquisition must not start before plan acknowledgement: {url}"
        )

    with pytest.raises(UCSCBundleAcquisitionPlanAcknowledgementRequired):
        _acquire_ucsc_resource_bundle(
            plan,
            tmp_path,
            transfer_plan_acknowledged=False,
            terms_acknowledged=True,
            refresh=False,
            acquire_resource=forbidden,
        )

    assert not list(tmp_path.iterdir())


def test_comparative_bundle_acquisition_returns_only_after_all_roles_complete(
    tmp_path: Path,
) -> None:
    plan = plan_ucsc_bundle_acquisition(_comparative_bundle())
    calls: list[tuple[str, bool, bool]] = []

    def acquire(
        url: str,
        cache_root: str | PathLike[str],
        *,
        terms_acknowledged: bool,
        refresh: bool = False,
    ) -> CachedResource:
        calls.append((url, terms_acknowledged, refresh))
        return _cached_for_url(url, Path(cache_root), cache_hit=False)

    result = _acquire_ucsc_resource_bundle(
        plan,
        tmp_path,
        transfer_plan_acknowledged=True,
        terms_acknowledged=True,
        refresh=True,
        acquire_resource=acquire,
    )

    assert isinstance(result, CachedUCSCResourceBundle)
    assert result.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
    assert result.chain.source_url == plan.items[0].url
    assert result.net is not None and result.net.source_url == plan.items[1].url
    assert result.syntenic_net is not None
    assert result.syntenic_net.source_url == plan.items[2].url
    assert result.reciprocal_best_chain is not None
    assert result.reciprocal_best_chain.source_url == plan.items[3].url
    assert result.reciprocal_best_net is not None
    assert result.reciprocal_best_net.source_url == plan.items[4].url
    assert calls == [(item.url, True, True) for item in plan.items]


def test_liftover_bundle_acquisition_returns_chain_only(tmp_path: Path) -> None:
    plan = plan_ucsc_bundle_acquisition(_liftover_bundle())

    def acquire(
        url: str,
        cache_root: str | PathLike[str],
        *,
        terms_acknowledged: bool,
        refresh: bool = False,
    ) -> CachedResource:
        return _cached_for_url(url, Path(cache_root), cache_hit=True)

    result = _acquire_ucsc_resource_bundle(
        plan,
        tmp_path,
        transfer_plan_acknowledged=True,
        terms_acknowledged=True,
        refresh=False,
        acquire_resource=acquire,
    )

    assert result.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY
    assert result.chain.cache_hit is True
    assert result.net is None
    assert result.syntenic_net is None
    assert result.reciprocal_best_chain is None
    assert result.reciprocal_best_net is None


def test_bundle_failure_propagates_without_returning_partial_bundle(tmp_path: Path) -> None:
    plan = plan_ucsc_bundle_acquisition(_comparative_bundle())
    acquired_urls: list[str] = []

    def acquire(
        url: str,
        cache_root: str | PathLike[str],
        *,
        terms_acknowledged: bool,
        refresh: bool = False,
    ) -> CachedResource:
        acquired_urls.append(url)
        if len(acquired_urls) == 3:
            raise UCSCResourceAcquisitionError("simulated third-resource failure")
        return _cached_for_url(url, Path(cache_root), cache_hit=False)

    with pytest.raises(UCSCResourceAcquisitionError, match="third-resource failure"):
        _acquire_ucsc_resource_bundle(
            plan,
            tmp_path,
            transfer_plan_acknowledged=True,
            terms_acknowledged=True,
            refresh=False,
            acquire_resource=acquire,
        )

    assert acquired_urls == [item.url for item in plan.items[:3]]


def test_public_bundle_acquisition_signature_requires_explicit_plan_acknowledgement(
    tmp_path: Path,
) -> None:
    plan = plan_ucsc_bundle_acquisition(_liftover_bundle())

    with pytest.raises(UCSCBundleAcquisitionPlanAcknowledgementRequired):
        acquire_ucsc_resource_bundle(
            plan,
            tmp_path,
            transfer_plan_acknowledged=False,
            terms_acknowledged=True,
        )



def test_bundle_plan_item_rejects_terms_that_do_not_match_url() -> None:
    comparative_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )
    restricted_terms = ucsc_resource_terms(
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
        "canFam3ToCanFam4.over.chain.gz"
    )

    with pytest.raises(ValueError, match="terms must match"):
        UCSCBundleAcquisitionItem(
            role=UCSCBundleResourceRole.NET,
            url=comparative_url,
            terms=restricted_terms,
        )


def test_bundle_terms_acknowledgement_still_precedes_provider_network(
    tmp_path: Path,
) -> None:
    plan = plan_ucsc_bundle_acquisition(_liftover_bundle())

    with pytest.raises(UCSCResourceTermsAcknowledgementRequired):
        acquire_ucsc_resource_bundle(
            plan,
            tmp_path,
            transfer_plan_acknowledged=True,
            terms_acknowledged=False,
        )



def test_comparative_bundle_plan_rejects_partial_role_set() -> None:
    item = plan_ucsc_bundle_acquisition(_comparative_bundle()).items[0]

    with pytest.raises(ValueError, match="exact ordered resource roles"):
        UCSCBundleAcquisitionPlan(
            source_db="canFam3",
            target_db="canFam4",
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
            items=(item,),
        )


def test_cached_comparative_bundle_rejects_partial_resource_state(
    tmp_path: Path,
) -> None:
    chain_url = _comparative_bundle().chain_url
    chain = _cached_for_url(chain_url, tmp_path, cache_hit=True)

    with pytest.raises(ValueError, match="COMPARATIVE cached bundle requires"):
        CachedUCSCResourceBundle(
            source_db="canFam3",
            target_db="canFam4",
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
            chain=chain,
        )


def test_bundle_plan_item_rejects_role_filename_mismatch() -> None:
    net_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
        "canFam3.canFam4.net.gz"
    )

    with pytest.raises(ValueError, match="role CHAIN does not match filename"):
        UCSCBundleAcquisitionItem(
            role=UCSCBundleResourceRole.CHAIN,
            url=net_url,
            terms=ucsc_resource_terms(net_url),
        )


def test_bundle_plan_rejects_right_role_from_wrong_directional_pair() -> None:
    valid = plan_ucsc_bundle_acquisition(_comparative_bundle())
    wrong_pair_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam5/"
        "canFam3.canFam5.net.gz"
    )
    items = list(valid.items)
    items[1] = UCSCBundleAcquisitionItem(
        role=UCSCBundleResourceRole.NET,
        url=wrong_pair_url,
        terms=ucsc_resource_terms(wrong_pair_url),
    )

    with pytest.raises(ValueError, match="must use directional filename"):
        UCSCBundleAcquisitionPlan(
            source_db=valid.source_db,
            target_db=valid.target_db,
            evidence_tier=valid.evidence_tier,
            items=tuple(items),
        )


def test_liftover_plan_rejects_chain_for_different_target() -> None:
    wrong_target_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
        "canFam3ToCanFam5.over.chain.gz"
    )
    item = UCSCBundleAcquisitionItem(
        role=UCSCBundleResourceRole.CHAIN,
        url=wrong_target_url,
        terms=ucsc_resource_terms(wrong_target_url),
    )

    with pytest.raises(ValueError, match="must use directional filename"):
        UCSCBundleAcquisitionPlan(
            source_db="canFam3",
            target_db="canFam4",
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
            items=(item,),
        )


def test_cached_bundle_rejects_swapped_chain_and_net_resources(
    tmp_path: Path,
) -> None:
    bundle = _comparative_bundle()
    assert bundle.net_url is not None
    assert bundle.syntenic_net_url is not None
    assert bundle.reciprocal_best_chain_url is not None
    assert bundle.reciprocal_best_net_url is not None

    with pytest.raises(ValueError, match="bundle resource role CHAIN"):
        CachedUCSCResourceBundle(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=_cached_for_url(bundle.net_url, tmp_path, cache_hit=True),
            net=_cached_for_url(bundle.chain_url, tmp_path, cache_hit=True),
            syntenic_net=_cached_for_url(
                bundle.syntenic_net_url, tmp_path, cache_hit=True
            ),
            reciprocal_best_chain=_cached_for_url(
                bundle.reciprocal_best_chain_url, tmp_path, cache_hit=True
            ),
            reciprocal_best_net=_cached_for_url(
                bundle.reciprocal_best_net_url, tmp_path, cache_hit=True
            ),
        )


def test_cached_bundle_rejects_resource_from_wrong_directional_pair(
    tmp_path: Path,
) -> None:
    bundle = _comparative_bundle()
    assert bundle.net_url is not None
    assert bundle.syntenic_net_url is not None
    assert bundle.reciprocal_best_chain_url is not None
    assert bundle.reciprocal_best_net_url is not None
    wrong_pair_net_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam5/"
        "canFam3.canFam5.net.gz"
    )

    with pytest.raises(ValueError, match="must use directional filename"):
        CachedUCSCResourceBundle(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=_cached_for_url(bundle.chain_url, tmp_path, cache_hit=True),
            net=_cached_for_url(wrong_pair_net_url, tmp_path, cache_hit=True),
            syntenic_net=_cached_for_url(
                bundle.syntenic_net_url, tmp_path, cache_hit=True
            ),
            reciprocal_best_chain=_cached_for_url(
                bundle.reciprocal_best_chain_url, tmp_path, cache_hit=True
            ),
            reciprocal_best_net=_cached_for_url(
                bundle.reciprocal_best_net_url, tmp_path, cache_hit=True
            ),
        )


class _HeadResponse(_Response):
    def __init__(self, headers: Mapping[str, str]) -> None:
        super().__init__(b"", headers=headers)
        if not any(key.casefold() == "content-length" for key in headers):
            self._headers.pop("content-length", None)

    def read(self, size: int | None = -1) -> bytes:
        raise AssertionError("HEAD metadata inspection must not read a response body")

    def getheader(self, name: str, default: str | None = None) -> str | None:
        return self._headers.get(name.casefold(), default)


def test_remote_metadata_inspection_uses_head_without_reading_body() -> None:
    url = _comparative_bundle().chain_url
    requests: list[Request] = []

    def open_url(request: Request) -> _HeadResponse:
        requests.append(request)
        return _HeadResponse(
            {
                "Content-Length": "2684354560",
                "Accept-Ranges": "bytes",
                "Last-Modified": "Tue, 12 May 2020 22:56:00 GMT",
                "ETag": '"example"',
            }
        )

    metadata = _inspect_ucsc_resource(
        url,
        terms_acknowledged=True,
        open_url=open_url,
    )

    assert len(requests) == 1
    assert requests[0].get_method() == "HEAD"
    assert requests[0].get_header("Accept-encoding") == "identity"
    assert metadata.url == url
    assert metadata.content_length_bytes == 2684354560
    assert metadata.accept_ranges == "bytes"
    assert metadata.last_modified == "Tue, 12 May 2020 22:56:00 GMT"
    assert metadata.etag == '"example"'
    assert metadata.content_encoding is None
    assert metadata.terms == ucsc_resource_terms(url)


def test_remote_metadata_inspection_requires_terms_before_network_access() -> None:
    url = _comparative_bundle().chain_url
    called = False

    def open_url(_: Request) -> _HeadResponse:
        nonlocal called
        called = True
        return _HeadResponse({})

    with pytest.raises(UCSCResourceTermsAcknowledgementRequired):
        _inspect_ucsc_resource(
            url,
            terms_acknowledged=False,
            open_url=open_url,
        )

    assert called is False


def test_remote_metadata_inspection_preserves_missing_optional_headers() -> None:
    url = _comparative_bundle().net_url
    assert url is not None

    metadata = _inspect_ucsc_resource(
        url,
        terms_acknowledged=True,
        open_url=lambda _: _HeadResponse({}),
    )

    assert metadata.content_length_bytes is None
    assert metadata.accept_ranges is None
    assert metadata.last_modified is None
    assert metadata.etag is None
    assert metadata.content_encoding is None


def test_remote_metadata_inspection_rejects_malformed_content_length() -> None:
    url = _comparative_bundle().chain_url

    with pytest.raises(UCSCResourceAcquisitionError, match="Content-Length"):
        _inspect_ucsc_resource(
            url,
            terms_acknowledged=True,
            open_url=lambda _: _HeadResponse({"Content-Length": "not-a-number"}),
        )


def test_remote_metadata_inspection_propagates_transport_failure() -> None:
    url = _comparative_bundle().chain_url

    def fail(_: Request) -> _HeadResponse:
        raise URLError("simulated metadata outage")

    with pytest.raises(UCSCResourceAcquisitionError, match="metadata outage"):
        _inspect_ucsc_resource(
            url,
            terms_acknowledged=True,
            open_url=fail,
        )


def test_bundle_transfer_inspection_requires_terms_before_any_item_inspection() -> None:
    plan = plan_ucsc_bundle_acquisition(_comparative_bundle())
    calls: list[str] = []

    def inspect(url: str) -> UCSCRemoteResourceMetadata:
        calls.append(url)
        raise AssertionError("bundle inspection must fail before contacting the provider")

    with pytest.raises(UCSCResourceTermsAcknowledgementRequired):
        _inspect_ucsc_bundle_transfer_plan(
            plan,
            terms_acknowledged=False,
            inspect_resource=inspect,
        )

    assert calls == []


def test_bundle_transfer_inspection_preserves_exact_plan_and_totals() -> None:
    plan = plan_ucsc_bundle_acquisition(_comparative_bundle())
    sizes = {item.url: (index + 1) * 100 for index, item in enumerate(plan.items)}

    def inspect(url: str) -> UCSCRemoteResourceMetadata:
        return UCSCRemoteResourceMetadata(
            url=url,
            terms=ucsc_resource_terms(url),
            content_length_bytes=sizes[url],
            accept_ranges="bytes",
            last_modified=None,
            etag=None,
            content_encoding=None,
        )

    result = _inspect_ucsc_bundle_transfer_plan(
        plan,
        terms_acknowledged=True,
        inspect_resource=inspect,
    )

    assert isinstance(result, UCSCBundleTransferInspection)
    assert result.source_db == plan.source_db
    assert result.target_db == plan.target_db
    assert result.evidence_tier is plan.evidence_tier
    assert tuple(item.role for item in result.items) == tuple(
        item.role for item in plan.items
    )
    assert tuple(item.metadata.url for item in result.items) == tuple(
        item.url for item in plan.items
    )
    assert result.known_content_length_bytes == 1500
    assert result.total_content_length_bytes == 1500


def test_bundle_transfer_inspection_excludes_non_identity_encoded_lengths() -> None:
    plan = plan_ucsc_bundle_acquisition(_liftover_bundle())

    def inspect(url: str) -> UCSCRemoteResourceMetadata:
        return UCSCRemoteResourceMetadata(
            url=url,
            terms=ucsc_resource_terms(url),
            content_length_bytes=1234,
            accept_ranges="bytes",
            last_modified=None,
            etag=None,
            content_encoding="gzip",
        )

    result = _inspect_ucsc_bundle_transfer_plan(
        plan,
        terms_acknowledged=True,
        inspect_resource=inspect,
    )

    assert result.items[0].metadata.content_length_bytes == 1234
    assert result.items[0].metadata.content_encoding == "gzip"
    assert result.known_content_length_bytes == 0
    assert result.total_content_length_bytes is None


def test_bundle_transfer_inspection_does_not_guess_unknown_total_size() -> None:
    plan = plan_ucsc_bundle_acquisition(_liftover_bundle())

    def inspect(url: str) -> UCSCRemoteResourceMetadata:
        return UCSCRemoteResourceMetadata(
            url=url,
            terms=ucsc_resource_terms(url),
            content_length_bytes=None,
            accept_ranges=None,
            last_modified=None,
            etag=None,
            content_encoding=None,
        )

    result = _inspect_ucsc_bundle_transfer_plan(
        plan,
        terms_acknowledged=True,
        inspect_resource=inspect,
    )

    assert result.known_content_length_bytes == 0
    assert result.total_content_length_bytes is None


def test_bundle_transfer_inspection_rejects_wrong_pair_metadata() -> None:
    plan = plan_ucsc_bundle_acquisition(_comparative_bundle())
    wrong_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam5/"
        "canFam3.canFam5.net.gz"
    )
    inspected = tuple(
        UCSCBundleTransferInspectionItem(
            role=item.role,
            metadata=UCSCRemoteResourceMetadata(
                url=(wrong_url if item.role is UCSCBundleResourceRole.NET else item.url),
                terms=ucsc_resource_terms(
                    wrong_url if item.role is UCSCBundleResourceRole.NET else item.url
                ),
                content_length_bytes=1,
                accept_ranges=None,
                last_modified=None,
                etag=None,
                content_encoding=None,
            ),
        )
        for item in plan.items
    )

    with pytest.raises(ValueError, match="must use directional filename"):
        UCSCBundleTransferInspection(
            source_db=plan.source_db,
            target_db=plan.target_db,
            evidence_tier=plan.evidence_tier,
            items=inspected,
        )


def test_public_remote_metadata_inspection_is_exported() -> None:
    assert callable(inspect_ucsc_resource)
    assert callable(inspect_ucsc_bundle_transfer_plan)
