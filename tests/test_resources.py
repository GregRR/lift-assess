from __future__ import annotations

from collections.abc import Callable, Mapping

import pytest

from liftassess import EvidenceAvailabilityTier
from liftassess.resources import (
    UCSCAssemblyMetadataResources,
    UCSCResourceBundle,
    UCSCResourceDiscoveryError,
    UCSCSegmentalDuplicationResource,
    _discover_ucsc_assembly_metadata,
    _discover_ucsc_resources,
    _discover_ucsc_segmental_duplication_resource,
    _ucsc_title_db,
)


def _reader(
    listings: Mapping[str, frozenset[str] | None],
) -> Callable[[str], frozenset[str] | None]:
    def read(url: str) -> frozenset[str] | None:
        return listings.get(url)

    return read


def test_discovers_ucsc_assembly_metadata_from_observed_database_tables() -> None:
    database = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/"

    result = _discover_ucsc_assembly_metadata(
        "canFam3",
        _reader(
            {
                database: frozenset(
                    {"chromInfo.txt.gz", "chromAlias.txt.gz", "README.txt"}
                )
            }
        ),
    )

    assert result == UCSCAssemblyMetadataResources(
        db="canFam3",
        chrom_info_url=f"{database}chromInfo.txt.gz",
        chrom_alias_url=f"{database}chromAlias.txt.gz",
    )


def test_discovers_segmental_duplication_table_only_when_observed() -> None:
    database = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/"

    result = _discover_ucsc_segmental_duplication_resource(
        "hg38",
        _reader({database: frozenset({"genomicSuperDups.txt.gz"})}),
    )

    assert result == UCSCSegmentalDuplicationResource(
        db="hg38",
        url=f"{database}genomicSuperDups.txt.gz",
    )

    absent = _discover_ucsc_segmental_duplication_resource(
        "hg38",
        _reader({database: frozenset({"chromInfo.txt.gz"})}),
    )
    assert absent is None


def test_assembly_metadata_discovery_requires_observed_chrom_info() -> None:
    database = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/"

    result = _discover_ucsc_assembly_metadata(
        "canFam3",
        _reader({database: frozenset({"chromAlias.txt.gz"})}),
    )

    assert result is None


def test_discovers_complete_comparative_resource_set() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    reciprocal = f"{comparative}reciprocalBest/"
    result = _discover_ucsc_resources(
        "canFam3",
        "canFam4",
        _reader(
            {
                comparative: frozenset(
                    {
                        "canFam3.canFam4.all.chain.gz",
                        "canFam3.canFam4.net.gz",
                        "canFam3.canFam4.syn.net.gz",
                        "reciprocalBest/",
                        "README.txt",
                    }
                ),
                reciprocal: frozenset(
                    {
                        "canFam3.canFam4.rbest.chain.gz",
                        "canFam3.canFam4.rbest.net.gz",
                    }
                ),
            }
        ),
    )

    assert result == UCSCResourceBundle(
        source_db="canFam3",
        target_db="canFam4",
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        chain_url=f"{comparative}canFam3.canFam4.all.chain.gz",
        net_url=f"{comparative}canFam3.canFam4.net.gz",
        syntenic_net_url=f"{comparative}canFam3.canFam4.syn.net.gz",
        reciprocal_best_chain_url=(f"{reciprocal}canFam3.canFam4.rbest.chain.gz"),
        reciprocal_best_net_url=f"{reciprocal}canFam3.canFam4.rbest.net.gz",
    )


def test_complete_forward_reciprocal_best_does_not_require_reverse_listing() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    reciprocal = f"{comparative}reciprocalBest/"
    reverse_comparative = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/"
    )

    listings = {
        comparative: frozenset(
            {
                "canFam3.canFam4.all.chain.gz",
                "canFam3.canFam4.net.gz",
                "canFam3.canFam4.syn.net.gz",
                "reciprocalBest/",
            }
        ),
        reciprocal: frozenset(
            {
                "canFam3.canFam4.rbest.chain.gz",
                "canFam3.canFam4.rbest.net.gz",
            }
        ),
    }

    def read(url: str) -> frozenset[str] | None:
        if url == reverse_comparative:
            raise AssertionError("reverse listing should not be checked")
        return listings.get(url)

    result = _discover_ucsc_resources("canFam3", "canFam4", read)

    assert result is not None
    assert result.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
    assert result.reciprocal_best_chain_url == (
        f"{reciprocal}canFam3.canFam4.rbest.chain.gz"
    )


def test_discovers_reciprocal_best_from_reverse_pair_directory() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    reverse_comparative = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/"
    )
    reverse_reciprocal = f"{reverse_comparative}reciprocalBest/"

    result = _discover_ucsc_resources(
        "canFam3",
        "canFam4",
        _reader(
            {
                comparative: frozenset(
                    {
                        "canFam3.canFam4.all.chain.gz",
                        "canFam3.canFam4.net.gz",
                        "canFam3.canFam4.syn.net.gz",
                    }
                ),
                reverse_comparative: frozenset({"reciprocalBest/"}),
                reverse_reciprocal: frozenset(
                    {
                        "canFam3.canFam4.rbest.chain.gz",
                        "canFam3.canFam4.rbest.net.gz",
                        "canFam4.canFam3.rbest.chain.gz",
                        "canFam4.canFam3.rbest.net.gz",
                    }
                ),
            }
        ),
    )

    assert result == UCSCResourceBundle(
        source_db="canFam3",
        target_db="canFam4",
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        chain_url=f"{comparative}canFam3.canFam4.all.chain.gz",
        net_url=f"{comparative}canFam3.canFam4.net.gz",
        syntenic_net_url=f"{comparative}canFam3.canFam4.syn.net.gz",
        reciprocal_best_chain_url=(
            f"{reverse_reciprocal}canFam3.canFam4.rbest.chain.gz"
        ),
        reciprocal_best_net_url=(f"{reverse_reciprocal}canFam3.canFam4.rbest.net.gz"),
    )


def test_reverse_pair_lookup_transport_failure_propagates() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    reverse_comparative = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/"
    )

    listings = {
        comparative: frozenset(
            {
                "canFam3.canFam4.all.chain.gz",
                "canFam3.canFam4.net.gz",
                "canFam3.canFam4.syn.net.gz",
            }
        ),
    }

    def read(url: str) -> frozenset[str] | None:
        if url == reverse_comparative:
            raise UCSCResourceDiscoveryError("simulated sibling transport failure")
        return listings.get(url)

    with pytest.raises(
        UCSCResourceDiscoveryError, match="simulated sibling transport failure"
    ):
        _discover_ucsc_resources("canFam3", "canFam4", read)


def test_reverse_pair_fallback_requires_exact_directional_rbest_files() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    reverse_comparative = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/"
    )
    reverse_reciprocal = f"{reverse_comparative}reciprocalBest/"
    liftover = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"

    result = _discover_ucsc_resources(
        "canFam3",
        "canFam4",
        _reader(
            {
                comparative: frozenset(
                    {
                        "canFam3.canFam4.all.chain.gz",
                        "canFam3.canFam4.net.gz",
                        "canFam3.canFam4.syn.net.gz",
                    }
                ),
                reverse_comparative: frozenset({"reciprocalBest/"}),
                # The sibling directory exists but has only the opposite direction.
                reverse_reciprocal: frozenset(
                    {
                        "canFam4.canFam3.rbest.chain.gz",
                        "canFam4.canFam3.rbest.net.gz",
                    }
                ),
                liftover: frozenset({"canFam3ToCanFam4.over.chain.gz"}),
            }
        ),
    )

    assert result is not None
    assert result.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY


def test_missing_forward_directional_rbest_files_checks_reverse_pair() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    reciprocal = f"{comparative}reciprocalBest/"
    reverse_comparative = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/"
    )
    reverse_reciprocal = f"{reverse_comparative}reciprocalBest/"

    result = _discover_ucsc_resources(
        "canFam3",
        "canFam4",
        _reader(
            {
                comparative: frozenset(
                    {
                        "canFam3.canFam4.all.chain.gz",
                        "canFam3.canFam4.net.gz",
                        "canFam3.canFam4.syn.net.gz",
                        "reciprocalBest/",
                    }
                ),
                reciprocal: frozenset({"README.txt"}),
                reverse_comparative: frozenset({"reciprocalBest/"}),
                reverse_reciprocal: frozenset(
                    {
                        "canFam3.canFam4.rbest.chain.gz",
                        "canFam3.canFam4.rbest.net.gz",
                    }
                ),
            }
        ),
    )

    assert result is not None
    assert result.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
    assert result.reciprocal_best_chain_url == (
        f"{reverse_reciprocal}canFam3.canFam4.rbest.chain.gz"
    )


def test_discovery_accepts_absolute_directory_hrefs() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    reciprocal = f"{comparative}reciprocalBest/"
    result = _discover_ucsc_resources(
        "canFam3",
        "canFam4",
        _reader(
            {
                comparative: frozenset(
                    {
                        f"{comparative}canFam3.canFam4.all.chain.gz",
                        f"{comparative}canFam3.canFam4.net.gz",
                        f"{comparative}canFam3.canFam4.syn.net.gz",
                        reciprocal,
                    }
                ),
                reciprocal: frozenset(
                    {
                        f"{reciprocal}canFam3.canFam4.rbest.chain.gz",
                        f"{reciprocal}canFam3.canFam4.rbest.net.gz",
                    }
                ),
            }
        ),
    )

    assert result is not None
    assert result.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE


def test_partial_comparative_publication_falls_back_to_liftover_only() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    liftover = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
    result = _discover_ucsc_resources(
        "canFam3",
        "canFam4",
        _reader(
            {
                comparative: frozenset(
                    {
                        "canFam3.canFam4.all.chain.gz",
                        "canFam3.canFam4.net.gz",
                        # syn.net and reciprocalBest are deliberately absent.
                    }
                ),
                liftover: frozenset({"canFam3ToCanFam4.over.chain.gz"}),
            }
        ),
    )

    assert result is not None
    assert result.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY
    assert result.chain_url == f"{liftover}canFam3ToCanFam4.over.chain.gz"
    assert result.net_url is None
    assert result.reciprocal_best_chain_url is None


def test_missing_reciprocal_best_files_prevent_comparative_tier() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    reciprocal = f"{comparative}reciprocalBest/"
    liftover = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
    result = _discover_ucsc_resources(
        "canFam3",
        "canFam4",
        _reader(
            {
                comparative: frozenset(
                    {
                        "canFam3.canFam4.all.chain.gz",
                        "canFam3.canFam4.net.gz",
                        "canFam3.canFam4.syn.net.gz",
                        "reciprocalBest/",
                    }
                ),
                reciprocal: frozenset({"README.txt"}),
                liftover: frozenset({"canFam3ToCanFam4.over.chain.gz"}),
            }
        ),
    )

    assert result is not None
    assert result.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY


def test_returns_none_when_no_verified_resource_exists() -> None:
    result = _discover_ucsc_resources(
        "canFam3",
        "canFam6",
        _reader(
            {
                "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam6/": None,
                "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/": (
                    frozenset({"canFam3ToCanFam4.over.chain.gz"})
                ),
            }
        ),
    )

    assert result is None


def test_missing_comparative_directory_still_checks_liftover_listing() -> None:
    liftover = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
    result = _discover_ucsc_resources(
        "canFam3",
        "canFam6",
        _reader(
            {
                "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam6/": None,
                liftover: frozenset({"canFam3ToCanFam6.over.chain.gz"}),
            }
        ),
    )

    assert result is not None
    assert result.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY
    assert result.chain_url == f"{liftover}canFam3ToCanFam6.over.chain.gz"


def test_comparative_bundle_rejects_missing_required_url() -> None:
    with pytest.raises(ValueError, match="COMPARATIVE resource bundle requires"):
        UCSCResourceBundle(
            source_db="canFam3",
            target_db="canFam4",
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
            chain_url="https://example.invalid/all.chain.gz",
        )


def test_liftover_only_bundle_rejects_comparative_url() -> None:
    with pytest.raises(ValueError, match="cannot carry comparative URLs"):
        UCSCResourceBundle(
            source_db="canFam3",
            target_db="canFam4",
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
            chain_url="https://example.invalid/over.chain.gz",
            net_url="https://example.invalid/net.gz",
        )


def test_ucsc_title_db_matches_kent_ucfirst_convention() -> None:
    assert _ucsc_title_db("canFam4") == "CanFam4"
    assert _ucsc_title_db("hg38") == "Hg38"
    assert _ucsc_title_db("GCF_009762305.2") == "GCF_009762305.2"


@pytest.mark.parametrize("db", ["", "../canFam4", "canFam4/x", "can Fam4"])
def test_rejects_unsafe_or_empty_ucsc_database_identifiers(db: str) -> None:
    with pytest.raises(ValueError, match="UCSC database identifier"):
        _ucsc_title_db(db)


def test_directory_listing_parser_preserves_exact_href_entries() -> None:
    from liftassess.resources import _parse_directory_links

    html = """
    <html><body>
      <a href="../">Parent Directory</a>
      <a href="canFam3.canFam4.all.chain.gz">chain</a>
      <a href="reciprocalBest/">reciprocalBest</a>
    </body></html>
    """

    assert _parse_directory_links(html) == frozenset(
        {"../", "canFam3.canFam4.all.chain.gz", "reciprocalBest/"}
    )


def test_directory_listing_404_means_candidate_directory_absent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from http.client import HTTPMessage
    from typing import NoReturn
    from urllib.error import HTTPError
    from urllib.request import Request

    import liftassess.resources as resource_module

    def raise_404(request: Request, *, timeout: int) -> NoReturn:
        del timeout
        raise HTTPError(request.full_url, 404, "not found", HTTPMessage(), None)

    monkeypatch.setattr(resource_module, "urlopen", raise_404)

    assert resource_module._read_directory_links("https://example.invalid/") is None


def test_directory_listing_server_error_is_not_resource_absence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from http.client import HTTPMessage
    from typing import NoReturn
    from urllib.error import HTTPError
    from urllib.request import Request

    import liftassess.resources as resource_module

    def raise_500(request: Request, *, timeout: int) -> NoReturn:
        del timeout
        raise HTTPError(request.full_url, 500, "server error", HTTPMessage(), None)

    monkeypatch.setattr(resource_module, "urlopen", raise_500)

    with pytest.raises(resource_module.UCSCResourceDiscoveryError, match="HTTP 500"):
        resource_module._read_directory_links("https://example.invalid/")


def test_directory_listing_transport_error_is_not_resource_absence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from typing import NoReturn
    from urllib.error import URLError
    from urllib.request import Request

    import liftassess.resources as resource_module

    def raise_transport_error(request: Request, *, timeout: int) -> NoReturn:
        del request, timeout
        raise URLError("temporary DNS failure")

    monkeypatch.setattr(resource_module, "urlopen", raise_transport_error)

    with pytest.raises(
        resource_module.UCSCResourceDiscoveryError, match="temporary DNS failure"
    ):
        resource_module._read_directory_links("https://example.invalid/")


def test_directory_listing_timeout_is_not_resource_absence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from typing import NoReturn
    from urllib.request import Request

    import liftassess.resources as resource_module

    def raise_timeout(request: Request, *, timeout: int) -> NoReturn:
        del request, timeout
        raise TimeoutError("timed out while reading response")

    monkeypatch.setattr(resource_module, "urlopen", raise_timeout)

    with pytest.raises(
        resource_module.UCSCResourceDiscoveryError, match="request timed out"
    ):
        resource_module._read_directory_links("https://example.invalid/")


def test_exact_liftover_discovery_skips_comparative_publication() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    liftover = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"

    def read(url: str) -> frozenset[str] | None:
        if url == comparative:
            raise AssertionError("exact LIFTOVER-ONLY discovery must skip comparative")
        if url == liftover:
            return frozenset({"canFam3ToCanFam4.over.chain.gz"})
        return None

    result = _discover_ucsc_resources(
        "canFam3",
        "canFam4",
        read,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert result == UCSCResourceBundle(
        source_db="canFam3",
        target_db="canFam4",
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain_url=f"{liftover}canFam3ToCanFam4.over.chain.gz",
    )


def test_exact_comparative_discovery_does_not_fall_back_to_liftover() -> None:
    comparative = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    liftover = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"

    def read(url: str) -> frozenset[str] | None:
        if url == comparative:
            return frozenset()
        if url == liftover:
            raise AssertionError("exact COMPARATIVE discovery must not fall back")
        return None

    result = _discover_ucsc_resources(
        "canFam3",
        "canFam4",
        read,
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert result is None
