"""Discover UCSC-hosted evidence resources without downloading them.

The UCSC Golden Path directory layout is useful for discovery, but liftAssess
must not treat a constructed URL as evidence that a resource exists.  This
module therefore reads the published directory listings and only returns URLs
corresponding to entries that were actually observed.  Directory links are
resolved against the listing URL before comparison so discovery does not depend
on whether UCSC renders an entry as a relative or absolute ``href``.

Implementation references (checked through 2026-08-12):

- UCSC Golden Path canFam3/canFam4 comparative directory and README:
  https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/
- UCSC canFam4/canFam3 reciprocal-best directory.  The live publication layout
  stores both directional reciprocal-best resources under this sibling/reverse
  comparison directory, including ``canFam3.canFam4.rbest.*``:
  https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/reciprocalBest/
- UCSC liftOver download README, which defines
  ``<db1>To<Db2>.over.chain.gz`` naming:
  https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/
- UCSC Kent ``doRecipBest.pl``.  The download step uses ``vs$QDb`` where
  ``$QDb = ucfirst($qDb)`` and publishes ``$tDb.$qDb.rbest.*`` files under
  ``reciprocalBest/``:
  https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/utils/automation/doRecipBest.pl

These are provider layout semantics, not a biological method.  Discovery is
kept separate from the later downloader/cache layer so network availability or
UCSC licensing acceptance cannot be confused with evidence interpretation.
In particular, discovering a ``liftOver/*.over.chain.gz`` URL is only an
availability result: it must not be treated as permission to download or use a
resource.  A future downloader must apply the terms published for the actual
resource class, surface UCSC's restricted liftOver-chain terms when applicable,
and obtain explicit user acknowledgement before retrieval.  Comparative
``vsTarget/`` resources must likewise follow their own published directory
terms rather than inheriting a license solely from their file format.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from html.parser import HTMLParser
from re import fullmatch
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen

from .models import EvidenceAvailabilityTier

_UCSC_GOLDEN_PATH = "https://hgdownload.soe.ucsc.edu/goldenPath/"
_UCSC_DB_PATTERN = r"[A-Za-z0-9_.-]+"
_USER_AGENT = "liftAssess/0.0 resource-discovery"


class UCSCResourceDiscoveryError(RuntimeError):
    """A UCSC listing could not be checked reliably.

    Missing directories are normal discovery results and do not raise this
    exception.  Transport/server failures do raise: treating a temporary
    network error as "resource absent" could silently downgrade the evidence
    tier, which would be scientifically misleading.
    """


@dataclass(frozen=True)
class UCSCResourceBundle:
    """Verified UCSC URLs for one source-to-target assembly comparison.

    ``COMPARATIVE`` means all v1 comparative resources represented here were
    observed in UCSC's directory listings.  ``LIFTOVER_ONLY`` carries only the
    liftOver chain.  This is an evidence-*availability* statement; it says
    nothing about the support for any candidate or eventual assessment verdict.
    """

    source_db: str
    target_db: str
    evidence_tier: EvidenceAvailabilityTier
    chain_url: str
    net_url: str | None = None
    syntenic_net_url: str | None = None
    reciprocal_best_chain_url: str | None = None
    reciprocal_best_net_url: str | None = None

    def __post_init__(self) -> None:
        _validate_ucsc_db(self.source_db)
        _validate_ucsc_db(self.target_db)
        if not self.chain_url:
            raise ValueError("UCSC resource bundle chain_url must not be empty")

        comparative_urls = (
            self.net_url,
            self.syntenic_net_url,
            self.reciprocal_best_chain_url,
            self.reciprocal_best_net_url,
        )
        if self.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
            if any(url is None for url in comparative_urls):
                raise ValueError(
                    "COMPARATIVE resource bundle requires net, syntenic net, "
                    "and reciprocal-best chain/net URLs"
                )
        elif any(url is not None for url in comparative_urls):
            raise ValueError(
                "LIFTOVER_ONLY resource bundle cannot carry comparative URLs"
            )


class _HrefParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: set[str] = set()

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.casefold() != "a":
            return
        for name, value in attrs:
            if name.casefold() == "href" and value is not None:
                self.hrefs.add(value)


ListingReader = Callable[[str], frozenset[str] | None]


def discover_ucsc_resources(
    source_db: str, target_db: str
) -> UCSCResourceBundle | None:
    """Discover verified UCSC resources for one assembly direction.

    Return a ``COMPARATIVE`` bundle when the full v1 comparative set is
    published, otherwise fall back to the UCSC liftOver chain and return
    ``LIFTOVER_ONLY``.  Return ``None`` when the checked candidate listings are
    reachable but contain no usable resource for this assembly pair.

    Network/server failures raise ``UCSCResourceDiscoveryError`` instead of
    being interpreted as resource absence.
    """

    return _discover_ucsc_resources(source_db, target_db, _read_directory_links)


def _discover_ucsc_resources(
    source_db: str,
    target_db: str,
    read_listing: ListingReader,
) -> UCSCResourceBundle | None:
    _validate_ucsc_db(source_db)
    _validate_ucsc_db(target_db)

    target_title = _ucsc_title_db(target_db)
    comparative_base = urljoin(_UCSC_GOLDEN_PATH, f"{source_db}/vs{target_title}/")
    comparative_links = read_listing(comparative_base)
    if comparative_links is not None:
        chain_name = f"{source_db}.{target_db}.all.chain.gz"
        net_name = f"{source_db}.{target_db}.net.gz"
        syn_net_name = f"{source_db}.{target_db}.syn.net.gz"

        required_parent_entries = (chain_name, net_name, syn_net_name)
        if all(
            _listing_contains(comparative_base, comparative_links, entry)
            for entry in required_parent_entries
        ):
            reciprocal_base = _discover_reciprocal_best_base(
                source_db,
                target_db,
                comparative_base=comparative_base,
                comparative_links=comparative_links,
                read_listing=read_listing,
            )
            if reciprocal_base is not None:
                rbest_chain_name = f"{source_db}.{target_db}.rbest.chain.gz"
                rbest_net_name = f"{source_db}.{target_db}.rbest.net.gz"
                return UCSCResourceBundle(
                    source_db=source_db,
                    target_db=target_db,
                    evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
                    chain_url=urljoin(comparative_base, chain_name),
                    net_url=urljoin(comparative_base, net_name),
                    syntenic_net_url=urljoin(comparative_base, syn_net_name),
                    reciprocal_best_chain_url=urljoin(
                        reciprocal_base, rbest_chain_name
                    ),
                    reciprocal_best_net_url=urljoin(reciprocal_base, rbest_net_name),
                )

    # A partial comparative directory must not be mislabeled COMPARATIVE.  The
    # lightweight liftOver chain remains useful on its own, so check that
    # independently rather than failing discovery just because comparative
    # publication is incomplete.
    liftover_base = urljoin(_UCSC_GOLDEN_PATH, f"{source_db}/liftOver/")
    liftover_links = read_listing(liftover_base)
    liftover_name = f"{source_db}To{target_title}.over.chain.gz"
    if liftover_links is not None and _listing_contains(
        liftover_base, liftover_links, liftover_name
    ):
        return UCSCResourceBundle(
            source_db=source_db,
            target_db=target_db,
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
            chain_url=urljoin(liftover_base, liftover_name),
        )

    return None


def _discover_reciprocal_best_base(
    source_db: str,
    target_db: str,
    *,
    comparative_base: str,
    comparative_links: frozenset[str],
    read_listing: ListingReader,
) -> str | None:
    """Return the verified directory containing directional rbest files.

    UCSC sometimes publishes both directions of reciprocal-best output under one
    member of an assembly pair rather than under both ``source/vsTarget/`` trees.
    Check the forward comparison first, then the sibling/reverse comparison.  The
    fallback is accepted only when the exact ``source.target.rbest`` chain *and* net
    filenames are observed in the directory listing.  Which pair directory hosts a
    file is publication layout only; chain/net coordinate semantics come from the
    directional filename and file contents consumed by the parsers.

    This two-location search is based on measured UCSC publication behavior, not a
    provider guarantee.  Every candidate location is still verified by reading the
    live directory listing rather than treating the constructed path as existence.
    """

    reciprocal_dir = "reciprocalBest/"
    rbest_chain_name = f"{source_db}.{target_db}.rbest.chain.gz"
    rbest_net_name = f"{source_db}.{target_db}.rbest.net.gz"

    if _listing_contains(comparative_base, comparative_links, reciprocal_dir):
        forward_reciprocal_base = urljoin(comparative_base, reciprocal_dir)
        if _reciprocal_listing_has_directional_files(
            forward_reciprocal_base,
            read_listing(forward_reciprocal_base),
            rbest_chain_name=rbest_chain_name,
            rbest_net_name=rbest_net_name,
        ):
            return forward_reciprocal_base

    source_title = _ucsc_title_db(source_db)
    reverse_comparative_base = urljoin(
        _UCSC_GOLDEN_PATH, f"{target_db}/vs{source_title}/"
    )
    if reverse_comparative_base == comparative_base:
        return None

    reverse_links = read_listing(reverse_comparative_base)
    if reverse_links is None or not _listing_contains(
        reverse_comparative_base, reverse_links, reciprocal_dir
    ):
        return None

    reverse_reciprocal_base = urljoin(reverse_comparative_base, reciprocal_dir)
    if _reciprocal_listing_has_directional_files(
        reverse_reciprocal_base,
        read_listing(reverse_reciprocal_base),
        rbest_chain_name=rbest_chain_name,
        rbest_net_name=rbest_net_name,
    ):
        return reverse_reciprocal_base

    return None


def _reciprocal_listing_has_directional_files(
    reciprocal_base: str,
    links: frozenset[str] | None,
    *,
    rbest_chain_name: str,
    rbest_net_name: str,
) -> bool:
    if links is None:
        return False
    return all(
        _listing_contains(reciprocal_base, links, entry)
        for entry in (rbest_chain_name, rbest_net_name)
    )


def _read_directory_links(url: str) -> frozenset[str] | None:
    """Return exact href entries from one UCSC directory listing.

    A 404 means the candidate directory does not exist and is represented by
    ``None``.  Other HTTP or transport failures are not evidence of absence and
    therefore raise ``UCSCResourceDiscoveryError``.
    """

    request = Request(url, headers={"User-Agent": _USER_AGENT})
    try:
        with urlopen(request, timeout=15) as response:
            html = response.read().decode("utf-8", errors="replace")
    except HTTPError as exc:
        if exc.code == 404:
            return None
        raise UCSCResourceDiscoveryError(
            f"UCSC resource listing request failed with HTTP {exc.code}: {url}"
        ) from exc
    except URLError as exc:
        raise UCSCResourceDiscoveryError(
            f"UCSC resource listing request failed: {url}: {exc.reason}"
        ) from exc
    except TimeoutError as exc:
        # urllib may surface a socket/read timeout directly as TimeoutError rather
        # than wrapping it in URLError.  It is still a transport failure, not
        # evidence that the requested UCSC directory or resource is absent.
        raise UCSCResourceDiscoveryError(
            f"UCSC resource listing request timed out: {url}"
        ) from exc

    return _parse_directory_links(html)


def _listing_contains(base_url: str, links: frozenset[str], entry: str) -> bool:
    """Return whether a listing contains ``entry`` independent of href style.

    Apache-style indexes normally emit relative hrefs, but that presentation is
    not part of the scientific resource contract.  Resolve both observed and
    expected links against the listing URL so a future switch to absolute hrefs
    cannot silently make an existing resource look absent and downgrade the
    evidence-availability tier.
    """

    expected_url = urljoin(base_url, entry)
    return any(urljoin(base_url, href) == expected_url for href in links)


def _parse_directory_links(html: str) -> frozenset[str]:
    """Parse exact href values from a UCSC/Apache directory index."""

    parser = _HrefParser()
    parser.feed(html)
    return frozenset(parser.hrefs)


def _ucsc_title_db(db: str) -> str:
    """Return UCSC's ``ucfirst`` form used in vsTarget/liftOver names.

    Kent's current ``doRecipBest.pl`` sets ``$QDb = ucfirst($qDb)`` before
    constructing the Golden Path ``vs$QDb`` directory.  Keep that provider
    convention isolated here and still verify every resulting path by reading
    the actual directory listing.
    """

    _validate_ucsc_db(db)
    return db[0].upper() + db[1:]


def _validate_ucsc_db(db: str) -> None:
    if not db:
        raise ValueError("UCSC database identifier must not be empty")
    if fullmatch(_UCSC_DB_PATTERN, db) is None:
        raise ValueError(
            "UCSC database identifier may contain only letters, digits, '.', '_', and '-'"
        )
