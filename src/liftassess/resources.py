"""Discover UCSC-hosted evidence resources without downloading them.

The UCSC Golden Path directory layout is useful for discovery, but liftAssess
must not treat a constructed URL as evidence that a resource exists.  This
module therefore reads the published directory listings and only returns URLs
corresponding to entries that were actually observed.  Directory links are
resolved against the listing URL before comparison so discovery does not depend
on whether UCSC renders an entry as a relative or absolute ``href``.

Implementation references (checked 2026-08-10):

- UCSC Golden Path canFam3/canFam4 comparative directory and README:
  https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/
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

    def handle_starttag(
        self, tag: str, attrs: list[tuple[str, str | None]]
    ) -> None:
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
    ``LIFTOVER_ONLY``.  Return ``None`` when both checked listings are reachable
    but contain no usable resource for this assembly pair.

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
    comparative_base = urljoin(
        _UCSC_GOLDEN_PATH, f"{source_db}/vs{target_title}/"
    )
    comparative_links = read_listing(comparative_base)
    if comparative_links is not None:
        chain_name = f"{source_db}.{target_db}.all.chain.gz"
        net_name = f"{source_db}.{target_db}.net.gz"
        syn_net_name = f"{source_db}.{target_db}.syn.net.gz"
        reciprocal_dir = "reciprocalBest/"

        required_parent_entries = (
            chain_name,
            net_name,
            syn_net_name,
            reciprocal_dir,
        )
        if all(
            _listing_contains(comparative_base, comparative_links, entry)
            for entry in required_parent_entries
        ):
            reciprocal_base = urljoin(comparative_base, reciprocal_dir)
            reciprocal_links = read_listing(reciprocal_base)
            if reciprocal_links is not None:
                rbest_chain_name = f"{source_db}.{target_db}.rbest.chain.gz"
                rbest_net_name = f"{source_db}.{target_db}.rbest.net.gz"
                if all(
                    _listing_contains(reciprocal_base, reciprocal_links, entry)
                    for entry in (rbest_chain_name, rbest_net_name)
                ):
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
                        reciprocal_best_net_url=urljoin(
                            reciprocal_base, rbest_net_name
                        ),
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
