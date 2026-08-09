"""Shared provenance-graph operations.

Provenance identity is structural evidence about dependence, not a declaration of
statistical independence. These helpers deliberately compare the source IDs that
callers supplied; they cannot detect a caller assigning the same ID to unrelated
real-world sources.
"""

from .models import ProvenanceSource


def shares_upstream_source(first: ProvenanceSource, second: ProvenanceSource) -> bool:
    """Return whether two provenance graphs contain at least one common source ID."""

    return bool(_provenance_source_ids(first) & _provenance_source_ids(second))


def _provenance_source_ids(source: ProvenanceSource) -> set[str]:
    """Collect unique source IDs without assuming the graph is a tree or acyclic."""

    source_ids: set[str] = set()
    pending = [source]

    while pending:
        current = pending.pop()
        if current.source_id in source_ids:
            continue
        source_ids.add(current.source_id)
        pending.extend(current.derived_from)

    return source_ids
