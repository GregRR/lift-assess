"""Verify synchronized liftAssess release metadata."""

from __future__ import annotations

import argparse
import re
import tomllib
from datetime import date
from pathlib import Path
from typing import cast

_PROJECT_NAME = "liftassess"
_RELEASE_HEADING = re.compile(r"^## (?P<version>\S+) - (?P<date>\d{4}-\d{2}-\d{2})$")


def _project_version(root: Path) -> str:
    payload = tomllib.loads((root / "pyproject.toml").read_text())
    project = cast(dict[str, object], payload["project"])
    version = project.get("version")
    if not isinstance(version, str):
        raise TypeError("pyproject.toml project.version is not a string")
    return version


def _lock_version(root: Path) -> str:
    payload = tomllib.loads((root / "uv.lock").read_text())
    packages = payload.get("package")
    if not isinstance(packages, list):
        raise TypeError("uv.lock package table is missing")

    matches: list[str] = []
    for raw_package in packages:
        if not isinstance(raw_package, dict):
            continue
        package = cast(dict[str, object], raw_package)
        if package.get("name") != _PROJECT_NAME:
            continue
        source = package.get("source")
        if not isinstance(source, dict):
            continue
        source_metadata = cast(dict[str, object], source)
        if source_metadata.get("editable") != ".":
            continue
        version = package.get("version")
        if not isinstance(version, str):
            raise TypeError("uv.lock liftassess version is not a string")
        matches.append(version)

    if len(matches) != 1:
        raise RuntimeError(
            "uv.lock must contain exactly one editable liftassess package entry"
        )
    return matches[0]


def _citation_scalar(root: Path, key: str) -> str:
    prefix = f"{key}:"
    matches: list[str] = []
    for line in (root / "CITATION.cff").read_text().splitlines():
        if not line.startswith(prefix):
            continue
        value = line.removeprefix(prefix).strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
            value = value[1:-1]
        matches.append(value)

    if len(matches) != 1 or not matches[0]:
        raise RuntimeError(f"CITATION.cff must contain exactly one top-level {key!r}")
    return matches[0]


def _changelog_release_date(root: Path, version: str) -> str:
    matches: list[str] = []
    for line in (root / "CHANGELOG.md").read_text().splitlines():
        match = _RELEASE_HEADING.fullmatch(line)
        if match is not None and match.group("version") == version:
            matches.append(match.group("date"))

    if len(matches) != 1:
        raise RuntimeError(
            f"CHANGELOG.md must contain exactly one release heading for {version}"
        )
    return matches[0]


def verify_release_metadata(
    root: Path,
    *,
    expected_version: str | None,
    expected_tag: str | None,
) -> None:
    project_version = _project_version(root)
    version = expected_version or project_version

    if project_version != version:
        raise RuntimeError(
            f"pyproject.toml version {project_version!r} != expected {version!r}"
        )

    lock_version = _lock_version(root)
    if lock_version != version:
        raise RuntimeError(f"uv.lock version {lock_version!r} != expected {version!r}")

    citation_version = _citation_scalar(root, "version")
    if citation_version != version:
        raise RuntimeError(
            f"CITATION.cff version {citation_version!r} != expected {version!r}"
        )

    citation_date = _citation_scalar(root, "date-released")
    date.fromisoformat(citation_date)
    changelog_date = _changelog_release_date(root, version)
    if changelog_date != citation_date:
        raise RuntimeError(
            "CHANGELOG.md release date does not match CITATION.cff: "
            f"{changelog_date!r} != {citation_date!r}"
        )

    if expected_tag is not None:
        required_tag = f"v{version}"
        if expected_tag != required_tag:
            raise RuntimeError(
                f"release tag {expected_tag!r} != expected {required_tag!r}"
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected-version")
    parser.add_argument("--expected-tag")
    args = parser.parse_args()

    expected_version = cast(str | None, args.expected_version)
    expected_tag = cast(str | None, args.expected_tag)
    verify_release_metadata(
        Path.cwd(),
        expected_version=expected_version,
        expected_tag=expected_tag,
    )

    version = expected_version or _project_version(Path.cwd())
    print(f"release metadata consistent for liftassess {version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
