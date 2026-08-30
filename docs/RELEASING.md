# Releasing liftAssess

This is the maintainer checklist for public liftAssess releases. The scientific and
validation requirements remain authoritative in [`DESIGN.md`](DESIGN.md) and
[`ROADMAP.md`](ROADMAP.md); this document makes the release mechanics reproducible.

## Before preparing a release

Do not prepare or tag a release while its roadmap validation gate is still open. For
`v0.2.0a1`, Milestone 23 must first be complete with any blocking outside-user/domain
feedback resolved.

Begin from a clean `main` branch that is up to date with `origin/main`.

## Synchronize release metadata

For the release-preparation commit:

1. set `[project].version` in `pyproject.toml`;
2. run `uv lock` so the editable project entry in `uv.lock` has the same version;
3. set `version` and `date-released` in `CITATION.cff`;
4. move the accumulated changelog entries from `Unreleased` under a versioned heading
   using the same release date, leaving a fresh `Unreleased` section for later work; and
5. update current-version/status wording in README and user documentation where needed,
   while preserving historical references to older releases.

Then verify the synchronized metadata:

```bash
uv run python scripts/check_release_metadata.py
```

The checker requires `pyproject.toml`, `uv.lock`, `CITATION.cff`, and the versioned
`CHANGELOG.md` heading to agree. During the tag-triggered release workflow it also
requires the Git tag to be exactly `v<package-version>`.

## Run the release-candidate gate

Run the complete native gate before building release artifacts:

```bash
uv run pytest

uv run ruff check \
  src tests scripts

uv run ruff format --check \
  src tests scripts

uv run mypy --strict \
  src tests scripts

git diff --check
```

## Build and review both distributions

Remove earlier local artifacts and build from the release-candidate tree:

```bash
rm -rf \
  dist \
  build \
  /tmp/liftassess-smoke-wheel \
  /tmp/liftassess-smoke-sdist

uv build
```

Smoke-test the wheel and source distribution independently. Derive the expected version
from the prepared package metadata rather than typing it again.

```bash
VERSION="$(uv run --no-project python -c 'import pathlib, tomllib; print(tomllib.loads(pathlib.Path("pyproject.toml").read_text())["project"]["version"])')"

uv venv \
  --python 3.13 \
  /tmp/liftassess-smoke-wheel

uv pip install \
  --python /tmp/liftassess-smoke-wheel/bin/python \
  dist/*.whl

/tmp/liftassess-smoke-wheel/bin/python \
  scripts/release_smoke.py \
  --expected-version "$VERSION"

uv venv \
  --python 3.13 \
  /tmp/liftassess-smoke-sdist

uv pip install \
  --python /tmp/liftassess-smoke-sdist/bin/python \
  dist/*.tar.gz

/tmp/liftassess-smoke-sdist/bin/python \
  scripts/release_smoke.py \
  --expected-version "$VERSION"
```

Review the contents and metadata of both artifacts before tagging. The tag workflow
publishes to PyPI after its release gate passes, so the tag should not be pushed until
these local artifacts and the release-preparation diff have been reviewed.

## Commit and tag

After the release-preparation commit is pushed and normal CI is green, create and push
an annotated version tag. Derive the tag version from the prepared package metadata.

```bash
VERSION="$(uv run --no-project python -c 'import pathlib, tomllib; print(tomllib.loads(pathlib.Path("pyproject.toml").read_text())["project"]["version"])')"

git tag \
  -a "v$VERSION" \
  -m "v$VERSION"

git push origin \
  "v$VERSION"
```

The tag starts the release workflow. That workflow re-runs the quality gate, verifies
release metadata, rebuilds and smoke-tests both distributions, publishes them to PyPI
through trusted publishing, and creates a draft GitHub release containing the artifacts.

Review the workflow in GitHub Actions. Do not create a second tag to work around a
failed release; diagnose the failure first.

## Smoke-test the published package

After PyPI publication, test the public package through the same pip-based installation
path documented for users. Derive the expected version from the release checkout.

```bash
rm -rf /tmp/liftassess-pypi-smoke

VERSION="$(uv run --no-project python -c 'import pathlib, tomllib; print(tomllib.loads(pathlib.Path("pyproject.toml").read_text())["project"]["version"])')"

uv venv \
  --python 3.13 \
  --seed \
  /tmp/liftassess-pypi-smoke

/tmp/liftassess-pypi-smoke/bin/python \
  -m pip install \
  --no-cache-dir \
  --pre \
  "liftassess==$VERSION"

/tmp/liftassess-pypi-smoke/bin/assess-liftover \
  --help

/tmp/liftassess-pypi-smoke/bin/python \
  scripts/release_smoke.py \
  --expected-version "$VERSION"
```

Confirm the expected version is visible on PyPI and review the automatically created
draft GitHub release. Use the finalized changelog section as the release-note basis,
then publish the draft through GitHub's web interface.

Finally, record the completed release state in the roadmap and verify that `main`, the
release tag, PyPI, and the published GitHub release all identify the same version.
