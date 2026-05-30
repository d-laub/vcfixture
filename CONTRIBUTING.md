# Contributing to vcforge

## Requirements

vcforge uses [**uv**](https://docs.astral.sh/uv/) for environment management and
as its build backend. **uv is required** — the project is not set up for bare
`pip`/`venv` workflows, and the dev toolchain is pinned in `uv.lock`.

Install uv (see the [official guide](https://docs.astral.sh/uv/getting-started/installation/)):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

uv manages the Python interpreter too; you do not need a system Python matching
`requires-python` (`>=3.10`).

## Setup

```bash
uv sync            # create .venv and install the project + dev dependencies
uv run prek install  # install the git pre-commit and commit-msg hooks
```

## Everyday commands

Run tools through `uv run` so they use the locked environment:

```bash
uv run pytest                 # test suite
uv run ruff check             # lint
uv run ruff format            # format
uv run pyrefly check          # type-check the library (src/)
uv run prek run --all-files   # run all hooks against the whole tree
```

## Code quality gates

`prek` runs these hooks (config in `.pre-commit-config.yaml`); they also run in
CI and on every commit:

- **ruff** — linting and formatting.
- **pyrefly** — type checking. vcforge ships `py.typed` and is a type-safe
  library; `src/` must type-check cleanly. New code needs accurate annotations.
- **commitizen** — commit messages must follow
  [Conventional Commits](https://www.conventionalcommits.org/) (e.g.
  `feat: ...`, `fix: ...`, `docs: ...`, `test: ...`, `chore: ...`). The
  `commit-msg` hook rejects non-conforming messages.

## Workflow

1. Branch off `main`.
2. Make changes with tests (the suite round-trips generated VCFs through the
   independent `cyvcf2` parser — keep that oracle intact).
3. Ensure `uv run pytest`, `uv run ruff check`, and `uv run pyrefly check` pass.
4. Commit with a Conventional Commit message (hooks enforce this).
