# Development

Master branch was renamed to main as is standard these days. Main branch is now protected, so submit PRs if needed.

within `tools/` subdir is a mini script to refresh your docker dev stack. cos any dev is probably lazy

```sh
uv run pytest                        # tests
uv run ruff check packages apps      # lint
uv run ruff format packages apps     # format
uv run mypy                          # types
(cd apps/web && npm ci)              # frontend dependencies (initial, and after they change)
(cd apps/web && npm run build)       # frontend type-check and build
uv run python packages/core/benchmarks/parse_accuracy.py
uv run python packages/core/benchmarks/genotype_simulated.py
uv run python packages/core/benchmarks/legacy_matrix.py
uv run python packages/core/benchmarks/stutter_robustness.py    # stutter unlike the priors
uv run python packages/core/benchmarks/posterior_calibration.py # is a 0.99 call right 99%?
uv run python packages/core/benchmarks/prior_crossval.py        # priors from the other half
```

`tools/check.sh` runs every check above in order (lint, format, types, frontend build, tests) and stops at the first failure. A git pre-push hook runs it before each push and
blocks the push if anything fails. Turn the hook on in your local clone:

```sh
git config core.hooksPath tools/git-hooks
```

`git push --no-verify` skips it for one push.  Maybe don't do that.

`tools/dev.sh` runs the server and web interface for development purposes i.e auto-reloading when pages/files are updated. Uses a separate dev db.


## Layout

```
packages/core/        scalehd: pure-Python library and CLI (no web or DB dependencies)
  src/scalehd/        structure, amplicon, parse, pairs, counts, calibration,
                      simulate, seqio, cli
    genotype/         model/ (the model-based caller), legacy/ (ScaleHD 1.x, not yet implemented)
  tests/              unit, property-based and end-to-end tests
  benchmarks/         accuracy on simulated data and the legacy labelled matrix, and how
                      the model holds up when its assumptions are off (see MODEL.md)
apps/server/          scalehd-server: FastAPI, job runner, SQLite database
apps/web/             web interface: React, TypeScript, Vite
tools/                check.sh, dev.sh, refresh.sh and the pre-push git hook
Dockerfile            one image with the server, the core and the built frontend
compose.yaml          runs the docker image with a database volume, the data (input) folder read-only and
                      the workspace (job results) read-write
legacy/               ScaleHD 1.x, for reference only
```
