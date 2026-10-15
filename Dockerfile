# syntax=docker/dockerfile:1
#
# ScaleHD: the API, job runner and web interface in one image.
#   docker compose up --build        (see compose.yaml)
#
# The new genotyping parses reads directly and needs no external tools. Legacy genotyping
# (ScaleHD 1.x) needs its own Python and the aligner and other programs it uses. Those
# are built in their own stage into /opt/legacy, kept off the server's $PATH.

ARG PYTHON=3.14-slim
ARG NODE=24-slim
ARG UV=0.12.21
ARG MICROMAMBA=1.5.10

# The web frontend, built to static files.
FROM node:${NODE} AS web
WORKDIR /web
COPY apps/web/package.json apps/web/package-lock.json ./
RUN --mount=type=cache,target=/root/.npm npm ci
COPY apps/web/ ./
RUN npm run build

FROM ghcr.io/astral-sh/uv:${UV} AS uv

# ScaleHD 1.x, as maintained upstream (legacy/), with Python 3.8 in /opt/legacy/python and
# the programs it uses in /opt/legacy/tools. Versions are pinned in apps/legacy-runner/.
FROM mambaorg/micromamba:${MICROMAMBA} AS legacy
USER root
COPY apps/legacy-runner/environment-tools.yml apps/legacy-runner/environment-python.yml /tmp/env/
RUN micromamba create -y -p /opt/legacy/tools -f /tmp/env/environment-tools.yml \
 && micromamba create -y -p /opt/legacy/python -f /tmp/env/environment-python.yml \
 && micromamba create -y -p /tmp/build -c conda-forge python=3.8.20 pip c-compiler \
 && micromamba clean -a -y
# cutadapt 1.18 has no wheel for Python 3.8, compile manually here
RUN micromamba run -p /tmp/build pip wheel --no-cache-dir --no-deps cutadapt==1.18 -w /tmp/wheels
COPY apps/legacy-runner/requirements.txt /tmp/env/
RUN /opt/legacy/python/bin/pip install --no-cache-dir --no-deps --find-links /tmp/wheels \
    -r /tmp/env/requirements.txt
COPY legacy /tmp/legacy
RUN /opt/legacy/python/bin/pip install --no-cache-dir --no-deps /tmp/legacy \
 && mkdir -p /opt/legacy/references \
 && cp /tmp/legacy/ScaleHD/config/*.fasta /opt/legacy/references/ \
 && cp /tmp/legacy/UPSTREAM /opt/legacy/UPSTREAM
COPY apps/legacy-runner/picard /opt/legacy/bin/picard

# The Python environment: scalehd and scalehd-server, installed into /app/.venv.
FROM python:${PYTHON} AS python
COPY --from=uv /uv /bin/uv
ENV UV_COMPILE_BYTECODE=1 UV_LINK_MODE=copy UV_PYTHON_DOWNLOADS=never
WORKDIR /app
# Third-party dependencies first, so code changes don't reinstall them.
COPY pyproject.toml uv.lock ./
COPY packages/core/pyproject.toml packages/core/README.md packages/core/
COPY apps/server/pyproject.toml apps/server/README.md apps/server/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --no-install-workspace --package scalehd-server
COPY packages/core/ packages/core/
COPY apps/server/ apps/server/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --no-editable --package scalehd-server

FROM python:${PYTHON}
RUN useradd --create-home --uid 1000 scalehd && mkdir /data && chown scalehd /data
COPY --from=python /app/.venv /app/.venv
COPY --from=web /web/dist /app/web
COPY --from=legacy /opt/legacy /opt/legacy
# wrapper to run scalehd 1.x per-sample and report back to new web frontend
COPY --chmod=755 apps/legacy-runner/scalehd_legacy_run.py /opt/legacy/bin/scalehd-legacy-run
ENV PATH=/app/.venv/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    SCALEHD_DATABASE_DIR=/data \
    SCALEHD_WEB_DIR=/app/web
USER scalehd
VOLUME /data
EXPOSE 8000
HEALTHCHECK CMD ["python", "-c", "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/api/health')"]
CMD ["scalehd-server", "--host", "0.0.0.0", "--port", "8000"]
