# syntax=docker/dockerfile:1
#
# ScaleHD: the API, job runner and web interface in one image.
#   docker compose up --build        (see compose.yaml)
#
# No external tools are needed: reads are parsed directly, without an aligner.
# Anything a later stage needs (e.g. for PDF reports) gets installed in the last stage.

ARG PYTHON=3.14-slim
ARG NODE=24-slim
ARG UV=0.12.21

# The web frontend, built to static files.
FROM node:${NODE} AS web
WORKDIR /web
COPY apps/web/package.json apps/web/package-lock.json ./
RUN --mount=type=cache,target=/root/.npm npm ci
COPY apps/web/ ./
RUN npm run build

FROM ghcr.io/astral-sh/uv:${UV} AS uv

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
RUN useradd --create-home --uid 1000 scalehd && mkdir /data /input && chown scalehd /data
COPY --from=python /app/.venv /app/.venv
COPY --from=web /web/dist /app/web
ENV PATH=/app/.venv/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    SCALEHD_DATA_DIR=/data \
    SCALEHD_INPUT_DIR=/input \
    SCALEHD_WEB_DIR=/app/web
USER scalehd
VOLUME /data
EXPOSE 8000
HEALTHCHECK CMD ["python", "-c", "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/api/health')"]
CMD ["scalehd-server", "--host", "0.0.0.0", "--port", "8000"]
