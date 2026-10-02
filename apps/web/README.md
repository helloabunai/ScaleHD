# scalehd-web

The ScaleHD web interface: React and TypeScript, built with Vite. Mostly a skeleton.
Logging in, registering (the first account on a new server is the admin), the welcome
page, the account page (change password) and default job settings (the genotyping
method) are present but mostly just placeholders.

Stylings etc are all placeholders. wip wip wip

```sh
npm install
npm run dev        # http://localhost:5173, forwards /api to the server on :8000
npm run build      # type-check, then write dist/ for the server to serve
```

Run the API alongside with `uv run scalehd-server --reload` from the repository root.
The types in `src/api.ts` mirror `apps/server/src/scalehd_server/schemas.py` by hand.

again wip wip wip