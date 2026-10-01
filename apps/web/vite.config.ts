import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";

// In development the API runs separately (`uv run scalehd-server --reload`) and Vite
// forwards /api to it. In production the server serves the built files itself.
export default defineConfig({
  plugins: [react()],
  server: { proxy: { "/api": "http://127.0.0.1:8000" } },
});
