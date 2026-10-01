import { NavLink, Route, Routes } from "react-router";
import { api } from "./api";
import { JobDetail } from "./pages/JobDetail";
import { Jobs } from "./pages/Jobs";
import { Login } from "./pages/Login";
import { NewJob } from "./pages/NewJob";
import { Settings } from "./pages/Settings";
import { useApi } from "./useApi";

// TODO: send the user to /login when /api/auth/me says 401, once accounts exist.
export function App() {
  return (
    <div className="layout">
      <header>
        <strong>ScaleHD</strong>
        <nav>
          <NavLink to="/" end>
            Jobs
          </NavLink>
          <NavLink to="/jobs/new">New job</NavLink>
          <NavLink to="/settings">Settings</NavLink>
        </nav>
      </header>
      <main>
        <Routes>
          <Route path="/" element={<Jobs />} />
          <Route path="/jobs/new" element={<NewJob />} />
          <Route path="/jobs/:jobId" element={<JobDetail />} />
          <Route path="/settings" element={<Settings />} />
          <Route path="/login" element={<Login />} />
          <Route path="*" element={<p>Page not found.</p>} />
        </Routes>
      </main>
      <ServerVersion />
    </div>
  );
}

function ServerVersion() {
  const health = useApi(api.health);
  return (
    <footer>
      {health.state === "done"
        ? `server ${health.data.version} · core ${health.data.core_version}`
        : health.state === "error"
          ? "server unreachable"
          : "connecting…"}
    </footer>
  );
}
