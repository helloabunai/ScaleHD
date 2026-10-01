import { Link, NavLink, Route, Routes } from "react-router";
import { api } from "./api";
import { AuthProvider, RequireLogin, useAuth } from "./auth";
import { Account } from "./pages/Account";
import { JobDetail } from "./pages/JobDetail";
import { Jobs } from "./pages/Jobs";
import { Login } from "./pages/Login";
import { NewJob } from "./pages/NewJob";
import { Register } from "./pages/Register";
import { Settings } from "./pages/Settings";
import { Welcome } from "./pages/Welcome";
import { useApi } from "./useApi";

export function App() {
  return (
    <AuthProvider>
      <div className="layout">
        <Header />
        <main>
          <Routes>
            <Route path="/login" element={<Login />} />
            <Route path="/register" element={<Register />} />
            <Route element={<RequireLogin />}>
              <Route path="/" element={<Welcome />} />
              <Route path="/jobs" element={<Jobs />} />
              <Route path="/jobs/new" element={<NewJob />} />
              <Route path="/jobs/:jobId" element={<JobDetail />} />
              <Route path="/settings" element={<Settings />} />
              <Route path="/account" element={<Account />} />
            </Route>
            <Route path="*" element={<p>Page not found.</p>} />
          </Routes>
        </main>
        <ServerVersion />
      </div>
    </AuthProvider>
  );
}

function Header() {
  const { user, logOut } = useAuth();
  return (
    <header>
      <Link to="/" className="brand">
        ScaleHD
      </Link>
      {user && (
        <>
          <nav>
            <NavLink to="/jobs" end>
              Jobs
            </NavLink>
            <NavLink to="/jobs/new">New job</NavLink>
            <NavLink to="/settings">Settings</NavLink>
          </nav>
          <div className="account">
            <NavLink to="/account">{user.username}</NavLink>
            <button type="button" className="link" onClick={logOut}>
              Log out
            </button>
          </div>
        </>
      )}
    </header>
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
