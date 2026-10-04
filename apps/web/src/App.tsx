import { Link, NavLink, Route, Routes } from "react-router";
import { api } from "./api";
import { AuthProvider, RequireLogin, useAuth } from "./auth";
import { Account } from "./pages/Account";
import { Admin } from "./pages/Admin";
import { JobDetail } from "./pages/JobDetail";
import { Jobs } from "./pages/Jobs";
import { Login } from "./pages/Login";
import { NewJob } from "./pages/NewJob";
import { Register } from "./pages/Register";
import { SampleResults } from "./pages/SampleResults";
import { Settings } from "./pages/Settings";
import { Welcome } from "./pages/Welcome";
import { ThemeProvider } from "./theme";
import { useApi } from "./useApi";

export function App() {
  return (
    <AuthProvider>
      <ThemeProvider>
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
                <Route path="/jobs/:jobId/samples/:sampleId" element={<SampleResults />} />
                <Route path="/settings" element={<Settings />} />
                <Route path="/account" element={<Account />} />
                <Route path="/admin" element={<Admin />} />
              </Route>
              <Route path="*" element={<p>Page not found.</p>} />
            </Routes>
          </main>
          <ServerVersion />
        </div>
      </ThemeProvider>
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
            {user.is_admin && <NavLink to="/admin">Admin</NavLink>}
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
  if (health.state !== "done") {
    return <footer>{health.state === "error" ? "server unreachable" : "connecting…"}</footer>;
  }
  const { version, core_version, python, platform, sqlite, libraries } = health.data;
  const fundamentals = [
    `ScaleHD server ${version}`,
    `ScaleHD core ${core_version}`
  ]
  const libraryinfo = [
    `Docker ${platform}`,
    `Python ${python}`,
    `SQLite ${sqlite}`,
    ...Object.entries(libraries).map(([name, v]) => `${name} ${v}`)
  ]
  const techFooter = (
    <footer>
      <span>{fundamentals.join(" · ")}</span>
      <span>{libraryinfo.join(" · ")}</span>
    </footer>
  );

  return techFooter;
}
