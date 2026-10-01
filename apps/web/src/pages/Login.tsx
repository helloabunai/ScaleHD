import { type SubmitEvent, useState } from "react";
import { Link, Navigate, useLocation } from "react-router";
import { api } from "../api";
import { useAuth } from "../auth";
import { useApi } from "../useApi";

export function Login() {
  const { user, loggedIn } = useAuth();
  const from = (useLocation().state as { from?: string } | null)?.from ?? "/";
  const registration = useApi(api.registration);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  if (user) return <Navigate to={from} replace />;
  if (registration.state === "loading") return null;
  // A new server.. nobody can log in until the admin is created.
  if (registration.state === "done" && registration.data.first_account) {
    return <Navigate to="/register" replace />;
  }

  async function submit(event: SubmitEvent<HTMLFormElement>) {
    event.preventDefault();
    const form = new FormData(event.currentTarget);
    setBusy(true);
    setError(null);
    try {
      loggedIn(await api.login(String(form.get("username")), String(form.get("password"))));
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }

  return (
    <form className="form" onSubmit={submit}>
      <h1>Log in</h1>
      <label>
        Username
        <input name="username" autoComplete="username" autoFocus required />
      </label>
      <label>
        Password
        <input name="password" type="password" autoComplete="current-password" required />
      </label>
      <button type="submit" disabled={busy}>
        Log in
      </button>
      {error && <p className="error">{error}</p>}
      {registration.state === "done" && registration.data.open && (
        <p className="muted">
          No account? <Link to="/register">Create one.</Link>
        </p>
      )}
    </form>
  );
}
