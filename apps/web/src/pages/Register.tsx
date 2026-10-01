import { type SubmitEvent, useState } from "react";
import { Link, Navigate } from "react-router";
import { api } from "../api";
import { useAuth } from "../auth";
import { useApi } from "../useApi";
import { Status } from "./Status";

export function Register() {
  const { user, loggedIn } = useAuth();
  const registration = useApi(api.registration);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  if (user) return <Navigate to="/" replace />;

  async function submit(event: SubmitEvent<HTMLFormElement>) {
    event.preventDefault();
    const form = new FormData(event.currentTarget);
    const password = String(form.get("password"));
    if (password !== form.get("confirm")) {
      setError("The passwords don't match.");
      return;
    }
    setBusy(true);
    setError(null);
    try {
      loggedIn(await api.register(String(form.get("username")), password));
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }

  return (
    <Status of={registration}>
      {({ open, first_account }) =>
        !open ? (
          <>
            <h1>Registration is closed</h1>
            <p>
              Ask the server's admin for an account, or <Link to="/login">log in</Link>.
            </p>
          </>
        ) : (
          <form className="form" onSubmit={submit}>
            <h1>{first_account ? "Welcome to ScaleHD" : "Create an account"}</h1>
            {first_account && (
              <p>This server has no accounts yet. The first one you create is the admin.</p>
            )}
            <label>
              Username
              <input
                name="username"
                autoComplete="username"
                autoFocus
                required
                maxLength={64}
                pattern="[A-Za-z0-9][A-Za-z0-9._\-]*"
                title="Letters, digits, dots, dashes and underscores, starting with a letter or digit"
              />
            </label>
            <label>
              <span>
                Password <span className="muted">(at least 8 characters)</span>
              </span>
              <input
                name="password"
                type="password"
                autoComplete="new-password"
                required
                minLength={8}
              />
            </label>
            <label>
              Password again
              <input name="confirm" type="password" autoComplete="new-password" required />
            </label>
            <button type="submit" disabled={busy}>
              Create account
            </button>
            {error && <p className="error">{error}</p>}
            {!first_account && (
              <p className="muted">
                Already have one? <Link to="/login">Log in.</Link>
              </p>
            )}
          </form>
        )
      }
    </Status>
  );
}
