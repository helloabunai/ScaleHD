import { type SubmitEvent, useState } from "react";
import { Link } from "react-router";
import { api } from "../api";
import { useUser } from "../auth";

type Outcome = { ok: boolean; message: string };

export function Account() {
  const user = useUser();
  const [outcome, setOutcome] = useState<Outcome | null>(null);
  const [busy, setBusy] = useState(false);

  async function submit(event: SubmitEvent<HTMLFormElement>) {
    event.preventDefault();
    const formElement = event.currentTarget;
    const form = new FormData(formElement);
    const password = String(form.get("password"));
    if (password !== form.get("confirm")) {
      setOutcome({ ok: false, message: "The new passwords don't match." });
      return;
    }
    setBusy(true);
    try {
      await api.changePassword(String(form.get("current")), password);
      formElement.reset();
      setOutcome({
        ok: true,
        message: "Password changed. Other browsers using this account have been logged out.",
      });
    } catch (e) {
      setOutcome({ ok: false, message: (e as Error).message });
    } finally {
      setBusy(false);
    }
  }

  return (
    <>
      <h1>{user.username}</h1>
      <p className="muted">
        {user.is_admin ? "Admin" : "User"} · account created{" "}
        {new Date(user.created_at).toLocaleDateString()}
      </p>
      <p>
        <Link to="/settings">Default job settings</Link>: the genotyping method and thresholds
        your new jobs start with.
      </p>
      <form className="form" onSubmit={submit}>
        <h2>Change password</h2>
        {/* For password managers: which account this password belongs to. */}
        <input type="hidden" name="username" autoComplete="username" defaultValue={user.username} />
        <label>
          Current password
          <input name="current" type="password" autoComplete="current-password" required />
        </label>
        <label>
          <span>
            New password <span className="muted">(at least 8 characters)</span>
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
          New password again
          <input name="confirm" type="password" autoComplete="new-password" required />
        </label>
        <button type="submit" disabled={busy}>
          Change password
        </button>
        {outcome && <p className={outcome.ok ? "success" : "error"}>{outcome.message}</p>}
      </form>
    </>
  );
}
