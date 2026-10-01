import { type FormEvent, useState } from "react";
import { useNavigate } from "react-router";
import { api } from "../api";

export function Login() {
  const navigate = useNavigate();
  const [error, setError] = useState<string | null>(null);

  async function submit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const form = new FormData(event.currentTarget);
    try {
      await api.login(String(form.get("username")), String(form.get("password")));
      navigate("/");
    } catch (e) {
      setError((e as Error).message);
    }
  }

  return (
    <form className="login" onSubmit={submit}>
      <h1>Log in</h1>
      <label>
        Username <input name="username" autoComplete="username" required />
      </label>
      <label>
        Password <input name="password" type="password" autoComplete="current-password" required />
      </label>
      <button type="submit">Log in</button>
      {error && <p className="error">{error}</p>}
    </form>
  );
}
