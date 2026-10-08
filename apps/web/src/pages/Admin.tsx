import { type KeyboardEvent, type SubmitEvent, useState } from "react";
import { api, MAX_TAG_LENGTH, type Tag, type User } from "../api";
import { useAuth, useUser } from "../auth";
import { useApi } from "../useApi";
import { formatTime } from "./jobDisplay";
import { Status } from "./Status";

export function Admin() {
  const user = useUser();
  if (!user.is_admin) return <p>Only an admin can see this page.</p>;
  return (
    <>
      <h1>Admin</h1>
      <Tags />
      <Users />
    </>
  );
}

function Tags() {
  // Bumped to load the list again after a change.
  const [version, setVersion] = useState(0);
  const tags = useApi(api.listTags, [version]);
  return (
    <section>
      <h2>Tags</h2>
      <p className="muted">
        Shared by every user. Anyone can make one when starting a job but only an admin can
        rename or delete one, or change a job's tags afterwards (to prevent user 'conflicts').
      </p>
      <Status of={tags}>
        {(list) =>
          list.length === 0 ? (
            <p className="muted">No tags yet.</p>
          ) : (
            <table className="admin-table">
              <thead>
                <tr>
                  <th>Tag</th>
                  <th>Jobs</th>
                  <th aria-label="Actions" />
                </tr>
              </thead>
              <tbody>
                {list.map((tag) => (
                  <TagRow key={tag.id} tag={tag} onChanged={() => setVersion((v) => v + 1)} />
                ))}
              </tbody>
            </table>
          )
        }
      </Status>
    </section>
  );
}

function TagRow({ tag, onChanged }: Readonly<{ tag: Tag; onChanged: () => void }>) {
  const [name, setName] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  async function rename() {
    if (name === null || !name.trim()) return;
    try {
      await api.renameTag(tag.id, name);
      setName(null);
      onChanged();
    } catch (e) {
      setError((e as Error).message);
    }
  }

  async function remove() {
    const jobs = tag.jobs === 1 ? "1 job" : `${tag.jobs} jobs`;
    const note = tag.jobs > 0 ? ` It comes off ${jobs}; the jobs themselves stay.` : "";
    if (!window.confirm(`Delete the tag "${tag.name}"?${note}`)) return;
    try {
      await api.deleteTag(tag.id);
      onChanged();
    } catch (e) {
      setError((e as Error).message);
    }
  }

  function onKeyDown(event: KeyboardEvent<HTMLInputElement>) {
    if (event.key === "Enter") void rename();
    if (event.key === "Escape") setName(null);
  }

  return (
    <tr>
      <td>
        {name === null ? (
          <span className="tag">{tag.name}</span>
        ) : (
          <input
            aria-label={`New name for ${tag.name}`}
            value={name}
            maxLength={MAX_TAG_LENGTH}
            onChange={(e) => setName(e.target.value)}
            onKeyDown={onKeyDown}
          />
        )}
        {error && <span className="error"> {error}</span>}
      </td>
      <td>{tag.jobs}</td>
      <td className="row-actions">
        {name === null ? (
          <>
            <button
              type="button"
              className="link"
              onClick={() => {
                setName(tag.name);
                setError(null);
              }}
            >
              Rename
            </button>
            <button type="button" className="link danger" onClick={() => void remove()}>
              Delete
            </button>
          </>
        ) : (
          <>
            <button type="button" onClick={() => void rename()} disabled={!name.trim()}>
              Save
            </button>
            <button type="button" className="link" onClick={() => setName(null)}>
              Cancel
            </button>
          </>
        )}
      </td>
    </tr>
  );
}

function Users() {
  const me = useUser();
  const { loggedIn } = useAuth();
  const [version, setVersion] = useState(0);
  const users = useApi(api.listUsers, [version]);
  const [error, setError] = useState<string | null>(null);
  const [notice, setNotice] = useState<string | null>(null);

  async function toggle(user: User) {
    const self = user.id === me.id;
    if (self && !window.confirm("Stop being an admin? You'll lose access to this page.")) return;
    try {
      const changed = await api.setAdmin(user.id, !user.is_admin);
      setError(null);
      if (self) loggedIn(changed);
      else setVersion((v) => v + 1);
    } catch (e) {
      setError((e as Error).message);
    }
  }

  return (
    <section>
      <h2>Users</h2>
      <p className="muted">
        Admins manage tags and users. You should make at least one other admin, 
        so the server never depends on a single account.
      </p>
      {error && <p className="error">{error}</p>}
      {notice && <p className="success">{notice}</p>}
      <Status of={users}>
        {(list) => (
          <table className="admin-table">
            <thead>
              <tr>
                <th>User</th>
                <th>Account made</th>
                <th>Role</th>
                <th aria-label="Actions" />
              </tr>
            </thead>
            <tbody>
              {list.map((user) => (
                <tr key={user.id}>
                  <td>
                    {user.username}
                    {user.id === me.id && <span className="muted"> (you)</span>}
                  </td>
                  <td>{formatTime(user.created_at)}</td>
                  <td>{user.is_admin ? "Admin" : "User"}</td>
                  <td className="row-actions">
                    <button type="button" className="link" onClick={() => void toggle(user)}>
                      {user.is_admin ? "Remove admin" : "Make admin"}
                    </button>
                    {user.id !== me.id && (
                      <SetPassword
                        user={user}
                        onSet={() => {
                          setError(null);
                          setNotice(
                            `Password reset for ${user.username}.`,
                          );
                        }}
                      />
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        )}
      </Status>
    </section>
  );
}

function SetPassword({ user, onSet }: Readonly<{ user: User; onSet: () => void }>) {
  const [password, setPassword] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function submit(event: SubmitEvent<HTMLFormElement>) {
    event.preventDefault();
    if (password === null) return;
    setBusy(true);
    try {
      await api.setPassword(user.id, password);
      setPassword(null);
      onSet();
    } catch (e) {
      setError((e as Error).message);
    }
    setBusy(false);
  }

  if (password === null) {
    return (
      <button
        type="button"
        className="link"
        onClick={() => {
          setPassword("");
          setError(null);
        }}
      >
        Set password
      </button>
    );
  }
  return (
    <form className="set-password" onSubmit={(e) => void submit(e)}>
      <input
        type="password"
        aria-label={`New password for ${user.username}`}
        placeholder="New password (8 or more)"
        autoComplete="new-password"
        minLength={8}
        maxLength={256}
        required
        value={password}
        onChange={(e) => setPassword(e.target.value)}
        onKeyDown={(e) => e.key === "Escape" && setPassword(null)}
      />
      <button type="submit" disabled={busy || password.length < 8}>
        Set
      </button>
      <button type="button" className="link" onClick={() => setPassword(null)}>
        Cancel
      </button>
      {error && <span className="error">{error}</span>}
    </form>
  );
}
