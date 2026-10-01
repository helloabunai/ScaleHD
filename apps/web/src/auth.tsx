import { createContext, type ReactNode, useContext, useEffect, useState } from "react";
import { Navigate, Outlet, useLocation } from "react-router";
import { api, onUnauthorized, type User } from "./api";

interface Auth {
  /** undefined while checking for an existing login, null when logged out. */
  user: User | null | undefined;
  loggedIn: (user: User) => void;
  logOut: () => Promise<void>;
}

const AuthContext = createContext<Auth | null>(null);

export function AuthProvider({ children }: { children: ReactNode }) {
  const [user, setUser] = useState<User | null | undefined>(undefined);
  useEffect(() => {
    onUnauthorized(() => setUser(null));
    api.me().then(setUser, () => setUser(null));
  }, []);

  async function logOut() {
    await api.logout().catch(() => undefined);
    setUser(null);
  }

  return <AuthContext value={{ user, loggedIn: setUser, logOut }}>{children}</AuthContext>;
}

export function useAuth(): Auth {
  const auth = useContext(AuthContext);
  if (!auth) throw new Error("useAuth needs an AuthProvider");
  return auth;
}

/** The logged-in user, in pages under RequireLogin. */
export function useUser(): User {
  const { user } = useAuth();
  if (!user) throw new Error("useUser needs RequireLogin");
  return user;
}

/** Routes for logged-in users. Anyone else logs in first, then comes back. */
export function RequireLogin() {
  const { user } = useAuth();
  const location = useLocation();
  if (user === undefined) return <p className="muted">Loading…</p>;
  if (user === null) {
    return <Navigate to="/login" replace state={{ from: location.pathname + location.search }} />;
  }
  return <Outlet />;
}
