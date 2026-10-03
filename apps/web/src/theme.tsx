import { createContext, type ReactNode, useContext, useEffect, useState } from "react";
import { api, type Theme } from "./api";
import { useAuth } from "./auth";

export type Shown = "light" | "dark";
const STORED = "scalehd.theme";
const DARK = "(prefers-color-scheme: dark)";

function resolve(theme: Theme): Shown {
  if (theme === "system") return matchMedia(DARK).matches ? "dark" : "light";
  return theme;
}

function apply(theme: Theme): Shown {
  const shown = resolve(theme);
  document.documentElement.dataset.theme = shown;
  try {
    localStorage.setItem(STORED, theme);
  } catch {
    // Storage blocked (e.g. a private window)
  }
  return shown;
}

function stored(): Theme {
  try {
    const value = localStorage.getItem(STORED);
    if (value === "system" || value === "light" || value === "dark") return value;
  } catch {
    // Storage blocked? fall back to system setting
  }
  return "system";
}

interface ThemeState {
  /** light, dark, or whatever the system is */
  theme: Theme;
  shown: Shown;
  choose: (theme: Theme) => Promise<void>;
}

const ThemeContext = createContext<ThemeState | null>(null);

export function ThemeProvider({ children }: { children: ReactNode }) {
  const { user, loggedIn } = useAuth();
  const [theme, setTheme] = useState<Theme>(stored);
  const [shown, setShown] = useState<Shown>(() => resolve(stored()));

  function show(next: Theme) {
    setShown(apply(next));
    setTheme(next);
  }

  // user choice setting
  useEffect(() => {
    if (user && user.theme !== theme) show(user.theme);
  }, [user]);

  // system setting
  useEffect(() => {
    if (theme !== "system") return;
    const media = matchMedia(DARK);
    const changed = () => setShown(apply("system"));
    media.addEventListener("change", changed);
    return () => media.removeEventListener("change", changed);
  }, [theme]);

  async function choose(next: Theme) {
    const before = theme;
    show(next);
    if (!user) return;
    try {
      loggedIn(await api.setTheme(next));
    } catch (e) {
      show(before);
      throw e;
    }
  }

  return <ThemeContext value={{ theme, shown, choose }}>{children}</ThemeContext>;
}

export function useTheme(): ThemeState {
  const theme = useContext(ThemeContext);
  if (!theme) throw new Error("useTheme needs a ThemeProvider");
  return theme;
}
