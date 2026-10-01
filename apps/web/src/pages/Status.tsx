import type { ReactNode } from "react";
import type { Loaded } from "../useApi";

/** Loading and error states for a page's data, or the page once it has loaded. */
export function Status<T>({ of, children }: { of: Loaded<T>; children: (data: T) => ReactNode }) {
  if (of.state === "loading") return <p className="muted">Loading…</p>;
  if (of.state === "error") return <p className="error">{of.error.message}</p>;
  return <>{children(of.data)}</>;
}
