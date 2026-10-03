import { useEffect, useState } from "react";

export type Loaded<T> =
  | { state: "loading" }
  | { state: "error"; error: Error }
  | { state: "done"; data: T };

export interface Refresh<T> {
  /** Milliseconds between refreshes. */
  every: number;
  /** Keep refreshing while this holds for the latest data. */
  while: (data: T) => boolean;
}

/** Load once per change of `deps`, then optionally refresh on a timer. */
export function useApi<T>(
  load: () => Promise<T>,
  deps: unknown[] = [],
  refresh?: Refresh<T>,
): Loaded<T> {
  const [result, setResult] = useState<Loaded<T>>({ state: "loading" });
  useEffect(() => {
    let current = true;
    let timer: ReturnType<typeof setTimeout> | undefined;
    let loadedOnce = false;
    const fetchOnce = () => {
      load().then(
        (data) => {
          if (!current) return;
          loadedOnce = true;
          setResult({ state: "done", data });
          if (refresh?.while(data)) timer = setTimeout(fetchOnce, refresh.every);
        },
        (error: Error) => {
          if (!current) return;
          // A failed refresh keeps the last data on screen and tries again.
          if (loadedOnce && refresh) {
            timer = setTimeout(fetchOnce, refresh.every);
            return;
          }
          setResult({ state: "error", error });
        },
      );
    };
    setResult({ state: "loading" });
    fetchOnce();
    return () => {
      current = false;
      clearTimeout(timer);
    };
  }, deps);
  return result;
}
