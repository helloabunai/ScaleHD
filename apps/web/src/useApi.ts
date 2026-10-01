import { useEffect, useState } from "react";

export type Loaded<T> =
  | { state: "loading" }
  | { state: "error"; error: Error }
  | { state: "done"; data: T };

/** Load once per change of `deps`. TODO: polling for running jobs. */
export function useApi<T>(load: () => Promise<T>, deps: unknown[] = []): Loaded<T> {
  const [result, setResult] = useState<Loaded<T>>({ state: "loading" });
  useEffect(() => {
    let current = true;
    setResult({ state: "loading" });
    load().then(
      (data) => current && setResult({ state: "done", data }),
      (error: Error) => current && setResult({ state: "error", error }),
    );
    return () => {
      current = false;
    };
  }, deps);
  return result;
}
