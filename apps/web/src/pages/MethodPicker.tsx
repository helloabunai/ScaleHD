import type { GenotypeMethod } from "../api";

interface Method {
  label: string;
  description: string;
  /** Shown as a badge next to the label. */
  badge?: string;
  /** A caveat under the description. */
  note?: string;
  /** Why a server might not offer it. */
  unavailable: string;
}

// In the order the picker lists them.
export const METHODS: Record<GenotypeMethod, Method> = {
  model: {
    label: "New (model-based)",
    description:
      "Read the repeat structure straight from each read, then the model-based caller.",
    badge: "Beta",
    note: "Work in progress: calls and confidence values may change between versions.",
    unavailable: "Not on this server.",
  },
  legacy: {
    label: "Legacy ScaleHD",
    description:
      "Legacy ScaleHD itself. Align reads to its reference library, then call genotypes.",
    unavailable: "Not on this server: legacy ScaleHD is only installed in the Docker image.",
  },
};

export function MethodPicker({
  value,
  onChange,
  defaultMethod,
  available,
}: Readonly<{
  value: GenotypeMethod;
  onChange: (method: GenotypeMethod) => void;
  /** Marked "your default", e.g. when choosing for one job. */
  defaultMethod?: GenotypeMethod;
  available: GenotypeMethod[];
}>) {
  return (
    <fieldset className="methods">
      <legend>Genotyping method</legend>
      {(Object.entries(METHODS) as [GenotypeMethod, Method][]).map(([method, info]) => {
        const offered = available.includes(method);
        return (
          <label key={method} className={offered ? undefined : "unavailable"}>
            <input
              type="radio"
              name="method"
              value={method}
              checked={value === method}
              disabled={!offered}
              onChange={() => onChange(method)}
            />
            <span>
              <span className="method-name">
                {info.label}
                {info.badge && <span className="badge">{info.badge}</span>}
                {method === defaultMethod && <span className="muted">your default</span>}
              </span>
              <span className="muted">{info.description}</span>
              {info.note && <span className="note">{info.note}</span>}
              {!offered && <span className="note">{info.unavailable}</span>}
            </span>
          </label>
        );
      })}
    </fieldset>
  );
}
