import type { GenotypeMethod } from "../api";

interface Method {
  label: string;
  description: string;
  /** Shown as a badge next to the label. */
  badge?: string;
  /** A caveat under the description. */
  note?: string;
}

export const METHODS: Record<GenotypeMethod, Method> = {
  legacy: {
    label: "Legacy (ScaleHD 1.x)",
    description: "Align reads to the ScaleHD 1.x reference library, then the original genotyper.",
    // TODO: remove once the 1.x genotyper (and alignment) is brought over; see
    // RUNNABLE_METHODS in apps/server/src/scalehd_server/runner.py.
    note: "Not available yet: still being brought over from ScaleHD 1.x.",
  },
  model: {
    label: "New (model-based)",
    description:
      "Read the repeat structure straight from each read, then the model-based caller.",
    badge: "Beta",
    note: "Work in progress: calls and confidence values may change between versions.",
  },
};

export function MethodPicker({
  value,
  onChange,
  defaultMethod,
}: {
  value: GenotypeMethod;
  onChange: (method: GenotypeMethod) => void;
  /** Marked "your default", e.g. when choosing for one job. */
  defaultMethod?: GenotypeMethod;
}) {
  return (
    <fieldset className="methods">
      <legend>Genotyping method</legend>
      {(Object.entries(METHODS) as [GenotypeMethod, Method][]).map(([method, info]) => (
        <label key={method}>
          <input
            type="radio"
            name="method"
            value={method}
            checked={value === method}
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
          </span>
        </label>
      ))}
    </fieldset>
  );
}
