import { type SubmitEvent, useState } from "react";
import { api, type JobSettings } from "../api";
import { useApi } from "../useApi";
import { MethodPicker } from "./MethodPicker";
import { Status } from "./Status";

type Outcome = { ok: boolean; message: string };

// The GeM-HD MiSeq primers (Cell 2019), the server's default.
// likely to be corrected by real science people
const GEM_HD = { forward: "ATGAAGGCCTTCGAGTCCC", reverse: "GGCTGAGGAAGCTGAGGA" };

export function Settings() {
  const settings = useApi(api.getSettings);
  return (
    <>
      <h1>Default job settings</h1>
      <p className="muted">What your new jobs start with. Each job can still change them.</p>
      <Status of={settings}>{(current) => <SettingsForm saved={current} />}</Status>
    </>
  );
}

// TODO: the flag thresholds too. Unset ones use the scalehd defaults, so show those as
// placeholders.
function SettingsForm({ saved }: { saved: JobSettings }) {
  const [method, setMethod] = useState(saved.method);
  const [forward, setForward] = useState(saved.forward_primer);
  const [reverse, setReverse] = useState(saved.reverse_primer);
  const [trim, setTrim] = useState(saved.trim_primers);
  const [outcome, setOutcome] = useState<Outcome | null>(null);
  const [busy, setBusy] = useState(false);

  function changed<T>(set: (value: T) => void) {
    return (value: T) => {
      set(value);
      setOutcome(null);
    };
  }

  async function submit(event: SubmitEvent<HTMLFormElement>) {
    event.preventDefault();
    setBusy(true);
    try {
      const settings = { ...saved, method, forward_primer: forward, reverse_primer: reverse };
      const stored = await api.saveSettings({ ...settings, trim_primers: trim });
      setForward(stored.forward_primer);
      setReverse(stored.reverse_primer);
      setOutcome({ ok: true, message: "Saved." });
    } catch (e) {
      setOutcome({ ok: false, message: (e as Error).message });
    } finally {
      setBusy(false);
    }
  }

  return (
    <form className="form wide" onSubmit={submit}>
      <MethodPicker
        value={method}
        onChange={changed(setMethod)}
      />
      <fieldset className="methods primers">
        <legend>PCR primers</legend>
        <p className="muted">
          5' to 3', as ordered. A read's bases from a primer show the primer, not the sample,
          so SNPs are only called between the primers.
        </p>
        <label>
          Forward primer
          <input
            name="forward_primer"
            value={forward}
            onChange={(e) => changed(setForward)(e.target.value)}
            spellCheck={false}
            autoComplete="off"
            required
          />
        </label>
        <label>
          Reverse primer
          <input
            name="reverse_primer"
            value={reverse}
            onChange={(e) => changed(setReverse)(e.target.value)}
            spellCheck={false}
            autoComplete="off"
            required
          />
        </label>
        <label className="tick">
          <input
            type="checkbox"
            name="trim_primers"
            checked={trim}
            onChange={(e) => changed(setTrim)(e.target.checked)}
          />
          <span>
            Trim primers
            <span className="muted">
              Leave the primers' bases out of SNP calling. Reads aren't cut, they're placed
              by the sequence next to the repeats.
            </span>
          </span>
        </label>
        {(forward !== GEM_HD.forward || reverse !== GEM_HD.reverse) && (
          <button
            type="button"
            className="link"
            onClick={() => {
              changed(setForward)(GEM_HD.forward);
              setReverse(GEM_HD.reverse);
            }}
          >
            Use the GeM-HD primers
          </button>
        )}
      </fieldset>
      <button type="submit" disabled={busy}>
        Save
      </button>
      {outcome && <p className={outcome.ok ? "success" : "error"}>{outcome.message}</p>}
    </form>
  );
}
