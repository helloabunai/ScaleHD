import { type SubmitEvent, useState } from "react";
import { api, type JobSettings } from "../api";
import { useApi } from "../useApi";
import { MethodPicker } from "./MethodPicker";
import { Status } from "./Status";

type Outcome = { ok: boolean; message: string };

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
  const [outcome, setOutcome] = useState<Outcome | null>(null);
  const [busy, setBusy] = useState(false);

  async function submit(event: SubmitEvent<HTMLFormElement>) {
    event.preventDefault();
    setBusy(true);
    try {
      await api.saveSettings({ ...saved, method });
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
        onChange={(m) => {
          setMethod(m);
          setOutcome(null);
        }}
      />
      <button type="submit" disabled={busy}>
        Save
      </button>
      {outcome && <p className={outcome.ok ? "success" : "error"}>{outcome.message}</p>}
    </form>
  );
}
