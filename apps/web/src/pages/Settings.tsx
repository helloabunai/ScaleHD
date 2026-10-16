import { type SubmitEvent, useState } from "react";
import { api, type JobSettings, type ServerOptions } from "../api";
import { useApi } from "../useApi";
import { LegacySettingsForm } from "./LegacySettingsForm";
import { MethodPicker } from "./MethodPicker";
import { Status } from "./Status";

type Outcome = { ok: boolean; message: string };

export function Settings() {
  const settings = useApi(api.getSettings);
  const options = useApi(api.getOptions);
  return (
    <>
      <h1>Default job settings</h1>
      <p className="muted">What your new jobs start with. Each job can still change them.</p>
      <Status of={settings}>
        {(current) => (
          <Status of={options}>
            {(server) => <SettingsForm saved={current} options={server} />}
          </Status>
        )}
      </Status>
    </>
  );
}

// TODO: the flag thresholds too. Unset ones use the scalehd defaults, so show those as
// placeholders.
function SettingsForm({ saved, options }: Readonly<{ saved: JobSettings; options: ServerOptions }>) {
  const [method, setMethod] = useState(saved.method);
  const [legacy, setLegacy] = useState(saved.legacy ?? options.legacy_defaults);
  const [outcome, setOutcome] = useState<Outcome | null>(null);
  const [busy, setBusy] = useState(false);

  async function submit(event: SubmitEvent<HTMLFormElement>) {
    event.preventDefault();
    setBusy(true);
    try {
      await api.saveSettings({ ...saved, method, legacy });
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
        available={options.methods}
        onChange={(m) => {
          setMethod(m);
          setOutcome(null);
        }}
      />
      {(method === "legacy" || options.methods.includes("legacy")) && (
        <>
          {method !== "legacy" && (
            <p className="muted">
              Jobs you run with the legacy method pass these settings to the ScaleHD binary.
            </p>
          )}
          <LegacySettingsForm
            value={legacy}
            cores={options.cores}
            revertTo={options.legacy_defaults}
            revertLabel="Restore to legacy ScaleHD defaults"
            onChange={(changed) => {
              setLegacy(changed);
              setOutcome(null);
            }}
          />
        </>
      )}
      <button type="submit" disabled={busy}>
        Save
      </button>
      {outcome && <p className={outcome.ok ? "success" : "error"}>{outcome.message}</p>}
    </form>
  );
}
