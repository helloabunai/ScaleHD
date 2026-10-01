import { api } from "../api";
import { useApi } from "../useApi";
import { Status } from "./Status";

// TODO: a form for these, saved with api.saveSettings. Unset thresholds use the
// scalehd defaults, so show those as placeholders.
export function Settings() {
  const settings = useApi(api.getSettings);
  return (
    <>
      <h1>Default job settings</h1>
      <Status of={settings}>{(current) => <pre>{JSON.stringify(current, null, 2)}</pre>}</Status>
    </>
  );
}
