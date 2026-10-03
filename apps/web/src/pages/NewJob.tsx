import { useState } from "react";
import { api, type JobSettings } from "../api";
import { useApi } from "../useApi";
import { MethodPicker } from "./MethodPicker";
import { Status } from "./Status";

export function NewJob() {
  const defaults = useApi(api.getSettings);
  return (
    <>
      <h1>New job</h1>
      <Status of={defaults}>{(settings) => <NewJobForm defaults={settings} />}</Status>
    </>
  );
}

// TODO: browse folders of the input directory, tick samples, name the job, submit with
// api.createJob({ name, samples, settings: { ...defaults, method } }), then go to the
// job's page.
function NewJobForm({ defaults }: { defaults: JobSettings }) {
  const [method, setMethod] = useState(defaults.method);
  const inputs = useApi(api.listInputs);
  return (
    <div className="form wide">
      <MethodPicker value={method} onChange={setMethod} defaultMethod={defaults.method} />
      <h2>Samples</h2>
      <Status of={inputs}>
        {(pairs) => (
          <ul>
            {pairs.map((pair) => (
              <li key={pair.r1}>
                {pair.name} <span className="muted">{pair.r2 ? "paired" : "R1 only"}</span>
              </li>
            ))}
          </ul>
        )}
      </Status>
    </div>
  );
}
