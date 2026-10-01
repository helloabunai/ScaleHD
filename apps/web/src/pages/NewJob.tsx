import { api } from "../api";
import { useApi } from "../useApi";
import { Status } from "./Status";

// TODO: browse folders of the input directory, tick samples, name the job, adjust
// settings (prefilled from /api/settings), submit, then go to the job's page.
export function NewJob() {
  const inputs = useApi(api.listInputs);
  return (
    <>
      <h1>New job</h1>
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
    </>
  );
}
