import { useParams } from "react-router";
import { api } from "../api";
import { useApi } from "../useApi";
import { Status } from "./Status";

// TODO: poll while the job runs, a cancel button, each sample's full call
// (alleles, stutter, alternatives) and its molecule-count plot.
export function JobDetail() {
  const id = Number(useParams().jobId);
  const job = useApi(() => api.getJob(id), [id]);
  return (
    <Status of={job}>
      {(job) => (
        <>
          <h1>{job.name}</h1>
          <p>
            {job.status} · {job.samples_done}/{job.sample_count} samples ·{" "}
            <a href={api.reportUrl(job.id)}>PDF report</a>
          </p>
          <table>
            <thead>
              <tr>
                <th>Sample</th>
                <th>Status</th>
                <th>Genotype</th>
                <th>Quality</th>
                <th>Flags</th>
              </tr>
            </thead>
            <tbody>
              {job.samples.map((sample) => (
                <tr key={sample.id}>
                  <td>{sample.name}</td>
                  <td>{sample.error ?? sample.status}</td>
                  <td>{sample.genotype ?? "–"}</td>
                  <td>{sample.quality?.toFixed(1) ?? "–"}</td>
                  <td>{sample.flags.join(", ")}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </>
      )}
    </Status>
  );
}
