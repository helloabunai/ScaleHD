import { useParams } from "react-router";
import { api, type Sample } from "../api";
import { useApi } from "../useApi";
import { DemoTag, formatTime, isActive } from "./jobDisplay";
import { METHODS } from "./MethodPicker";
import { ServerFolders } from "./ServerFolders";
import { Status } from "./Status";

// TODO: a cancel button, and a page per sample with its full call (alleles, stutter,
// alternatives) and molecule-count plot, reachable as soon as that sample finishes.
// TODO: consider large samples may take time to generate results even if the job
// is finished processing/genotyping.
export function JobDetail() {
  const id = Number(useParams().jobId);
  const job = useApi(() => api.getJob(id), [id], {
    every: 1000,
    while: (loaded) => isActive(loaded.status),
  });
  return (
    <Status of={job}>
      {(job) => (
        <>
          <h1>
            {job.name} {job.demo && <DemoTag />}
          </h1>
          <p>
            {job.status} · {METHODS[job.method].label} · {job.samples_done}/{job.sample_count}{" "}
            samples
          </p>
          <p className="muted">
            {job.started_at && <>started {formatTime(job.started_at)}</>}
            {job.finished_at && <> · finished {formatTime(job.finished_at)}</>}
          </p>
          {job.output_dir && (
            <p className="muted">
              Results in <code>{job.output_dir}</code>
            </p>
          )}
          <table>
            <thead>
              <tr>
                <th>Sample</th>
                <th>Status</th>
                <th>Genotype</th>
                <th>Quality</th>
                <th>Flags</th>
                {job.demo && <th>Truth</th>}
                {job.demo && <th>Match</th>}
              </tr>
            </thead>
            <tbody>
              {job.samples.map((sample) => (
                <tr key={sample.id}>
                  <td>{sample.name}</td>
                  <td>
                    <SampleStatus sample={sample} />
                  </td>
                  <td>{sample.genotype ?? "–"}</td>
                  <td>{sample.quality?.toFixed(1) ?? "–"}</td>
                  <td>{sample.flags.join(", ")}</td>
                  {job.demo && <td>{sample.truth ?? "–"}</td>}
                  {job.demo && (
                    <td>
                      <Match matches={sample.matches_truth} />
                    </td>
                  )}
                </tr>
              ))}
            </tbody>
          </table>
          <ServerFolders />
        </>
      )}
    </Status>
  );
}

function SampleStatus({ sample }: { sample: Sample }) {
  if (sample.status === "failed") return <span className="error">{sample.error ?? "failed"}</span>;
  if (sample.status === "running") return <>running…</>;
  return <>{sample.status}</>;
}

function Match({ matches }: { matches: boolean | null }) {
  if (matches === null) return <>–</>;
  return matches ? <span className="success">✓</span> : <span className="error">✗</span>;
}
