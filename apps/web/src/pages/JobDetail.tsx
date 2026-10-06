import { useState } from "react";
import { Link, useNavigate, useParams } from "react-router";
import { api, type Job, type Sample, type SampleStatus as Stage } from "../api";
import { useUser } from "../auth";
import { useApi } from "../useApi";
import { DemoTag, formatTime, isActive, JobTags } from "./jobDisplay";
import { METHODS } from "./MethodPicker";
import { ServerFolders } from "./ServerFolders";
import { Status } from "./Status";
import { type PickedTag, TagPicker, tagIds } from "./TagPicker";

export function JobDetail() {
  const id = Number(useParams().jobId);
  const user = useUser();
  // Bumped to load the job again, e.g. after its tags change.
  const [version, setVersion] = useState(0);
  const job = useApi(() => api.getJob(id), [id, version], {
    every: 1000,
    while: (loaded) => isActive(loaded.status),
  });
  return (
    <Status of={job}>
      {(job) => (
        <>
          <h1>
            {job.name} {job.demo && <DemoTag />}
            <JobTags tags={job.tags} />
          </h1>
          <dl className="job-facts">
            <dt>Status</dt>
            <dd>{job.status}</dd>
            <dt>Method</dt>
            <dd>{METHODS[job.method].label}</dd>
            <dt>Tags</dt>
            <dd>
              {job.tags.length > 0 ? <JobTags tags={job.tags} /> : "none"}
              {user.is_admin && (
                <EditJobTags job={job} onSaved={() => setVersion((v) => v + 1)} />
              )}
            </dd>
            <dt>Samples</dt>
            <dd>{samplesDone(job.samples)}</dd>
            <dt>Started</dt>
            <dd>{job.started_at ? formatTime(job.started_at) : "–"}</dd>
            <dt>Finished</dt>
            <dd>{job.finished_at ? formatTime(job.finished_at) : "–"}</dd>
            {job.output_dir && (
              <>
                <dt>Results</dt>
                <dd>
                  <code>{job.output_dir}</code>
                </dd>
              </>
            )}
          </dl>
          <JobProgress samples={job.samples} />
          <div className="job-actions">
            <button type="button" disabled title="Coming soon">
              Export job results
            </button>
            <span className="muted">coming soon</span>
            {(job.status === "queued" || job.status === "running") && (
              <CancelJob job={job} onCancelled={() => setVersion((v) => v + 1)} />
            )}
            {job.status === "cancelling" && (
              <span className="muted">Cancelling.. waiting for samples already processing to finish.</span>
            )}
            {!isActive(job.status) && <DeleteJob job={job} />}
          </div>
          <table>
            <thead>
              <tr>
                <th>Sample</th>
                <th>Status</th>
                <th>Progress</th>
                <th>Genotype</th>
                <th>Confidence</th>
                <th>Flags</th>
                {hasTruth(job) && <th>Truth</th>}
                {hasTruth(job) && <th>Match</th>}
              </tr>
            </thead>
            <tbody>
              {job.samples.map((sample) => (
                <tr key={sample.id}>
                  <td>
                    {sample.status === "finished" || sample.status === "failed" ? (
                      <Link to={`/jobs/${job.id}/samples/${sample.id}`}>{sample.name}</Link>
                    ) : (
                      sample.name
                    )}
                  </td>
                  <td>
                    <SampleStatus sample={sample} />
                  </td>
                  <td>
                    <SampleProgress stage={sample.status} />
                  </td>
                  <td>{sample.genotype ?? "–"}</td>
                  <td>{sample.confidence?.toFixed(1) ?? "–"}</td>
                  <td>{sample.flags.join(", ")}</td>
                  {hasTruth(job) && <td>{sample.truth ?? "–"}</td>}
                  {hasTruth(job) && (
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

/** Simulated samples (the demo, or a generated run) know the genotype they were made from. 
 *  Real data jobs obviously won't have this element.
*/
function hasTruth(job: Job): boolean {
  return job.samples.some((sample) => sample.truth !== null);
}

/** Only an admin can change a job's tags once it's been made. */
function EditJobTags({ job, onSaved }: Readonly<{ job: Job; onSaved: () => void }>) {
  const [editing, setEditing] = useState(false);
  const [tags, setTags] = useState<PickedTag[]>(job.tags);
  const [error, setError] = useState<string | null>(null);

  if (!editing) {
    return (
      <button
        type="button"
        className="link edit-tags"
        onClick={() => {
          setTags(job.tags);
          setError(null);
          setEditing(true);
        }}
      >
        Edit tags
      </button>
    );
  }

  async function save() {
    try {
      await api.setJobTags(job.id, await tagIds(tags));
      setEditing(false);
      onSaved();
    } catch (e) {
      setError((e as Error).message);
    }
  }

  return (
    <div className="tag-editor">
      <TagPicker value={tags} onChange={setTags} />
      <div className="job-actions">
        <button type="button" onClick={() => void save()}>
          Save tags
        </button>
        <button type="button" className="link" onClick={() => setEditing(false)}>
          Cancel
        </button>
        {error && <span className="error">{error}</span>}
      </div>
    </div>
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

function CancelJob({ job, onCancelled }: Readonly<{ job: Job; onCancelled: () => void }>) {
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function cancel() {
    const running = count(job.samples, "running");
    const note =
      running > 0
        ? `\n\nSamples not started yet won't run. The ${running} running will finish gracefully.`
        : "\n\nNone of its samples have started, so none will run.";
    if (!window.confirm(`Cancel "${job.name}"?${note}`)) return;
    setBusy(true);
    try {
      await api.cancelJob(job.id);
      onCancelled();
    } catch (e) {
      setError((e as Error).message);
    }
    setBusy(false);
  }

  return (
    <>
      <button type="button" className="danger" onClick={() => void cancel()} disabled={busy}>
        Cancel job
      </button>
      {error && <span className="error">{error}</span>}
    </>
  );
}

function DeleteJob({ job }: { job: Job }) {
  const navigate = useNavigate();
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function remove() {
    const folder = job.output_dir ? `\n\nThis also deletes the folder:\n${job.output_dir}` : "";
    if (!window.confirm(`Delete "${job.name}"?${folder}`)) return;
    setBusy(true);
    try {
      await api.deleteJob(job.id);
      navigate("/jobs");
    } catch (e) {
      setError((e as Error).message);
      setBusy(false);
    }
  }

  return (
    <>
      <button type="button" className="danger" onClick={remove} disabled={busy}>
        Delete job
      </button>
      {error && <span className="error">{error}</span>}
    </>
  );
}

function count(samples: Sample[], stage: Stage): number {
  return samples.filter((sample) => sample.status === stage).length;
}

/** e.g. "4 of 9 done, 1 failed, 2 running, 2 cancelled". Failures are technically 'done'. */
function samplesDone(samples: Sample[]): string {
  const failed = count(samples, "failed");
  const running = count(samples, "running");
  const cancelled = count(samples, "cancelled");
  const done = count(samples, "finished") + failed;
  return [
    `${done} of ${samples.length} done`,
    ...(failed ? [`${failed} failed`] : []),
    ...(running ? [`${running} running`] : []),
    ...(cancelled ? [`${cancelled} cancelled`] : []),
  ].join(", ");
}

/** overall progress bar for job */
function JobProgress({ samples }: Readonly<{ samples: Sample[] }>) {
  const total = Math.max(samples.length, 1);
  const parts: Stage[] = ["finished", "failed", "running", "cancelled"];
  return (
    <div
      className="progress"
      role="progressbar"
      aria-label="samples done"
      aria-valuemin={0}
      aria-valuemax={samples.length}
      aria-valuenow={count(samples, "finished") + count(samples, "failed")}
    >
      {parts.map((stage) => {
        const n = count(samples, stage);
        return n ? (
          <span key={stage} className={`segment ${stage}`} style={{ width: `${(100 * n) / total}%` }} />
        ) : null;
      })}
    </div>
  );
}

// TODO: real progress for running samples (stage and share of reads counted) once
// jobs read real FASTQ files; until then a running sample shows a moving bar.
function SampleProgress({ stage }: Readonly<{ stage: Stage }>) {
  return (
    <div className="progress small" aria-label={stage}>
      {stage !== "queued" && <span className={`segment ${stage}`} style={{ width: "100%" }} />}
    </div>
  );
}
