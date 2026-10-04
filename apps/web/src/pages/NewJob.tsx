import { type SubmitEvent, useEffect, useId, useState } from "react";
import { useLocation, useNavigate } from "react-router";
import { api, type InputFolder, type InputSample, type JobSettings } from "../api";
import { useApi } from "../useApi";
import { MethodPicker } from "./MethodPicker";
import { Status } from "./Status";
import { type PickedTag, TagPicker, tagIds } from "./TagPicker";

export function NewJob() {
  const defaults = useApi(api.getSettings);
  // Each visit is a new location, so the New job link starts the page afresh even from here.
  const { key } = useLocation();
  return (
    <>
      <h1>New job</h1>
      <p className="muted">
        Pick a data subdirectory from the server's data root, then the samples to genotype.
      </p>
      <Status of={defaults}>{(settings) => <NewJobForm key={key} defaults={settings} />}</Status>
    </>
  );
}

function NewJobForm({ defaults }: Readonly<{ defaults: JobSettings }>) {
  const [folder, setFolder] = useState("");
  const listing = useApi(() => api.listInputs(folder), [folder]);
  return (
    <div className="new-job-layout">
      <FolderTree open={folder} onOpen={setFolder} />
      <div>
        <Status of={listing}>
          {(found) =>
            found.samples.length > 0 ? (
              <JobForm key={found.folder} found={found} defaults={defaults} />
            ) : (
              <p className="muted">
                No FASTQ files in {found.path}. <br/>
                Pick a run's folder from the folder tree.
              </p>
            )
          }
        </Status>
      </div>
    </div>
  );
}

/** The data folder's folders as a tree. */
function FolderTree({ open, onOpen }: Readonly<{ open: string; onOpen: (f: string) => void }>) {
  return (
    <nav className="folder-tree" aria-label="Data folders">
      <ul>
        <FolderNode
          folder=""
          name="ScaleHD Data Root"
          samples={0}
          hasFolders
          open={open}
          onOpen={onOpen}
        />
      </ul>
    </nav>
  );
}

function FolderNode({
  folder,
  name,
  samples,
  hasFolders,
  open,
  onOpen,
}: Readonly<{
  folder: string;
  name: string;
  samples: number;
  hasFolders: boolean;
  open: string;
  onOpen: (f: string) => void;
}>) {
  // The data folder starts with its folders showing, the rest start closed.
  const [expanded, setExpanded] = useState(folder === "");
  const inside = useApi<InputFolder | null>(
    () => (expanded && hasFolders ? api.listInputs(folder) : Promise.resolve(null)),
    [expanded, hasFolders, folder],
  );
  const isOpen = folder === open;
  return (
    <li>
      <div className={isOpen ? "folder-row open" : "folder-row"}>
        {hasFolders ? (
          <button
            type="button"
            className="link folder-toggle"
            aria-expanded={expanded}
            aria-label={`Folders in ${name}`}
            onClick={() => setExpanded(!expanded)}
          >
            {expanded ? "▾" : "▸"}
          </button>
        ) : (
          <span className="folder-toggle" />
        )}
        <button
          type="button"
          className="link folder-name"
          aria-current={isOpen ? "true" : undefined}
          onClick={() => {
            onOpen(folder);
            setExpanded(true);
          }}
        >
          {name}
        </button>
        {samples > 0 && (
          <span className="folder-count muted">
            {samples} sample{samples === 1 ? "" : "s"}
          </span>
        )}
      </div>
      {inside.state === "error" && <p className="error">{inside.error.message}</p>}
      {inside.state === "done" && inside.data && inside.data.folders.length > 0 && (
        <ul>
          {inside.data.folders.map((sub) => (
            <FolderNode
              key={sub.name}
              folder={folder ? `${folder}/${sub.name}` : sub.name}
              name={sub.name}
              samples={sub.samples}
              hasFolders={sub.folders > 0}
              open={open}
              onOpen={onOpen}
            />
          ))}
        </ul>
      )}
    </li>
  );
}

/** A sample that can be run (i.e. at least R1 found. R2 only = nono) */
type Runnable = InputSample & { r1: string };

function JobForm({ found, defaults }: Readonly<{ found: InputFolder; defaults: JobSettings }>) {
  const navigate = useNavigate();
  const ids = useId();
  const folderName = found.folder.split("/").at(-1) ?? "";
  const runnable = found.samples.filter((s): s is Runnable => s.r1 !== null);
  const [name, setName] = useState(folderName || "New job");
  const [method, setMethod] = useState(defaults.method);
  const [tags, setTags] = useState<PickedTag[]>([]);
  // Every sample but the Undetermined reads, to start with.
  const [chosen, setChosen] = useState(
    () => new Set(runnable.filter((s) => !s.undetermined).map((s) => s.r1)),
  );
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  useEffect(() => setError(null), [chosen, name, method]);

  function toggle(r1: string) {
    const next = new Set(chosen);
    if (next.has(r1)) next.delete(r1);
    else next.add(r1);
    setChosen(next);
  }

  async function submit(event: SubmitEvent<HTMLFormElement>) {
    event.preventDefault();
    setBusy(true);
    try {
      const job = await api.createJob({
        name,
        samples: runnable
          .filter((s) => chosen.has(s.r1))
          .map(({ name, r1, r2 }) => ({ name, r1, r2 })),
        settings: { ...defaults, method },
        tags: await tagIds(tags),
      });
      navigate(`/jobs/${job.id}`);
    } catch (e) {
      setError((e as Error).message);
      setBusy(false);
    }
  }

  const problems = found.samples.filter((s) => s.r1 === null);
  const plural = (n: number, word: string) => `${n} ${word}${n === 1 ? "" : "s"}`;
  return (
    <form className="new-job" onSubmit={submit}>
      <h2>{folderName || "ScaleHD Data Root"}</h2>
      <p className="muted">{found.path}</p>

      <details className="file-group" open>
        <summary>
          Samples{" "}
          <span className="muted">
            {chosen.size} of {runnable.length} ticked
          </span>
        </summary>
        {runnable.length === 0 ? (
          <p className="muted">None of the FASTQ files here can be run.</p>
        ) : (
          <>
            <p className="sample-tools">
              <button
                type="button"
                className="link"
                onClick={() => setChosen(new Set(runnable.map((s) => s.r1)))}
              >
                Select all
              </button>
              <button type="button" className="link" onClick={() => setChosen(new Set())}>
                Select none
              </button>
            </p>
            <table className="sample-choice">
              <thead>
                <tr>
                  <th aria-label="Chosen" />
                  <th>Sample and files</th>
                  <th>Reads</th>
                  <th>Size</th>
                </tr>
              </thead>
              <tbody>
                {runnable.map((sample, i) => {
                  const id = `${ids}-${i}`;
                  return (
                    <tr key={sample.r1}>
                      <td>
                        <input
                          id={id}
                          type="checkbox"
                          aria-label={sample.name}
                          checked={chosen.has(sample.r1)}
                          onChange={() => toggle(sample.r1)}
                        />
                      </td>
                      <td>
                        <SampleFiles sample={sample} labelFor={id} />
                      </td>
                      <td>{sample.r2 ? "R1 + R2" : "R1 only"}</td>
                      <td>{fileSize(sample.size)}</td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </>
        )}
      </details>

      {problems.length > 0 && (
        <details className="file-group problems" open>
          <summary>
            Can't be run <span className="muted">{plural(problems.length, "sample")}</span>
          </summary>
          <table>
            <thead>
              <tr>
                <th>Sample and files</th>
                <th>Size</th>
              </tr>
            </thead>
            <tbody>
              {problems.map((sample) => (
                <tr key={sample.files[0]}>
                  <td>
                    <SampleFiles sample={sample} />
                  </td>
                  <td>{fileSize(sample.size)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </details>
      )}

      {found.other_files.length > 0 && (
        <details className="file-group">
          <summary>
            Other files{" "}
            <span className="muted">{plural(found.other_files.length, "non-FASTQ file")}</span>
          </summary>
          <ul className="other-files">
            {found.other_files.map((file) => (
              <li key={file}>{file}</li>
            ))}
          </ul>
        </details>
      )}

      <div className="form wide">
        <label>
          Job name
          <input value={name} onChange={(e) => setName(e.target.value)} required maxLength={200} />
        </label>
        <TagPicker value={tags} onChange={setTags} />
        <MethodPicker value={method} onChange={setMethod} defaultMethod={defaults.method} />
        <button type="submit" disabled={busy || chosen.size === 0}>
          Run {chosen.size} sample{chosen.size === 1 ? "" : "s"}
        </button>
        {error && <p className="error">{error}</p>}
      </div>
    </form>
  );
}

/** A sample's name, the files it was paired from, and why it can't be run if it can't. */
function SampleFiles({ sample, labelFor }: Readonly<{ sample: InputSample; labelFor?: string }>) {
  return (
    <>
      {labelFor ? (
        <label htmlFor={labelFor} className="sample-name">
          {sample.name}
        </label>
      ) : (
        <span className="sample-name">{sample.name}</span>
      )}
      {sample.undetermined && <span className="muted"> reads that matched no sample</span>}
      <ul className="sample-files">
        {sample.files.map((file) => (
          <li key={file}>{file}</li>
        ))}
      </ul>
      {sample.skipped && <p className="sample-skipped">{sample.skipped}</p>}
    </>
  );
}

function fileSize(bytes: number): string {
  return bytes < 1e5 ? `${Math.max(1, Math.round(bytes / 1e3))} KB` : `${(bytes / 1e6).toFixed(1)} MB`;
}
