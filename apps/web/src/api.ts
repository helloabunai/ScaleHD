// Client for the ScaleHD HTTP API. The types mirror
// apps/server/src/scalehd_server/schemas.py by hand. TODO: generate them from
// /api/openapi.json (e.g. openapi-typescript) once the API settles.

export type JobStatus =
  | "queued"
  | "running"
  | "cancelling"
  | "finished"
  | "failed"
  | "cancelled";
/** Light or dark pages, or "system" */
export type Theme = "system" | "light" | "dark";
export type SampleStatus = "queued" | "running" | "finished" | "failed" | "cancelled";

export interface Health {
  status: "ok";
  version: string;
  core_version: string;
  python: string;
  platform: string;
  sqlite: string;
  libraries: Record<string, string>;
}

export interface User {
  id: number;
  username: string;
  is_admin: boolean;
  created_at: string;
  theme: Theme;
}

export interface Folders {
  /** Sequencing data to pick from, read-only. Null when the server has none set. */
  data_root: string | null;
  /** Where every user's results go. */
  workspace: string;
  /** This user's folder in the workspace: their jobs are saved here. */
  your_folder: string;
}

export interface Registration {
  open: boolean;
  first_account: boolean;
}

export type GenotypeMethod = "legacy" | "model";

export interface JobSettings {
  method: GenotypeMethod;
  call: boolean;
  discordant: "drop" | "prefer";
  min_molecules: number | null;
  min_posterior: number | null;
  max_background: number | null;
  max_dropped: number | null;
}

/** A sample's FASTQ files, as paths relative to the server's data folder. */
export interface InputPair {
  name: string;
  r1: string;
  r2: string | null;
}

/** A sample found in the data folder & its FASTQ files, paired by name. */
export interface InputSample {
  name: string;
  files: string[];
  r1: string | null;/** Null if not runnable for whatever reason */
  r2: string | null; /** Null if not runnable for whatever reason */
  size: number;
  undetermined: boolean;
  skipped: string | null;
}

/** A folder inside the open one, with what it contains (for tree) */
export interface InputSubfolder {
  name: string;
  folders: number;
  samples: number;
}

/** Subdir of the server data root, and the samples its FASTQ files make. */
export interface InputFolder {
  folder: string;
  path: string;
  folders: InputSubfolder[];
  samples: InputSample[];
  other_files: string[];
}

export interface JobCreate {
  name: string;
  samples: InputPair[];
  /** Left out, the job uses the user's default settings. */
  settings?: JobSettings;
  tags?: number[];
}

export const MAX_TAGS = 5;
export const MAX_TAG_LENGTH = 15;

export interface JobTag {
  id: number;
  name: string;
}

export interface Tag extends JobTag {
  jobs: number;
}

export interface Sample {
  id: number;
  name: string;
  r1: string | null;
  r2: string | null;
  status: SampleStatus;
  genotype: string | null;
  /** Phred-scaled confidence of the call (the call's quality). */
  confidence: number | null;
  flags: string[];
  /** Simulated samples: the true genotype, and whether the call matched it. */
  truth: string | null;
  matches_truth: boolean | null;
  error: string | null;
}

export interface JobSummary {
  id: number;
  name: string;
  demo: boolean;
  method: GenotypeMethod;
  status: JobStatus;
  created_at: string;
  started_at: string | null;
  finished_at: string | null;
  /** The job's folder in the workspace, a path on the server machine. */
  output_dir: string | null;
  sample_count: number;
  samples_done: number;
  tags: JobTag[];
}

export interface Job extends JobSummary {
  settings: JobSettings;
  samples: Sample[];
  error: string | null;
}

export interface Stutter {
  contraction: number;
  contraction_step: number;
  contraction_tail: number;
  expansion: number;
  expansion_step: number;
  expansion_tail: number;
}

/** One allele of a saved call (scalehd.call/2). */
export interface CalledAllele {
  structure: string;
  beyond_read_length: boolean;
  cag: number;
  caacag: number;
  ccgcca: number;
  ccg: number;
  cct: number;
  typical: boolean;
  polyglutamine_length: number | null;
  fraction: number;
  molecules: number;
  stutter: Stutter;
  backward_slippage: number | null;
  somatic_mosaicism: number | null;
  expansion_index: number | null;
  contraction_index: number | null;
  /** For an allele beyond read length: [best, low, high] estimate of its CAG. */
  cag_estimate: [number, number, number] | null;
}

/** A sample's saved genotype call (scalehd.call/2). */
export interface GenotypeCall {
  genotype: string;
  posterior: number;
  quality: number;
  flags: string[];
  alleles: CalledAllele[];
  alternatives: { genotype: string; posterior: number }[];
  molecules: number;
  background: number;
  ccg_slippage: number;
  misread: number;
  unexplained: { structure: string; molecules: number }[];
}

export interface CagBar {
  cag: number;
  molecules: number;
  /** Molecules only known to be at least this long. */
  lower_bound: number;
}

export interface CagChart {
  caacag: number;
  ccgcca: number;
  ccg: number;
  cct: number;
  alleles: string[];
  /** The called alleles' CAG lengths, ascending (lower bounds for alleles beyond read length). */
  called: number[];
  bars: CagBar[];
}

export interface CcgBar {
  ccg: number;
  molecules: number;
}

export interface Cell {
  cag: number;
  ccg: number;
  molecules: number;
}

export interface Reads {
  molecules: number;
  complete: number;
  partial: number;
  dropped: number;
  unusable: number;
  read_outcomes: Record<string, number>;
  discordant: Record<string, number>;
}

export type SampleFile = "call" | "counts" | "r1" | "r2";

export interface SampleDetail {
  sample: Sample;
  job_id: number;
  job_name: string;
  demo: boolean;
  tags: JobTag[];
  folder: string | null;
  call: GenotypeCall | null;
  cag_charts: CagChart[];
  ccg: CcgBar[];
  cells: Cell[];
  reads: Reads | null;
  files: SampleFile[];
  previous_id: number | null;
  next_id: number | null;
}

export class ApiError extends Error {
  readonly status: number;

  constructor(status: number, message: string) {
    super(message);
    this.status = status;
  }
}

let unauthorized = () => {};

/** Called whenever the server says 401, e.g. when a login has expired. */
export function onUnauthorized(handler: () => void) {
  unauthorized = handler;
}

async function request<T>(path: string, init: RequestInit = {}): Promise<T> {
  const response = await fetch(`/api${path}`, {
    credentials: "same-origin",
    headers: init.body ? { "Content-Type": "application/json" } : undefined,
    ...init,
  });
  if (!response.ok) {
    if (response.status === 401) unauthorized();
    const body = await response.json().catch(() => null);
    throw new ApiError(response.status, describe(body?.detail) ?? response.statusText);
  }
  return (response.status === 204 ? undefined : await response.json()) as T;
}

// FastAPI's detail is a string, or for invalid input a list of {loc, msg} objects.
function describe(detail: unknown): string | undefined {
  if (typeof detail === "string") return detail;
  if (Array.isArray(detail)) {
    return detail
      .map((d: { loc?: unknown[]; msg?: string }) => `${d.loc?.at(-1)}: ${d.msg}`)
      .join("; ");
  }
  return undefined;
}

const post = (body?: unknown): RequestInit => ({
  method: "POST",
  body: body === undefined ? undefined : JSON.stringify(body),
});

export const api = {
  health: () => request<Health>("/health"),
  folders: () => request<Folders>("/folders"),

  registration: () => request<Registration>("/auth/registration"),
  register: (username: string, password: string) =>
    request<User>("/auth/register", post({ username, password })),
  me: () => request<User>("/auth/me"),
  login: (username: string, password: string) =>
    request<User>("/auth/login", post({ username, password })),
  logout: () => request<void>("/auth/logout", post()),
  changePassword: (current_password: string, new_password: string) =>
    request<void>("/auth/password", {
      method: "PUT",
      body: JSON.stringify({ current_password, new_password }),
    }),
  setTheme: (theme: Theme) =>
    request<User>("/auth/theme", { method: "PUT", body: JSON.stringify({ theme }) }),

  listInputs: (folder = "") =>
    request<InputFolder>(`/inputs?folder=${encodeURIComponent(folder)}`),

  listJobs: () => request<JobSummary[]>("/jobs"),
  getJob: (id: number) => request<Job>(`/jobs/${id}`),
  deleteJob: (id: number) => request<void>(`/jobs/${id}`, { method: "DELETE" }),
  getSample: (jobId: number, sampleId: number) =>
    request<SampleDetail>(`/jobs/${jobId}/samples/${sampleId}`),
  sampleFileUrl: (jobId: number, sampleId: number, file: SampleFile) =>
    `/api/jobs/${jobId}/samples/${sampleId}/files/${file}`,
  createJob: (job: JobCreate) => request<Job>("/jobs", post(job)),
  createDemoJob: () => request<Job>("/jobs/demo", post()),
  cancelJob: (id: number) => request<Job>(`/jobs/${id}/cancel`, post()),
  reportUrl: (id: number) => `/api/jobs/${id}/report`,

  listTags: () => request<Tag[]>("/tags"),
  createTag: (name: string) => request<Tag>("/tags", post({ name })),
  renameTag: (id: number, name: string) =>
    request<Tag>(`/tags/${id}`, { method: "PATCH", body: JSON.stringify({ name }) }),
  deleteTag: (id: number) => request<void>(`/tags/${id}`, { method: "DELETE" }),
  setJobTags: (jobId: number, tags: number[]) =>
    request<Job>(`/jobs/${jobId}/tags`, { method: "PUT", body: JSON.stringify({ tags }) }),

  listUsers: () => request<User[]>("/admin/users"),
  setAdmin: (userId: number, isAdmin: boolean) =>
    request<User>(`/admin/users/${userId}/admin`, {
      method: "PUT",
      body: JSON.stringify({ is_admin: isAdmin }),
    }),

  getSettings: () => request<JobSettings>("/settings"),
  saveSettings: (settings: JobSettings) =>
    request<JobSettings>("/settings", { method: "PUT", body: JSON.stringify(settings) }),
};
