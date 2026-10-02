// Client for the ScaleHD HTTP API. The types mirror
// apps/server/src/scalehd_server/schemas.py by hand. TODO: generate them from
// /api/openapi.json (e.g. openapi-typescript) once the API settles.

export type JobStatus = "queued" | "running" | "finished" | "failed" | "cancelled";
export type SampleStatus = "queued" | "running" | "finished" | "failed";

export interface Health {
  status: "ok";
  version: string;
  core_version: string;
}

export interface User {
  id: number;
  username: string;
  is_admin: boolean;
  created_at: string;
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

export interface InputPair {
  name: string;
  r1: string;
  r2: string | null;
}

export interface JobCreate {
  name: string;
  samples: InputPair[];
  /** Left out, the job uses the user's default settings. */
  settings?: JobSettings;
}

export interface Sample {
  id: number;
  name: string;
  r1: string | null;
  r2: string | null;
  status: SampleStatus;
  genotype: string | null;
  quality: number | null;
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
}

export interface Job extends JobSummary {
  settings: JobSettings;
  samples: Sample[];
  error: string | null;
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

  listInputs: (folder = "") => request<InputPair[]>(`/inputs?folder=${encodeURIComponent(folder)}`),

  listJobs: () => request<JobSummary[]>("/jobs"),
  getJob: (id: number) => request<Job>(`/jobs/${id}`),
  createJob: (job: JobCreate) => request<Job>("/jobs", post(job)),
  createDemoJob: () => request<Job>("/jobs/demo", post()),
  cancelJob: (id: number) => request<Job>(`/jobs/${id}/cancel`, post()),
  reportUrl: (id: number) => `/api/jobs/${id}/report`,

  getSettings: () => request<JobSettings>("/settings"),
  saveSettings: (settings: JobSettings) =>
    request<JobSettings>("/settings", { method: "PUT", body: JSON.stringify(settings) }),
};
