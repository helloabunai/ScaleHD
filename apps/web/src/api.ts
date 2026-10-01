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

export interface JobSettings {
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
  settings?: Partial<JobSettings>;
}

export interface Sample extends InputPair {
  id: number;
  status: SampleStatus;
  genotype: string | null;
  quality: number | null;
  flags: string[];
  error: string | null;
}

export interface JobSummary {
  id: number;
  name: string;
  status: JobStatus;
  created_at: string;
  started_at: string | null;
  finished_at: string | null;
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

async function request<T>(path: string, init: RequestInit = {}): Promise<T> {
  const response = await fetch(`/api${path}`, {
    credentials: "same-origin",
    headers: init.body ? { "Content-Type": "application/json" } : undefined,
    ...init,
  });
  if (!response.ok) {
    const body = await response.json().catch(() => null);
    throw new ApiError(response.status, body?.detail ?? response.statusText);
  }
  return (response.status === 204 ? undefined : await response.json()) as T;
}

const post = (body?: unknown): RequestInit => ({
  method: "POST",
  body: body === undefined ? undefined : JSON.stringify(body),
});

export const api = {
  health: () => request<Health>("/health"),

  me: () => request<User>("/auth/me"),
  login: (username: string, password: string) =>
    request<User>("/auth/login", post({ username, password })),
  logout: () => request<void>("/auth/logout", post()),

  listInputs: (folder = "") => request<InputPair[]>(`/inputs?folder=${encodeURIComponent(folder)}`),

  listJobs: () => request<JobSummary[]>("/jobs"),
  getJob: (id: number) => request<Job>(`/jobs/${id}`),
  createJob: (job: JobCreate) => request<Job>("/jobs", post(job)),
  cancelJob: (id: number) => request<Job>(`/jobs/${id}/cancel`, post()),
  reportUrl: (id: number) => `/api/jobs/${id}/report`,

  getSettings: () => request<JobSettings>("/settings"),
  saveSettings: (settings: JobSettings) =>
    request<JobSettings>("/settings", { method: "PUT", body: JSON.stringify(settings) }),
};
