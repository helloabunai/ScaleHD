// Client for the ScaleHD HTTP API. The types gen by
// tools/gen-api-types.sh from the server's OpenAPI desc.
// Rerun after any spec changes.

import type { components } from "./api.gen";

type Schemas = components["schemas"];

export type JobStatus = Schemas["JobStatus"];
/** Light or dark pages, or "system" */
export type Theme = Schemas["Theme"];
export type SampleStatus = Schemas["SampleStatus"];
export type Health = Schemas["Health"];
export type User = Schemas["UserOut"];
export type Folders = Schemas["Folders"];
export type Registration = Schemas["Registration"];
export type GenotypeMethod = Schemas["GenotypeMethod"];
export type JobSettings = Schemas["JobSettings"];
export type InputPair = Schemas["InputPair"];
export type InputSample = Schemas["InputSample"];
export type InputSubfolder = Schemas["InputSubfolder"];
export type InputFolder = Schemas["InputFolder"];
export type JobCreate = Schemas["JobCreate"];

// Also in schemas.py, as limits on the request bodies. OpenAPI has no constants.
export const MAX_TAGS = 5;
export const MAX_TAG_LENGTH = 15;

export type JobTag = Schemas["JobTag"];
export type Tag = Schemas["TagOut"];
export type Sample = Schemas["SampleOut"];
export type SampleStage = Schemas["SampleStage"];
export type SampleProgress = Schemas["SampleProgress"];
export type JobSummary = Schemas["JobSummary"];
export type Job = Schemas["JobOut"];
export type Flag = Schemas["Flag"];
export type Stutter = Schemas["StutterJson"];
/** One allele of a saved call (scalehd.call/2). */
export type CalledAllele = Schemas["AlleleCallJson"];
/** A sample's saved genotype call (scalehd.call/2). */
export type GenotypeCall = Schemas["GenotypeCallJson"];
export type CagBar = Schemas["CagBar"];
export type CagChart = Schemas["CagChart"];
export type CcgBar = Schemas["CcgBar"];
export type Cell = Schemas["Cell"];
export type Reads = Schemas["Reads"];
export type SampleDetail = Schemas["SampleDetail"];
export type SampleFile = SampleDetail["files"][number];

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
  jobResultsUrl: (jobId: number) => `/api/jobs/${jobId}/results.csv`,
  jobValidationUrl: (jobId: number) => `/api/jobs/${jobId}/validation.csv`,
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
  setPassword: (userId: number, password: string) =>
    request<void>(`/admin/users/${userId}/password`, {
      method: "PUT",
      body: JSON.stringify({ password }),
    }),

  getSettings: () => request<JobSettings>("/settings"),
  saveSettings: (settings: JobSettings) =>
    request<JobSettings>("/settings", { method: "PUT", body: JSON.stringify(settings) }),
};
