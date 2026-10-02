import type { JobStatus } from "../api";

/** Whether a job is not entirely done, so its page should keep polling for updates. */
export function isActive(status: JobStatus): boolean {
  return status === "queued" || status === "running";
}

export function DemoTag() {
  return <span className="tag">demo</span>;
}

export function formatTime(iso: string): string {
  return new Date(iso).toLocaleString();
}
