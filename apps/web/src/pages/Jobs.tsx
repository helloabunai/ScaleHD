import { Link } from "react-router";
import { api } from "../api";
import { useApi } from "../useApi";
import { DemoTag, formatTime, isActive } from "./jobDisplay";
import { METHODS } from "./MethodPicker";
import { ServerFolders } from "./ServerFolders";
import { Status } from "./Status";

export function Jobs() {
  const jobs = useApi(api.listJobs, [], {
    every: 2000,
    while: (list) => list.some((job) => isActive(job.status)),
  });
  return (
    <>
      <h1>Jobs</h1>
      <ServerFolders />
      <Status of={jobs}>
        {(list) =>
          list.length === 0 ? (
            <p>
              No jobs yet. <Link to="/">Run the demo</Link> or{" "}
              <Link to="/jobs/new">start one</Link>.
            </p>
          ) : (
            <table>
              <thead>
                <tr>
                  <th>Name</th>
                  <th>Method</th>
                  <th>Status</th>
                  <th>Samples</th>
                  <th>Created</th>
                </tr>
              </thead>
              <tbody>
                {list.map((job) => (
                  <tr key={job.id}>
                    <td>
                      <Link to={`/jobs/${job.id}`}>{job.name}</Link> {job.demo && <DemoTag />}
                    </td>
                    <td>{METHODS[job.method].label}</td>
                    <td>{job.status}</td>
                    <td>
                      {job.samples_done}/{job.sample_count}
                    </td>
                    <td>{formatTime(job.created_at)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          )
        }
      </Status>
    </>
  );
}
