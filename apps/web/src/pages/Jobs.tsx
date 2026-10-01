import { Link } from "react-router";
import { api } from "../api";
import { useApi } from "../useApi";
import { Status } from "./Status";

export function Jobs() {
  const jobs = useApi(api.listJobs);
  return (
    <>
      <h1>Jobs</h1>
      <Status of={jobs}>
        {(list) =>
          list.length === 0 ? (
            <p>
              No jobs yet. <Link to="/jobs/new">Start one.</Link>
            </p>
          ) : (
            <table>
              <thead>
                <tr>
                  <th>Name</th>
                  <th>Status</th>
                  <th>Samples</th>
                  <th>Created</th>
                </tr>
              </thead>
              <tbody>
                {list.map((job) => (
                  <tr key={job.id}>
                    <td>
                      <Link to={`/jobs/${job.id}`}>{job.name}</Link>
                    </td>
                    <td>{job.status}</td>
                    <td>
                      {job.samples_done}/{job.sample_count}
                    </td>
                    <td>{new Date(job.created_at).toLocaleString()}</td>
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
