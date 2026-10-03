import { api } from "../api";
import { useApi } from "../useApi";
import { Status } from "./Status";

/** Where the server reads sequencing data from and saves results to. */
export function ServerFolders() {
  const folders = useApi(api.folders);
  return (
    <section className="folders">
      <h2>Server folders</h2>
      <Status of={folders}>
        {(folders) => (
          <>
            <h3>Data</h3>
            <p>
              {folders.data_root ? (
                <code>{folders.data_root}</code>
              ) : (
                <>
                  No data folder is set (<code>SCALEHD_DATA_ROOT</code>).
                </>
              )}
            </p>
            <p className="muted">Read-only directory containing your sequencing data. 
              <br/>The server admin must define this directory before starting the ScaleHD server.</p>
            <h3>Workspace</h3>
            <p>
              <code>{folders.workspace}</code>
            </p>
            <p className="muted">
              Your jobs are saved to <code>{folders.your_folder}/&lt;job id&gt;-&lt;job name&gt;/</code>,
              with a folder per sample inside.
            </p>
          </>
        )}
      </Status>
    </section>
  );
}
