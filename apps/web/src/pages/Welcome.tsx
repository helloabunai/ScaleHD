import { Link } from "react-router";
import { useUser } from "../auth";

export function Welcome() {
  const user = useUser();
  return (
    <>
      <h1>Welcome, {user.username}</h1>
      <p>
        ScaleHD genotypes the HTT CAG/CCG repeat from paired-end amplicon sequencing: the
        repeat structure of each allele, how confident the call is, and anything that deserves
        a second look.

        This website is served from a docker container which talks to the command line package
        via an API. In theory. Not yet .. :) 
      </p>
      <div className="cards">
        <Link className="card" to="/jobs/new">
          <h2>Start a job</h2>
          <p>Pick FASTQ files from the server's input folder and genotype them.</p>
        </Link>
        <Link className="card" to="/jobs">
          <h2>Jobs</h2>
          <p>Follow running jobs and look back at past results.</p>
        </Link>
        <Link className="card" to="/settings">
          <h2>Default settings</h2>
          <p>The thresholds new jobs start with.</p>
        </Link>
      </div>
      {user.is_admin && <p className="muted">You're this server's admin. Go nuts.</p>}
    </>
  );
}
