# Web interface and API

The web interface is a completely separate application running within the same docker container/machine.
This allows people to choose how they want to interact with the software. Jobs can be launched via the web
interface, or via the command line. If the front end breaks (for whatever reason), then jobs running won't
be taken down with it.

It's a fairly simple interface at the moment, where users can run a demo job (eventually to be removed) to
demonstrate the Frontend <-> API <-> Backend communication is functional. Real jobs can be submitted from the
jobs page, where users can also browse previous jobs that were analysed on the same ScaleHD server.

An admin user can see jobs from every user account. Individual regular users can only see their own jobs.

What has been implemented so far:

- Accounts: register and log in (the first account is made admin), change password,
  and a light/dark setting per account. Any admin can promote regular users to have admin rights,
  and set a new password for a user who has forgotten theirs (which logs them out everywhere).
- "Run demo" on the home page genotypes thirteen simulated samples with known
  genotypes, end to end, with the model-based caller.
- Jobs: launch genotyping jobs. Browse running jobs (and progress), and previous jobs
  with their genotyping results / analysis flags. Cancel a job which you started. Samples not
  yet started don't continue after cancellation, those running finish. Jobs cancelled can also be deleted.
  Users only see their own jobs.
- A results page per sample with each allele's repeat structure rendered,
  the CAG and CCG distributions, instability figures, alternative genotypes,
  unexplained peaks and how its reads fared, with its counts and call JSON and its FASTQ
  files to download. No results exporting yet (e.g. PDF file or something).
- Settings: the default genotyping method and thresholds new jobs start with. Only basic settings
  present for now.

Not yet: exporting reports. API docs are at `/api/docs`. See
[`apps/server/README.md`](apps/server/README.md) for what is real and what's a stub.

Things will move around a lot as development proceeds.
