# From scratch, running the server with Docker

Everything (server, web interface, genotyping algorithm(s)) runs in one container, so the host machine
only needs 3 things. Docker, Docker Compose and git.

1. Install Docker and Compose.

Example from scratch with ubuntu, except installing docker direct
from docker's repository as the OS repo is typically out of date.

```sh
sudo apt-get update
sudo apt-get install -y ca-certificates curl git
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc
echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] \
https://download.docker.com/linux/ubuntu $(. /etc/os-release && echo "$VERSION_CODENAME") stable" \
    | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

Then let your user run Docker without `sudo`, and log out and back in for it to take
effect (or run `newgrp docker` in the current shell):

```sh
sudo usermod -aG docker "$USER"
```

Membership of the `docker` group is root-equivalent on that machine, so only add
users who should have it.

Check both work:

```sh
docker run --rm hello-world
docker compose version
```

2. Get this code

```sh
git clone https://github.com/helloabunai/ScaleHD.git
cd ScaleHD
```

3. Pick the folders and write `.env`.

The server needs two folders on the host. The container sees each at the same path as
the host does, so paths shown in the web interface are real host paths.

- A data folder with the sequencing runs (FASTQ files), mounted read-only. A network
  share works if it is mounted on the host before the container starts.
- A workspace folder for results, mounted read-write. ScaleHD will create a folder per user,
  and a subfolder per job. The root workspace/results folder must exist before starting. On Linux it must be writable by the container's user (the following paths are made up examples):

  ```sh
  sudo mkdir -p /path/to/desired/workspace
  sudo chown 1000 /path/to/desired/workspace
  ```

Specify such folder locations in the env file (copy the example template).

```sh
cp .env.example .env
```

Then edit `.env`. Only the first two settings are required at the moment:

| setting | default | meaning |
|---|---|---|
| `SCALEHD_DATA_ROOT` | required | sequencing data users can browse, read-only (absolute path) |
| `SCALEHD_WORKSPACE` | required | where results go (absolute path, must exist) |
| `SCALEHD_WORKERS` | every core | samples processed at once, one process each |
| `SCALEHD_ALLOW_REGISTRATION` | `true` | whether anyone who can reach the server may create an account (the first account, the admin, can always be created) |
| `SCALEHD_SESSION_DAYS` | `30` | how long a login lasts |

The database lives in a Docker volume (`scalehd_scalehd-data`), separate from the
workspace. `SCALEHD_DATABASE_DIR` and `SCALEHD_WEB_DIR` are set by the image; leave
them alone.

Most settings are work in progress and/or missing.

- Of the genotyping methods (Settings page, and per job), Legacy (ScaleHD 1.x) is not
  available yet and New (model-based) is a beta. The demo always uses New.
- The flag thresholds and the other job settings are not implemented or don't exist yet so ignore those.

4. Start it

```sh
docker compose up -d --build
docker compose ps        # STATUS should become "healthy"
docker compose logs -f   # server log (Ctrl+C stops following, not the server)
```

The container restarts on its own after a reboot (`restart: unless-stopped`).
`docker compose down` stops it. The database volume and the folders are kept.

5. Open the web interface.

The docker image (i.e. "production") will launch by default at <http://localhost:8000>.
The development environment will launch by default at <http://localhost:5173>.

The first account created is automatically an administrator. Other user accounts can be
promoted to have admin rights by any other admin on the server.

From there, "Run demo" on the home page runs the simulated samples
end to end. You can also provide real input FastQ/gz files and launch genotyping from the
jobs page.

By default the server only listens on the host machine itself. To reach it from other
computers on the same local network (<http://whatever_server_dedicated_ip_address:8000>), drop the
`127.0.0.1:` from the `ports` line in `compose.yaml` (the first `8000` there is the port
to use, if you want a different one).

Because this is just a docker image based setup with the intention of it running on an offline server
that is within the same local network as your client(s), the server doesn't serve with HTTPS.

If you want secure remote access then you'd need to do a reverse proxy in front of the server on your own.

## Updating

```sh
git pull
docker compose up -d --build
```

Databases are managed with alembic migrations. If one doesn't exist (first boot) then it's created.
If it exists, updates are attempted if required. Failures can be reverted with backup copies of db
(auto created at time of migration attempt). To go back to an older ScaleHD DB, stop the server, 
put the named copy back as main `scalehd.db`, and run the older version of the server.
Still subject to change. Development work in progress baby !!!
