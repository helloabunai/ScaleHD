-- A ScaleHD database as main (848b161) made it before migrations. Tables from
-- create_all, unnamed constraints, etc
BEGIN TRANSACTION;
CREATE TABLE job_tags (
	job_id INTEGER NOT NULL, 
	tag_id INTEGER NOT NULL, 
	PRIMARY KEY (job_id, tag_id), 
	FOREIGN KEY(job_id) REFERENCES jobs (id) ON DELETE CASCADE, 
	FOREIGN KEY(tag_id) REFERENCES tags (id) ON DELETE CASCADE
);
INSERT INTO "job_tags" VALUES(1,1);
INSERT INTO "job_tags" VALUES(1,2);
CREATE TABLE jobs (
	id INTEGER NOT NULL, 
	owner_id INTEGER NOT NULL, 
	name VARCHAR(200) NOT NULL, 
	status VARCHAR(10) NOT NULL, 
	settings JSON NOT NULL, 
	created_at DATETIME NOT NULL, 
	started_at DATETIME, 
	finished_at DATETIME, 
	error VARCHAR, 
	demo BOOLEAN NOT NULL, 
	output_dir VARCHAR, 
	PRIMARY KEY (id), 
	FOREIGN KEY(owner_id) REFERENCES users (id)
);
INSERT INTO "jobs" VALUES(1,1,'run-01','FINISHED','{"method": "model"}','2026-10-06 09:00:00.000000','2026-10-06 09:00:00.000000','2026-10-06 09:02:00.000000',NULL,0,'/srv/scalehd/workspace/first-admin/1-run-01');
INSERT INTO "jobs" VALUES(2,2,'cancelled run','CANCELLED','{"method": "model"}','2026-10-06 09:00:00.000000',NULL,NULL,NULL,1,NULL);
CREATE TABLE login_sessions (
	id INTEGER NOT NULL, 
	user_id INTEGER NOT NULL, 
	token_hash VARCHAR(64) NOT NULL, 
	created_at DATETIME NOT NULL, 
	expires_at DATETIME NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(user_id) REFERENCES users (id) ON DELETE CASCADE, 
	UNIQUE (token_hash)
);
INSERT INTO "login_sessions" VALUES(1,1,'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa','2026-10-06 09:00:00.000000','2026-11-05 09:00:00.000000');
CREATE TABLE samples (
	id INTEGER NOT NULL, 
	job_id INTEGER NOT NULL, 
	name VARCHAR(200) NOT NULL, 
	r1 VARCHAR, 
	r2 VARCHAR, 
	simulation JSON, 
	truth VARCHAR, 
	matches_truth BOOLEAN, 
	status VARCHAR(9) NOT NULL, 
	genotype VARCHAR, 
	confidence DOUBLE, 
	flags JSON NOT NULL, 
	call JSON, 
	error VARCHAR, 
	PRIMARY KEY (id), 
	FOREIGN KEY(job_id) REFERENCES jobs (id)
);
INSERT INTO "samples" VALUES(1,1,'normal-17-21','/srv/data/run-01/a_R1.fastq.gz','/srv/data/run-01/a_R2.fastq.gz',NULL,'17_1_1_7_2/21_1_1_7_2',1,'FINISHED','17_1_1_7_2/21_1_1_7_2',99.0,'[]','{"schema": "scalehd.call/2", "genotype": "17_1_1_7_2/21_1_1_7_2"}',NULL);
INSERT INTO "samples" VALUES(2,1,'orphan','/srv/data/run-01/b_R1.fastq.gz',NULL,NULL,NULL,NULL,'FAILED',NULL,NULL,'[]',NULL,'ValueError: no molecules');
INSERT INTO "samples" VALUES(3,2,'never-ran',NULL,NULL,'{"alleles": ["21_1_1_7_2"], "pairs": 100, "seed": 1}',NULL,NULL,'CANCELLED',NULL,NULL,'[]',NULL,NULL);
CREATE TABLE tags (
	id INTEGER NOT NULL, 
	name VARCHAR(15) NOT NULL, 
	created_by_id INTEGER, 
	created_at DATETIME NOT NULL, 
	PRIMARY KEY (id), 
	UNIQUE (name), 
	FOREIGN KEY(created_by_id) REFERENCES users (id)
);
INSERT INTO "tags" VALUES(1,'Paper XYZ',1,'2026-10-06 09:00:00.000000');
INSERT INTO "tags" VALUES(2,'cohort 2',2,'2026-10-06 09:00:00.000000');
CREATE TABLE users (
	id INTEGER NOT NULL, 
	username VARCHAR(64) NOT NULL, 
	password_hash VARCHAR NOT NULL, 
	is_admin BOOLEAN NOT NULL, 
	created_at DATETIME NOT NULL, 
	default_settings JSON NOT NULL, 
	theme VARCHAR(6) NOT NULL, 
	PRIMARY KEY (id), 
	UNIQUE (username)
);
INSERT INTO "users" VALUES(1,'first-admin','$argon2id$placeholder',1,'2026-10-06 09:00:00.000000','{"method": "legacy"}','SYSTEM');
INSERT INTO "users" VALUES(2,'second-user','$argon2id$placeholder',0,'2026-10-06 09:00:00.000000','{}','SYSTEM');
COMMIT;
