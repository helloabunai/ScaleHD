import { Fragment, useMemo } from "react";
import type {
  LegacyAlignedReads,
  LegacyAllele,
  LegacyCagChart,
  LegacyCall,
  LegacyCcgChart,
  LegacyFastqc,
  LegacyQualityControl,
  LegacySample,
  LegacyTrimming,
  SampleDetail,
} from "../api";
import { type Bar, BarChart, type Scale } from "../charts/BarChart";
import { BaseKey, ReadBlocks } from "../charts/ReadBlocks";
import { LEGACY_FAILED_STAGES, LEGACY_FLAGS, legacyConfidenceInWords } from "../legacyFlags";
import { ScaleToggle } from "./SampleResults";

export function LegacyResults({
  detail,
  legacy,
  scale,
  setScale,
}: Readonly<{
  detail: SampleDetail;
  legacy: LegacySample;
  scale: Scale;
  setScale: (s: Scale) => void;
}>) {
  return (
    <>
      {legacy.status === "failed" && <p className="error">{failure(legacy)}</p>}
      {legacy.call && (
        <>
          <Summary detail={detail} call={legacy.call} />
          <Alleles alleles={distinct(legacy.call.alleles)} />
          <section>
            <div className="section-head">
              <h2>CAG distribution</h2>
              <ScaleToggle scale={scale} onChange={setScale} />
            </div>
            <p className="muted">
              Forward reads aligned to each CAG length at the allele's CCG.
              Called lengths are highlighted. Hover a column for details.
            </p>
            {legacy.cag_charts.map((chart) => (
              <CagChartView key={chart.alleles.join()} chart={chart} scale={scale} />
            ))}
          </section>
          <section>
            <h2>CCG distribution</h2>
            <p className="muted">Reverse reads aligned to each CCG length.</p>
            {legacy.ccg_charts.map((chart) => (
              <CcgChartView
                key={chart.alleles.join()}
                chart={chart}
                call={legacy.call as LegacyCall}
                scale={scale}
                caption={legacy.ccg_charts.length > 1}
              />
            ))}
          </section>
          <Mapping alleles={distinct(legacy.call.alleles)} />
          <Flags call={legacy.call} />
        </>
      )}
      {legacy.aligned_reads.length > 0 && <AlignedReads shown={legacy.aligned_reads} />}
      {legacy.quality_control && (
        <QualityControl checks={legacy.quality_control} scale={scale} />
      )}
      <Run legacy={legacy} />
    </>
  );
}

function failure(legacy: LegacySample): string {
  const stage = legacy.failed_stage;
  const where = stage ? ` in ${LEGACY_FAILED_STAGES[stage] ?? stage} (${stage})` : "";
  return `Legacy ScaleHD failed${where}: ${legacy.error ?? "it gave no reason"}`;
}

function Summary({ detail, call }: Readonly<{ detail: SampleDetail; call: LegacyCall }>) {
  const { truth, matches_truth } = detail.sample;
  return (
    <section>
      <p className="genotype">{call.genotype.replace("/", " / ")}</p>
      <dl className="job-facts">
        <dt>Confidence</dt>
        <dd>
          {distinct(call.alleles).map((allele) => (
            <div key={allele.structure}>
              {allele.structure}: {allele.confidence}%{" "}
              <span className="muted">{legacyConfidenceInWords(allele.confidence)}</span>
            </div>
          ))}
        </dd>
        {truth && (
          <>
            <dt>Truth</dt>
            <dd>
              {truth}{" "}
              {matches_truth === null ? null : matches_truth ? (
                <span className="success">✓ matches</span>
              ) : (
                <span className="error">✗ doesn't match</span>
              )}
            </dd>
          </>
        )}
        <dt>Flags</dt>
        <dd>
          {call.flags.length === 0 ? (
            "none"
          ) : (
            <ul className="flags">
              {call.flags.map((flag) => (
                <li key={flag}>
                  <span className="badge">{flag}</span> {LEGACY_FLAGS[flag]?.meaning ?? ""}
                </li>
              ))}
            </ul>
          )}
        </dd>
      </dl>
      <p className="muted">
        Legacy ScaleHD lists its primary allele (the one with most reads) first. Its confidence is a
        linear percentage based on iterative development. Over 80% is good, over 60% believable,
        and less may justify manual inspection.
      </p>
    </section>
  );
}

const ALLELE_ROWS: { name: string; value: (allele: LegacyAllele) => string }[] = [
  { name: "Status", value: (a) => a.status },
  { name: "CAG", value: (a) => String(a.cag ?? "–") },
  { name: "CCG", value: (a) => String(a.ccg ?? "–") },
  { name: "Intervening sequence", value: (a) => a.intervening },
  { name: "Confidence", value: (a) => `${a.confidence}%` },
  { name: "Confidence interval", value: (a) => a.confidence_interval || "–" },
  { name: "Backward slippage", value: (a) => number(a.backward_slippage) },
  { name: "Somatic mosaicism", value: (a) => number(a.somatic_mosaicism) },
  { name: "SNPs", value: (a) => a.variant_call.trim().split(/\s+/).join(", ") },
  { name: "SNP score", value: (a) => String(a.variant_score) },
];

function Alleles({ alleles }: Readonly<{ alleles: LegacyAllele[] }>) {
  return (
    <section>
      <h2>Alleles</h2>
      <AlleleTable alleles={alleles} rows={ALLELE_ROWS} />
      <p className="muted">
        Backward slippage: reads at N-1 and N-2 over those at N. Somatic mosaicism: reads at N+1
        to N+10 over those at N.
      </p>
    </section>
  );
}

function Mapping({ alleles }: Readonly<{ alleles: LegacyAllele[] }>) {
  const rows: typeof ALLELE_ROWS = [
    { name: "Forward reads mapped", value: (a) => mapped(a.map_fw, a.map_percent_fw) },
    { name: "Forward reads removed", value: (a) => a.purged_fw.toLocaleString() },
    { name: "Reverse reads mapped", value: (a) => mapped(a.map_rv, a.map_percent_rv) },
    { name: "Reverse reads removed", value: (a) => a.purged_rv.toLocaleString() },
  ];
  return (
    <section>
      <h2>Alignment</h2>
      <AlleleTable alleles={alleles} rows={rows} />
      <p className="muted">
        Reads removed: those mapping to more than one reference, unless the job kept them. Legacy
        ScaleHD considers under 90% of forward reads mapped significant.
      </p>
    </section>
  );
}

function AlleleTable({
  alleles,
  rows,
}: Readonly<{ alleles: LegacyAllele[]; rows: typeof ALLELE_ROWS }>) {
  return (
    <table>
      <thead>
        <tr>
          <th />
          {alleles.map((allele) => (
            <th key={allele.structure}>{allele.structure}</th>
          ))}
        </tr>
      </thead>
      <tbody>
        {rows.map((row) => (
          <tr key={row.name}>
            <td>{row.name}</td>
            {alleles.map((allele) => (
              <td key={allele.structure}>{row.value(allele)}</td>
            ))}
          </tr>
        ))}
      </tbody>
    </table>
  );
}

function CagChartView({ chart, scale }: Readonly<{ chart: LegacyCagChart; scale: Scale }>) {
  const bars = useMemo<Bar[]>(() => {
    const total = chart.reads.reduce((sum, n) => sum + n, 0);
    const [from, to] = shown(chart.reads, chart.called);
    return chart.reads.slice(from - 1, to).map((reads, i) => {
      const cag = from + i;
      const called = chart.called.includes(cag);
      return {
        label: String(cag),
        value: reads,
        highlight: called,
        tooltip: [
          `CAG ${cag}${called ? " (called allele)" : ""}`,
          `${aligned(reads)}, ${percent(reads / (total || 1))} at CCG ${chart.ccg}`,
        ],
      };
    });
  }, [chart]);
  return (
    <figure>
      <figcaption>
        {chart.alleles.join(" and ")} <span className="muted">(CCG {chart.ccg})</span>
      </figcaption>
      <BarChart
        bars={bars}
        scale={scale}
        xTitle="CAG length"
        yTitle="reads"
        label={`CAG distribution for ${chart.alleles.join(" and ")}`}
      />
    </figure>
  );
}

function CcgChartView({
  chart,
  call,
  scale,
  caption,
}: Readonly<{ chart: LegacyCcgChart; call: LegacyCall; scale: Scale; caption: boolean }>) {
  const bars = useMemo<Bar[]>(() => {
    const total = chart.reads.reduce((sum, n) => sum + n, 0);
    const called = new Set(call.alleles.map((allele) => allele.ccg));
    return chart.reads.map((reads, i) => ({
      label: String(i + 1),
      value: reads,
      highlight: called.has(i + 1),
      tooltip: [
        `CCG ${i + 1}${called.has(i + 1) ? " (called)" : ""}`,
        `${aligned(reads)}, ${percent(reads / (total || 1))}`,
      ],
    }));
  }, [chart, call]);
  return (
    <figure>
      {caption && <figcaption>{chart.alleles.join(" and ")}</figcaption>}
      <BarChart
        bars={bars}
        scale={scale}
        xTitle="CCG length"
        yTitle="reads"
        label={`CCG distribution for ${chart.alleles.join(" and ")}`}
      />
    </figure>
  );
}

function Flags({ call }: Readonly<{ call: LegacyCall }>) {
  return (
    <section>
      <h2>Legacy flags</h2>
      <table>
        <thead>
          <tr>
            <th>Flag</th>
            <th>Raised</th>
            <th>Significance</th>
            <th>Meaning</th>
          </tr>
        </thead>
        <tbody>
          {Object.entries(call.flag_values).map(([flag, raised]) => (
            <tr key={flag}>
              <td>{flag}</td>
              <td>{raised ? <strong>yes</strong> : "no"}</td>
              <td>{LEGACY_FLAGS[flag]?.significance || "–"}</td>
              <td className="muted">{LEGACY_FLAGS[flag]?.meaning ?? ""}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </section>
  );
}

function Run({ legacy }: Readonly<{ legacy: LegacySample }>) {
  const upstream = legacy.versions?.upstream;
  return (
    <section>
      <h2>The run</h2>
      <dl className="job-facts">
        <dt>Stages</dt>
        <dd>{legacy.stages.join(", ").replaceAll("_", " ") || "–"}</dd>
        <dt>Took</dt>
        <dd>
          {legacy.seconds === null ? "–" : `${Math.round(legacy.seconds)} s`}
          {legacy.threads !== null && (
            <span className="muted">
              {" "}
              on {legacy.threads} core{legacy.threads === 1 ? "" : "s"}
            </span>
          )}
        </dd>
        <dt>Legacy ScaleHD</dt>
        <dd>
          {legacy.versions?.scalehd ?? "–"}
          {upstream && (
            <span className="muted">
              {" "}
              ({upstream.url.replace("https://github.com/", "")} at{" "}
              <code>{upstream.commit.slice(0, 7)}</code>, {upstream.date}), Python{" "}
              {legacy.versions?.python}
            </span>
          )}
        </dd>
      </dl>
    </section>
  );
}

/** Subsection visualisation of alignments in-app */
function AlignedReads({ shown }: Readonly<{ shown: LegacyAlignedReads[] }>) {
  const alleles = shown.filter((a, i) => shown.findIndex((b) => b.allele === a.allele) === i);
  return (
    <section>
      <h2>Aligned reads</h2>
      <p className="muted">
        Forward reads aligned to each allele's reference, one read to a row and a block for each
        base. Scroll over the reads to zoom in and out, drag to move about once zoomed in.
      </p>
      <BaseKey />
      {alleles.map((allele) => (
        <figure key={allele.allele}>
          <figcaption>
            {allele.allele}{" "}
            {allele.reads !== null && (
              <span className="muted">
                A subsection of reads (of {allele.reads.toLocaleString()} total) reads aligned to{" "}
                {allele.contig}
              </span>
            )}
          </figcaption>
          {allele.notes.map((note) => (
            <p key={note} className="note">
              {note}
            </p>
          ))}
          {allele.sequences.length > 0 && (
            <ReadBlocks
              sequences={allele.sequences}
              label={`${allele.sequences.length} aligned reads for ${allele.allele}`}
            />
          )}
        </figure>
      ))}
    </section>
  );
}

/** Trimming (cutadapt) and FastQC, as legacy ScaleHD ran. */
function QualityControl({
  checks,
  scale,
}: Readonly<{ checks: LegacyQualityControl; scale: Scale }>) {
  return (
    <section>
      <h2>Quality control</h2>
      {checks.trimming.length > 0 && <Trimming trimming={checks.trimming} />}
      {checks.fastqc && <Fastqc fastqc={checks.fastqc} scale={scale} />}
    </section>
  );
}

function Trimming({ trimming }: Readonly<{ trimming: LegacyTrimming[] }>) {
  const byRead = trimming.map(
    (read) =>
      new Map([
        ...read.summary,
        ...read.adapter.map(([what, value]): [string, string] => [
          `Adapter ${what.toLowerCase()}`,
          value,
        ]),
      ]),
  );
  const labels = [...new Set(byRead.flatMap((values) => [...values.keys()]))];
  return (
    <>
      <h3>Trimming</h3>
      <p className="muted">cutadapt's summary for each read file.</p>
      <table>
        <thead>
          <tr>
            <th />
            {trimming.map((read) => (
              <th key={read.read}>
                {read.read === "R1" ? "Forward reads (R1)" : "Reverse reads (R2)"}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {labels.map((label) => (
            <tr key={label}>
              <td>{label}</td>
              {byRead.map((values, i) => (
                <td
                  key={trimming[i].read}
                  className={label === "Adapter sequence" ? "bases" : undefined}
                >
                  {values.get(label) ?? "–"}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </>
  );
}

function Fastqc({ fastqc, scale }: Readonly<{ fastqc: LegacyFastqc; scale: Scale }>) {
  const quality = useMemo<Bar[]>(
    () =>
      fastqc.per_base_quality.map((bin) => ({
        label: bin.base,
        value: bin.mean ?? 0,
        tooltip: [
          `Base ${bin.base}`,
          `mean ${figure(bin.mean)}, median ${figure(bin.median)}`,
          `middle half ${figure(bin.lower_quartile)}-${figure(bin.upper_quartile)}`,
          `10th-90th percentile ${figure(bin.low)}-${figure(bin.high)}`,
        ],
      })),
    [fastqc],
  );
  const lengths = useMemo<Bar[]>(
    () =>
      fastqc.lengths.map((bin) => ({
        label: bin.length,
        value: bin.reads ?? 0,
        tooltip: [`${bin.length} bases`, aligned(bin.reads ?? 0)],
      })),
    [fastqc],
  );
  return (
    <>
      <h3>FastQC</h3>
      <p className="muted">On the forward reads after trimming, as per legacy ScaleHD.</p>
      <ul className="fastqc-modules">
        {fastqc.modules.map((module) => (
          <li key={module.module}>
            <span className={`fastqc-status ${module.status}`}>{module.status || "–"}</span>{" "}
            {module.module}
          </li>
        ))}
      </ul>
      <dl className="job-facts">
        {fastqc.basic_statistics.map(([what, value]) => (
          <Fragment key={what}>
            <dt>{what}</dt>
            <dd>{value}</dd>
          </Fragment>
        ))}
      </dl>
      <figure>
        <figcaption>Quality along the reads</figcaption>
        <BarChart
          bars={quality}
          scale="linear"
          xTitle="position in read (bases)"
          yTitle="mean Phred quality"
          label="FastQC quality along the reads"
        />
      </figure>
      <figure>
        <figcaption>Read lengths</figcaption>
        <BarChart
          bars={lengths}
          scale={scale}
          xTitle="read length (bases)"
          yTitle="reads"
          label="FastQC read lengths"
        />
      </figure>
    </>
  );
}

function figure(value: number | null): string {
  return value === null ? "–" : String(Math.round(value * 10) / 10);
}

function shown(reads: number[], called: number[]): [number, number] {
  const present = reads.flatMap((n, i) => (n > 0 ? [i + 1] : []));
  const lengths = [...present, ...called];
  const from = Math.max(1, Math.min(...lengths) - 2);
  const to = Math.min(reads.length, Math.max(...lengths) + 2);
  return [from, to];
}

function distinct(alleles: LegacyAllele[]): LegacyAllele[] {
  return alleles.filter((a, i) => alleles.findIndex((b) => b.structure === a.structure) === i);
}

function mapped(count: number, share: number): string {
  return `${count.toLocaleString()} (${share.toFixed(1)}%)`;
}

function aligned(n: number): string {
  return `${n.toLocaleString()} read${n === 1 ? "" : "s"}`;
}

function percent(fraction: number): string {
  return `${(100 * fraction).toFixed(fraction < 0.01 ? 2 : 1)}%`;
}

function number(value: number): string {
  return typeof value === "number" ? value.toFixed(3) : String(value);
}
