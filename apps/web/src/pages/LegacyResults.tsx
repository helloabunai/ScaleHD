import { useMemo } from "react";
import type {
  LegacyAllele,
  LegacyCagChart,
  LegacyCall,
  LegacyCcgChart,
  LegacySample,
  SampleDetail,
} from "../api";
import { type Bar, BarChart, type Scale } from "../charts/BarChart";
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
        percentage: over 80% is good, over 60% believable, and less may justify manual
        inspection.
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
