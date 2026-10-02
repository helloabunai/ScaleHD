import { useMemo, useState } from "react";
import { Link, useParams } from "react-router";
import {
  api,
  type CagChart,
  type CalledAllele,
  type CcgBar,
  type GenotypeCall,
  type Reads,
  type SampleDetail,
  type SampleFile,
} from "../api";
import { type Bar, BarChart, type Scale } from "../charts/BarChart";
import { Heatmap } from "../charts/Heatmap";
import { FLAGS } from "../flags";
import { useApi } from "../useApi";
import { DemoTag } from "./jobDisplay";
import { Status } from "./Status";
import { StructureDiagrams } from "./StructureDiagram";

const FILE_NAMES: Record<SampleFile, string> = {
  call: "Call (JSON)",
  counts: "Molecule counts (JSON)",
  r1: "R1 FASTQ",
  r2: "R2 FASTQ",
};

export function SampleResults() {
  const params = useParams();
  const jobId = Number(params.jobId);
  const sampleId = Number(params.sampleId);
  const detail = useApi(() => api.getSample(jobId, sampleId), [jobId, sampleId], {
    every: 1000,
    while: (loaded) => ["queued", "running"].includes(loaded.sample.status),
  });
  const [scale, setScale] = useState<Scale>("linear");

  return (
    <Status of={detail}>
      {(detail) => (
        <>
          <SampleNav detail={detail} />
          <h1>
            {detail.sample.name} {detail.demo && <DemoTag />}
          </h1>
          {detail.call === null ? (
            <NotCalled detail={detail} />
          ) : (
            <>
              <Summary detail={detail} call={detail.call} />
              <ScaleToggle scale={scale} onChange={setScale} />
              <CagDistributions charts={detail.cag_charts} scale={scale} />
              <Structures alleles={detail.call.alleles} />
              <CcgDistribution bars={detail.ccg} call={detail.call} scale={scale} />
              <section>
                <h2>CAG and CCG</h2>
                <Heatmap
                  cells={detail.cells}
                  called={detail.call.alleles.map((a) => ({ cag: a.cag, ccg: a.ccg }))}
                  scale={scale}
                />
              </section>
              <Instability alleles={distinct(detail.call.alleles)} demo={detail.demo} />
              <Alternatives call={detail.call} />
            </>
          )}
          {detail.reads && <ReadSummary reads={detail.reads} />}
          <Downloads detail={detail} />
        </>
      )}
    </Status>
  );
}

function SampleNav({ detail }: Readonly<{ detail: SampleDetail }>) {
  const at = (id: number) => `/jobs/${detail.job_id}/samples/${id}`;
  return (
    <nav className="sample-nav">
      <Link to={`/jobs/${detail.job_id}`}>← {detail.job_name}</Link>
      <span>
        {detail.previous_id !== null && <Link to={at(detail.previous_id)}>Previous sample</Link>}
        {detail.next_id !== null && <Link to={at(detail.next_id)}>Next sample</Link>}
      </span>
    </nav>
  );
}

function NotCalled({ detail }: Readonly<{ detail: SampleDetail }>) {
  const { status, error } = detail.sample;
  if (status === "failed") return <p className="error">This sample failed: {error}</p>;
  if (status === "finished") return <p>This sample was counted but not called.</p>;
  return <p className="muted">This sample is {status}. Its results appear here when it finishes.</p>;
}

function Summary({ detail, call }: Readonly<{ detail: SampleDetail; call: GenotypeCall }>) {
  const { truth, matches_truth } = detail.sample;
  return (
    <section>
      <p className="genotype">{call.genotype.replace("/", " / ")}</p>
      <dl className="job-facts">
        <dt>Confidence</dt>
        <dd>
          Q {call.quality.toFixed(1)}, posterior {call.posterior.toFixed(4)}:{" "}
          <span className="muted">{inWords(call.quality)}</span>
        </dd>
        <dt>Molecules</dt>
        <dd>{call.molecules.toLocaleString()} used</dd>
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
                  <span className="badge">{flag}</span> {FLAGS[flag] ?? ""}
                </li>
              ))}
            </ul>
          )}
        </dd>
      </dl>
      <table>
        <thead>
          <tr>
            <th>Allele</th>
            <th>CAG</th>
            <th>Polyglutamine</th>
            <th>Share</th>
            <th>Molecules</th>
          </tr>
        </thead>
        <tbody>
          {call.alleles.map((allele, i) => (
            <tr key={i}>
              <td>{allele.structure}</td>
              <td>{cagLength(allele)}</td>
              <td>{polyglutamine(allele)}</td>
              <td>{percent(allele.fraction)}</td>
              <td>{Math.round(allele.molecules).toLocaleString()}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <p className="muted">
        Polyglutamine counts the glutamine codons: every CAG, plus two for each CAACAG (CAA also codes
        glutamine). A typical allele's polyglutamine length is its CAG + 2.
      </p>
    </section>
  );
}

function ScaleToggle({ scale, onChange }: Readonly<{ scale: Scale; onChange: (s: Scale) => void }>) {
  return (
    <div className="scale-toggle" role="group" aria-label="Chart scale">
      <span className="muted">Chart scale</span>
      {(["linear", "log"] as const).map((option) => (
        <button
          key={option}
          type="button"
          className={option === scale ? "selected" : ""}
          aria-pressed={option === scale}
          onClick={() => onChange(option)}
        >
          {option === "linear" ? "Linear" : "Log"}
        </button>
      ))}
    </div>
  );
}

function CagDistributions({ charts, scale }: Readonly<{ charts: CagChart[]; scale: Scale }>) {
  return (
    <section>
      <h2>CAG distribution</h2>
      <p className="muted">
        Molecules at each CAG length within each called allele's structure. Called lengths are
        highlighted. Striped bars are reads that ended inside the CAG tract, so their CAG is at
        least that long. Hover a column for details.
      </p>
      {charts.map((chart) => (
        <CagChartView key={`${chart.caacag}-${chart.ccgcca}-${chart.ccg}-${chart.cct}`} chart={chart} scale={scale} />
      ))}
    </section>
  );
}

function CagChartView({ chart, scale }: Readonly<{ chart: CagChart; scale: Scale }>) {
  const bars = useMemo<Bar[]>(() => {
    const total = chart.bars.reduce((sum, bar) => sum + bar.molecules + bar.lower_bound, 0);
    return chart.bars.map((bar) => {
      const called = chart.called.includes(bar.cag);
      return {
        label: String(bar.cag),
        value: bar.molecules,
        lowerBound: bar.lower_bound,
        highlight: called,
        tooltip: [
          `CAG ${bar.cag}${called ? " (called allele)" : ""}`,
          `${molecules(bar.molecules)}, ${percent(bar.molecules / (total || 1))} of this structure`,
          ...(bar.lower_bound
            ? [`${bar.lower_bound.toLocaleString()} more only known to be at least ${bar.cag}`]
            : []),
        ],
      };
    });
  }, [chart]);
  // An allele beyond read length is labelled like "83+_1_1_7_2": its CAG is a lower bound.
  const lowerBounds = new Set(
    chart.alleles.filter((label) => label.split("_")[0].endsWith("+")).map((label) => parseInt(label)),
  );
  const cag = chart.called.map((n) => (lowerBounds.has(n) ? `>=${n}` : String(n))).join(" and ");
  const name = `CAG ${cag}, CAACAG ${chart.caacag}, CCGCCA ${chart.ccgcca}, CCG ${chart.ccg}, CCT ${chart.cct}`;
  return (
    <figure>
      <figcaption>
        {chart.alleles.join(" and ")} <span className="muted">({name})</span>
      </figcaption>
      <BarChart
        bars={bars}
        scale={scale}
        xTitle="CAG length"
        yTitle="molecules"
        label={`CAG distribution for ${chart.alleles.join(" and ")}`}
      />
    </figure>
  );
}

function Structures({ alleles }: Readonly<{ alleles: CalledAllele[] }>) {
  return (
    <section>
      <h2>Repeat structure</h2>
      <p className="muted">
        Every tract drawn to one scale in bases. Hover a block for details; atypical tracts are
        outlined.
      </p>
      <StructureDiagrams alleles={distinct(alleles)} />
    </section>
  );
}

function CcgDistribution({ bars, call, scale }: Readonly<{ bars: CcgBar[]; call: GenotypeCall; scale: Scale }>) {
  const shown = useMemo<Bar[]>(() => {
    const total = bars.reduce((sum, bar) => sum + bar.molecules, 0);
    const called = new Set(call.alleles.map((allele) => allele.ccg));
    return bars.map((bar) => ({
      label: String(bar.ccg),
      value: bar.molecules,
      highlight: called.has(bar.ccg),
      tooltip: [
        `CCG ${bar.ccg}${called.has(bar.ccg) ? " (called)" : ""}`,
        `${molecules(bar.molecules)}, ${percent(bar.molecules / (total || 1))}`,
      ],
    }));
  }, [bars, call]);
  return (
    <section>
      <h2>CCG distribution</h2>
      <p className="muted">Complete molecules at each CCG length, across every structure.</p>
      <BarChart bars={shown} scale={scale} xTitle="CCG length" yTitle="molecules" label="CCG distribution" />
    </section>
  );
}

const INSTABILITY: { key: keyof CalledAllele; name: string; meaning: string }[] = [
  {
    key: "backward_slippage",
    name: "Backward slippage",
    meaning: "Molecules at N-1 and N-2, relative to those at N (as in ScaleHD 1.x).",
  },
  {
    key: "somatic_mosaicism",
    name: "Somatic mosaicism",
    meaning: "Molecules at N+1 to N+10, relative to those at N (as in ScaleHD 1.x).",
  },
  {
    key: "expansion_index",
    name: "Expansion index",
    meaning: "Average CAG gained above N, weighted by molecules.",
  },
  {
    key: "contraction_index",
    name: "Contraction index",
    meaning: "Average CAG lost below N, weighted by molecules (negative).",
  },
];

function Instability({ alleles, demo }: Readonly<{ alleles: CalledAllele[]; demo: boolean }>) {
  return (
    <section>
      <h2>Instability</h2>
      {demo && (
        <p className="note">
          Simulated sample: these numbers describe the simulator's stutter model, not biology.
        </p>
      )}
      <table>
        <thead>
          <tr>
            <th>Measure</th>
            {alleles.map((allele) => (
              <th key={allele.structure}>{allele.structure}</th>
            ))}
            <th>Meaning</th>
          </tr>
        </thead>
        <tbody>
          {INSTABILITY.map((row) => (
            <tr key={row.key}>
              <td>{row.name}</td>
              {alleles.map((allele) => (
                <td key={allele.structure}>{number(allele[row.key] as number | null)}</td>
              ))}
              <td className="muted">{row.meaning}</td>
            </tr>
          ))}
          {(["contraction", "expansion"] as const).map((side) => (
            <tr key={side}>
              <td>Fitted {side} stutter</td>
              {alleles.map((allele) => (
                <td key={allele.structure}>
                  {number(allele.stutter[side])}, {number(allele.stutter[`${side}_step`])},{" "}
                  {number(allele.stutter[`${side}_tail`])}
                </td>
              ))}
              <td className="muted">
                First step, next step and tail ratios of the model's stutter kernel.
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </section>
  );
}

function Alternatives({ call }: Readonly<{ call: GenotypeCall }>) {
  return (
    <section>
      <h2>Alternatives and unexplained peaks</h2>
      {call.alternatives.length === 0 ? (
        <p className="muted">No other genotype came close.</p>
      ) : (
        <table>
          <thead>
            <tr>
              <th>Next most likely genotype</th>
              <th>Probability</th>
            </tr>
          </thead>
          <tbody>
            {call.alternatives.slice(0, 5).map((alternative) => (
              <tr key={alternative.genotype}>
                <td>{alternative.genotype}</td>
                <td>{alternative.posterior > 0 ? alternative.posterior.toExponential(1) : "too small to show"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
      {call.unexplained.length > 0 && (
        <ul>
          {call.unexplained.map((peak) => (
            <li key={peak.structure}>
              Unexplained peak <strong>{peak.structure}</strong>:{" "}
              {molecules(peak.molecules)}
            </li>
          ))}
        </ul>
      )}
    </section>
  );
}

function ReadSummary({ reads }: Readonly<{ reads: Reads }>) {
  const outcomes = Object.entries(reads.read_outcomes).sort(([, a], [, b]) => b - a);
  const discordant = Object.entries(reads.discordant);
  return (
    <section>
      <h2>Reads</h2>
      <dl className="job-facts">
        <dt>Molecules</dt>
        <dd>
          {reads.molecules.toLocaleString()} ({reads.complete.toLocaleString()} complete,{" "}
          {reads.partial.toLocaleString()} partial)
        </dd>
        <dt>Dropped</dt>
        <dd>
          {reads.dropped.toLocaleString()}{" "}
          <span className="muted">read base pairings disagreed</span>
        </dd>
        <dt>Unusable</dt>
        <dd>{reads.unusable.toLocaleString()}</dd>
        <dt>Read outcomes</dt>
        <dd>
          {outcomes.map(([outcome, n]) => `${outcome} ${n.toLocaleString()}`).join(", ") || "–"}
        </dd>
        {discordant.length > 0 && (
          <>
            <dt>Disagreeing tracts</dt>
            <dd>{discordant.map(([tract, n]) => `${tract} ${n.toLocaleString()}`).join(", ")}</dd>
          </>
        )}
      </dl>
    </section>
  );
}

function Downloads({ detail }: Readonly<{ detail: SampleDetail }>) {
  return (
    <section>
      <h2>Files</h2>
      {detail.files.length === 0 ? (
        <p className="muted">No files yet.</p>
      ) : (
        <ul>
          {detail.files.map((file) => (
            <li key={file}>
              <a href={api.sampleFileUrl(detail.job_id, detail.sample.id, file)} download>
                {FILE_NAMES[file]}
              </a>
            </li>
          ))}
        </ul>
      )}
      {detail.folder && (
        <p className="muted">
          On the server: <code>{detail.folder}</code>
        </p>
      )}
    </section>
  );
}

/** Alleles without repeats: a homozygote's allele once. */
function distinct(alleles: CalledAllele[]): CalledAllele[] {
  return alleles.filter((allele, i) => alleles.findIndex((a) => a.structure === allele.structure) === i);
}

/** e.g. "43", or for an allele beyond read length ">=83 (roughly 95, 88-104)". */
function cagLength(allele: CalledAllele): string {
  if (!allele.beyond_read_length) return String(allele.cag);
  const estimate = allele.cag_estimate;
  return estimate
    ? `>=${allele.cag} (roughly ${estimate[0]}, ${estimate[1]}-${estimate[2]})`
    : `>=${allele.cag}`;
}

/** Beyond read length, the readable part gives a lower bound. */
function polyglutamine(allele: CalledAllele): string {
  if (allele.polyglutamine_length !== null) return String(allele.polyglutamine_length);
  return `>=${allele.cag + 2 * allele.caacag}`;
}

function inWords(confidence: number): string {
  if (confidence >= 99) return "the highest confidence shown (about 1 in 8 billion or better)";
  return `about 1 in ${Math.round(10 ** (confidence / 10)).toLocaleString()} chance it's wrong`;
}

function percent(fraction: number): string {
  return `${(100 * fraction).toFixed(fraction < 0.01 ? 2 : 1)}%`;
}

function number(value: number | null): string {
  return value === null ? "–" : value.toFixed(3);
}

export function molecules(n: number): string {
  return `${n.toLocaleString()} molecule${n === 1 ? "" : "s"}`;
}
