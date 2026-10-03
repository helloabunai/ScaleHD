import { Fragment } from "react";
import type { CalledAllele } from "../api";

type TractKey = "cag" | "caacag" | "ccgcca" | "ccg" | "cct";

interface Tract {
  unit: string;
  key: TractKey;
  /** The count in a typical allele, or null where any count is typical (CAG, CCG). */
  typical: number | null;
  /** Part of the intervening sequence between the CAG and CCG tracts. */
  intervening?: boolean;
}

// (CAG)n (CAACAG)a (CCGCCA)b (CCG)m (CCT)k, in read order.
const TRACTS: Tract[] = [
  { unit: "CAG", key: "cag", typical: null },
  { unit: "CAACAG", key: "caacag", typical: 1, intervening: true },
  { unit: "CCGCCA", key: "ccgcca", typical: 1, intervening: true },
  { unit: "CCG", key: "ccg", typical: null },
  { unit: "CCT", key: "cct", typical: 2 },
];

/** Units a tract takes on the diagram i.e. its own, or the typical count where it has fewer. */
function drawnUnits(allele: CalledAllele, tract: Tract): number {
  return Math.max(allele[tract.key], tract.typical ?? 0);
}

/** Drawn bases from the start of the intervening sequence to the end of the CCT tract. */
function rightOfJunction(allele: CalledAllele): number {
  return TRACTS.slice(1).reduce(
    (sum, tract) => sum + drawnUnits(allele, tract) * tract.unit.length,
    0,
  );
}

/** Whether any allele is drawn with extra/missing units. */
function markings(alleles: CalledAllele[]): { inserted: boolean; deleted: boolean } {
  const differs = (by: (count: number, typical: number) => boolean) =>
    TRACTS.some((tract) => {
      const typical = tract.typical;
      return typical !== null && alleles.some((allele) => by(allele[tract.key], typical));
    });
  return {
    inserted: differs((count, typical) => count > typical),
    deleted: differs((count, typical) => count < typical),
  };
}

/**
 * Every allele's repeat structure to one scale, lined up where the intervening
 * sequence starts: the CAG tract grows left from there, CCG and CCT grow right.
 */
export function StructureDiagrams({ alleles }: Readonly<{ alleles: CalledAllele[] }>) {
  const left = Math.max(...alleles.map((allele) => allele.cag * 3));
  const total = left + Math.max(...alleles.map(rightOfJunction));
  return (
    <>
      <Legend {...markings(alleles)} />
      {alleles.map((allele) => (
        <figure key={allele.structure} className="structure">
          <figcaption>
            <strong>{allele.structure}</strong>
            {allele.polyglutamine_length !== null && (
              <span className="muted"> polyglutamine {allele.polyglutamine_length}</span>
            )}
            {!allele.typical && <span className="badge">atypical</span>}
          </figcaption>
          <Row allele={allele} left={left} total={total} />
        </figure>
      ))}
    </>
  );
}

function Row({ allele, left, total }: Readonly<{ allele: CalledAllele; left: number; total: number }>) {
  const at = (bases: number) => `${(100 * bases) / total}%`;
  let x = left - allele.cag * 3;
  const blocks = TRACTS.map((tract) => {
    const count = allele[tract.key];
    const length = tract.unit.length;
    const start = x;
    x += drawnUnits(allele, tract) * length;
    const missing = Math.max((tract.typical ?? 0) - count, 0);
    const extra = tract.typical === null ? 0 : Math.max(count - tract.typical, 0);
    const open = tract.key === "cag" && allele.beyond_read_length;
    const times = `${open ? ">=" : "×"}${count}`;
    // Extra CCT units get a block of their own, so the count label stays clear of the glow.
    const whole = extra > 0 && tract.intervening === true;
    const own = whole ? count : count - extra;
    return (
      <Fragment key={tract.key}>
        {own > 0 && (
          <div
            className={[
              "tract",
              `tract-${tract.key}`,
              open ? "open" : "",
              whole ? "tract-extra" : "",
            ].join(" ")}
            style={{ left: at(start), width: at(own * length) }}
            title={`${tract.unit} repeated ${count} time${count === 1 ? "" : "s"}${
              open ? " or more (reads ended inside this tract)" : ""
            }${whole ? `, ${extra} more than the typical structure (which has ${tract.typical})` : ""}: ${
              count * length
            } bases`}
          >
            <span className="full">
              {tract.unit} {times}
            </span>
            <span className="short">{times}</span>
          </div>
        )}
        {extra > 0 && !whole && (
          <div
            className={`tract tract-${tract.key} tract-extra`}
            style={{ left: at(start + own * length), width: at(extra * length) }}
            title={`${extra} more ${tract.unit} than the typical structure (which has ${tract.typical})`}
          >
            <span className="full">+{extra}</span>
            <span className="short">+{extra}</span>
          </div>
        )}
        {missing > 0 && (
          <div
            className="tract-missing"
            style={{ left: at(start + count * length), width: at(missing * length) }}
            title={`${missing} ${tract.unit} fewer than the typical structure (which has ${tract.typical})`}
          />
        )}
      </Fragment>
    );
  });
  return (
    <div className="structure-row" role="img" aria-label={describe(allele)}>
      {blocks}
      <div className="junction" style={{ left: at(left) }} aria-hidden />
    </div>
  );
}

/**
 * HTT structure key. atypical keys only drawn if present in allele(s)
 */
function Legend({ inserted, deleted }: Readonly<{ inserted: boolean; deleted: boolean }>) {
  return (
    <div className="structure-legend" aria-hidden>
      <span>HTT Repeat Units:</span>
      <span className="legend-group">
        {TRACTS.map((tract) => (
          <span key={tract.key}>
            <span className={`swatch tract-${tract.key}`} />
            {tract.unit}
          </span>
        ))}
      </span>
      {(deleted || inserted) && (
        <>
          <span>Intervening sequence structure:</span>
          <span className="legend-group">
            {deleted && (
              <span>
                <span className="swatch tract-missing" />
                atypical deletion
              </span>
            )}
            {inserted && (
              <span>
                <span className="swatch tract-extra" />
                atypical insertion
              </span>
            )}
          </span>
        </>
      )}
    </div>
  );
}

function describe(allele: CalledAllele): string {
  return TRACTS.map((tract) => {
    const open = tract.key === "cag" && allele.beyond_read_length;
    return `${tract.unit} ${open ? "at least " : "×"}${allele[tract.key]}`;
  }).join(", ");
}
