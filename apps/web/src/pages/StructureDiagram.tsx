import type { CalledAllele } from "../api";

type TractKey = "cag" | "caacag" | "ccgcca" | "ccg" | "cct";

interface Tract {
  unit: string;
  key: TractKey;
  /** The count in a typical allele, or null where any count is typical (CAG, CCG). */
  typical: number | null;
}

// (CAG)n (CAACAG)a (CCGCCA)b (CCG)m (CCT)k, in read order.
const TRACTS: Tract[] = [
  { unit: "CAG", key: "cag", typical: null },
  { unit: "CAACAG", key: "caacag", typical: 1 },
  { unit: "CCGCCA", key: "ccgcca", typical: 1 },
  { unit: "CCG", key: "ccg", typical: null },
  { unit: "CCT", key: "cct", typical: 2 },
];

/** Bases from the start of the intervening sequence to the end of the CCT tract. */
function rightOfJunction(allele: CalledAllele): number {
  return TRACTS.slice(1).reduce((sum, tract) => sum + allele[tract.key] * tract.unit.length, 0);
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
      <Legend />
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
          {!allele.typical && (
            <p className="muted structure-differences">Differs from typical: {differences(allele)}</p>
          )}
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
    const bases = count * tract.unit.length;
    const start = x;
    x += bases;
    if (count === 0) {
      return (
        <div
          key={tract.key}
          className="tract-missing"
          style={{ left: at(start) }}
          title={`no ${tract.unit} (a typical allele has ${tract.typical})`}
        />
      );
    }
    const open = tract.key === "cag" && allele.beyond_read_length;
    const atypical = tract.typical !== null && count !== tract.typical;
    const times = `${open ? ">=" : "×"}${count}`;
    return (
      <div
        key={tract.key}
        className={["tract", `tract-${tract.key}`, atypical ? " atypical" : "", open ? "open" : ""].join(" ")}
        style={{ left: at(start), width: at(bases) }}
        title={`${tract.unit} repeated ${count} time${count === 1 ? "" : "s"}${
          open ? " or more (reads ended inside this tract)" : ""
        }: ${bases} bases`}
      >
        <span className="full">
          {tract.unit} {times}
        </span>
        <span className="short">{times}</span>
      </div>
    );
  });
  return (
    <div className="structure-row" role="img" aria-label={describe(allele)}>
      {blocks}
      <div className="junction" style={{ left: at(left) }} aria-hidden />
    </div>
  );
}

function Legend() {
  return (
    <div className="structure-legend" aria-hidden>
      {TRACTS.map((tract) => (
        <span key={tract.key}>
          <span className={`swatch tract-${tract.key}`} />
          {tract.unit}
        </span>
      ))}
      <span>
        <span className="swatch strand" />
        rest of the strand, outside the repeat
      </span>
      <span className="muted">expected start of the intervening sequence,  visualisation aligned here</span>
    </div>
  );
}

function describe(allele: CalledAllele): string {
  return TRACTS.map((tract) => {
    const open = tract.key === "cag" && allele.beyond_read_length;
    return `${tract.unit} ${open ? "at least " : "×"}${allele[tract.key]}`;
  }).join(", ");
}

function differences(allele: CalledAllele): string {
  return TRACTS.filter((tract) => tract.typical !== null && allele[tract.key] !== tract.typical)
    .map((tract) => `${tract.unit} ×${allele[tract.key]} (typical ${tract.typical})`)
    .join(", ");
}
