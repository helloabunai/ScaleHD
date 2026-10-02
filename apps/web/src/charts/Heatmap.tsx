import { useId, useState } from "react";
import type { Cell } from "../api";
import type { Scale } from "./BarChart";

const WIDTH = 10;
const HEIGHT = 16;
const LEFT = 34;
const BOTTOM = 22;

/** Complete molecules by CAG (across) and CCG (down). Hover a cell for its count. */
export function Heatmap({
  cells,
  called,
  scale,
}: Readonly<{
  cells: Cell[];
  /** genotyped alleles' (CAG, CCG) are outlined. */
  called: { cag: number; ccg: number }[];
  scale: Scale;
}>) {
  const [hover, setHover] = useState<Cell | null>(null);
  const titleId = useId();
  if (cells.length === 0) return <p className="muted">No complete molecules.</p>;

  const cags = cells.map((cell) => cell.cag);
  const ccgs = cells.map((cell) => cell.ccg);
  const [cagFrom, cagTo] = [Math.min(...cags), Math.max(...cags)];
  const [ccgFrom, ccgTo] = [Math.min(...ccgs), Math.max(...ccgs)];
  const most = Math.max(...cells.map((cell) => cell.molecules));
  const shade = (n: number) =>
    scale === "log" ? Math.log1p(n) / Math.log1p(most) : n / most;
  const x = (cag: number) => LEFT + (cag - cagFrom) * WIDTH;
  const y = (ccg: number) => (ccg - ccgFrom) * HEIGHT;
  const width = LEFT + (cagTo - cagFrom + 1) * WIDTH;
  const height = (ccgTo - ccgFrom + 1) * HEIGHT + BOTTOM;
  const ccgRows = range(ccgFrom, ccgTo);
  const cagTicks = range(cagFrom, cagTo).filter((cag) => cag % 5 === 0);

  return (
    <div className="heatmap">
      <svg
        viewBox={`0 0 ${width} ${height}`}
        width={width * 1.4}
        height={height * 1.4}
        role="img"
        aria-labelledby={titleId}
        onMouseLeave={() => setHover(null)}
      >
        {/* The SVG's text alternative (as <svg> has no alt attribute). */}
        <title id={titleId}>
          Heatmap of molecules by CAG length (across, {cagFrom} to {cagTo}) and CCG length
          (down, {ccgFrom} to {ccgTo})
        </title>
        <rect
          x={LEFT}
          y={0}
          width={width - LEFT}
          height={height - BOTTOM}
          className="heatmap-empty"
        />
        {cells.map((cell) => (
          <rect
            key={`${cell.cag}-${cell.ccg}`}
            x={x(cell.cag)}
            y={y(cell.ccg)}
            width={WIDTH}
            height={HEIGHT}
            className="heatmap-cell"
            style={{ fillOpacity: Math.max(shade(cell.molecules), 0.06) }}
            onMouseEnter={() => setHover(cell)}
          >
            <title>{describe(cell)}</title>
          </rect>
        ))}
        {called.map((allele) => (
          <rect
            key={`called-${allele.cag}-${allele.ccg}`}
            x={x(allele.cag)}
            y={y(allele.ccg)}
            width={WIDTH}
            height={HEIGHT}
            className="heatmap-called"
          />
        ))}
        {ccgRows.map((ccg) => (
          <text key={ccg} x={LEFT - 4} y={y(ccg) + HEIGHT * 0.7} className="heatmap-label end">
            {ccg}
          </text>
        ))}
        {cagTicks.map((cag) => (
          <text key={cag} x={x(cag) + WIDTH / 2} y={height - 6} className="heatmap-label middle">
            {cag}
          </text>
        ))}
      </svg>
      <p className="muted heatmap-readout">
        {hover
          ? describe(hover)
          : "CAG across, CCG down. Hover a cell for its count; called alleles are outlined."}
      </p>
    </div>
  );
}

function range(from: number, to: number): number[] {
  return Array.from({ length: to - from + 1 }, (_, i) => from + i);
}

function describe(cell: Cell): string {
  const s = cell.molecules === 1 ? "" : "s";
  return `CAG ${cell.cag}, CCG ${cell.ccg}: ${cell.molecules.toLocaleString()} molecule${s}`;
}
