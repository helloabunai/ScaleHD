import {
  BarController,
  BarElement,
  CategoryScale,
  Chart,
  type ChartConfiguration,
  LinearScale,
  LogarithmicScale,
  type TooltipModel,
  Tooltip,
} from "chart.js";
import { useEffect, useRef } from "react";
import { useTheme } from "../theme";

Chart.register(BarController, BarElement, CategoryScale, LinearScale, LogarithmicScale, Tooltip);

export type Scale = "linear" | "log";

export interface Bar {
  /** The x-axis label, e.g. a CAG length. */
  label: string;
  /** Molecules seen at exactly this length. */
  value: number;
  /** Molecules only known to be at least this long (reads ended inside the tract). */
  lowerBound?: number;
  /** A called allele's length. */
  highlight?: boolean;
  /** Lines shown when the mouse is over this bar's column. */
  tooltip: string[];
}

/** An interactive bar chart: hover a column for its numbers; linear or log y axis. */
export function BarChart({
  bars,
  scale,
  xTitle,
  yTitle,
  label,
}: Readonly<{
  bars: Bar[];
  scale: Scale;
  xTitle: string;
  yTitle: string;
  /** What the chart shows, for screen readers. */
  label: string;
}>) {
  const canvas = useRef<HTMLCanvasElement>(null);
  const tip = useRef<HTMLDivElement>(null);
  const { shown } = useTheme();

  useEffect(() => {
    if (!canvas.current) return;
    const colors = themeColors();
    // A log axis can't show zero, so empty bars are simply left out.
    const shown = (n: number | undefined) => (scale === "log" && !n ? null : (n ?? 0));
    const config: ChartConfiguration<"bar", (number | null)[], string> = {
      type: "bar",
      data: {
        labels: bars.map((bar) => bar.label),
        datasets: [
          {
            label: "molecules",
            data: bars.map((bar) => shown(bar.value)),
            backgroundColor: bars.map((bar) => (bar.highlight ? colors.accent : colors.muted)),
            maxBarThickness: 48,
            grouped: false,
          },
          {
            label: "at least this long",
            data: bars.map((bar) => shown(bar.lowerBound)),
            backgroundColor: hatch(colors.warning),
            borderColor: colors.warning,
            borderWidth: 1,
            maxBarThickness: 48,
            grouped: false,
          },
        ],
      },
      options: {
        animation: false,
        responsive: true,
        maintainAspectRatio: false,
        interaction: { mode: "index", intersect: false },
        scales: {
          x: {
            title: { display: true, text: xTitle, color: colors.text },
            ticks: { color: colors.text, autoSkip: true, maxRotation: 0 },
            grid: { display: false },
          },
          y: {
            type: scale === "log" ? "logarithmic" : "linear",
            title: { display: true, text: yTitle, color: colors.text },
            ticks: { color: colors.text },
            grid: { color: colors.line },
            ...(scale === "linear" ? { beginAtZero: true } : {}),
          },
        },
        plugins: {
          legend: { display: false },
          tooltip: {
            enabled: false,
            external: ({ tooltip }: { tooltip: TooltipModel<"bar"> }) =>
              showTip(tip.current, tooltip, bars),
          },
        },
      },
    };
    const chart = new Chart(canvas.current, config);
    return () => chart.destroy();
  }, [bars, scale, xTitle, yTitle, shown]);

  return (
    <div className="chart" role="img" aria-label={label}>
      <canvas ref={canvas} />
      <div ref={tip} className="chart-tip" hidden />
    </div>
  );
}

function showTip(element: HTMLDivElement | null, tooltip: TooltipModel<"bar">, bars: Bar[]) {
  if (!element) return;
  const index = tooltip.dataPoints?.[0]?.dataIndex;
  if (tooltip.opacity === 0 || index === undefined) {
    element.hidden = true;
    return;
  }
  element.replaceChildren(
    ...bars[index].tooltip.map((line, i) => {
      const row = document.createElement("div");
      row.textContent = line;
      if (i === 0) row.className = "chart-tip-title";
      return row;
    }),
  );
  element.hidden = false;
  element.style.left = `${tooltip.caretX}px`;
  element.style.top = `${tooltip.caretY}px`;
}

/** Diagonal stripes, for molecules whose length is only a lower bound. */
function hatch(color: string): CanvasPattern | string {
  const tile = document.createElement("canvas");
  tile.width = tile.height = 8;
  const context = tile.getContext("2d");
  if (!context) return color;
  context.strokeStyle = color;
  context.lineWidth = 2;
  context.beginPath();
  context.moveTo(0, 8);
  context.lineTo(8, 0);
  context.stroke();
  return context.createPattern(tile, "repeat") ?? color;
}

function themeColors() {
  const style = getComputedStyle(document.documentElement);
  const read = (name: string, fallback: string) => style.getPropertyValue(name).trim() || fallback;
  return {
    accent: read("--accent", "#36c"),
    muted: read("--muted", "#777"),
    warning: read("--warning", "#8a5300"),
    line: read("--line", "#8884"),
    text: getComputedStyle(document.body).color || "#000",
  };
}
