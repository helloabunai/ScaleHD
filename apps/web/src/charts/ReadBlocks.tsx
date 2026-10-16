import {
  type PointerEvent,
  useCallback,
  useEffect,
  useLayoutEffect,
  useRef,
  useState,
} from "react";
import { useTheme } from "../theme";

const CELL = 3;
const ROW = 8;
const GAP = 1;
const MOST_ZOOM = 8;
const WHEEL_ZOOM = 0.002;
const BUTTON_ZOOM = 1.5;
const LETTER_WIDTH = 10;

export const BASES: { base: string; name: string }[] = [
  { base: "A", name: "adenine" },
  { base: "C", name: "cytosine" },
  { base: "G", name: "guanine" },
  { base: "T", name: "thymine" },
  { base: "N", name: "unknown" },
];

interface View {
  zoom: number;
  x: number;
  y: number;
}

const ALL: View = { zoom: 0, x: 0, y: 0 };

/** Reads drawn as rows of coloured blocks, one block per base.
 *
 * Wheel to zoom, drag to pan, and hover a block for info.
 */
export function ReadBlocks({
  sequences,
  label,
}: Readonly<{
  sequences: string[];
  label: string;
}>) {
  const frame = useRef<HTMLDivElement>(null);
  const canvas = useRef<HTMLCanvasElement>(null);
  const tip = useRef<HTMLDivElement>(null);
  const view = useRef<View>(ALL);
  const drag = useRef<{ x: number; y: number; from: View } | null>(null);
  const { shown } = useTheme();
  const [space, setSpace] = useState(0);
  const [zoom, setZoom] = useState(0);
  const [dragging, setDragging] = useState(false);
  const width = Math.max(0, ...sequences.map((sequence) => sequence.length)) * CELL;
  const height = sequences.length * (ROW + GAP);
  const fit = space > 0 && width > 0 ? Math.min(1, space / width) : 1;
  const bounds: Bounds = { fit, width, height, frameWidth: space, frameHeight: height * fit };

  useLayoutEffect(() => {
    const element = frame.current;
    if (!element) return;
    const measure = () => setSpace(element.clientWidth);
    measure();
    const observer = new ResizeObserver(measure);
    observer.observe(element);
    return () => observer.disconnect();
  }, []);

  const draw = useCallback(() => {
    const element = canvas.current;
    const context = element?.getContext("2d");
    if (!element || !context || bounds.frameWidth === 0) return;
    const at = clamp(view.current, bounds);
    const ratio = window.devicePixelRatio || 1;
    element.width = Math.round(bounds.frameWidth * ratio);
    element.height = Math.round(bounds.frameHeight * ratio);
    context.scale(ratio, ratio);
    // Edges on whole screen pixels, so one run of blocks meets the next without a seam.
    const snap = (value: number) => Math.round(value * ratio) / ratio;
    const colors = baseColors();
    const cell = CELL * at.zoom;
    const pitch = (ROW + GAP) * at.zoom;
    const letters = cell >= LETTER_WIDTH;
    if (letters) {
      const size = Math.round(Math.min(cell * 0.9, ROW * at.zoom * 0.7));
      context.font = `${size}px ui-monospace, monospace`;
      context.textAlign = "center";
      context.textBaseline = "middle";
    }
    // Only what's in view.
    const firstRow = Math.max(0, Math.floor(-at.y / pitch));
    const lastRow = Math.min(sequences.length, Math.ceil((bounds.frameHeight - at.y) / pitch));
    const firstBase = Math.max(0, Math.floor(-at.x / cell));
    for (let row = firstRow; row < lastRow; row++) {
      const sequence = sequences[row];
      const lastBase = Math.min(sequence.length, Math.ceil((bounds.frameWidth - at.x) / cell));
      const top = snap(at.y + row * pitch);
      const bottom = snap(at.y + row * pitch + ROW * at.zoom);

      let start = firstBase;
      for (let i = start + 1; i <= lastBase; i++) {
        if (i < lastBase && sequence[i] === sequence[start]) continue;
        const left = snap(at.x + start * cell);
        context.fillStyle = colors.get(sequence[start]) ?? colors.get("N") ?? "#999";
        context.fillRect(left, top, snap(at.x + i * cell) - left, bottom - top);
        start = i;
      }
      if (!letters) continue;
      context.fillStyle = colors.get("letter") ?? "#111";
      for (let i = firstBase; i < lastBase; i++) {
        context.fillText(sequence[i], at.x + (i + 0.5) * cell, (top + bottom) / 2);
      }
    }
    setZoom(at.zoom);
  }, [sequences, bounds.fit, bounds.frameWidth, bounds.frameHeight, width, height, shown]);

  useEffect(draw, [draw]);

  /** Show this view, as far as it goes (it's clamped when drawn). */
  function show(next: View) {
    view.current = next;
    draw();
  }

  /** Zoom by this much about a point in the box (its middle if none). */
  function zoomBy(factor: number, x = bounds.frameWidth / 2, y = bounds.frameHeight / 2): boolean {
    const at = clamp(view.current, bounds);
    const zoom = Math.min(MOST_ZOOM, Math.max(bounds.fit, at.zoom * factor));
    if (zoom === at.zoom) return false;
    const scale = zoom / at.zoom;

    show(zoom <= bounds.fit ? ALL : { zoom, x: x - (x - at.x) * scale, y: y - (y - at.y) * scale });
    return true;
  }

  // Not a React onWheel as it's passive, so couldn't stop the page scrolling instead.
  useEffect(() => {
    const element = canvas.current;
    if (!element) return;
    function wheel(event: WheelEvent) {
      const box = (event.currentTarget as HTMLCanvasElement).getBoundingClientRect();
      const lines = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? box.height : 1;
      const factor = Math.exp(-event.deltaY * lines * WHEEL_ZOOM);
      // Already zoomed as far as it goes = the page scrolls as usual.
      if (zoomBy(factor, event.clientX - box.left, event.clientY - box.top)) {
        event.preventDefault();
        hideTip();
      }
    }
    element.addEventListener("wheel", wheel, { passive: false });
    return () => element.removeEventListener("wheel", wheel);
  });

  function press(event: PointerEvent<HTMLCanvasElement>) {
    if (event.button !== 0 || clamp(view.current, bounds).zoom <= bounds.fit) return;
    event.currentTarget.setPointerCapture(event.pointerId);
    drag.current = { x: event.clientX, y: event.clientY, from: clamp(view.current, bounds) };
    setDragging(true);
    hideTip();
  }

  function move(event: PointerEvent<HTMLCanvasElement>) {
    const from = drag.current;
    if (!from) {
      hover(event);
      return;
    }
    show({
      ...from.from,
      x: from.from.x + event.clientX - from.x,
      y: from.from.y + event.clientY - from.y,
    });
  }

  function release() {
    drag.current = null;
    setDragging(false);
  }

  function hover(event: PointerEvent<HTMLCanvasElement>) {
    const element = tip.current;
    if (!element) return;
    const at = clamp(view.current, bounds);
    const box = event.currentTarget.getBoundingClientRect();
    const pitch = (ROW + GAP) * at.zoom;
    const row = Math.floor((event.clientY - box.top - at.y) / pitch);
    const position = Math.floor((event.clientX - box.left - at.x) / (CELL * at.zoom));
    const base = sequences[row]?.[position];
    if (base === undefined) {
      element.hidden = true;
      return;
    }
    const name = BASES.find((known) => known.base === base)?.name ?? "unknown";
    const title = document.createElement("div");
    title.className = "chart-tip-title";
    title.textContent = `${base} (${name})`;
    const where = document.createElement("div");
    where.textContent = `read ${row + 1}, base ${position + 1} of ${sequences[row].length}`;
    element.replaceChildren(title, where);
    element.hidden = false;

    element.style.left = `${event.clientX}px`;
    element.style.top = `${box.top + Math.max(0, at.y + row * pitch)}px`;
  }

  function hideTip() {
    if (tip.current) tip.current.hidden = true;
  }

  const zoomed = zoom > fit;
  return (
    <div className="read-blocks">
      <div className="scale-toggle read-blocks-zoom" role="group" aria-label="Zoom">
        <span className="muted">Zoom</span>
        <button
          type="button"
          aria-label="Zoom in"
          disabled={zoom >= MOST_ZOOM}
          onClick={() => zoomBy(BUTTON_ZOOM)}
        >
          +
        </button>
        <button
          type="button"
          aria-label="Zoom out"
          disabled={!zoomed}
          onClick={() => zoomBy(1 / BUTTON_ZOOM)}
        >
          -
        </button>
        <button type="button" disabled={!zoomed} onClick={() => show(ALL)}>
          All reads
        </button>
      </div>
      <div className="read-blocks-frame">
        <div ref={frame}>
          <canvas
            ref={canvas}
            role="img"
            aria-label={label}
            className={dragging ? "dragging" : zoomed ? "zoomed" : undefined}
            style={{ width: `${bounds.frameWidth}px`, height: `${bounds.frameHeight}px` }}
            onPointerDown={press}
            onPointerMove={move}
            onPointerUp={release}
            onPointerCancel={release}
            onPointerLeave={hideTip}
          />
        </div>
        <div ref={tip} className="chart-tip read-tip" hidden />
      </div>
    </div>
  );
}

export function BaseKey() {
  return (
    <ul className="base-key" aria-label="Key">
      {BASES.map(({ base, name }) => (
        <li key={base}>
          <span className={`swatch base-${base}`} /> {base} <span className="muted">{name}</span>
        </li>
      ))}
    </ul>
  );
}

interface Bounds {
  fit: number;
  width: number;
  height: number;
  frameWidth: number;
  frameHeight: number;
}

function clamp(view: View, bounds: Bounds): View {
  const zoom = Math.min(MOST_ZOOM, Math.max(bounds.fit, view.zoom));
  return {
    zoom,
    x: Math.min(0, Math.max(bounds.frameWidth - bounds.width * zoom, view.x)),
    y: Math.min(0, Math.max(bounds.frameHeight - bounds.height * zoom, view.y)),
  };
}

function baseColors(): Map<string, string> {
  const style = getComputedStyle(document.documentElement);
  return new Map(
    [...BASES.map(({ base }) => base), "letter"].map((name) => [
      name,
      style.getPropertyValue(`--base-${name.toLowerCase()}`).trim() || "#999",
    ]),
  );
}
