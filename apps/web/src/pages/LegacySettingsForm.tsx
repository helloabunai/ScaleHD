import { type ReactNode, useId, useState } from "react";
import type { AdapterPosition, LegacySettings } from "../api";

type Flags<K extends keyof LegacySettings> = LegacySettings[K];

// Error tolerances 1.x refuses although they're between 0 and 1 (REFUSED_ERROR_TOLERANCES
// in apps/server/src/scalehd_server/schemas.py).
const REFUSED_TOLERANCES = [0.35, 0.41, 0.47, 0.57, 0.69, 0.7, 0.82, 0.83, 0.94, 0.95];

/** ScaleHD 1.x's settings for a legacy job taken from XML config file. */
export function LegacySettingsForm({
  value,
  onChange,
  cores,
  revertTo,
  revertLabel,
}: Readonly<{
  value: LegacySettings;
  onChange: (settings: LegacySettings) => void;
  cores: number;
  /** Revert to defaults (either default defaults or user account defaults). */
  revertTo: LegacySettings;
  revertLabel: string;
}>) {
  function set<K extends keyof LegacySettings>(element: K, change: Partial<Flags<K>>) {
    const current = value[element];
    onChange({
      ...value,
      [element]: typeof current === "object" ? { ...current, ...change } : change,
    });
  }
  const stages = value.instance_flags;
  const demultiplexing = value.demultiplex_flags;
  const trimming = value.trim_flags;
  const alignment = value.alignment_flags;
  const snps = value.prediction_flags;
  return (
    <fieldset className="legacy-settings">
      <legend>Legacy ScaleHD settings</legend>
      <p className="muted">
        Written into a legacy ScaleHD configuration XML.
        Defaults are taken from ScaleHD documentation.
      </p>

      <CoresSlider
        value={value.cores_per_sample}
        cores={cores}
        onChange={(n) => onChange({ ...value, cores_per_sample: n })}
      />

      <Group title="Stages" open>
        <Toggle
          label="Demultiplex first"
          name="demultiplex"
          checked={stages.demultiplex}
          onChange={(on) => set("instance_flags", { demultiplex: on })}
          help="Keep only reads with the demultiplexing adapters below (batchadapt)."
        />
        <Toggle
          label="Quality control"
          name="quality_control"
          checked={stages.quality_control}
          onChange={(on) => set("instance_flags", { quality_control: on })}
          help="Trim reads (cutadapt, below) and run FastQC. Off, reads are aligned as they are."
        />
        <Toggle
          label="Alignment"
          name="sequence_alignment"
          checked
          disabled
          help="Align reads to the legacy reference library (bwa mem). Always on."
        />
        <Toggle
          label="Atypical realignment"
          name="atypical_realignment"
          checked={stages.atypical_realignment}
          onChange={(on) => set("instance_flags", { atypical_realignment: on })}
          help="Re-align an atypical allele's reads to a reference made for its structure. Off, legacy ScaleHD genotypes it from the typical alignment."
        />
        <Toggle
          label="Genotyping"
          name="genotype_prediction"
          checked
          disabled
          help="Always on."
        />
        <Toggle
          label="SNP calling"
          name="snp_calling"
          checked={stages.snp_calling}
          onChange={(on) => set("instance_flags", { snp_calling: on })}
          help="Call SNPs (freebayes, below)."
        />
      </Group>

      <Group title="Demultiplexing" tool="batchadapt" off={!stages.demultiplex}>
        <p className="muted">
          Keeps the reads in each read file with its adapter (cutadapt, each file on its own).
        </p>
        <div className="field-pair">
          <Bases
            label="Forward adapter"
            name="forward_adapter"
            value={demultiplexing.forward_adapter}
            onChange={(bases) => set("demultiplex_flags", { forward_adapter: bases })}
          />
          <Position
            label="Forward position"
            name="forward_position"
            value={demultiplexing.forward_position}
            onChange={(at) => set("demultiplex_flags", { forward_position: at })}
          />
        </div>
        <div className="field-pair">
          <Bases
            label="Reverse adapter"
            name="reverse_adapter"
            value={demultiplexing.reverse_adapter}
            onChange={(bases) => set("demultiplex_flags", { reverse_adapter: bases })}
          />
          <Position
            label="Reverse position"
            name="reverse_position"
            value={demultiplexing.reverse_position}
            onChange={(at) => set("demultiplex_flags", { reverse_position: at })}
          />
        </div>
        <div className="field-pair">
          <NumberField
            label="Error rate"
            name="error_rate"
            value={demultiplexing.error_rate}
            min={0}
            onChange={(n) => set("demultiplex_flags", { error_rate: n })}
            help="Whole numbers only in legacy ScaleHD: 0 means an exact match."
          />
          <NumberField
            label="Minimum overlap"
            name="min_overlap"
            value={demultiplexing.min_overlap}
            min={0}
            onChange={(n) => set("demultiplex_flags", { min_overlap: n })}
          />
        </div>
        <p className="muted">
          Either end ("AP") and minimum and maximum read lengths aren't offered: with any of
          them, legacy demultiplexing keeps no reads.
        </p>
      </Group>

      <Group title="Trimming" tool="cutadapt" off={!stages.quality_control}>
        <div className="field-pair">
          <Choice
            label="Trim"
            name="trim_type"
            value={trimming.trim_type}
            options={[
              ["Adapter", "Adapters"],
              ["Quality", "Low quality ends"],
              ["Both", "Both"],
            ]}
            onChange={(type) => set("trim_flags", { trim_type: type })}
          />
          <NumberField
            label="Quality threshold"
            name="quality_threshold"
            value={trimming.quality_threshold}
            min={0}
            max={38}
            onChange={(n) => set("trim_flags", { quality_threshold: n })}
          />
        </div>
        <Choice
          label="Adapter"
          name="adapter_flag"
          value={trimming.adapter_flag}
          options={[
            ["-a", "3' adapter (-a)"],
            ["-a$", "3' adapter, anchored (-a ...$)"],
            ["-g", "5' adapter (-g)"],
            ["-g^", "5' adapter, anchored (-g ^...)"],
            ["-b", "Either end (-b)"],
          ]}
          onChange={(flag) => set("trim_flags", { adapter_flag: flag })}
        />
        <Bases
          label="Forward adapter"
          name="forward_adapter"
          value={trimming.forward_adapter}
          onChange={(bases) => set("trim_flags", { forward_adapter: bases })}
        />
        <Bases
          label="Reverse adapter"
          name="reverse_adapter"
          value={trimming.reverse_adapter}
          onChange={(bases) => set("trim_flags", { reverse_adapter: bases })}
        />
        <NumberField
          label="Error tolerance"
          name="error_tolerance"
          value={trimming.error_tolerance}
          min={0}
          max={1}
          step={0.01}
          onChange={(n) => set("trim_flags", { error_tolerance: n })}
          warning={
            REFUSED_TOLERANCES.includes(trimming.error_tolerance)
              ? `Legacy ScaleHD refuses ${trimming.error_tolerance.toFixed(2)} (a rounding quirk in its check): try one either side.`
              : undefined
          }
        />
      </Group>

      <Group title="Alignment" tool="bwa mem">
        <div className="field-grid">
          {(
            [
              ["min_seed_length", "Minimum seed length", "-k", 1],
              ["band_width", "Band width", "-w", 1],
              ["seed_length_extension", "Re-seed longer than", "-r", 0.1],
              ["skip_seed_with_occurrence", "Skip seeds occurring more than", "-c", 1],
              ["chain_drop", "Drop chains shorter than", "-D", 0.01],
              ["seeded_chain_drop", "Discard chains seeded less than", "-W", 1],
              ["seq_match_score", "Match score", "-A", 1],
              ["mismatch_penalty", "Mismatch penalty", "-B", 1],
              ["unpaired_pairing_penalty", "Unpaired penalty", "-U", 1],
            ] as const
          ).map(([name, label, option, step]) => (
            <NumberField
              key={name}
              label={`${label} (${option})`}
              name={name}
              value={alignment[name]}
              min={0}
              step={step}
              onChange={(n) => set("alignment_flags", { [name]: n } as Partial<Flags<"alignment_flags">>)}
            />
          ))}
          {(
            [
              ["indel_penalty", "Gap open penalties (-O)", "insertions", "deletions"],
              ["gap_extend_penalty", "Gap extension penalties (-E)", "insertions", "deletions"],
              ["prime_clipping_penalty", "Clipping penalties (-L)", "5'", "3'"],
            ] as const
          ).map(([name, label, first, second]) => (
            <Pair
              key={name}
              label={label}
              name={name}
              value={alignment[name]}
              parts={[first, second]}
              onChange={(pair) =>
                set("alignment_flags", { [name]: pair } as Partial<Flags<"alignment_flags">>)
              }
            />
          ))}
        </div>
      </Group>

      <Group title="SNP calling" tool="freebayes" off={!stages.snp_calling}>
        <div className="field-pair">
          <NumberField
            label="Observation threshold"
            name="snp_observation_threshold"
            value={snps.snp_observation_threshold}
            min={1}
            max={5}
            onChange={(n) => set("prediction_flags", { snp_observation_threshold: n })}
            help="Reads that must show a variant (freebayes -C): 1 most lenient, 5 harshest."
          />
          <NumberField
            label="Quality cutoff"
            name="quality_cutoff"
            value={snps.quality_cutoff}
            min={0}
            onChange={(n) => set("prediction_flags", { quality_cutoff: n })}
            help="The lowest variant quality (VCF QUAL) kept."
          />
        </div>
      </Group>

      <Group title="Reads">
        <Toggle
          label="Keep reads that map to more than one reference"
          name="-e"
          checked={value.enshrine}
          onChange={(on) => onChange({ ...value, enshrine: on })}
          help="Legacy ScaleHD's -e (enshrine). Off, it removes them."
        />
        <Toggle
          label="Use every read"
          name="-b"
          checked={value.broadscope}
          onChange={(on) => onChange({ ...value, broadscope: on })}
          help="Legacy ScaleHD's -b (broadscope). Off, it subsamples a sample with 25,000 reads or more."
        />
      </Group>

      <p>
        <button
          type="button"
          className="link"
          disabled={sameSettings(value, revertTo)}
          onClick={() => onChange(revertTo)}
        >
          {revertLabel}
        </button>
      </p>
    </fieldset>
  );
}

export function sameSettings(a: LegacySettings, b: LegacySettings): boolean {
  return canonical(a) === canonical(b);
}

function canonical(value: unknown): string {
  return JSON.stringify(value, (_, v: unknown) =>
    v && typeof v === "object" && !Array.isArray(v)
      ? Object.fromEntries(Object.entries(v).sort(([a], [b]) => a.localeCompare(b)))
      : v,
  );
}

/** e.g. 24 cores = slider from 1 core 24 jobs. to 24 cores 1 job. */
function CoresSlider({
  value,
  cores,
  onChange,
}: Readonly<{ value: number; cores: number; onChange: (n: number) => void }>) {
  const id = useId();
  const used = Math.min(value, cores);
  const atOnce = Math.floor(cores / used);
  const plural = (n: number, word: string) => `${n} ${word}${n === 1 ? "" : "s"}`;
  return (
    <div className="cores">
      <label htmlFor={id}>Cores per sample</label>
      <input
        id={id}
        type="range"
        min={1}
        max={cores}
        step={1}
        value={used}
        onChange={(e) => onChange(e.target.valueAsNumber)}
      />
      <p>
        <strong>{plural(used, "core")}</strong> per sample, so up to{" "}
        <strong>{plural(atOnce, "sample")}</strong> at once{" "}
        <span className="muted">
          ({plural(cores, "core")} on this server
          {cores % used ? `, ${plural(cores % used, "core")} left idle` : ""})
        </span>
      </p>
      {value > cores && (
        <p className="note">Set to {value}, more than this server has: samples get {cores}.</p>
      )}
      <p className="muted">
        We recommend using fewer cores per sample for large batch jobs. 
        For small jobs, use more cores per sample.
      </p>
    </div>
  );
}

function Group({
  title,
  tool,
  off,
  open,
  children,
}: Readonly<{
  title: string;
  tool?: string;
  off?: boolean;
  open?: boolean;
  children: ReactNode;
}>) {
  const [opened, setOpened] = useState(open ?? false);
  return (
    <details
      className={off ? "file-group disabled" : "file-group"}
      open={opened && !off}
      onToggle={(e) => setOpened(e.currentTarget.open)}
    >
      <summary
        aria-disabled={off || undefined}
        onClick={(e) => {
          if (off) e.preventDefault();
        }}
      >
        {title}
        {tool && <span className="muted">{tool}</span>}
        {off && <span className="muted">(editable if selected)</span>}
      </summary>
      <div className="settings-group">{children}</div>
    </details>
  );
}

function Label({ label, name }: Readonly<{ label: string; name?: string }>) {
  return (
    <span className="setting-label">
      {label} {name && <code className="muted">{name}</code>}
    </span>
  );
}

function Toggle({
  label,
  name,
  checked,
  onChange,
  disabled,
  help,
}: Readonly<{
  label: string;
  name: string;
  checked: boolean;
  onChange?: (on: boolean) => void;
  disabled?: boolean;
  help?: string;
}>) {
  return (
    <label className="toggle">
      <input
        type="checkbox"
        checked={checked}
        disabled={disabled}
        onChange={(e) => onChange?.(e.target.checked)}
      />
      <span>
        <Label label={label} name={name} />
        {help && <span className="muted">{help}</span>}
      </span>
    </label>
  );
}

function NumberField({
  label,
  name,
  value,
  onChange,
  min,
  max,
  step = 1,
  help,
  warning,
}: Readonly<{
  label: string;
  name?: string;
  value: number;
  onChange: (n: number) => void;
  min?: number;
  max?: number;
  step?: number;
  help?: string;
  warning?: string;
}>) {
  const [text, setText] = useState(String(value));
  const shown = text === "" || Number(text) === value ? text : String(value);
  return (
    <label>
      <Label label={label} name={name} />
      <input
        type="number"
        value={shown}
        min={min}
        max={max}
        step={step}
        required
        onChange={(e) => {
          setText(e.target.value);
          const n = e.target.valueAsNumber;
          if (!Number.isNaN(n)) onChange(step === 1 ? Math.round(n) : n);
        }}
      />
      {help && <span className="muted">{help}</span>}
      {warning && <span className="note">{warning}</span>}
    </label>
  );
}

function Pair({
  label,
  name,
  value,
  parts,
  onChange,
}: Readonly<{
  label: string;
  name: string;
  value: [number, number];
  parts: readonly [string, string];
  onChange: (pair: [number, number]) => void;
}>) {
  const id = useId();
  return (
    <div className="pair" role="group" aria-labelledby={id}>
      <span id={id}>
        <Label label={label} name={name} />
      </span>
      {([0, 1] as const).map((i) => (
        <NumberField
          key={i}
          label={parts[i]}
          value={value[i]}
          min={0}
          onChange={(n) => onChange(i === 0 ? [n, value[1]] : [value[0], n])}
        />
      ))}
    </div>
  );
}

/** A, C, G, T . */
function Bases({
  label,
  name,
  value,
  onChange,
}: Readonly<{ label: string; name: string; value: string; onChange: (bases: string) => void }>) {
  return (
    <label>
      <Label label={label} name={name} />
      <input
        className="bases"
        value={value}
        required
        pattern="[ACGTacgt]+"
        title="A, C, G and T only"
        spellCheck={false}
        onChange={(e) => onChange(e.target.value.toUpperCase().trim())}
      />
    </label>
  );
}

function Position({
  label,
  name,
  value,
  onChange,
}: Readonly<{
  label: string;
  name: string;
  value: AdapterPosition;
  onChange: (at: AdapterPosition) => void;
}>) {
  return (
    <Choice
      label={label}
      name={name}
      value={value}
      options={[
        ["5P", "5' end"],
        ["3P", "3' end"],
      ]}
      onChange={onChange}
    />
  );
}

function Choice<T extends string>({
  label,
  name,
  value,
  options,
  onChange,
}: Readonly<{
  label: string;
  name: string;
  value: T;
  options: [T, string][];
  onChange: (value: T) => void;
}>) {
  return (
    <label>
      <Label label={label} name={name} />
      <select value={value} onChange={(e) => onChange(e.target.value as T)}>
        {options.map(([option, text]) => (
          <option key={option} value={option}>
            {text}
          </option>
        ))}
      </select>
    </label>
  );
}
