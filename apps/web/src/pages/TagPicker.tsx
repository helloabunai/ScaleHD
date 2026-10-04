import { type KeyboardEvent, useId, useState } from "react";
import { api, MAX_TAG_LENGTH, MAX_TAGS } from "../api";
import { useApi } from "../useApi";

export interface PickedTag {
  id: number | null;
  name: string;
}

export async function tagIds(tags: PickedTag[]): Promise<number[]> {
  const ids: number[] = [];
  for (const tag of tags) ids.push(tag.id ?? (await api.createTag(tag.name)).id);
  return ids;
}

export function TagPicker({
  value,
  onChange,
}: Readonly<{ value: PickedTag[]; onChange: (tags: PickedTag[]) => void }>) {
  const known = useApi(api.listTags);
  const [text, setText] = useState("");
  const listId = useId();
  const full = value.length >= MAX_TAGS;
  const existing = known.state === "done" ? known.data : [];
  const same = (a: string, b: string) => a.toLowerCase() === b.toLowerCase();

  function add() {
    const name = text.trim();
    if (!name || full) return;
    if (!value.some((tag) => same(tag.name, name))) {
      const tag = existing.find((t) => same(t.name, name));
      onChange([...value, tag ? { id: tag.id, name: tag.name } : { id: null, name }]);
    }
    setText("");
  }

  function onKeyDown(event: KeyboardEvent<HTMLInputElement>) {
    // Enter adds the tag rather than submitting the form around it.
    if (event.key === "Enter") {
      event.preventDefault();
      add();
    }
  }

  return (
    <fieldset className="tag-picker">
      <legend>
        Tags <span className="muted">up to {MAX_TAGS}, e.g. the paper or cohort</span>
      </legend>
      {value.length > 0 && (
        <ul className="chosen-tags">
          {value.map((tag) => (
            <li
              key={tag.name}
              className={tag.id === null ? "tag new" : "tag"}
              title={tag.id === null ? "New: made when this is saved" : undefined}
            >
              {tag.name}
              <button
                type="button"
                className="link"
                aria-label={`Remove ${tag.name}`}
                onClick={() => onChange(value.filter((t) => t.name !== tag.name))}
              >
                ×
              </button>
            </li>
          ))}
        </ul>
      )}
      <div className="tag-add">
        <input
          list={listId}
          aria-label="Tag"
          value={text}
          maxLength={MAX_TAG_LENGTH}
          disabled={full}
          placeholder={full ? `${MAX_TAGS} tags is the most` : "Pick a tag or type a new one"}
          onChange={(e) => setText(e.target.value)}
          onKeyDown={onKeyDown}
        />
        <span className="muted">
          {text.length}/{MAX_TAG_LENGTH}
        </span>
        <button type="button" onClick={add} disabled={full || !text.trim()}>
          Add tag
        </button>
      </div>
      <datalist id={listId}>
        {existing
          .filter((tag) => !value.some((t) => same(t.name, tag.name)))
          .map((tag) => (
            <option key={tag.id} value={tag.name} />
          ))}
      </datalist>
    </fieldset>
  );
}
