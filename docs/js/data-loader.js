/* data-loader.js, fetch + parse TSV/JSON with caching (no backend). */
const DataLoader = (() => {
  const cache = new Map();

  function parseTSV(text) {
    const lines = text.replace(/\r/g, "").split("\n").filter(l => l.length);
    if (!lines.length) return { columns: [], rows: [] };
    const columns = lines[0].split("\t");
    // sample row 2 to detect numeric columns
    const sample = (lines[1] || "").split("\t");
    const numeric = columns.map((_, i) => i > 0 && sample[i] !== undefined && sample[i] !== "" && !isNaN(+sample[i]));
    const rows = new Array(lines.length - 1);
    for (let r = 1; r < lines.length; r++) {
      const cells = lines[r].split("\t"), o = {};
      for (let c = 0; c < columns.length; c++) {
        const v = cells[c];
        o[columns[c]] = (numeric[c] && v !== "" && v !== undefined) ? +v : v;
      }
      rows[r - 1] = o;
    }
    return { columns, rows };
  }

  async function loadTSV(path) {
    if (cache.has(path)) return cache.get(path);
    const p = fetch(path, { cache: "no-cache" }).then(async res => {
      if (!res.ok) throw new Error(`${path}: HTTP ${res.status}`);
      return parseTSV(await res.text());
    });
    cache.set(path, p);
    return p;
  }
  async function loadJSON(path) {
    if (cache.has(path)) return cache.get(path);
    const p = fetch(path, { cache: "no-cache" }).then(res => {
      if (!res.ok) throw new Error(`${path}: HTTP ${res.status}`);
      return res.json();
    });
    cache.set(path, p);
    return p;
  }
  return { loadTSV, loadJSON, parseTSV };
})();
