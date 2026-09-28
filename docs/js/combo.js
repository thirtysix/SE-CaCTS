/* combo.js — a searchable dropdown (combobox): a visible, clickable list that also filters as you type.
   Adapted from pyCaCTS. Options carry a `search` string of aliases, and matching ignores case and
   punctuation, so "ovcar-3", "OVCAR3" and "NIH:OVCAR-3" all find NIH:OVCAR-3 (DepMap: NIHOVCAR3).
   Usage: const c = Combo.make(inputEl, key => onSelect(key)); c.setOptions([{key,label,search}]); c.setValue(key). */
const Combo = (() => {
  const CAP = 300;                                   // max rows rendered at once (the line level has 386)
  const norm = s => String(s || "").toLowerCase().replace(/[^a-z0-9]/g, "");

  function make(input, onSelect) {
    const wrap = input.closest(".combo"), pop = wrap.querySelector(".combo-pop");
    let opts = [], view = [], hi = -1, curKey = null, curLabel = "";

    const isOpen = () => wrap.classList.contains("open");
    function render(q) {
      const s = (q || "").trim().toLowerCase(), n = norm(q);
      view = s ? opts.filter(o => o._s.includes(s) || (n && o._n.includes(n))) : opts;
      const shown = view.slice(0, CAP);
      pop.innerHTML = (shown.length ? shown.map((o, i) =>
        `<div class="combo-opt${i === hi ? " hi" : ""}${o.key === curKey ? " on" : ""}" data-i="${i}" role="option" aria-selected="${o.key === curKey}">${U.esc(o.label)}</div>`).join("")
        : `<div class="combo-empty">no match</div>`)
        + (view.length > CAP ? `<div class="combo-empty">…and ${view.length - CAP} more, keep typing to narrow</div>` : "");
    }
    function open() { if (!isOpen()) { wrap.classList.add("open"); input.setAttribute("aria-expanded", "true"); } hi = -1; render(""); scrollToSel(); }
    function close() { wrap.classList.remove("open"); input.setAttribute("aria-expanded", "false"); input.value = curLabel; hi = -1; }
    function scrollToSel() { const el = pop.querySelector(".combo-opt.on") || pop.querySelector(".combo-opt.hi"); if (el) el.scrollIntoView({ block: "nearest" }); }
    function choose(o) { if (!o) return; curKey = o.key; curLabel = o.label; input.value = o.label; close(); onSelect(o.key); }

    input.addEventListener("focus", () => { input.select(); open(); });
    input.addEventListener("mousedown", () => { if (!isOpen()) setTimeout(open, 0); });   // click reopens after a close
    input.addEventListener("input", () => { wrap.classList.add("open"); input.setAttribute("aria-expanded", "true"); hi = 0; render(input.value); });
    input.addEventListener("keydown", e => {
      if (e.key === "ArrowDown") { e.preventDefault(); if (!isOpen()) return open(); hi = Math.min(hi + 1, view.length - 1); render(input.value); scrollToSel(); }
      else if (e.key === "ArrowUp") { e.preventDefault(); hi = Math.max(hi - 1, 0); render(input.value); scrollToSel(); }
      else if (e.key === "Enter") { e.preventDefault(); if (isOpen() && view[hi]) choose(view[hi]); }
      else if (e.key === "Escape") { close(); }
    });
    // mousedown (not click) so selecting an option fires before the input blurs
    pop.addEventListener("mousedown", e => { const d = e.target.closest(".combo-opt"); if (d) { e.preventDefault(); choose(view[+d.dataset.i]); } });
    document.addEventListener("click", e => { if (isOpen() && !wrap.contains(e.target)) close(); });

    return {
      setOptions(o) {
        opts = (o || []).map(x => ({ ...x, _s: String(x.search || x.label).toLowerCase(), _n: norm(x.search || x.label) }));
      },
      setValue(key) { const o = opts.find(x => x.key === key); curKey = o ? o.key : null; curLabel = o ? o.label : ""; input.value = curLabel; },
      getKey() { return curKey; },
    };
  }
  return { make };
})();
