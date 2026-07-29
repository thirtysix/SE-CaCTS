/* tab-manager.js — sidebar/tab nav, hashchange routing, lazy per-tab init. */
const TabManager = (() => {
  const inits = new Map(), initialized = new Set();
  let current = null;
  const TABS = ["overview", "atlas", "cn", "concordance", "finder", "about"];

  function registerInit(id, fn) { inits.set(id, fn); }

  function switchTab(id) {
    if (!TABS.includes(id)) id = "overview";
    if (current === id) return;
    document.querySelectorAll(".panel").forEach(p => p.classList.toggle("active", p.dataset.tab === id));
    document.querySelectorAll("[data-tab]").forEach(b => {
      if (b.tagName === "BUTTON") b.classList.toggle("active", b.dataset.tab === id);
    });
    current = id;
    if (history?.replaceState) history.replaceState(null, "", `#${id}`);
    if (!initialized.has(id) && inits.has(id)) {
      initialized.add(id);
      Promise.resolve(inits.get(id)()).catch(e => console.error(`tab ${id} init failed`, e));
    }
  }

  function init() {
    const bar = document.getElementById("tabbar");
    if (bar) bar.innerHTML = [...document.querySelectorAll(".nav button")]
      .map(b => `<button data-tab="${b.dataset.tab}">${b.textContent.trim()}</button>`).join("");
    // Only nav / tabbar BUTTONS get the intercepting handler. The panels are also [data-tab] sections;
    // binding them would let a click inside a panel bubble up and preventDefault() the Overview guide's
    // <a href="#tab"> links (they navigate via the hashchange listener below instead).
    document.querySelectorAll("button[data-tab]").forEach(el =>
      el.addEventListener("click", e => { e.preventDefault(); switchTab(el.dataset.tab); }));
    const fromHash = () => (location.hash || "#overview").slice(1) || "overview";
    switchTab(fromHash());
    window.addEventListener("hashchange", () => switchTab(fromHash()));
  }
  return { init, switchTab, registerInit };
})();
