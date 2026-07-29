/* theme.js — light/dark toggle, persisted. (No Plotly; charts are CSS.) */
const Theme = (() => {
  const KEY = "secacts-theme";
  const current = () => document.documentElement.getAttribute("data-theme") || "light";
  function apply(t) {
    document.documentElement.setAttribute("data-theme", t);
    try { localStorage.setItem(KEY, t); } catch (_) {}
    const lab = document.querySelector(".theme-toggle-label");
    if (lab) lab.textContent = t === "dark" ? "Light mode" : "Dark mode";
  }
  function toggle() { apply(current() === "dark" ? "light" : "dark"); }
  function init() {
    apply(current());
    document.getElementById("theme-toggle")?.addEventListener("click", toggle);
  }
  return { init, toggle, current };
})();
