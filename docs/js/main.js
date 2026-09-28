/* main.js, register per-tab inits and boot. */
document.addEventListener("DOMContentLoaded", () => {
  Theme.init();
  Sidebar.init().catch(e => console.error("sidebar init failed", e));   // every tab, not just Overview
  TabManager.registerInit("overview", Overview.init);
  TabManager.registerInit("atlas", Atlas.init);
  TabManager.registerInit("cn", CN.init);
  TabManager.registerInit("concordance", Concordance.init);
  TabManager.registerInit("finder", Finder.init);
  TabManager.registerInit("about", About.init);
  TabManager.init();
});
