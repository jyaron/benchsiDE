const { defineConfig } = require("@playwright/test");
// Timing benchmark: one browser at a time (workers: 1) so the runs do not compete for the CPU.
// Results: bench-results/<browser>_<date>.json. Run on the machine whose timings you want to report,
// plugged in, with other applications closed.
module.exports = defineConfig({
  testDir: "./tests/bench",
  workers: 1,
  fullyParallel: false,
  timeout: 3 * 60 * 60 * 1000,
  reporter: "list",
  projects: [
    { name: "chromium", use: { browserName: "chromium" } },
    { name: "firefox",  use: { browserName: "firefox" } },
    { name: "webkit",   use: { browserName: "webkit" } },
  ],
});
