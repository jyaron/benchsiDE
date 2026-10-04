const { defineConfig } = require("@playwright/test");
// Correctness tests in Chromium, Firefox and WebKit. The timing benchmark is separate
// (playwright.bench.config.js, run with: npm run bench).
module.exports = defineConfig({
  testDir: "./tests",
  testIgnore: "**/bench/**",
  timeout: 120000,
  retries: process.env.CI ? 1 : 0,
  reporter: process.env.CI ? [["list"], ["html", { open: "never" }]] : "list",
  projects: [
    { name: "chromium", use: { browserName: "chromium" } },
    { name: "firefox",  use: { browserName: "firefox" } },
    { name: "webkit",   use: { browserName: "webkit" } },
  ],
});
