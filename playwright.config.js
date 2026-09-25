const { defineConfig } = require("@playwright/test");
module.exports = defineConfig({
  testDir: "./tests",
  timeout: 60000,
  projects: [
    { name: "chromium", use: { browserName: "chromium" } },
    { name: "firefox",  use: { browserName: "firefox" } },
    { name: "webkit",   use: { browserName: "webkit" } },
  ],
});
