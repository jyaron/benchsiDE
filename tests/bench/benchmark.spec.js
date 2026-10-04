// benchsiDE timing benchmark in real browser engines (see playwright.bench.config.js).
// Uses the offline build, so no network request affects the timings.
const { test } = require("@playwright/test");
const path = require("path");
const fs = require("fs");

const ROOT = path.resolve(__dirname, "..", "..");
const OFFLINE = "file://" + path.join(ROOT, "dist", "benchside-offline.html");
const REPS = parseInt(process.env.BENCH_REPS || "5", 10);
const SIZES = (process.env.BENCH_SIZES || "6,24,96,192").split(",").map(Number);

test("timing benchmark", async ({ page, browserName }) => {
  test.skip(!fs.existsSync(path.join(ROOT, "dist", "benchside-offline.html")), "run: python3 scripts/build_offline.py");
  await page.goto(OFFLINE);
  const bver = page.context().browser().version();
  const { res, csv } = await page.evaluate(async ([sizes, reps, name, ver]) => {
    const r = await runBenchmark({ sizes, reps, warmup: 1, gseaPerm: 1000 });
    r.env.playwrightBrowser = name; r.env.playwrightBrowserVersion = ver; r.env.headless = true;
    return { res: r, csv: benchCSV(r) };
  }, [SIZES, REPS, browserName, bver]);
  const dir = path.join(ROOT, "bench-results");
  fs.mkdirSync(dir, { recursive: true });
  const file = path.join(dir, `${browserName}_${res.env.run_at.slice(0, 10)}.json`);
  fs.writeFileSync(file, JSON.stringify(res, null, 1));
  fs.writeFileSync(file.replace(/\.json$/, "_runs.csv"), csv);
  console.log(`${browserName} ${bver}: wrote ${file} and its _runs.csv`);
  for (const d of res.datasets) console.log(`  ${d.nSamples} samples: total ${(d.totalMedian / 1000).toFixed(2)} s`);
});
