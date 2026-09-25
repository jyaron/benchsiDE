// benchsiDE CI: numerical validation and demo pipeline in real browser engines.
// The in-app self-test regenerates a deterministic dataset and compares 16 statistical
// components against embedded edgeR, limma, scipy and statsmodels reference values.
const { test, expect } = require("@playwright/test");
const path = require("path");
const fs = require("fs");

const PAGE = "file://" + path.resolve(__dirname, "..", "index.html");
const OFFLINE = "file://" + path.resolve(__dirname, "..", "dist", "benchside-offline.html");
const N_COMPONENTS = 16;

async function runSelfTestOn(page, url) {
  await page.goto(url);
  const res = await page.evaluate(() => runSelfTest());
  expect(res.rows.length).toBe(N_COMPONENTS);
  for (const row of res.rows) {
    expect(row.pass, `${row.name}: |delta| ${row.delta} > tol ${row.tol}`).toBe(true);
  }
  return res;
}

test("statistical self-test passes (CDN build)", async ({ page }) => {
  await runSelfTestOn(page, PAGE);
});

test("statistical self-test passes (offline build)", async ({ page }) => {
  test.skip(!fs.existsSync(path.resolve(__dirname, "..", "dist", "benchside-offline.html")),
            "offline build not present; run scripts/build_offline.py");
  await runSelfTestOn(page, OFFLINE);
});

// Pinned values were reproduced independently in R (edgeR calcNormFactors, filterByExpr or
// CPM filter, limma voom + eBayes on the Basal and LP samples; BH FDR <= 0.05, |log2FC| >= 1).
const DEMO_PINS = {
  cpm: { nG: 14490, up: 2410, dn: 2835, df0: 5.2046 },
  fbe: { nG: 16624, up: 2810, dn: 3338, df0: 5.4806 },
};

for (const mode of Object.keys(DEMO_PINS)) {
  test(`demo dataset reproduces pinned results (filter: ${mode})`, async ({ page }) => {
    await page.goto(PAGE);
    const matrix = fs.readFileSync(path.resolve(__dirname, "..", "demo", "GSE63310_counts.tsv"), "utf-8");
    const design = fs.readFileSync(path.resolve(__dirname, "..", "demo", "GSE63310_design.tsv"), "utf-8");
    const out = await page.evaluate(([m, d, fm]) => {
      document.getElementById("normsel").value = "tmm";
      document.getElementById("fmode").value = fm;
      loadFromText(m, d);
      GROUPTYPE = "auto"; analyzeNow();
      document.getElementById("selA").value = "Basal";
      document.getElementById("selB").value = "LP";
      document.getElementById("selFDR").value = "0.05";
      document.getElementById("fcthr").value = "1";
      deMethod = "voom"; deKey = "";
      const r = computeDE();
      let up = 0, dn = 0;
      for (let g = 0; g < nG; g++)
        if (r.q[g] <= 0.05 && Math.abs(r.lfc[g]) >= 1) { if (r.lfc[g] > 0) up++; else dn++; }
      return { nG, nS, groups: GROUPS.join(","), up, dn, df0: r.prior.df0 };
    }, [matrix, design, mode]);
    const pin = DEMO_PINS[mode];
    expect(out.nS).toBe(9);
    expect(out.groups).toBe("Basal,LP,ML");
    expect(out.nG).toBe(pin.nG);
    expect(out.up).toBe(pin.up);
    expect(out.dn).toBe(pin.dn);
    expect(Math.abs(out.df0 - pin.df0)).toBeLessThan(1e-3);
  });
}

test("UI smoke: self-test button renders a passing evidence panel", async ({ page }) => {
  await page.goto(PAGE);
  await page.click("#selftestbtn");
  await expect(page.locator("#selftestout")).toContainText(
    `All ${N_COMPONENTS} components reproduce their references`, { timeout: 30000 });
});
