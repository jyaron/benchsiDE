// benchsiDE CI: external validation on the public demo dataset (GEO GSE63310).
// The application's computations are compared with R reference values produced by
// validation/reference.R (edgeR, limma, RRHO, RNASeqPower). See VALIDATION.md.
const { test, expect } = require("@playwright/test");
const path = require("path");
const fs = require("fs");

const ROOT = path.resolve(__dirname, "..");
const PAGE = "file://" + path.join(ROOT, "index.html");
const read = (p) => fs.readFileSync(path.join(ROOT, p), "utf8");

test("GSE63310: application matches R reference values", async ({ page }) => {
  test.setTimeout(120000);
  await page.goto(PAGE);
  await page.addScriptTag({ content: read("validation/compare.js") });
  const res = await page.evaluate(([m, d, ref]) => runExternalValidation(m, d, JSON.parse(ref)),
    [read("demo/GSE63310_counts.tsv"), read("demo/GSE63310_design.tsv"), read("validation/gse63310_reference.json")]);
  expect(res.fixture.nG).toBe(16624);
  expect(res.rows.length).toBe(33);
  for (const row of res.rows) {
    expect(row.pass, `${row.name}: max deviation ${row.delta} exceeds ${row.tol}`).toBe(true);
  }
});
