// benchsiDE CI: user-interface checks for the features added at the design freeze, in real browser engines.
const { test, expect } = require("@playwright/test");
const path = require("path");
const fs = require("fs");
const ROOT = path.resolve(__dirname, "..");
const PAGE = "file://" + path.join(ROOT, "index.html");
const read = (p) => fs.readFileSync(path.join(ROOT, p), "utf8");

test("meta-analysis controls default to REML with the Hartung-Knapp test", async ({ page }) => {
  await page.goto(PAGE);
  expect(await page.locator("#metamethod").inputValue()).toBe("REML");
  expect(await page.locator("#metahk").inputValue()).toBe("hk");
});

test("interaction card: runs, draws, lists genes and exports", async ({ page }) => {
  await page.goto(PAGE);
  // a 2 x 2 layout built from the demo: cell type x an artificial batch (first vs second half of each cell type)
  const m = read("demo/GSE63310_counts.tsv");
  const d = read("demo/GSE63310_design.tsv").trim().split("\n");
  const seen = {}; const rows = ["sample\tcell"];
  d.slice(1).forEach((l) => { const [s, ct] = l.split("\t"); seen[ct] = (seen[ct] || 0) + 1; if (ct !== "ML") rows.push(`${s}\t${ct}_${seen[ct] <= 2 ? "b1" : "b2"}`); else rows.push(`${s}\t${ct}`); });
  await page.evaluate(([mm, dd]) => {
    document.getElementById("normsel").value = "tmm"; document.getElementById("fmode").value = "fbe";
    loadFromText(mm, dd); GROUPTYPE = "auto"; analyzeNow(); switchTab("de");
  }, [m, rows.join("\n") + "\n"]);
  await page.selectOption("#ixA1", "Basal_b1"); await page.selectOption("#ixB1", "LP_b1");
  await page.selectOption("#ixA2", "Basal_b2"); await page.selectOption("#ixB2", "LP_b2");
  await page.selectOption("#ixfit", "all");
  await page.click("#ixrun");
  await expect(page.locator("#ixinfo")).toContainText("A positive interaction log");
  await expect(page.locator("#p_ixvolcano .main-svg").first()).toBeVisible();
  expect(await page.locator("#ixtable tr.clickable").count()).toBeGreaterThan(10);
  const [dl] = await Promise.all([page.waitForEvent("download"), page.click("#ixexp")]);
  expect(dl.suggestedFilename()).toBe("interaction_contrast.csv");
  // compare with the engine directly
  const n = await page.evaluate(() => IXRES.lfc.length);
  expect(n).toBe(await page.evaluate(() => nG));
});
