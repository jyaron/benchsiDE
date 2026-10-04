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

test("angled category labels stay inside the figure after the window is narrowed", async ({ page }) => {
  await page.setViewportSize({ width: 1400, height: 900 });
  await page.goto(PAGE);
  await page.evaluate(() => {
    const d = document.createElement("div"); d.id = "refittest"; d.style.width = "100%"; d.style.height = "420px";
    document.body.prepend(d);
    const cats = Array.from({ length: 15 }, (_, i) => "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_" + i);
    Plotly.react(d, [{ type: "heatmap", x: cats, y: cats, z: cats.map(() => cats.map(() => 1)) }],
      { margin: { t: 30, b: 130, l: 220, r: 10 }, height: 420, font: { size: 9 }, xaxis: { tickangle: 45 } }, { responsive: true });
  });
  await page.setViewportSize({ width: 700, height: 900 });
  await page.waitForTimeout(900);
  const over = await page.evaluate(() => {
    const d = document.getElementById("refittest"), box = d.querySelector(".main-svg").getBoundingClientRect();
    let worst = { px: -1e9, edge: "", label: "" };
    d.querySelectorAll(".xtick text").forEach((t) => {
      const r = t.getBoundingClientRect();
      for (const [edge, px] of [["right", r.right - box.right], ["left", box.left - r.left], ["bottom", r.bottom - box.bottom]])
        if (px > worst.px) worst = { px, edge, label: t.textContent };
    });
    return worst;
  });
  expect(over.px, `tick label "${over.label}" crosses the ${over.edge} edge by ${over.px} px`).toBeLessThanOrEqual(1);
});

test("FRY results replace the panels and notes of an earlier GSEA or ORA run", async ({ page }) => {
  await page.goto(PAGE);
  await page.evaluate(([m, d]) => {
    document.getElementById("species").value = "mouse"; useBuiltinAnno();
    document.getElementById("normsel").value = "tmm"; document.getElementById("fmode").value = "fbe";
    loadFromText(m, d); GROUPTYPE = "auto"; analyzeNow(); switchTab("enrich"); useBuiltin("hallmark_mm");
  }, [read("demo/GSE63310_counts.tsv"), read("demo/GSE63310_design.tsv")]);
  await page.click('#enrmethod button[data-e="gsea"]');
  await expect(page.locator("#gseapermwrap")).toBeVisible();
  await page.click("#enrrun");
  await expect(page.locator("#enrinfo")).toContainText("GSEA", { timeout: 60000 });
  await page.click('#enrmethod button[data-e="fry"]');
  await expect(page.locator("#gseapermwrap")).toBeHidden();
  await expect(page.locator("#enrmethodnote")).toContainText("FRY (limma)");
  await page.click("#enrrun");
  await expect(page.locator("#enrinfo")).toContainText("FRY self-contained test", { timeout: 60000 });
  expect(await page.locator("#p_leadedge").innerHTML()).toBe("");
  await expect(page.locator("#enrfootnote")).toBeHidden();
  await expect(page.locator("#enrtable")).toContainText("p (mixed)");
  // GSEA again, then ORA: the GSEA running-score plot must not remain under the ORA result
  await page.click('#enrmethod button[data-e="gsea"]'); await page.click("#enrrun");
  await expect(page.locator("#esplotwrap")).toBeVisible({ timeout: 60000 });
  await page.click('#enrmethod button[data-e="ora"]'); await page.click("#enrrun");
  await expect(page.locator("#enrinfo")).toContainText("Query:", { timeout: 60000 });
  await expect(page.locator("#esplotwrap")).toBeHidden();
  await expect(page.locator("#enrfootnote")).toBeVisible();
  // printing: the scrolling result table is printed in full
  await page.emulateMedia({ media: "print" });
  expect(await page.evaluate(() => getComputedStyle(document.getElementById("enrtable")).maxHeight)).toBe("none");
  await page.emulateMedia({ media: "screen" });
});
