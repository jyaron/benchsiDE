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

test("voom with a covariate and TREAT reproduce the limma workflow on the demo", async ({ page }) => {
  const errors = []; page.on("pageerror", (e) => errors.push(String(e)));
  await page.goto(PAGE);
  // sequencing lane as in Law et al. (F1000Research 5:1408, 2016)
  const lane = { GSM1545535: "L004", GSM1545536: "L004", GSM1545538: "L004", GSM1545539: "L006", GSM1545540: "L006", GSM1545541: "L006", GSM1545542: "L006", GSM1545544: "L008", GSM1545545: "L008" };
  const m = read("demo/GSE63310_counts.tsv");
  const d = read("demo/GSE63310_design.tsv").trim().split("\n");
  const rows = [d[0] + "\tlane"].concat(d.slice(1).map((l) => l + "\t" + lane[l.split("\t")[0]]));
  await page.evaluate(([mm, dd]) => {
    document.getElementById("species").value = "mouse"; useBuiltinAnno();
    document.getElementById("normsel").value = "tmm"; document.getElementById("fmode").value = "fbe";
    loadFromText(mm, dd); applyFactor("celltype"); GROUPTYPE = "auto"; analyzeNow(); switchTab("de");
  }, [m, rows.join("\n") + "\n"]);
  await page.click('#demethod button[data-m="voom"]');
  await page.evaluate(() => { DECOVSEL = ["lane"]; deKey = ""; });
  await page.selectOption("#selA", "Basal"); await page.selectOption("#selB", "LP");
  await page.selectOption("#selFDR", "0.05"); await page.fill("#fcthr", "1");
  await page.check("#detreat");
  await expect(page.locator("#decounts")).toContainText("TREAT");
  const res = await page.evaluate(() => { drawDE(); const r = computeDE(); let n = 0; for (let g = 0; g < nG; g++) if (r.q[g] <= 0.05) n++; return { n, treat: !!r.treat, cov: r.prior.cov }; });
  expect(res.treat).toBe(true);
  expect(res.cov).toContain("lane");
  expect(res.n).toBe(3647);   // limma 3.66: voom, ~0 + group + lane, contrasts.fit, treat(lfc = 1)
  // record what the export handler does, so that a missing download can be traced
  await page.evaluate(() => {
    window.__saved = []; window.__errs = [];
    window.addEventListener("error", (e) => window.__errs.push(String(e.message)));
    const orig = saveBlob;
    saveBlob = (b, f) => { window.__saved.push({ f, size: b.size, type: b.type }); return orig(b, f); };
  });
  await page.locator("#deexport").scrollIntoViewIfNeeded();
  const dlP = page.waitForEvent("download", { timeout: 20000 }).catch(() => null);
  await page.click("#deexport", { timeout: 15000 });   // fails here if the button cannot be clicked
  const dl = await dlP;
  const diag = await page.evaluate(() => { const r = computeDE(); return { saved: window.__saved, errs: window.__errs, deErr: r.err || null, method: deMethod, treat: DETREAT, cov: DECOVSEL }; });
  expect(errors, "page errors: " + errors.join("; ")).toEqual([]);
  expect(dl, "no download event; page state: " + JSON.stringify(diag)).not.toBeNull();
  const csv = fs.readFileSync(await dl.path(), "utf8");
  expect(csv).toContain("t_treat,p_treat,FDR_treat");
  await page.click('#demethod button[data-m="welch"]');
  expect(await page.locator("#detreat").isDisabled()).toBe(true);
});
