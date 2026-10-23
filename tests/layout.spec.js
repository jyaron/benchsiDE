// benchsiDE CI: every checkbox, radio button, select and text box stays on the same line as the text that names it,
// on every tab and at three window widths (the line may wrap between controls, not inside one).
const { test, expect } = require("@playwright/test");
const path = require("path");
const fs = require("fs");
const ROOT = path.resolve(__dirname, "..");
const PAGE = "file://" + path.join(ROOT, "index.html");
const read = (p) => fs.readFileSync(path.join(ROOT, p), "utf8");

for (const width of [1400, 900, 560]) {
  test(`controls stay with their text at ${width} px`, async ({ page }) => {
    await page.setViewportSize({ width, height: 900 });
    await page.goto(PAGE);
    await page.evaluate(([m, d]) => {
      document.getElementById("normsel").value = "tmm"; document.getElementById("fmode").value = "fbe";
      loadFromText(m, d); GROUPTYPE = "auto"; analyzeNow();
    }, [read("demo/GSE63310_counts.tsv"), read("demo/GSE63310_design.tsv")]);
    const bad = [];
    for (const tab of ["overview", "gene", "de", "coexp", "heat", "enrich", "compare", "patterns", "disco", "traits"]) {
      await page.evaluate((t) => switchTab(t), tab);
      await page.waitForTimeout(250);
      const b = await page.evaluate(() => {
        const out = [], vis = (e) => { const r = e.getBoundingClientRect(); return r.width > 0 && r.height > 0; };
        const textRect = (n) => { const rg = document.createRange(); rg.selectNodeContents(n); const rs = rg.getClientRects(); return rs.length ? rs[0] : null; };
        const apart = (a, b) => Math.abs((a.top + a.bottom) / 2 - (b.top + b.bottom) / 2) > Math.max(a.height, b.height) * 0.75;
        for (const c of document.querySelectorAll("input[type=checkbox],input[type=radio]")) {
          if (!vis(c)) continue;
          const lab = c.closest("label"), n = lab ? [...lab.childNodes].find((x) => x.nodeType === 3 && x.nodeValue.trim()) : (c.nextSibling && c.nextSibling.nodeType === 3 && c.nextSibling.nodeValue.trim() ? c.nextSibling : null);
          if (!n) continue; const t = textRect(n);
          if (t && apart(c.getBoundingClientRect(), t)) out.push("checkbox " + (c.id || c.className) + ": " + n.nodeValue.trim().slice(0, 40));
        }
        for (const c of document.querySelectorAll(".ctl > select, .ctl > input[type=text]")) {
          if (!vis(c)) continue; const n = [...c.parentNode.childNodes].find((x) => x.nodeType === 3 && x.nodeValue.trim());
          const t = n && textRect(n); if (t && apart(c.getBoundingClientRect(), t)) out.push("control " + (c.id || c.className) + ": " + n.nodeValue.trim().slice(0, 40));
        }
        return out;
      });
      b.forEach((x) => bad.push(tab + " / " + x));
    }
    expect(bad).toEqual([]);
  });
}
