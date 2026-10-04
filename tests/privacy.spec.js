// benchsiDE CI: the page makes no network requests while it analyses data, and its
// Content-Security-Policy blocks any attempt to open one. Both builds, all three engines.
const { test, expect } = require("@playwright/test");
const path = require("path");
const fs = require("fs");

const ROOT = path.resolve(__dirname, "..");
const BUILDS = {
  online: { url: "file://" + path.join(ROOT, "index.html"), allowed: [/^https:\/\/cdn\.plot\.ly\/plotly-2\.32\.0\.min\.js$/] },
  offline: { url: "file://" + path.join(ROOT, "dist", "benchside-offline.html"), allowed: [] },
};
const read = (p) => fs.readFileSync(path.join(ROOT, p), "utf8");

for (const [name, b] of Object.entries(BUILDS)) {
  test(`no network requests during a full analysis (${name} build)`, async ({ page }) => {
    test.skip(name === "offline" && !fs.existsSync(path.join(ROOT, "dist", "benchside-offline.html")), "offline build not present");
    const requests = [];
    page.on("request", (r) => { const u = r.url(); if (!u.startsWith("file:") && !u.startsWith("data:") && !u.startsWith("blob:")) requests.push(u); });
    await page.addInitScript(() => {
      window.__cspViolations = [];
      document.addEventListener("securitypolicyviolation", (e) => window.__cspViolations.push(e.violatedDirective + " " + e.blockedURI));
    });
    await page.goto(b.url);
    const policy = await page.evaluate(() => { const m = document.querySelector('meta[http-equiv="Content-Security-Policy"]'); return m ? m.content : ""; });
    expect(policy).toContain("connect-src 'none'");
    expect(policy).toContain("default-src 'none'");
    if (name === "offline") expect(policy).not.toContain("http");
    // a full analysis: load, DE, plots (including WebGL scatter), enrichment, meta-analysis UI, figure export, session
    await page.evaluate(([m, d]) => {
      document.getElementById("normsel").value = "tmm"; document.getElementById("fmode").value = "fbe";
      loadFromText(m, d); GROUPTYPE = "auto"; analyzeNow();
      document.getElementById("selA").value = "Basal"; document.getElementById("selB").value = "LP"; deKey = "";
      for (const t of ["overview", "gene", "de", "coexp", "heat", "enrich", "compare", "patterns", "disco", "traits"]) switchTab(t);
      useBuiltin("hallmark_mm"); switchTab("enrich"); runEnrichment();
      sessionState();
    }, [read("demo/GSE63310_counts.tsv"), read("demo/GSE63310_design.tsv")]);
    await page.waitForTimeout(1500);
    const img = await page.evaluate(async () => { const d = document.querySelector(".js-plotly-plot"); return d ? (await Plotly.toImage(d, { format: "png", width: 600, height: 400 })).slice(0, 22) : "none"; });
    expect(img).toBe("data:image/png;base64,");
    // explicit attempts to reach another server must be blocked by the policy
    const attempts = await page.evaluate(async () => {
      const out = {};
      try { await fetch("https://example.org/x?d=1"); out.fetch = "sent"; } catch (e) { out.fetch = "blocked"; }
      await new Promise((res) => { const i = new Image(); i.onload = () => { out.img = "loaded"; res(); }; i.onerror = () => { out.img = "blocked"; res(); }; i.src = "https://example.org/p.png?d=1"; });
      await new Promise((res) => { try { const w = new WebSocket("wss://example.org/"); w.onerror = () => { out.ws = "blocked"; res(); }; w.onopen = () => { out.ws = "open"; res(); }; } catch (e) { out.ws = "blocked"; res(); } setTimeout(() => { out.ws = out.ws || "blocked"; res(); }, 3000); });
      return out;
    });
    expect(attempts.fetch).toBe("blocked");
    expect(attempts.img).toBe("blocked");
    expect(attempts.ws).toBe("blocked");
    const outside = requests.filter((u) => !b.allowed.some((re) => re.test(u)) && !u.startsWith("https://example.org"));
    expect(outside, "unexpected requests: " + outside.join(", ")).toEqual([]);
    const violations = await page.evaluate(() => window.__cspViolations);
    // the only violations must be the deliberate attempts above
    expect(violations.filter((v) => !/example\.org/.test(v)), "policy violations during normal use: " + violations.join("; ")).toEqual([]);
  });
}
