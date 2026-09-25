#!/usr/bin/env python3
"""Build dist/benchside-offline.html: index.html with Plotly.js inlined (no CDN dependency).

Reads the pinned Plotly version from the script tag in index.html so the pin lives in
one place. Verifies the downloaded file against a checksum and confirms the app script
is byte-identical between builds. Note: Plotly's source contains a default topojsonURL
string pointing at its CDN; it is used only by map traces, which benchsiDE does not
create. The offline build therefore performs no network requests.
"""
import hashlib, re, sys, urllib.request, pathlib
PLOTLY_SHA256 = "0a17719a72751704861215da0e5c5cdb3f9a8d50eff5cb84cb6f8b80786682b0"  # 2.32.0
root = pathlib.Path(__file__).resolve().parents[1]
src = (root / "index.html").read_text(encoding="utf-8")
m = re.search(r'<script src="(https://cdn\.plot\.ly/plotly-[\d.]+\.min\.js)"></script>', src)
if not m:
    sys.exit("CDN plotly tag not found in index.html")
js = urllib.request.urlopen(m.group(1)).read().decode("utf-8")
digest = hashlib.sha256(js.encode()).hexdigest()
if digest != PLOTLY_SHA256:
    sys.exit(f"Plotly checksum mismatch: {digest} (update PLOTLY_SHA256 when bumping the pin)")
out = src.replace(m.group(0), "<script>\n/* Plotly.js (MIT) inlined for offline use */\n" + js + "\n</script>")
assert '<script src="https://cdn.plot.ly' not in out
app_on  = re.findall(r"<script>\n(.*?)\n</script>", src, re.S)[-1]
app_off = re.findall(r"<script>\n(.*?)\n</script>", out, re.S)[-1]
assert app_on == app_off, "app script changed during inlining"
dist = root / "dist"; dist.mkdir(exist_ok=True)
(dist / "benchside-offline.html").write_text(out, encoding="utf-8")
print(f"wrote {dist/'benchside-offline.html'} ({len(out)/1e6:.1f} MB)")
