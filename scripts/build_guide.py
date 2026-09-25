#!/usr/bin/env python3
"""Compile docs/user-guide.md into dist/benchside-guide.html.

Standard library only. Supports the Markdown subset used by the guide: ATX headings
(with GitHub-style anchors), paragraphs, ordered and unordered lists, task lists,
tables, fenced code blocks, horizontal rules, inline code, bold, italic and links.
Links to knowledgebase pages (../wiki/Page.md) are rewritten to the repository wiki.
"""
import html
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "docs" / "user-guide.md"
OUT = ROOT / "dist" / "benchside-guide.html"
REPO = "https://github.com/jyaron/benchsiDE"


def slug(text):
    s = re.sub(r"<[^>]+>", "", text).strip().lower()
    s = re.sub(r"[^\w\- ]", "", s, flags=re.UNICODE)
    return s.replace(" ", "-")


def rewrite_link(url):
    m = re.match(r"^\.\./wiki/([A-Za-z0-9\-]+)\.md(#.*)?$", url)
    if m:
        return f"{REPO}/wiki/{m.group(1)}{m.group(2) or ''}"
    m = re.match(r"^\.\./([A-Za-z0-9_\-./]+)$", url)
    if m:
        return f"{REPO}/blob/main/{m.group(1)}"
    return url


def inline(s):
    codes = []

    def stash(m):
        codes.append(html.escape(m.group(1), quote=False))
        return f"\x00{len(codes) - 1}\x00"

    s = re.sub(r"`([^`]+)`", stash, s)
    s = html.escape(s, quote=False)
    s = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)",
               lambda m: f'<a href="{html.escape(rewrite_link(m.group(2)))}">{m.group(1)}</a>', s)
    s = re.sub(r"&lt;(https?://[^&\s]+)&gt;", r'<a href="\1">\1</a>', s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"(?<![\w*])\*([^*\s][^*]*)\*(?![\w*])", r"<i>\1</i>", s)
    s = re.sub(r"\x00(\d+)\x00", lambda m: f"<code>{codes[int(m.group(1))]}</code>", s)
    return s


def convert(md):
    out = []
    lines = md.split("\n")
    i = 0
    para = []

    def flush_para():
        if para:
            out.append("<p>" + inline(" ".join(para)) + "</p>")
            para.clear()

    while i < len(lines):
        line = lines[i]
        if line.startswith("```"):
            flush_para()
            j = i + 1
            block = []
            while j < len(lines) and not lines[j].startswith("```"):
                block.append(lines[j])
                j += 1
            out.append("<pre><code>" + html.escape("\n".join(block), quote=False) + "</code></pre>")
            i = j + 1
            continue
        m = re.match(r"^(#{1,6}) (.*)$", line)
        if m:
            flush_para()
            lvl, txt = len(m.group(1)), m.group(2).strip()
            out.append(f'<h{lvl} id="{slug(txt)}">{inline(txt)}</h{lvl}>')
            i += 1
            continue
        if re.match(r"^---+\s*$", line):
            flush_para()
            out.append("<hr>")
            i += 1
            continue
        if line.startswith("|"):
            flush_para()
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            head, body = rows[0], [r for r in rows[1:] if not all(re.match(r"^:?-+:?$", c) for c in r)]
            t = ["<table><thead><tr>" + "".join(f"<th>{inline(c)}</th>" for c in head) + "</tr></thead><tbody>"]
            for r in body:
                t.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>")
            out.append("".join(t) + "</tbody></table>")
            continue
        if re.match(r"^\s*(?:[-*]|\d+\.) ", line):
            flush_para()
            ordered = bool(re.match(r"^\s*\d+\. ", line))
            tag = "ol" if ordered else "ul"
            items = []
            while i < len(lines) and re.match(r"^\s*(?:[-*]|\d+\.) ", lines[i]):
                items.append(re.sub(r"^\s*(?:[-*]|\d+\.) ", "", lines[i]))
                i += 1
                while i < len(lines) and lines[i].startswith("   ") and lines[i].strip():
                    items[-1] += " " + lines[i].strip()
                    i += 1
            li = []
            for it in items:
                box = re.match(r"^\[( |x)\] (.*)$", it)
                if box:
                    li.append(f'<li class="task"><input type="checkbox" disabled{" checked" if box.group(1) == "x" else ""}> {inline(box.group(2))}</li>')
                else:
                    li.append(f"<li>{inline(it)}</li>")
            out.append(f"<{tag}>" + "".join(li) + f"</{tag}>")
            continue
        if not line.strip():
            flush_para()
            i += 1
            continue
        para.append(line.strip())
        i += 1
    flush_para()
    return "\n".join(out)


CSS = """
body{margin:0;font-family:-apple-system,'Segoe UI',system-ui,sans-serif;color:#1e293b;background:#f8fafc;line-height:1.55}
header{background:#0f172a;color:#fff;padding:14px 28px;display:flex;gap:18px;align-items:baseline;flex-wrap:wrap}
header h1{margin:0;font-size:19px;font-weight:400}header h1 b{color:#60a5fa;font-weight:800}
header a{color:#94a3b8;font-size:13px;text-decoration:none}
main{max-width:900px;margin:0 auto;padding:18px 28px 60px;background:#fff;border-left:1px solid #e2e8f0;border-right:1px solid #e2e8f0}
h1,h2,h3{color:#0f172a;line-height:1.25}h2{margin-top:2.2em;padding-bottom:4px;border-bottom:1px solid #e2e8f0}
h3{margin-top:1.6em}code{background:#f1f5f9;padding:1px 5px;border-radius:4px;font-size:.92em}
pre{background:#f1f5f9;padding:10px 14px;border-radius:6px;overflow-x:auto}pre code{background:none;padding:0}
table{border-collapse:collapse;margin:10px 0;font-size:14px;width:100%}th,td{border:1px solid #e2e8f0;padding:6px 10px;text-align:left;vertical-align:top;overflow-wrap:anywhere}
th{background:#f8fafc}hr{border:0;border-top:1px solid #e2e8f0;margin:2em 0}li{margin:3px 0}li.task{list-style:none;margin-left:-1.2em}
a{color:#2563eb}
"""


def main():
    body = convert(SRC.read_text(encoding="utf-8"))
    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>benchsiDE: User Guide</title>
<style>{CSS}</style>
</head>
<body>
<header><h1>benchsi<b>DE</b> User Guide</h1>
<a href="https://www.benchside.org">Application</a>
<a href="{REPO}/wiki">Knowledgebase</a>
<a href="{REPO}">Source</a></header>
<main>
{body}
</main>
</body>
</html>
"""
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT} ({len(page) / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
