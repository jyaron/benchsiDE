"""Download the public GEO count matrices used by the engine tests and verify their checksums.
Usage: python3 tests/engine/fetch_data.py   (writes tests/engine/data/<GSE>.tsv)"""
import gzip, hashlib, io, json, os, sys, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "data"); os.makedirs(OUT, exist_ok=True)
DS = json.load(open(os.path.join(HERE, "datasets.json")))
for acc, d in DS.items():
    dest = os.path.join(OUT, acc + ".tsv")
    if os.path.exists(dest): continue
    raw = gzip.decompress(urllib.request.urlopen(d["source"], timeout=600).read()).decode()
    if "subset" in d:   # GSE143688: wild-type and Il36r fl/fl columns only
        lines = raw.rstrip("\n").split("\n"); hdr = lines[0].split("\t")
        keep = [0] + [i for i, c in enumerate(hdr) if i and c.split("_")[1] in ("WT", "FL")]
        rows = []
        for l in lines:
            f = l.split("\t"); rows.append("\t".join(f[i].strip() for i in keep))
        raw = "\n".join(rows) + "\n"
        want = d["sha256_of_subset_tsv"]
    else:
        want = d["sha256_of_decompressed"]
    got = hashlib.sha256(raw.encode()).hexdigest()
    status = "ok" if got == want else "CHECKSUM DIFFERS (GEO may have regenerated the file)"
    open(dest, "w").write(raw); print(acc, status)
