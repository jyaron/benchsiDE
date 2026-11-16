#!/usr/bin/env python3
"""Rebuild the pig annotation (ANNO_BUILTIN.pig), the one-to-one ortholog tables for pairs other than mouse-human
(ORTHO_PAIRS) and the rat and pig members of the homology families (ORTHO_FAM) embedded in index.html.

Sources (downloaded when run):
  MGI HOM_AllOrganism.rpt                       human, mouse, rat homology classes
  NCBI Sus_scrofa.gene_info                     pig Entrez gene IDs and symbols
  Ensembl BioMart (current release)             pig genes, Entrez cross-references, Compara orthologs of pig with
                                                human, mouse and rat; the 12 pig breed assemblies and their
                                                one-to-one Sscrofa11.1 counterparts
Rules:
  one-to-one   a gene with exactly one partner in the ortholog graph of the two species (and that partner with
               exactly one partner)
  pig symbol   NCBI Gene symbol of the linked Entrez gene (non-LOC symbol, unique), else the Ensembl gene name,
               else the LOC symbol (a gene linked to several non-LOC symbols and without an Ensembl name is left unnamed); breed genes take the symbol of their one-to-one Sscrofa11.1 counterpart,
               else their Ensembl gene name
  families     existing human-mouse families keep their members; a rat or pig gene joins the family its orthologs
               belong to (genes linked to several families are left out); remaining genes form new families from
               connected components of the ortholog graph, kept when not one to one for some pair involving rat or
               pig and when no species has more members than the largest existing family
  libraries    rat and pig gene-set libraries: each human MSigDB set (from BUILTIN_ENC *_hs) with every member
               replaced by all its rat orthologs (MGI) or pig orthologs (Ensembl Compara); human symbols renamed since
               MSigDB v2024.1 are updated through NCBI Gene synonyms (unique match only); sets left empty are removed
Usage: python3 scripts/build_orthology.py   (writes build/orthology_tables.json; the human, mouse and rat
       annotation and ORTHO_MH / ORTHO_FAM human and mouse members are read from index.html)
"""
import gzip, io, json, os, re, sys
from collections import defaultdict
import requests
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
MART = "https://www.ensembl.org/biomart/martservice"
BREEDS = ["ssusmarc", "ssjinhua", "ssmeishan", "ssrongchang", "sslandrace", "sshampshire", "sstibetan", "ssberkshire",
          "ssbamei", "sspietrain", "sslargewhite", "sswuzhishan"]


def mart(ds, attrs):
    xml = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" '
           'header="1" uniqueRows="1" datasetConfigVersion="0.6"><Dataset name="%s" interface="default">' % ds
           + "".join('<Attribute name="%s"/>' % a for a in attrs) + "</Dataset></Query>")
    r = requests.get(MART, params={"query": xml}, timeout=900); r.raise_for_status()
    assert not r.text.startswith("Query ERROR"), r.text[:300]
    return pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)


def app_const(src, name):
    k = src.find("const %s=" % name); e = src.find("\n", k)
    return json.loads(src[k + len(name) + 7:e].rstrip(";"))


def decode(tbl, prefix, pad):
    syms = tbl["sym"].split("\t"); acc = 0; out = {}
    for i, p in enumerate(tbl["ids"].split(".")):
        acc += int(p, 36); out[prefix + str(acc).zfill(pad)] = syms[i]
    return out


def strict(E):
    da = defaultdict(int); db = defaultdict(int)
    for x, y in E: da[x] += 1; db[y] += 1
    return {x: y for x, y in E if da[x] == 1 and db[y] == 1}


def main():
    src = open(os.path.join(ROOT, "index.html"), encoding="utf-8").read()
    AB = app_const(src, "ANNO_BUILTIN"); FAM0 = [{"h": f["h"], "m": f["m"]} for f in app_const(src, "ORTHO_FAM")
                                                  if "hm" in f.get("x", "hm").split(",")]
    BM = {sp: {"ens": decode(AB[sp], AB[sp]["prefix"], AB[sp]["pad"]), "ent": decode(AB[sp]["ent"], "", 0)}
          for sp in ("human", "mouse", "rat")}
    # MGI homology classes, symbols as in the built-in tables (through the Entrez ID)
    H = pd.read_csv("https://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt", sep="\t", dtype=str)
    sp_of = {"human": "h", "mouse, laboratory": "m", "rat": "r"}; ent = {"h": "human", "m": "mouse", "r": "rat"}
    H = H[H["Common Organism Name"].isin(sp_of)].copy(); H["k"] = H["Common Organism Name"].map(sp_of)
    H["sym"] = [BM[ent[k]]["ent"].get(e, s) if isinstance(e, str) else s for k, e, s in zip(H.k, H["EntrezGene ID"], H.Symbol)]
    EM = {}
    for a, b in (("r", "h"), ("r", "m"), ("m", "h")):
        E = set()
        for _, g in H.groupby("DB Class Key"):
            for x in g[g.k == a].sym.unique():
                for y in g[g.k == b].sym.unique(): E.add((x, y))
        EM[a + "|" + b] = E
    # pig annotation
    GI = pd.read_csv("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Sus_scrofa.gene_info.gz", sep="\t", dtype=str)
    ent_sym = dict(zip(GI.GeneID, GI.Symbol))
    EP = mart("sscrofa_gene_ensembl", ["ensembl_gene_id", "external_gene_name", "entrezgene_id", "gene_biotype"])
    EP.columns = ["ens", "name", "entrez", "biotype"]
    links = defaultdict(set)
    for e, g in zip(EP.ens, EP.entrez):
        if isinstance(g, str): links[e].add(g)
    for g, x in zip(GI.GeneID, GI.dbXrefs):
        for e in re.findall(r"Ensembl:(ENSSSCG\d+)", x): links[e].add(g)
    ensname = EP.dropna(subset=["name"]).drop_duplicates("ens").set_index("ens").name.to_dict()
    pig_ens = {}
    for e in EP.ens.unique():
        gs = [g for g in links.get(e, ()) if g in ent_sym]
        nonloc = sorted({ent_sym[g] for g in gs if not ent_sym[g].startswith("LOC")})
        if len(nonloc) == 1: pig_ens[e] = nonloc[0]
        elif e in ensname: pig_ens[e] = ensname[e]
        elif gs and not nonloc: pig_ens[e] = sorted(ent_sym[g] for g in gs)[0]   # several NCBI symbols: left unnamed
    alt = {}
    for b in BREEDS:
        nm = mart(b + "_gene_ensembl", ["ensembl_gene_id", "external_gene_name"]); nm.columns = ["ens", "name"]
        rf = mart(b + "_gene_ensembl", ["ensembl_gene_id", "sscrofa_homolog_ensembl_gene", "sscrofa_homolog_orthology_type"])
        rf.columns = ["ens", "ref", "type"]
        o2o = rf[rf.type == "ortholog_one2one"].drop_duplicates("ens").set_index("ens").ref.to_dict()
        nmd = nm.dropna(subset=["name"]).set_index("ens").name.to_dict()
        for e in nm.ens:
            if o2o.get(e) in pig_ens: alt[e] = pig_ens[o2o[e]]
            elif e in nmd: alt[e] = nmd[e]
    # pig orthologs
    EE = {}
    for k, ds, mp in (("p|h", "hsapiens", BM["human"]["ens"]), ("p|m", "mmusculus", BM["mouse"]["ens"]), ("p|r", "rnorvegicus", BM["rat"]["ens"])):
        d = mart("sscrofa_gene_ensembl", ["ensembl_gene_id", ds + "_homolog_ensembl_gene", ds + "_homolog_associated_gene_name"])
        d.columns = ["ens", "o_ens", "o_name"]; E = set()
        for e, oe, on in zip(d.ens, d.o_ens, d.o_name):
            if not isinstance(oe, str) or e not in pig_ens: continue
            y = mp.get(oe) or (on if isinstance(on, str) else None)
            if y: E.add((pig_ens[e], y))
        EE[k] = E
    PAIRS = {"rat|human": strict(EM["r|h"]), "rat|mouse": strict(EM["r|m"]), "pig|human": strict(EE["p|h"]),
             "pig|mouse": strict(EE["p|m"]), "pig|rat": strict(EE["p|r"])}
    # families
    adj = defaultdict(set)
    for key, E in list(EM.items()) + list(EE.items()):
        a, b = key.split("|")
        for x, y in E: adj[(a, x)].add((b, y)); adj[(b, y)].add((a, x))
    owner = {(k, x): i for i, f in enumerate(FAM0) for k in ("h", "m") for x in f[k]}
    FAM = [{"h": list(f["h"]), "m": list(f["m"]), "r": [], "p": []} for f in FAM0]; conflict = set()
    for node in list(adj):
        if node[0] in ("r", "p"):
            fams = {owner[n] for n in adj[node] if n in owner}
            if len(fams) == 1: FAM[fams.pop()][node[0]].append(node[1])
            elif len(fams) > 1: conflict.add(node)
    MAXSIDE = max(max(len(f["h"]), len(f["m"])) for f in FAM0)
    def pairs_of(f, no_hm):
        sps = [k for k in "hmrp" if f[k]]
        P = [a + b for i, a in enumerate(sps) for b in sps[i + 1:] if len(f[a]) + len(f[b]) > 2 and max(len(f[a]), len(f[b])) <= MAXSIDE]
        return [p for p in P if not (no_hm and p == "hm")]
    OUT = [dict(f, x=",".join(pairs_of(f, False))) for f in FAM]
    assigned = set(owner) | {(k, x) for f in FAM for k in ("r", "p") for x in f[k]} | conflict; seen = set()
    for node in adj:
        if node in assigned or node in seen: continue
        comp, stack = [], [node]; seen.add(node)
        while stack:
            n = stack.pop(); comp.append(n)
            for m in adj[n]:
                if m not in assigned and m not in seen: seen.add(m); stack.append(m)
        by = defaultdict(list)
        for k, x in comp: by[k].append(x)
        if len(by) >= 2 and max(len(v) for v in by.values()) >= 2:
            f = {k: sorted(by.get(k, [])) for k in "hmrp"}; p = pairs_of(f, True)
            if p: OUT.append(dict(f, x=",".join(p)))
    # rat and pig gene-set libraries
    BE = app_const(src, "BUILTIN_ENC")
    hsyms = set(BM["human"]["ens"].values()) | set(BM["human"]["ent"].values())
    HG = pd.read_csv("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz", sep="\t", dtype=str)
    syn = defaultdict(set)
    for sy, ss in zip(HG.Symbol, HG.Synonyms):
        if isinstance(ss, str) and ss != "-":
            for x in ss.split("|"): syn[x].add(sy)
    def dec(lib):
        V = lib["vocab"].split("\t"); out = []
        for st in lib["sets"]:
            b = st.index("|"); acc = 0; g = []
            for d in st[b + 1:].split("."): acc += int(d, 36); g.append(V[acc])
            out.append((st[:b], g))
        return out
    HSL = {k[:-3]: dec(BE[k]) for k in BE if k.endswith("_hs")}
    allg = {g for v in HSL.values() for _, gs in v for g in gs}; cur = set(HG.Symbol); RENAME = {}
    for g in allg - hsyms:
        if g in cur: continue
        c = [x for x in syn.get(g, ()) if x in hsyms]
        if len(c) == 1 and c[0] not in allg: RENAME[g] = c[0]
    h2 = {"rn": defaultdict(set), "ss": defaultdict(set)}
    for x, y in EM["r|h"]: h2["rn"][y].add(x)
    for x, y in EE["p|h"]: h2["ss"][y].add(x)
    LIBS = {}
    for base, sets in HSL.items():
        for sx in ("rn", "ss"):
            out = []
            for name, gs in sets:
                mem = set()
                for g in gs: mem |= h2[sx].get(RENAME.get(g, g), set())
                if mem: out.append([name, sorted(mem)])
            LIBS[base + "_" + sx] = out
    os.makedirs(os.path.join(ROOT, "build"), exist_ok=True)
    json.dump({"pig_ens": pig_ens, "pig_ent": ent_sym, "pig_alt": alt, "ORTHO_PAIRS": PAIRS, "ORTHO_FAM": OUT, "LIBS": LIBS, "RENAME": RENAME},
              open(os.path.join(ROOT, "build", "orthology_tables.json"), "w"))
    print({k: len(v) for k, v in PAIRS.items()}, "families:", len(OUT), "pig IDs:", len(pig_ens), len(ent_sym), len(alt))


if __name__ == "__main__":
    main()
