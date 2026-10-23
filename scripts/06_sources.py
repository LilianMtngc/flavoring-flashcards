"""Étape 6 — remarques sourcées pour chaque molécule.

- PubChem (PUG-View) : numéro FEMA (avec le lien vers la fiche FEMA),
  descriptions d'odeur et de goût (HSDB…) avec leur source.
- Fiche FEMA Flavor Library : profil aromatique officiel (« Flavor Profile »)
  et numéro JECFA.

Cache : data/sources_cache.json (relancer ne réinterroge que les manquants).
"""
import html
import json
import re
import time
import urllib.parse
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
MOLS = ROOT / "data" / "molecules_verified.json"
CACHE = ROOT / "data" / "sources_cache.json"
VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading={h}"
UA = {"User-Agent": "flavoring-flashcards (https://github.com/LilianMtngc/flavoring-flashcards)"}


def get(url, timeout=30):
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=timeout, headers=UA)
            if r.status_code == 503:
                time.sleep(2 + 2 * attempt)
                continue
            return r
        except requests.RequestException:
            time.sleep(2)
    return None


def pubchem_heading(cid, heading):
    r = get(VIEW.format(cid=cid, h=urllib.parse.quote(heading)))
    time.sleep(0.25)
    if r is None or r.status_code != 200:
        return []
    rec = r.json()["Record"]
    refs = {x["ReferenceNumber"]: x for x in rec.get("Reference", [])}
    out = []

    def walk(sec):
        for info in sec.get("Information", []):
            v = info.get("Value", {})
            txt = " ".join(s.get("String", "") for s in v.get("StringWithMarkup", []))
            if not txt and "Number" in v:
                txt = str(v["Number"][0])
            ref = refs.get(info.get("ReferenceNumber"), {})
            if txt.strip():
                out.append({"texte": txt.strip(), "source": ref.get("SourceName", ""), "url": ref.get("URL", "")})
        for sub in sec.get("Section", []):
            walk(sub)

    walk(rec)
    return out


def fema_page(url):
    r = get(url)
    time.sleep(0.5)
    if r is None or r.status_code != 200:
        return {}
    t = re.sub(r"<script.*?</script>|<style.*?</style>", "", r.text, flags=re.S)
    lines = [l.strip() for l in html.unescape(re.sub(r"<[^>]+>", "\n", t)).splitlines() if l.strip()]
    res = {}
    for i, l in enumerate(lines):
        if l == "Flavor Profile" and i + 2 < len(lines):
            res["profil"] = lines[i + 2] if lines[i + 1] == "info" else lines[i + 1]
        if l == "JECFA Number" and i + 1 < len(lines) and lines[i + 1].isdigit():
            res["jecfa"] = lines[i + 1]
    return res


def main():
    mols = json.loads(MOLS.read_text(encoding="utf-8"))
    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    todo = [k for k, m in mols.items() if m.get("cid") and cache.get(k, {}).get("cid") != m["cid"]]
    print(f"{len(todo)} molécules à interroger")
    for i, k in enumerate(todo, 1):
        cid = mols[k]["cid"]
        entry = {"cid": cid,
                 "odeur": pubchem_heading(cid, "Odor"),
                 "gout": pubchem_heading(cid, "Taste"),
                 "fema": None}
        fema = pubchem_heading(cid, "FEMA Number")
        if fema:
            f = fema[0]
            entry["fema"] = {"numero": f["texte"], "url": f["url"], **(fema_page(f["url"]) if f["url"] else {})}
        cache[k] = entry
        if i % 20 == 0:
            print(f"  {i}/{len(todo)}", flush=True)
            CACHE.write_text(json.dumps(cache, indent=1, ensure_ascii=False), encoding="utf-8")
    CACHE.write_text(json.dumps(cache, indent=1, ensure_ascii=False), encoding="utf-8")
    print("terminé", flush=True)


if __name__ == "__main__":
    main()
