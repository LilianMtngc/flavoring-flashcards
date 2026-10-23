"""Étape 8 — classifications de danger « CMR » de chaque molécule (PubChem PUG-View).

- Classification harmonisée européenne (règlement CLP (CE) n° 1272/2008, annexe VI) :
  mentions de danger H340 / H341 (mutagène), H350 / H351 (cancérogène),
  H360 / H361 (toxique pour la reproduction). Seule la source « Regulation (EC)
  No 1272/2008 » est retenue : les autoclassifications d'autres pays ne le sont pas.
- Centre international de recherche sur le cancer (CIRC / IARC) : groupes 1, 2A, 2B.

Il s'agit de dangers de la substance pure, pas d'une évaluation du risque alimentaire.
Cache : data/securite_cache.json (relancer ne réinterroge que les manquants).
"""
import json
import re
import time
import urllib.parse
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
MOLS = ROOT / "data" / "molecules_verified.json"
CACHE = ROOT / "data" / "securite_cache.json"
VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading={h}"
UA = {"User-Agent": "flavoring-flashcards (https://github.com/LilianMtngc/flavoring-flashcards)"}
CLP = "Regulation (EC) No 1272/2008"
CMR_CODE = re.compile(r"\bH3(?:40|41|50|51|60|61)[A-Za-z]*\b")


def infos(cid, heading):
    for attempt in range(3):
        try:
            r = requests.get(VIEW.format(cid=cid, h=urllib.parse.quote(heading)), timeout=30, headers=UA)
            time.sleep(0.25)
            if r.status_code == 503:
                time.sleep(2 + 2 * attempt)
                continue
            break
        except requests.RequestException:
            time.sleep(2)
    else:
        return None
    if r.status_code == 404:
        return []
    if r.status_code != 200:
        return None
    rec = r.json()["Record"]
    refs = {x["ReferenceNumber"]: x for x in rec.get("Reference", [])}
    out = []

    def walk(sec):
        for i in sec.get("Information", []):
            txt = " | ".join(s.get("String", "") for s in i.get("Value", {}).get("StringWithMarkup", []))
            ref = refs.get(i.get("ReferenceNumber"), {})
            out.append({"nom": i.get("Name", ""), "source": ref.get("SourceName", ""), "url": ref.get("URL", ""), "texte": txt})
        for s in sec.get("Section", []):
            walk(s)

    walk(rec)
    return out


def classify(cid):
    ghs, carc = infos(cid, "GHS Classification"), infos(cid, "Carcinogen Classification")
    if ghs is None or carc is None:
        return None
    clp = sorted({c for i in ghs if i["source"].startswith(CLP) and i["nom"] == "GHS Hazard Statements"
                  for c in CMR_CODE.findall(i["texte"])})
    iarc = sorted({m.group(1) for i in carc if i["nom"] == "IARC Carcinogenic Classes"
                   for m in [re.search(r"Group (1|2A|2B)\b", i["texte"])] if m})
    iarc_url = next((i["url"] for i in carc if i["nom"] == "IARC Carcinogenic Classes" and i["url"]), "")
    return {"cid": cid, "clp": clp, "iarc": iarc, "iarc_url": iarc_url}


def main():
    mols = json.loads(MOLS.read_text(encoding="utf-8"))
    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    todo = [k for k, m in mols.items() if m.get("cid") and cache.get(k, {}).get("cid") != m["cid"]]
    print(f"{len(todo)} molécules à interroger", flush=True)
    for n, k in enumerate(todo, 1):
        res = classify(mols[k]["cid"])
        if res is not None:
            cache[k] = res
        if n % 25 == 0:
            print(f"  {n}/{len(todo)}", flush=True)
            CACHE.write_text(json.dumps(cache, indent=1, ensure_ascii=False), encoding="utf-8")
    CACHE.write_text(json.dumps(cache, indent=1, ensure_ascii=False), encoding="utf-8")
    flagged = {k: v for k, v in cache.items() if v["clp"] or v["iarc"]}
    print(f"terminé : {len(flagged)} molécules avec une classification CMR (CLP) ou CIRC")
    for k, v in sorted(flagged.items()):
        print(f"  {k}: CLP {v['clp']} CIRC {v['iarc']}")


if __name__ == "__main__":
    main()
