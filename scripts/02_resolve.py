"""Étape 2 — résout chaque molécule de data/molecules_ref.csv auprès de 3 sources.

- PubChem par nom et PubChem par CAS (CID, SMILES, InChIKey, formule, masse, nom IUPAC)
- NCI CACTUS par nom et par CAS (InChIKey standard)
- OPSIN par nom (InChIKey standard ; ne marche que pour les noms systématiques)

Les réponses sont mises en cache dans data/resolve_cache.json : relancer le
script ne réinterroge que les entrées manquantes ou modifiées.
"""
import csv
import json
import time
import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "data" / "molecules_ref.csv"
CACHE = ROOT / "data" / "resolve_cache.json"

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
PROPS = "SMILES,ConnectivitySMILES,InChIKey,MolecularFormula,MolecularWeight,IUPACName,Title"
session = requests.Session()


def get(url, timeout=30, retries=3):
    for attempt in range(retries):
        try:
            r = session.get(url, timeout=timeout)
            if r.status_code == 503:
                time.sleep(2 + attempt * 2)
                continue
            return r
        except requests.RequestException:
            time.sleep(2)
    return None


def pubchem(identifier):
    q = urllib.parse.quote(identifier, safe="")
    r = get(f"{PUG}/name/{q}/cids/JSON")
    time.sleep(0.25)
    if r is None or r.status_code != 200:
        return None
    cids = r.json().get("IdentifierList", {}).get("CID", [])
    if not cids:
        return None
    r = get(f"{PUG}/cid/{cids[0]}/property/{PROPS}/JSON")
    time.sleep(0.25)
    if r is None or r.status_code != 200:
        return None
    p = r.json()["PropertyTable"]["Properties"][0]
    p["n_cids"] = len(cids)
    return p


def cactus(identifier):
    q = urllib.parse.quote(identifier, safe="")
    r = get(f"https://cactus.nci.nih.gov/chemical/structure/{q}/stdinchikey", timeout=40, retries=2)
    if r is None or r.status_code != 200:
        return None
    lines = [l for l in r.text.splitlines() if l.startswith("InChIKey=")]
    return lines[0].split("=", 1)[1].strip() if lines else None


def opsin(name):
    q = urllib.parse.quote(name, safe="")
    r = get(f"https://opsin.ch.cam.ac.uk/opsin/{q}.json", timeout=30, retries=2)
    if r is None or r.status_code != 200:
        return None
    data = r.json()
    return data.get("stdinchikey") if data.get("status") == "SUCCESS" else None


def resolve_independent(entry):
    return {
        "cactus_name": cactus(entry["query"]),
        "cactus_cas": cactus(entry["cas"]) if entry["cas"] else None,
        "opsin": opsin(entry["query"]),
    }


def main():
    refs = list(csv.DictReader(REF.open(encoding="utf-8"), delimiter=";"))
    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    todo = [e for e in refs
            if cache.get(e["key"], {}).get("query") != e["query"]
            or cache.get(e["key"], {}).get("cas") != e["cas"]]
    print(f"{len(refs)} molécules, {len(todo)} à résoudre")

    with ThreadPoolExecutor(max_workers=6) as pool:
        independent = dict(zip([e["key"] for e in todo], pool.map(resolve_independent, todo)))

    for i, e in enumerate(todo, 1):
        res = {
            "query": e["query"],
            "cas": e["cas"],
            "pubchem_name": pubchem(e["query"]),
            "pubchem_cas": pubchem(e["cas"]) if e["cas"] else None,
            **independent[e["key"]],
        }
        cache[e["key"]] = res
        if i % 25 == 0:
            print(f"  {i}/{len(todo)}")
            CACHE.write_text(json.dumps(cache, indent=1, ensure_ascii=False), encoding="utf-8")
    overrides = ROOT / "data" / "overrides.csv"
    if overrides.exists():
        for o in csv.DictReader(overrides.open(encoding="utf-8"), delimiter=";"):
            current = cache.get(o["key"], {}).get("override") or {}
            if str(current.get("CID")) != o["cid"]:
                r = get(f"{PUG}/cid/{o['cid']}/property/{PROPS}/JSON")
                cache[o["key"]]["override"] = r.json()["PropertyTable"]["Properties"][0]
                time.sleep(0.25)
    CACHE.write_text(json.dumps(cache, indent=1, ensure_ascii=False), encoding="utf-8")
    print("terminé")


if __name__ == "__main__":
    main()
