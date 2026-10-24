"""Étape 7 — source Wikipédia pour chaque couple produit naturel → molécule marqueur.

Un marqueur est retenu si la page Wikipédia (anglais) du produit cite la
molécule, ou si la page de la molécule cite le produit. Les couples non
confirmés sont écartés ; un produit sans aucun marqueur confirmé n'a pas de carte.

Entrées : data/naturels_ref.csv, data/naturels_wiki.csv, data/molecules_ref.csv
Sortie  : data/marqueurs_sources.json
"""
import csv
import json
import re
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "data"
API = "https://en.wikipedia.org/w/api.php"
UA = {"User-Agent": "AromaCards (https://github.com/LilianMtngc/aromacards)"}

# noms usuels en anglais qui ne se déduisent pas du nom de requête
EXTRA = {
    "cineole": ["eucalyptol", "cineole", "cineol"],
    "furfurylthiol": ["furfurylthiol", "furan-2-ylmethanethiol", "furfuryl mercaptan", "2-furanmethanethiol"],
    "vinylguaiacol_4": ["vinylguaiacol", "4-vinylguaiacol"],
    "whisky_lactone": ["whisky lactone", "whiskey lactone", "oak lactone", "methyl-octalactone"],
    "egcg": ["epigallocatechin gallate", "egcg", "catechin"],
    "cis_3_hexenol": ["cis-3-hexen-1-ol", "3-hexen-1-ol", "leaf alcohol", "hexenol"],
    "msg": ["glutamate", "msg"],
    "tetramethylpyrazine": ["tetramethylpyrazine", "ligustrazine"],
    "isovaleraldehyde": ["3-methylbutanal", "isovaleraldehyde"],
    "mcp": ["cyclotene", "methylcyclopentenolone", "maple lactone"],
    "sotolon": ["sotolon", "sotolone"],
    "glycyrrhizic_acid": ["glycyrrhizin", "glycyrrhizic acid"],
    "beta_damascenone": ["damascenone"],
    "furaneol": ["furaneol", "furanone"],
    "acetylpyrroline": ["2-acetyl-1-pyrroline", "acetylpyrroline"],
    "rebaudioside_a": ["rebaudioside"],
    "ligustilide": ["ligustilide"],
    "butylidenephthalide": ["butylidenephthalide", "phthalide"],
    "thiomenthone": ["mercapto-p-menthan-3-one", "8-mercaptomenthone", "thiomenthone"],
    "diosphenol": ["diosphenol", "buchu camphor"],
    "isobutyl_angelate": ["angelate"],
    "himachalene": ["himachalene"],
    "patchoulol": ["patchoulol", "patchouli alcohol"],
    "santalol": ["santalol"],
    "linalyl_acetate": ["linalyl acetate"],
    "methyl_anthranilate": ["methyl anthranilate"],
    "ethyl_dimethylpyrazine": ["2-ethyl-3,5-dimethylpyrazine", "pyrazine"],
    "delta_decalactone": ["decalactone", "δ-decalactone", "delta-decalactone"],
    "delta_octalactone": ["octalactone", "δ-octalactone", "delta-octalactone"],
    "terpinyl_acetate": ["terpinyl acetate"],
    "bornyl_acetate": ["bornyl acetate"],
    "neryl_acetate": ["neryl acetate"],
}
# couples écartés après relecture de la phrase trouvée (trop indirecte)
REJECTED = {
    ("mountain_maple", "sotolon"): "la source parle d'arôme d'érable artificiel, pas de l'érable",
    ("matcha", "egcg"): "la phrase dit seulement que l'ombrage réduit les catéchines",
}
PREFIX = re.compile(r"^(\(?[+\-−]?\)?-?|\([rszeE, 0-9]+\)-|alpha-|beta-|gamma-|delta-|trans-|cis-|[dl]-|[0-9,]+-)+", re.I)


def extract(title):
    for attempt in range(4):
        r = requests.get(API, params={"action": "query", "prop": "extracts", "explaintext": 1, "redirects": 1,
                                      "titles": title, "format": "json"}, headers=UA, timeout=30)
        time.sleep(0.3)
        try:
            data = r.json()
            break
        except ValueError:
            time.sleep(3 + 3 * attempt)
    else:
        return None, None
    page = next(iter(data["query"]["pages"].values()))
    if "missing" in page:
        return None, None
    return page["title"], page.get("extract", "")


# sur la page de la molécule, la phrase qui cite le produit doit dire que la molécule s'y trouve
# (« smells like lime » ne suffit pas, « found in lime oil » oui)
OCCURRENCE = re.compile(r"found|occur|present|component|constituent|contain|isolated|responsible|extract|derived|"
                        r"produced|abundant|dominant|major|main|principal|primary|characteristic|key|source|"
                        r"identified|detected|essential oil", re.I)


def sentence_with(text, terms, need_occurrence=False):
    """Phrase où l'un des termes apparaît comme mot entier (pluriel toléré : « rum » ne matche pas « spectrum »)."""
    for sent in re.split(r"(?<=[.!?])\s+|\n+", text):
        low = sent.lower()
        if any(re.search(r"(?<![a-z0-9])" + re.escape(t.lower()) + r"(?:s|es)?(?![a-z0-9])", low) for t in terms):
            if not need_occurrence or OCCURRENCE.search(sent):
                return sent.strip()
    return None


def synonyms(key, ref, title):
    names = {ref["query"].lower(), title.lower()} | set(EXTRA.get(key, []))
    names |= {PREFIX.sub("", n) for n in list(names)}
    return {n for n in names if len(n) >= 4}


def main():
    nats = list(csv.DictReader((D / "naturels_ref.csv").open(encoding="utf-8"), delimiter=";"))
    wiki = {r["key"]: r for r in csv.DictReader((D / "naturels_wiki.csv").open(encoding="utf-8"), delimiter=";")}
    mref = {r["key"]: r for r in csv.DictReader((D / "molecules_ref.csv").open(encoding="utf-8"), delimiter=";")}
    mols = json.loads((D / "molecules_verified.json").read_text(encoding="utf-8"))
    pub = json.loads((D / "resolve_cache.json").read_text(encoding="utf-8"))

    pages = {}

    def page(title):
        if title not in pages:
            pages[title] = extract(title)
        return pages[title]

    def marker_page(key):
        for t in [mref[key]["query"], (pub[key].get("pubchem_cas") or pub[key].get("pubchem_name") or {}).get("Title", ""),
                  *EXTRA.get(key, [])[:1]]:
            if t:
                title, text = page(t[0].upper() + t[1:])
                if text:
                    return title, text
        return None, None

    out = {}
    for n in nats:
        w = wiki[n["key"]]
        p_title, p_text = page(w["wiki"])
        words = [x.lower() for x in w["mots"].split("|")]
        res = []
        for mk in n["marqueurs"].split("|"):
            pc_title = (pub[mk].get("pubchem_cas") or pub[mk].get("pubchem_name") or {}).get("Title", "")
            syn = synonyms(mk, mref[mk], pc_title)
            source = phrase = None
            if p_text and (phrase := sentence_with(p_text, syn)):
                source = p_title
            else:
                m_title, m_text = marker_page(mk)
                if m_text and (phrase := sentence_with(m_text, words, need_occurrence=True)):
                    source = m_title
            if (n["key"], mk) in REJECTED:
                source, phrase = None, "écarté : " + REJECTED[(n["key"], mk)]
            res.append({"marqueur": mk, "source": source, "phrase": phrase,
                        "url": f"https://en.wikipedia.org/wiki/{source.replace(' ', '_')}" if source else None})
            print(f"{'OK ' if source else '-- '} {n['key']:22} {mk:24} {source or ''}", flush=True)
        out[n["key"]] = res
    (D / "marqueurs_sources.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    ok = sum(1 for v in out.values() if any(x["source"] for x in v))
    print(f"{ok}/{len(out)} produits avec au moins un marqueur sourcé")


if __name__ == "__main__":
    main()
