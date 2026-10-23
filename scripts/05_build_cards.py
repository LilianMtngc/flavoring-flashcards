"""Étape 5 — assemble les cartes et génère les sorties.

Entrées : data/descripteurs.csv, data/molecules_ref.csv, data/molecules_verified.json,
          data/naturels_ref.csv, data/marqueurs_sources.json, data/sources_cache.json,
          data/usages.csv, structures/*.svg
Sorties : data/cartes.json, index.html (flash cards), flavoring-flashcards.apkg (Anki)
"""
import csv
import html
import json
import re
import warnings
from collections import defaultdict
from pathlib import Path

import genanki
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

# genanki ne connaît pas la balise <svg> : Anki l'affiche pourtant très bien
warnings.filterwarnings("ignore", message="Field contained the following invalid HTML tags")

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "data"
STRUCT = ROOT / "structures"

BROWN = {
    "CARAMEL", "CARAMEL AMBRÉ", "CARAMEL BRÛLÉ", "CACAO", "POUDRÉ CACAO", "CHOCOLAT", "CHOCOLAT NOIR",
    "CAFÉ", "MOUTURE", "GRILLÉ", "NOISETTE", "NUTTY", "PRALINÉ", "CACAHUÈTE", "ÉRABLE", "POUDRÉ VANILLÉ",
    "VANILLE", "FUMÉ", "BACON", "CENDRE", "MALT", "BISCUITÉ", "PAIN", "CÉRÉALE", "POP-CORN", "RIZ BASMATI",
    "BRÛLÉ", "SUCRE BRUN", "SUCRE CUIT", "MÉLASSE", "CASSONADE", "TABAC", "NOIX", "FENUGREC", "SÉSAME",
}
NIVEAUX = ["Facile", "Moyen", "Difficile"]


def read_csv(name):
    rows = list(csv.DictReader((D / name).open(encoding="utf-8"), delimiter=";"))
    bad = [r for r in rows if None in r]
    if bad:  # un « ; » non protégé par des guillemets décale les colonnes
        raise SystemExit(f"{name} : {len(bad)} ligne(s) avec un champ en trop, ex. {list(bad[0].values())[0]}")
    return rows


def is_brown(desc_lists):
    top2 = {d for l in desc_lists for d in l[:2]}
    every = {d for l in desc_lists for d in l}
    return bool(top2 & BROWN) or len(every & BROWN) >= 2


def svg_for(key):
    svg = (STRUCT / f"{key}.svg").read_text(encoding="utf-8")
    m = re.search(r"viewBox='0 0 ([\d.]+) ([\d.]+)'", svg)
    w, h = float(m.group(1)), float(m.group(2))
    # même échelle pour toutes les molécules ; les grosses sont réduites par max-width / max-height
    return svg.replace("<svg ", f"<svg class='mol' role='img' style='width:{w * 3:.0f}px;aspect-ratio:{w:.1f}/{h:.1f}' ", 1)


def bertz(smiles):
    return GraphDescriptors.BertzCT(Chem.MolFromSmiles(smiles))


def sourced_notes(src):
    """Remarques sourcées : fiche FEMA (n°, profil, JECFA) et descriptions PubChem (odeur, goût)."""
    if not src:
        return {}
    out = {}
    if src.get("fema"):
        out["fema"] = src["fema"]
    for kind in ("odeur", "gout"):
        seen, items = set(), []
        for it in src.get(kind, []):
            t = it["texte"].rstrip(".")
            if t.isupper():  # certaines fiches HSDB sont en capitales
                t = t.capitalize()
            if t.lower() in seen:
                continue
            seen.add(t.lower())
            items.append({**it, "texte": t if len(t) <= 180 else t[:177] + "…"})
        if items:
            out[kind] = items[:2]
    return out


def build():
    descr_rows = read_csv("descripteurs.csv")
    mref = {r["key"]: r for r in read_csv("molecules_ref.csv")}
    mols = json.loads((D / "molecules_verified.json").read_text(encoding="utf-8"))
    nats = {r["key"]: r for r in read_csv("naturels_ref.csv")}
    marker_src = json.loads((D / "marqueurs_sources.json").read_text(encoding="utf-8"))
    sources = json.loads((D / "sources_cache.json").read_text(encoding="utf-8"))
    usages = {r["key"]: r for r in read_csv("usages.csv")}

    groups = {}
    for r in descr_rows:
        groups.setdefault((r["type"], r["key"]), []).append(
            (r["forme"], [d for d in r["descripteurs"].split("|") if d]))

    # marqueurs confirmés par une source, et produits qui les partagent
    confirmed = {k: [x for x in v if x["source"]] for k, v in marker_src.items()}
    marker_users = defaultdict(list)
    for (typ, key) in groups:
        if typ == "nat":
            for x in confirmed.get(key, []):
                marker_users[x["marqueur"]].append(nats[key]["nom_fr"])

    cards, skipped = [], []
    for (typ, key), forms in groups.items():
        desc_lists = [l for _, l in forms]
        descripteurs = list(dict.fromkeys(d for l in desc_lists for d in l))
        par_forme = None
        if len({tuple(sorted(l)) for l in desc_lists if l}) > 1:
            par_forme = [{"forme": f, "descripteurs": l} for f, l in forms if l]
        u = usages[key]
        card = {
            "id": f"{typ}:{key}", "type": typ,
            "formes": [f for f, _ in forms if typ == "mol"],
            "descripteurs": descripteurs, "descripteurs_par_forme": par_forme,
            "brune": is_brown(desc_lists), "usages": u["usages"], "origine": u["origine"],
        }
        if typ == "mol":
            m = mols[key]
            note = mref[key]["note"]
            card.update({
                "titre": m["nom_fr"], "iupac": m["iupac"], "famille": m["famille"],
                "formule": m["formule"], "masse": m["masse"], "cas": m["cas"], "cid": m["cid"],
                "verif": {"status": m["status"], "sources": m["sources_ok"]},
                "structures": [{"key": key, "nom": m["nom_fr"], "svg": svg_for(key)}],
                "notes": [note] if note and note != "marqueur" else [],
                "sources": sourced_notes(sources.get(key)),
                "bertz": bertz(m["smiles"]),
            })
        else:
            markers = confirmed.get(key, [])
            if not markers:
                skipped.append(nats[key]["nom_fr"])
                continue
            n = nats[key]
            card.update({
                "titre": n["nom_fr"], "famille": "produit naturel",
                "structures": [{"key": x["marqueur"], "nom": mols[x["marqueur"]]["nom_fr"], "svg": svg_for(x["marqueur"]),
                                "formule": mols[x["marqueur"]]["formule"], "cas": mols[x["marqueur"]]["cas"],
                                "cid": mols[x["marqueur"]]["cid"], "source": x["source"], "url": x["url"]}
                               for x in markers],
                "verif": {"status": "ok", "sources": min(mols[x["marqueur"]]["sources_ok"] for x in markers)},
                "partage": sorted(set(p for x in markers for p in marker_users[x["marqueur"]]) - {n["nom_fr"]}),
                "notes": [], "sources": {},
                "bertz": max(bertz(mols[x["marqueur"]]["smiles"]) for x in markers),
            })
        cards.append(card)

    # niveau de difficulté = tercile de l'indice de complexité de Bertz (calculé sur les molécules)
    vals = sorted(c["bertz"] for c in cards if c["type"] == "mol")
    t1, t2 = vals[len(vals) // 3], vals[2 * len(vals) // 3]
    for c in cards:
        c["niveau"] = NIVEAUX[0] if c["bertz"] < t1 else NIVEAUX[1] if c["bertz"] < t2 else NIVEAUX[2]
        c["bertz"] = round(c["bertz"])
    (D / "cartes.json").write_text(json.dumps(cards, indent=1, ensure_ascii=False), encoding="utf-8")
    return cards, skipped, (round(t1), round(t2))


# ---------------------------------------------------------------- HTML

def build_html(cards, seuils):
    template = (ROOT / "scripts" / "flashcards_template.html").read_text(encoding="utf-8")
    payload = json.dumps({"cards": cards, "seuils": seuils}, ensure_ascii=False).replace("</", "<\\/")
    (ROOT / "index.html").write_text(template.replace("/*__DATA__*/null", payload), encoding="utf-8")


# ---------------------------------------------------------------- Anki

ANKI_CSS = """
.card{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;font-size:16px;color:#1d1d1b;background:#fbfaf7;text-align:left;padding:8px}
.nightMode.card,.night_mode .card{color:#ecebe6;background:#1c1c1a}
.kind{font-size:12px;letter-spacing:.04em;text-transform:uppercase;color:#8a8780;text-align:center}
.structs{display:flex;flex-wrap:wrap;gap:16px;justify-content:center;align-items:center;margin:12px 0}
.structs figure{margin:0;text-align:center}
svg.mol{max-width:100%;max-height:260px;height:auto}
.nom{font-size:22px;font-weight:600;margin:4px 0}
.sub{color:#8a8780;font-size:13px}
.chip{display:inline-block;border:1px solid #d8d5cc;border-radius:6px;padding:1px 7px;margin:2px;font-size:13px}
.chip.b{background:#faeeda;color:#633806;border-color:#efc98a}
.badge{display:inline-block;background:#faeeda;color:#633806;border-radius:6px;padding:2px 8px;font-size:12px}
.sec{margin-top:10px}.sec b{font-size:13px;color:#8a8780;font-weight:600}
a{color:inherit}
"""

FRONT = """<div class="kind">{{Type}}</div><div class="structs">{{Structure}}</div>"""
BACK = """{{FrontSide}}<hr id="answer">
<div class="nom">{{Nom}}</div><div class="sub">{{Sous_titre}}</div>
<div class="sec">{{Infos}}</div>
<div class="sec"><b>Descripteurs</b><br>{{Descripteurs}}</div>
<div class="sec"><b>Usages</b><br>{{Usages}}</div>
{{#Sources}}<div class="sec"><b>Sources</b><br>{{Sources}}</div>{{/Sources}}"""


def esc(s):
    return html.escape(str(s))


def link(text, url):
    return f"<a href='{esc(url)}'>{esc(text)}</a>" if url else esc(text)


def sources_html(c):
    parts = []
    if c["type"] == "nat":
        for s in c["structures"]:
            parts.append(f"{esc(s['nom'])} : {link('Wikipédia — ' + s['source'], s['url'])}")
    src = c.get("sources", {})
    if src.get("fema"):
        f = src["fema"]
        txt = f"FEMA n° {f['numero']}" + (f" · profil : {f['profil']}" if f.get("profil") else "") \
            + (f" · JECFA n° {f['jecfa']}" if f.get("jecfa") else "")
        parts.append(link(txt, f.get("url")))
    for kind, label in (("odeur", "Odeur"), ("gout", "Goût")):
        for it in src.get(kind, []):
            parts.append(f"{label} : « {esc(it['texte'])} » — {link(it['source'] or 'PubChem', it['url'])}")
    if c.get("cid"):
        parts.append(link(f"PubChem CID {c['cid']}", f"https://pubchem.ncbi.nlm.nih.gov/compound/{c['cid']}"))
    return "<br>".join(parts)


def anki_fields(c):
    structs = "".join(f"<figure>{s['svg']}</figure>" for s in c["structures"])
    if c["type"] == "mol":
        kind, sub = "Molécule", c.get("iupac", "")
        infos = (f"<span class='chip'>{esc(c['famille'])}</span><span class='chip'>{esc(c['formule'])}</span>"
                 f"<span class='chip'>{c['masse']} g/mol</span><span class='chip'>CAS {esc(c['cas'])}</span>")
    else:
        kind = "Produit naturel · molécule(s) marqueur(s)"
        sub = "Marqueur(s) : " + ", ".join(s["nom"] for s in c["structures"])
        infos = ""
    infos = (f"<span class='badge'>Note brune</span> " if c["brune"] else "") \
        + f"<span class='chip'>Complexité : {esc(c['niveau'])}</span>" + infos
    desc = "".join(f"<span class='chip{' b' if d in BROWN else ''}'>{esc(d.lower())}</span>" for d in c["descripteurs"]) \
        or "<i>descripteurs non renseignés</i>"
    usage = esc(c["usages"]) + (f"<br><i>Origine : {esc(c['origine'])}</i>" if c["origine"] else "")
    return [kind, structs, esc(c["titre"]), esc(sub), infos, desc, usage, sources_html(c)]


def build_anki(cards):
    model = genanki.Model(
        1607392320, "Flavoring flashcards — structure → nom",
        fields=[{"name": n} for n in ["Type", "Structure", "Nom", "Sous_titre", "Infos", "Descripteurs", "Usages", "Sources"]],
        templates=[{"name": "Structure → nom", "qfmt": FRONT, "afmt": BACK}],
        css=ANKI_CSS,
    )
    deck = genanki.Deck(2059400111, "Flavoring flashcards::Structures (notes brunes)")
    for c in cards:
        tags = [f"complexite_{c['niveau'].lower()}",
                "famille_" + re.sub(r"\W+", "_", c["famille"].lower()).strip("_"),
                "naturel" if c["type"] == "nat" else "molecule"]
        if c["brune"]:
            tags.append("note_brune")
        deck.add_note(genanki.Note(model=model, fields=anki_fields(c), tags=tags, guid=genanki.guid_for(c["id"])))
    genanki.Package(deck).write_to_file(str(ROOT / "flavoring-flashcards.apkg"))


def main():
    cards, skipped, seuils = build()
    build_html(cards, seuils)
    build_anki(cards)
    n_mol = sum(c["type"] == "mol" for c in cards)
    print(f"{len(cards)} cartes ({n_mol} molécules, {len(cards) - n_mol} produits naturels), "
          f"{sum(c['brune'] for c in cards)} notes brunes, seuils de Bertz {seuils}")
    print(f"{len(skipped)} produits naturels sans marqueur sourcé : {', '.join(skipped)}")


if __name__ == "__main__":
    main()
