"""Étape 5 — assemble les cartes et génère les sorties.

Entrées : data/descripteurs.csv, data/molecules_ref.csv, data/molecules_verified.json,
          data/naturels_ref.csv, data/marqueurs_sources.json, data/sources_cache.json,
          data/usages.csv, data/glossaire_fema.csv, data/theoriques_wiki.csv, structures/*.svg
Sorties : data/cartes.json, index.html (flash cards), flavoring-flashcards.apkg (Anki)

Trois catégories de cartes :
- « mol » : molécule évaluée (descripteurs et dosage du panel) ;
- « nat » : produit naturel, représenté par sa ou ses molécules marqueurs ;
- « theo » : molécule de référence ajoutée pour l'apprentissage (note « theorique »
  dans molecules_ref.csv) ; ses descripteurs viennent de la page Wikipédia de la
  molécule quand elle figure dans data/theoriques_wiki.csv (citation à l'appui),
  sinon du profil FEMA traduit (data/glossaire_fema.csv).
Les cartes dont le recto serait identique (même structure, stéréo comprise) sont
fusionnées : le nom de la molécule est affiché, les produits naturels sont cités.
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
    "BRÛLÉ", "CROÛTE DE PAIN", "SUCRE BRUN", "SUCRE CUIT", "MÉLASSE", "CASSONADE", "TABAC", "NOIX", "FENUGREC", "SÉSAME",
}
NIVEAUX = ["Facile", "Moyen", "Difficile"]
TYPES = {"mol": "Molécule", "theo": "Molécule · théorique", "nat": "Produit naturel"}


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


def fema_descriptors(src, glossaire):
    """Descripteurs d'une molécule théorique : profil FEMA traduit terme à terme (data/glossaire_fema.csv)."""
    profil = ((src or {}).get("fema") or {}).get("profil")
    if not profil:
        return []
    terms = [t.strip() for t in profil.split(",") if t.strip()]
    missing = [t for t in terms if t.lower() not in glossaire]
    if missing:
        raise SystemExit(f"glossaire_fema.csv : termes FEMA sans traduction : {', '.join(missing)}")
    return list(dict.fromkeys(glossaire[t.lower()] for t in terms))


def build():
    descr_rows = read_csv("descripteurs.csv")
    mref = {r["key"]: r for r in read_csv("molecules_ref.csv")}
    mols = json.loads((D / "molecules_verified.json").read_text(encoding="utf-8"))
    nats = {r["key"]: r for r in read_csv("naturels_ref.csv")}
    marker_src = json.loads((D / "marqueurs_sources.json").read_text(encoding="utf-8"))
    sources = json.loads((D / "sources_cache.json").read_text(encoding="utf-8"))
    usages = {r["key"]: r for r in read_csv("usages.csv")}
    glossaire = {r["en"].lower(): r["fr"] for r in read_csv("glossaire_fema.csv")}
    wiki = {r["key"]: r for r in read_csv("theoriques_wiki.csv")}

    groups = {}
    for r in descr_rows:
        groups.setdefault((r["type"], r["key"]), []).append({
            "forme": r["forme"],
            "descripteurs": [d for d in r["descripteurs"].split("|") if d],
            "ppm": [p for p in r["ppm"].split("|") if p]})
    confirmed = {k: [x for x in v if x["source"]] for k, v in marker_src.items()}

    # 1. entrées brutes : molécules évaluées, molécules théoriques, produits naturels avec marqueur sourcé
    entries, skipped = [], []
    for (typ, key), forms in groups.items():
        if typ == "mol":
            entries.append({"type": "mol", "key": key, "nom": mols[key]["nom_fr"], "formes": forms, "markers": [key]})
        elif confirmed.get(key):
            entries.append({"type": "nat", "key": key, "nom": nats[key]["nom_fr"], "formes": forms,
                            "markers": [x["marqueur"] for x in confirmed[key]], "liens": confirmed[key]})
        else:
            skipped.append(nats[key]["nom_fr"])
    for key, r in mref.items():
        if r["note"] == "theorique" and ("mol", key) not in groups:
            e = {"type": "theo", "key": key, "nom": mols[key]["nom_fr"], "markers": [key]}
            fema = fema_descriptors(sources.get(key), glossaire)
            if key in wiki:  # choisi explicitement (pas de profil FEMA, ou profil FEMA commun aux deux énantiomères)
                w = wiki[key]
                e["formes"] = [{"forme": "Wikipédia", "ppm": [], "descripteurs": w["descripteurs"].split("|")}]
                e["wiki"] = {"page": w["wiki"], "citation": w["citation"],
                             "url": "https://en.wikipedia.org/wiki/" + w["wiki"].replace(" ", "_")}
            elif fema:
                e["formes"] = [{"forme": "profil FEMA", "ppm": [], "descripteurs": fema}]
            else:
                raise SystemExit(f"molécule théorique sans descripteur sourcé : {key}")
            entries.append(e)

    # 2. un seul recto par structure : regroupement sur l'InChIKey complète (stéréo comprise)
    by_struct = defaultdict(list)
    for e in entries:
        by_struct[tuple(sorted(mols[k]["inchikey"] for k in e["markers"]))].append(e)

    cards = []
    for members in by_struct.values():
        members.sort(key=lambda e: ["mol", "theo", "nat"].index(e["type"]))
        main, prods = members[0], [e for e in members if e["type"] == "nat"]
        if main["type"] == "nat" and len(members) == 1:
            cards.append(nat_card(main, mols, usages[main["key"]]))
            continue
        assert len(main["markers"]) == 1, f"structures multiples partagées : {[e['nom'] for e in members]}"
        key = main["markers"][0]
        mol_forms = [f for e in members if e["type"] != "nat" for f in e["formes"]]
        prod_forms = [{**f, "forme": e["nom"]} for e in prods for f in e["formes"]]
        u = usages[key if main["type"] == "nat" else main["key"]]
        all_lists = [f["descripteurs"] for f in mol_forms + prod_forms]
        main_lists = [f["descripteurs"] for f in main["formes"]] if main["type"] != "nat" else all_lists
        m = mols[key]
        card = {
            # identifiant de la carte d'origine conservé : la progression enregistrée reste valable
            "id": f"{main['type']}:{main['key']}",
            "type": "theo" if all(e["type"] == "theo" for e in members) else "mol",
            "titre": m["nom_fr"], "iupac": m["iupac"], "famille": m["famille"],
            "formule": m["formule"], "masse": m["masse"], "cas": m["cas"], "cid": m["cid"],
            "verif": {"status": m["status"], "sources": m["sources_ok"]},
            "structures": [{"key": key, "nom": m["nom_fr"], "svg": svg_for(key)}],
            "descripteurs": list(dict.fromkeys(d for l in main_lists for d in l)),
            "formes": mol_forms + prod_forms,
            "marqueur_de": [{"nom": e["nom"], "source": x["source"], "url": x["url"]}
                            for e in prods for x in e["liens"]],
            "brune": is_brown(all_lists), "usages": u["usages"], "origine": u["origine"],
            "notes": [n for n in [mref[key]["note"]] if n and n not in ("marqueur", "theorique")],
            "sources": sourced_notes(sources.get(key)),
            "bertz": bertz(m["smiles"]),
        }
        for e in members:
            if e.get("wiki"):
                card["sources"]["wiki"] = e["wiki"]
        cards.append(card)

    # produits naturels qui partagent un marqueur (sans avoir exactement les mêmes)
    users = defaultdict(set)
    for c in cards:
        if c["type"] == "nat":
            for s in c["structures"]:
                users[s["key"]].add(c["titre"])
        for p in c["marqueur_de"]:
            users[c["structures"][0]["key"]].add(p["nom"])
    for c in cards:
        if c["type"] == "nat":
            c["partage"] = sorted(set().union(*(users[s["key"]] for s in c["structures"])) - {c["titre"]})

    # niveau de difficulté = tercile de l'indice de complexité de Bertz (calculé sur les molécules)
    vals = sorted(c["bertz"] for c in cards if c["type"] != "nat")
    t1, t2 = vals[len(vals) // 3], vals[2 * len(vals) // 3]
    for c in cards:
        c["niveau"] = NIVEAUX[0] if c["bertz"] < t1 else NIVEAUX[1] if c["bertz"] < t2 else NIVEAUX[2]
        c["bertz"] = round(c["bertz"])
    cards.sort(key=lambda c: (["mol", "nat", "theo"].index(c["type"]), c["id"]))
    (D / "cartes.json").write_text(json.dumps(cards, indent=1, ensure_ascii=False), encoding="utf-8")
    return cards, skipped, (round(t1), round(t2))


def nat_card(e, mols, u):
    structs = [mols[k] for k in e["markers"]]
    return {
        "id": f"nat:{e['key']}", "type": "nat", "titre": e["nom"], "famille": "produit naturel",
        "structures": [{"key": x["marqueur"], "nom": m["nom_fr"], "svg": svg_for(x["marqueur"]),
                        "formule": m["formule"], "cas": m["cas"], "cid": m["cid"],
                        "source": x["source"], "url": x["url"]} for x, m in zip(e["liens"], structs)],
        "verif": {"status": "ok", "sources": min(m["sources_ok"] for m in structs)},
        "descripteurs": list(dict.fromkeys(d for f in e["formes"] for d in f["descripteurs"])),
        "formes": e["formes"], "marqueur_de": [],
        "brune": is_brown([f["descripteurs"] for f in e["formes"]]),
        "usages": u["usages"], "origine": u["origine"], "notes": [], "sources": {},
        "bertz": max(bertz(m["smiles"]) for m in structs),
    }


def fmt_ppm(values):
    return " / ".join(v.replace(".", ",") for v in values) + " ppm" if values else ""


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
<div class="sec"><b>Usages (synthèse IA)</b><br>{{Usages}}</div>
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
    for p in c.get("marqueur_de", []):
        parts.append(f"Marqueur de {esc(p['nom'])} : {link('Wikipédia — ' + p['source'], p['url'])}")
    src = c.get("sources", {})
    if src.get("fema"):
        f = src["fema"]
        txt = f"FEMA n° {f['numero']}" + (f" · profil : {f['profil']}" if f.get("profil") else "") \
            + (f" · JECFA n° {f['jecfa']}" if f.get("jecfa") else "")
        parts.append(link(txt, f.get("url")))
    if src.get("wiki"):
        w = src["wiki"]
        parts.append(f"Descripteurs : « {esc(w['citation'])} » — {link('Wikipédia — ' + w['page'], w['url'])}")
    for kind, label in (("odeur", "Odeur"), ("gout", "Goût")):
        for it in src.get(kind, []):
            parts.append(f"{label} : « {esc(it['texte'])} » — {link(it['source'] or 'PubChem', it['url'])}")
    if c.get("cid"):
        parts.append(link(f"PubChem CID {c['cid']}", f"https://pubchem.ncbi.nlm.nih.gov/compound/{c['cid']}"))
    return "<br>".join(parts)


def chips(descr):
    return "".join(f"<span class='chip{' b' if d in BROWN else ''}'>{esc(d.lower())}</span>" for d in descr)


def anki_fields(c):
    structs = "".join(f"<figure>{s['svg']}</figure>" for s in c["structures"])
    if c["type"] != "nat":
        sub = c.get("iupac", "")
        if c["marqueur_de"]:
            sub += "<br>Marqueur de : " + esc(", ".join(dict.fromkeys(p["nom"] for p in c["marqueur_de"])))
        infos = (f"<span class='chip'>{esc(c['famille'])}</span><span class='chip'>{esc(c['formule'])}</span>"
                 f"<span class='chip'>{c['masse']} g/mol</span><span class='chip'>CAS {esc(c['cas'])}</span>")
    else:
        sub = "Marqueur(s) : " + esc(", ".join(s["nom"] for s in c["structures"]))
        infos = ""
    infos = (f"<span class='badge'>Note brune</span> " if c["brune"] else "") \
        + f"<span class='chip'>Complexité : {esc(c['niveau'])}</span>" + infos
    desc = chips(c["descripteurs"]) or "<i>descripteurs non renseignés</i>"
    if len(c["formes"]) > 1 or any(f["ppm"] for f in c["formes"]):
        desc += "".join(f"<br><small><b>{esc(f['forme'])}</b>{' · ' + fmt_ppm(f['ppm']) if f['ppm'] else ''} : "
                        f"{esc(', '.join(f['descripteurs']).lower())}</small>" for f in c["formes"])
    usage = esc(c["usages"]) + (f"<br><i>Origine : {esc(c['origine'])}</i>" if c["origine"] else "")
    kind = TYPES[c["type"]] + (" · molécule(s) marqueur(s)" if c["type"] == "nat" else "")
    return [kind, structs, esc(c["titre"]), sub, infos, desc, usage, sources_html(c)]


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
                {"nat": "naturel", "mol": "molecule", "theo": "theorique"}[c["type"]]]
        if c["brune"]:
            tags.append("note_brune")
        deck.add_note(genanki.Note(model=model, fields=anki_fields(c), tags=tags, guid=genanki.guid_for(c["id"])))
    genanki.Package(deck).write_to_file(str(ROOT / "flavoring-flashcards.apkg"))


def main():
    cards, skipped, seuils = build()
    build_html(cards, seuils)
    build_anki(cards)
    n = {t: sum(c["type"] == t for c in cards) for t in TYPES}
    print(f"{len(cards)} cartes ({n['mol']} molécules, {n['nat']} produits naturels, {n['theo']} théoriques), "
          f"{sum(c['brune'] for c in cards)} notes brunes, {sum(bool(c['marqueur_de']) for c in cards)} fusions, "
          f"seuils de Bertz {seuils}")
    print(f"{len(skipped)} produits naturels sans marqueur sourcé : {', '.join(skipped)}")


if __name__ == "__main__":
    main()
