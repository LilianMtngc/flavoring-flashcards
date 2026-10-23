"""Étape 3 — vérifie chaque structure et choisit celle qui sera dessinée.

Contrôles :
  A. PubChem par nom vs PubChem par CAS (même InChIKey ?)
  B. sources indépendantes (NCI CACTUS par nom et par CAS, OPSIN) : on vote
     sur le 1er bloc de l'InChIKey (connectivité, hors stéréo)
  C. cohérence chimique : groupes fonctionnels attendus d'après le nom,
     hétéroatomes, stéréo demandée par le nom (E/Z, énantiomère)

Sorties : data/molecules_verified.json, verification.html
Corrections manuelles : data/overrides.csv (key;cid;raison) force un CID PubChem ;
data/revue_manuelle.csv (key;verdict;commentaire) valide après examen une
structure signalée (verdict « ok ») — le motif reste tracé dans le rapport.
"""
import csv
import html
import json
import re
from collections import Counter
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "data" / "molecules_ref.csv"
CACHE = ROOT / "data" / "resolve_cache.json"
OVERRIDES = ROOT / "data" / "overrides.csv"
REVIEW = ROOT / "data" / "revue_manuelle.csv"
OUT = ROOT / "data" / "molecules_verified.json"
REPORT = ROOT / "verification.html"

S = Chem.MolFromSmarts
RULES = [
    # (condition sur le nom, SMARTS attendu, libellé)
    (lambda n: "pyrazine" in n, S("c1cnccn1"), "cycle pyrazine"),
    (lambda n: "quinoxaline" in n, S("c1cnccn1"), "cycle pyrazine"),
    (lambda n: "thiazol" in n and "thiazolin" not in n, S("s1cncc1"), "cycle thiazole"),
    (lambda n: "thiazolin" in n, S("[#16;R]1~[#6;R]~[#7;R]~[#6;R]~[#6;R]1"), "cycle thiazoline"),
    (lambda n: "pyridine" in n and "hydropyridine" not in n, S("n1ccccc1"), "cycle pyridine"),
    (lambda n: "pyrrol" in n, S("[#7]1~[#6]~[#6]~[#6]~[#6]~1"), "cycle pyrrole"),
    (lambda n: re.search(r"furan(?!one|eol)|furfur|furoate", n) and "tetrahydro" not in n and "furanone" not in n,
     S("c1ccoc1"), "cycle furane"),
    (lambda n: "thiophen" in n and "hydro" not in n, S("c1ccsc1"), "cycle thiophène"),
    (lambda n: re.search(r"gamma-\w*lactone|butyrolactone", n), S("[#6;R]1(=O)-[#8;R]-[#6;R]~[#6;R]~[#6;R]-1"), "γ-lactone (cycle à 5)"),
    (lambda n: re.search(r"delta-\w*lactone", n), S("[#6;R]1(=O)-[#8;R]-[#6;R]~[#6;R]~[#6;R]~[#6;R]-1"), "δ-lactone (cycle à 6)"),
    (lambda n: re.search(r"\bacid\b", n), S("[CX3](=O)[OX2H1,OX1-]"), "acide carboxylique"),
    (lambda n: re.search(r"\bacetate\b", n), S("[CH3]C(=O)O[#6]"), "ester acétique"),
    (lambda n: re.search(r"(?<!iso)(?<!thio)butyrate", n), S("[CH3][CH2,CH1][#6]C(=O)O[#6]"), "ester butyrique"),
    (lambda n: "isobutyrate" in n, S("[CH3][CH1]([CH3])C(=O)O[#6]"), "ester isobutyrique"),
    (lambda n: "formate" in n, S("[CH1](=O)O[#6]"), "ester formique"),
    (lambda n: re.search(r"(?<!thio)propionate", n) and not re.search(r"\d-\(", n), S("[CH3][CH2]C(=O)O[#6]"), "ester propionique"),
    (lambda n: re.search(r"aldehyde|al$|al\b(?!c)", n) and "alcohol" not in n, S("[CX3H1](=O)[#6]"), "aldéhyde"),
    (lambda n: (re.search(r"one$", n) or "ketone" in n) and not re.search(r"furanone|pyranone|pyran-2-one|trithioacetone|lactone", n),
     S("[#6][CX3](=O)[#6]"), "cétone"),
    (lambda n: (re.search(r"ol$|diol|alcohol", n)) and not re.search(r"eucalyptol|thiol$", n), S("[OX2H]"), "fonction OH"),
    (lambda n: re.search(r"thiol|mercapt", n), S("[SX2H]"), "thiol"),
    (lambda n: "trisulfide" in n, S("[#16][#16][#16]"), "trisulfure"),
    (lambda n: "disulfide" in n, S("[#6][#16][#16][#6]"), "disulfure"),
    (lambda n: re.search(r"(?<!di)(?<!tri)sulfide", n), S("[#6][SX2][#6]"), "sulfure"),
    (lambda n: "methoxy" in n, S("[CH3]O[#6]"), "méthoxy"),
    (lambda n: re.search(r"phenyl|benz|anis|cinnam|vanill|guaia", n), S("c1ccccc1"), "cycle benzénique"),
]
N_HINT = re.compile(r"azin|pyri|thiazol|indol|anthranil|amin|quinoxalin|pyrrol|caffeine|capsaicin|glutamate|thialdine|nitr|cyanat|piperin")
S_HINT = re.compile(r"thi|sulf|mercapt|methional")
STEREO_EZ = re.compile(r"\btrans\b|\bcis\b|\(e\)|\(z\)|\(\d?e,|\(2e|\(e,|,\d?z\)|\be-|\bz-")
STEREO_CHIRAL = re.compile(r"\bl-|\bd-|\(\+\)|\(-\)|\(r\)|\(s\)|\(r\)-")

FAMILY_RULES = [
    ("pyrazine", S("c1cnccn1")),
    ("thiazole", S("s1cncc1")),
    ("hétérocycle azoté", S("[n,N;R]")),
    ("thiol", S("[SX2H]")),
    ("composé soufré", S("[#16]")),
]


def family(key, name, mol):
    for label, patt in FAMILY_RULES:
        if mol.HasSubstructMatch(patt):
            return label
    if "furanone" in name or key in {"furaneol", "homofuraneol", "sotolon", "maple_furanone"}:
        return "furanone"
    if key in {"maltol", "ethyl_maltol"}:  # énol d'une 4-pyranone, pas un phénol
        return "pyranone"
    if mol.HasSubstructMatch(S("[#6;R](=O)[#8;R]")):
        return "lactone"
    if mol.HasSubstructMatch(S("c1ccoc1")):
        return "furane"
    if mol.HasSubstructMatch(S("c[OX2H]")):
        return "phénol"
    if mol.HasSubstructMatch(S("[CX3](=O)[OX2H1,OX1-]")):
        return "acide"
    if mol.HasSubstructMatch(S("[CX3](=O)O[#6]")):
        return "ester"
    if mol.HasSubstructMatch(S("[CX3H1](=O)")):
        return "aldéhyde"
    if mol.HasSubstructMatch(S("[#6][CX3](=O)[#6]")):
        return "cétone"
    if mol.HasSubstructMatch(S("[OX2H][CX4,CX3]")):
        return "alcool"
    if all(a.GetSymbol() in "CH" for a in mol.GetAtoms()):
        return "hydrocarbure terpénique" if mol.GetNumAtoms() in (10, 15) else "hydrocarbure"
    if mol.HasSubstructMatch(S("[#6][OX2][#6]")):
        return "éther"
    return "autre"


def block(ik):
    return ik.split("-")[0] if ik else None


def chem_checks(query, mol, smiles):
    n = query.lower()
    problems, passed = [], []
    for cond, patt, label in RULES:
        if cond(n):
            (passed if mol.HasSubstructMatch(patt) else problems).append(label)
    elements = {a.GetSymbol() for a in mol.GetAtoms()}
    if N_HINT.search(n) and "N" not in elements:
        problems.append("azote attendu d'après le nom")
    if "N" in elements and not N_HINT.search(n) and "msg" not in n:
        problems.append("azote présent mais absent du nom")
    if S_HINT.search(n) and "S" not in elements:
        problems.append("soufre attendu d'après le nom")
    if "S" in elements and not S_HINT.search(n):
        problems.append("soufre présent mais absent du nom")
    if STEREO_EZ.search(n) and "/" not in smiles and "\\" not in smiles:
        problems.append("stéréo E/Z demandée mais absente")
    if STEREO_CHIRAL.search(n) and "@" not in smiles:
        problems.append("énantiomère demandé mais stéréo absente")
    return passed, problems


def main():
    refs = list(csv.DictReader(REF.open(encoding="utf-8"), delimiter=";"))
    cache = json.loads(CACHE.read_text(encoding="utf-8"))
    overrides = {}
    if OVERRIDES.exists():
        overrides = {r["key"]: r for r in csv.DictReader(OVERRIDES.open(encoding="utf-8"), delimiter=";")}
    review = {}
    if REVIEW.exists():
        review = {r["key"]: r for r in csv.DictReader(REVIEW.open(encoding="utf-8"), delimiter=";")}
    out, rows = {}, []
    for e in refs:
        key, c = e["key"], cache.get(e["key"], {})
        pn, pc = c.get("pubchem_name"), c.get("pubchem_cas")
        votes = {
            "PubChem (nom)": block(pn["InChIKey"]) if pn else None,
            "PubChem (CAS)": block(pc["InChIKey"]) if pc else None,
            "CACTUS (nom)": block(c.get("cactus_name")),
            "CACTUS (CAS)": block(c.get("cactus_cas")),
            "OPSIN": block(c.get("opsin")),
        }
        counts = Counter(v for v in votes.values() if v)
        notes = []
        if key in overrides:
            rec = c.get("override")
            winner = block(rec["InChIKey"]) if rec else None
            notes.append(f"CID imposé manuellement : {overrides[key]['raison']}")
        else:
            winner = counts.most_common(1)[0][0] if counts else None
            candidates = [r for r in (pc, pn) if r and block(r["InChIKey"]) == winner]
            rec = candidates[0] if candidates else None
        status, problems = "ok", []
        if rec is None:
            status = "fail"
            problems.append("aucune structure PubChem cohérente")
            out[key] = {"key": key, "status": status, "problems": problems}
            rows.append((e, votes, None, status, problems, notes, winner))
            continue
        smiles = rec["SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        agree = [s for s, v in votes.items() if v == winner]
        dissent = [s for s, v in votes.items() if v and v != winner]
        independent = [s for s in agree if not s.startswith("PubChem")]
        if dissent:
            problems.append("désaccord : " + ", ".join(dissent))
        if key not in overrides and (len(agree) < 3 or not independent):
            problems.append(f"seulement {len(agree)} source(s) concordante(s)")
        if pn and pc and pn["InChIKey"] != pc["InChIKey"] and block(pn["InChIKey"]) == block(pc["InChIKey"]):
            notes.append("stéréo différente entre fiche nom et fiche CAS : fiche CAS retenue")
        passed, chem_problems = chem_checks(e["query"], mol, smiles)
        problems += chem_problems
        if problems:
            status = "warn"
            if review.get(key, {}).get("verdict") == "ok":
                status = "ok"
                notes.append(f"revue manuelle : {review[key]['commentaire']} (signalé : {'; '.join(problems)})")
                problems = []
        fam = family(key, e["query"].lower(), mol)
        out[key] = {
            "key": key,
            "nom_fr": e["nom_fr"],
            "cas": e["cas"],
            "cid": rec["CID"],
            "smiles": smiles,
            "inchikey": rec["InChIKey"],
            "formule": rdMolDescriptors.CalcMolFormula(mol),
            "masse": round(Descriptors.MolWt(mol), 2),
            "iupac": rec.get("IUPACName", ""),
            "famille": fam,
            "sources_ok": len(agree),
            "controles_chimiques": passed,
            "status": status,
            "problems": problems,
            "notes": notes + ([e["note"]] if e["note"] and e["note"] not in ("marqueur", "theorique") else []),
            "marqueur_seul": e["note"] == "marqueur",
        }
        rows.append((e, votes, rec, status, problems, notes, winner))
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    write_report(rows)
    st = Counter(v["status"] for v in out.values())
    print(f"ok={st['ok']} warn={st['warn']} fail={st['fail']}")
    for e, votes, rec, status, problems, notes, winner in rows:
        if status != "ok":
            print(f"[{status}] {e['key']}: {'; '.join(problems)}")


def write_report(rows):
    icon = {"ok": "✅", "warn": "⚠️", "fail": "❌"}
    order = {"fail": 0, "warn": 1, "ok": 2}
    rows = sorted(rows, key=lambda r: (order[r[3]], r[0]["key"]))
    st = Counter(r[3] for r in rows)
    trs = []
    for e, votes, rec, status, problems, notes, winner in rows:
        vote_cells = "".join(
            f"<td class='{'y' if v and v == winner else ('n' if v else 'e')}'>{'✓' if v and v == winner else ('✗' if v else '–')}</td>"
            for v in votes.values())
        trs.append(
            f"<tr><td>{icon[status]}</td><td><b>{html.escape(e['nom_fr'])}</b><br><small>{html.escape(e['query'])}</small></td>"
            f"<td>{e['cas']}</td><td>{rec['CID'] if rec else ''}</td>"
            f"<td><code>{html.escape(rec['SMILES']) if rec else ''}</code></td>{vote_cells}"
            f"<td>{html.escape('; '.join(problems + notes))}</td></tr>")
    REPORT.parent.mkdir(exist_ok=True)
    REPORT.write_text(f"""<!doctype html><html lang="fr"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Vérification des structures</title>
<style>
:root{{--bg:#fff;--fg:#1d1d1b;--mut:#6b6a66;--line:#e4e2dc;--y:#e6f4ea;--n:#fdecea}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1c1c1a;--fg:#ecebe6;--mut:#a3a29c;--line:#3a3936;--y:#1f3a27;--n:#4a2020}}}}
body{{font-family:system-ui,sans-serif;background:var(--bg);color:var(--fg);margin:16px;font-size:14px}}
table{{border-collapse:collapse;width:100%}}td,th{{border-bottom:1px solid var(--line);padding:6px;vertical-align:top;text-align:left}}
small{{color:var(--mut)}}code{{font-size:12px;word-break:break-all}}.y{{background:var(--y);text-align:center}}.n{{background:var(--n);text-align:center}}.e{{text-align:center;color:var(--mut)}}
</style></head><body>
<h1>Vérification des structures</h1>
<p>{len(rows)} molécules · ✅ {st['ok']} vérifiées · ⚠️ {st['warn']} à revoir · ❌ {st['fail']} non résolues.</p>
<p><small>✓ = la source donne la même connectivité (1er bloc InChIKey) que la structure retenue · ✗ = désaccord · – = source sans réponse.</small></p>
<table><tr><th></th><th>Molécule</th><th>CAS</th><th>CID</th><th>SMILES retenu</th><th>PubChem nom</th><th>PubChem CAS</th><th>CACTUS nom</th><th>CACTUS CAS</th><th>OPSIN</th><th>Remarques</th></tr>
{''.join(trs)}</table></body></html>""", encoding="utf-8")


if __name__ == "__main__":
    main()
