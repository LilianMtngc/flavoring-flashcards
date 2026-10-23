"""Détection des « sosies » : molécules à structures très proches (utilisé par 05_build_cards.py).

Deux molécules sont sosies si :
- elles ont la même formule brute, ou ne diffèrent que d'une double liaison C=C
  (deux hydrogènes) en gardant les mêmes fonctions chimiques ;
- leur plus grande sous-structure commune (MCS, ordres de liaison ignorés, cycles
  appariés aux cycles) couvre tous leurs atomes lourds sauf un au plus, ou sauf
  deux pour des isomères à cycles non aromatiques en même nombre (sotolon / furanéol) ;
- elles sont toutes deux des acides carboxyliques, ou aucune des deux.
Les stéréo-isomères (même connectivité) forment toujours un groupe à part.
Cela regroupe les isomères de position (eugénol / isoeugénol), les doubles liaisons
déplacées ou en plus (β-ionone / β-damascone / β-damascénone), les isomères E/Z et
les énantiomères (carvones).

La différence est décrite automatiquement et surlignée sur les dessins, alignés sur
le squelette commun : atomes hors du squelette commun, liaisons dont l'ordre diffère
d'une molécule à l'autre, ou centres stéréogènes / doubles liaisons E-Z quand la
connectivité est identique.
"""
import re
from collections import Counter
from itertools import combinations

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D

MIN_ATOMES = 6
TAILLE_MAX = 4  # au-delà, le groupe est découpé en sous-groupes plus serrés
HIGHLIGHT = (0.93, 0.55, 0.13, 0.55)
ATOM_COLORS = {8: (0.82, 0.19, 0.18), 7: (0.18, 0.44, 0.86), 16: (0.75, 0.54, 0.0), 11: (0.55, 0.35, 0.8)}


def heavy_formula(mol):
    return tuple(sorted(Counter(a.GetSymbol() for a in mol.GetAtoms()).items()))


def mcs(mols, exact=False, timeout=5):
    p = rdFMCS.MCSParameters()
    p.AtomTyper = rdFMCS.AtomCompare.CompareElements
    p.BondTyper = rdFMCS.BondCompare.CompareOrderExact if exact else rdFMCS.BondCompare.CompareAny
    p.BondCompareParameters.RingMatchesRingOnly = True
    p.AtomCompareParameters.RingMatchesRingOnly = True
    p.Timeout = timeout
    return rdFMCS.FindMCS(list(mols), p)


ACIDE = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
FONCTIONS = [Chem.MolFromSmarts(x) for x in (
    "[CX3]=O", "[OX2H]", "[#6][OX2][#6]", "[CX3](=O)[OX2][#6]", "[CX3H1]=O", "[SX2H]", "[#16][#16]", "a")]


def n_h(mol):
    return sum(a.GetTotalNumHs() for a in mol.GetAtoms())


def fonctions(mol):
    return tuple(len(mol.GetSubstructMatches(p)) for p in FONCTIONS)


def similarity(a, b):
    """Score dans [0, 1] (part des atomes dans la MCS), 0 si les molécules ne sont pas comparables."""
    if heavy_formula(a) != heavy_formula(b) or a.GetNumAtoms() < MIN_ATOMES:
        return 0.0
    dh = abs(n_h(a) - n_h(b))
    if dh > 2 or (dh == 2 and fonctions(a) != fonctions(b)):
        return 0.0
    if a.HasSubstructMatch(ACIDE) != b.HasSubstructMatch(ACIDE):
        return 0.0
    rings = a.GetRingInfo().NumRings()
    aromatic = any(x.GetIsAromatic() for x in a.GetAtoms()) or any(x.GetIsAromatic() for x in b.GetAtoms())
    ecart = 2 if dh == 0 and rings and rings == b.GetRingInfo().NumRings() and not aromatic else 1
    k = mcs([a, b]).numAtoms
    return k / a.GetNumAtoms() if a.GetNumAtoms() - k <= ecart else 0.0


def stereo_labels(mol):
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    atoms = {a.GetIdx(): a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")}
    bonds = {b.GetIdx(): str(b.GetStereo()) for b in mol.GetBonds() if b.GetStereo() != Chem.BondStereo.STEREONONE}
    return atoms, bonds


def describe(mols):
    """Nature de la différence entre les membres d'un groupe (texte calculé)."""
    keys = {Chem.MolToInchiKey(m).split("-")[0] for m in mols}
    if len(keys) == 1:
        if any(stereo_labels(m)[1] for m in mols):
            return "Même enchaînement d'atomes, géométrie différente d'une double liaison : isomères E / Z"
        return "Même enchaînement d'atomes, configuration différente : énantiomères (images dans un miroir)"
    same_skeleton = mcs(mols).numAtoms == mols[0].GetNumAtoms()
    if len({n_h(m) for m in mols}) > 1:
        return ("Même squelette, une double liaison C=C en plus ou en moins" if same_skeleton
                else "Atomes enchaînés autrement, et une double liaison C=C en plus ou en moins")
    if same_skeleton:
        return "Même squelette et même formule, double liaison ou fonction déplacée"
    return "Même formule, atomes enchaînés autrement : isomères de position"


def groups(items):
    """items : [(id, mol)] -> liste de groupes (listes d'id), et sosies de chaque id avec leur score."""
    n = len(items)
    edges = {}
    by_formula = {}
    for i, (_, m) in enumerate(items):
        by_formula.setdefault(heavy_formula(m), []).append(i)
    for idx in by_formula.values():
        for i, j in combinations(idx, 2):
            s = similarity(items[i][1], items[j][1])
            if s > 0:
                edges[(i, j)] = s
    neighbours = {i: [] for i in range(n)}
    for (i, j), s in edges.items():
        neighbours[i].append((s, j))
        neighbours[j].append((s, i))
    # regroupement par liens décroissants ; un groupe ne dépasse pas TAILLE_MAX et tous ses membres sont sosies deux à deux
    parent = list(range(n))
    size = [1] * n
    conn = [Chem.MolToInchiKey(m).split("-")[0] for _, m in items]
    fixed = set()  # stéréo-isomères : groupe fermé

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (i, j), s in edges.items():
        if conn[i] == conn[j]:
            a, b = find(i), find(j)
            if a != b:
                parent[b] = a
                size[a] += size[b]
            fixed.add(find(i))
    for (i, j), s in sorted(edges.items(), key=lambda e: (-e[1], e[0])):
        a, b = find(i), find(j)
        if a != b and a not in fixed and b not in fixed and size[a] + size[b] <= TAILLE_MAX:
            # toutes les paires du groupe fusionné doivent être sosies (groupe homogène)
            ma = [k for k in range(n) if find(k) == a]
            mb = [k for k in range(n) if find(k) == b]
            if all((min(x, y), max(x, y)) in edges for x in ma for y in mb):
                parent[b] = a
                size[a] += size[b]
    comp = {}
    for i in range(n):
        if neighbours[i]:
            comp.setdefault(find(i), []).append(i)
    out = [[items[i][0] for i in sorted(c)] for c in comp.values() if len(c) > 1]
    sos = {items[i][0]: [(items[j][0], round(s, 2)) for s, j in sorted(neighbours[i], reverse=True)] for i in range(n)}
    return out, sos


def highlight(mols):
    """Ce qui distingue les membres d'un groupe, pour chaque membre : (atomes, liaisons) à surligner.

    On superpose les molécules sur leur squelette commun (MCS, ordres de liaison ignorés) :
    - atomes hors du squelette commun (substituant ou hétéroatome placé ailleurs) ;
    - liaisons du squelette dont l'ordre n'est pas le même partout (double liaison déplacée ou en plus) ;
    - à défaut (même connectivité), centres stéréogènes et doubles liaisons E / Z.
    """
    patt = Chem.MolFromSmarts(mcs(mols).smartsString)
    matches = [m.GetSubstructMatch(patt) for m in mols]
    orders = {}  # liaison du motif -> ordres rencontrés
    for m, match in zip(mols, matches):
        for b in patt.GetBonds():
            i, j = match[b.GetBeginAtomIdx()], match[b.GetEndAtomIdx()]
            orders.setdefault(b.GetIdx(), set()).add(m.GetBondBetweenAtoms(i, j).GetBondType())
    differing = [k for k, v in orders.items() if len(v) > 1]
    out = []
    for m, match in zip(mols, matches):
        atoms = {a.GetIdx() for a in m.GetAtoms() if a.GetIdx() not in match}
        bonds = {b.GetIdx() for b in m.GetBonds() if b.GetBeginAtomIdx() in atoms or b.GetEndAtomIdx() in atoms}
        for k in differing:
            pb = patt.GetBondWithIdx(k)
            bonds.add(m.GetBondBetweenAtoms(match[pb.GetBeginAtomIdx()], match[pb.GetEndAtomIdx()]).GetIdx())
        out.append((atoms, bonds))
    if not any(a or b for a, b in out):  # même connectivité et mêmes ordres : la différence est stéréochimique
        out = []
        for m in mols:
            a_lab, b_lab = stereo_labels(m)
            out.append((set(a_lab), set(b_lab)))
    return [(sorted(a), sorted(b)) for a, b in out]


def draw(mols, hl):
    """SVG alignés sur le squelette commun (même orientation), différences surlignées."""
    rdDepictor.SetPreferCoordGen(True)
    ref = mols[0]
    rdDepictor.Compute2DCoords(ref)
    patt = Chem.MolFromSmarts(mcs(mols).smartsString)
    ref_match = ref.GetSubstructMatch(patt)
    svgs = []
    for k, (m, (atoms, bonds)) in enumerate(zip(mols, hl)):
        if k:
            try:  # forceRDKit : CoordGen ne respecte pas la correspondance d'atomes imposée
                rdDepictor.GenerateDepictionMatching2DStructure(m, ref, list(zip(ref_match, m.GetSubstructMatch(patt))), -1, True)
            except ValueError:
                rdDepictor.Compute2DCoords(m)
        d = rdMolDraw2D.MolDraw2DSVG(-1, -1)
        o = d.drawOptions()
        o.fixedBondLength = 30
        o.padding = 0.08
        o.bondLineWidth = 2
        o.clearBackground = False
        o.minFontSize = 13
        o.addStereoAnnotation = False
        o.updateAtomPalette(ATOM_COLORS)
        o.setSymbolColour((0, 0, 0))
        o.setAnnotationColour((0, 0, 0))
        o.setHighlightColour(HIGHLIGHT)
        o.highlightBondWidthMultiplier = 10
        rdMolDraw2D.PrepareAndDrawMolecule(d, m, highlightAtoms=atoms, highlightBonds=bonds,
                                          highlightAtomColors={a: HIGHLIGHT for a in atoms},
                                          highlightBondColors={b: HIGHLIGHT for b in bonds},
                                          highlightAtomRadii={a: 0.35 for a in atoms})
        d.FinishDrawing()
        svg = d.GetDrawingText()
        svg = re.sub(r"<\?xml[^>]*\?>\s*", "", svg)
        svg = re.sub(r"<!-- .*? -->\s*", "", svg)
        svg = svg.replace("#000000", "currentColor")
        svg = re.sub(r"(<svg[^>]*?)\s+width='[^']*'\s+height='[^']*'", r"\1", svg, count=1)
        svgs.append(svg)
    return svgs
