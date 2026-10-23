"""Étape 4 — dessine chaque structure retenue en SVG 2D (RDKit).

Style homogène : même longueur de liaison pour toutes les molécules,
traits en `currentColor` (le dessin suit la couleur du texte, donc le mode
sombre), hétéroatomes colorés, stéréo (coins / E-Z) visible.
"""
import json
import re
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "molecules_verified.json"
OUT = ROOT / "structures"

ATOM_COLORS = {8: (0.82, 0.19, 0.18), 7: (0.18, 0.44, 0.86), 16: (0.75, 0.54, 0.0), 11: (0.55, 0.35, 0.8)}


def render(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    opts = drawer.drawOptions()
    opts.fixedBondLength = 30
    opts.padding = 0.08
    opts.bondLineWidth = 2
    opts.clearBackground = False
    opts.minFontSize = 13
    opts.addStereoAnnotation = False
    opts.updateAtomPalette(ATOM_COLORS)
    opts.setSymbolColour((0, 0, 0))
    opts.setAnnotationColour((0, 0, 0))
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = re.sub(r"<\?xml[^>]*\?>\s*", "", svg)
    svg = re.sub(r"<!-- .*? -->\s*", "", svg)
    svg = svg.replace("#000000", "currentColor")
    # taille gérée par le CSS : on garde seulement le viewBox
    svg = re.sub(r"(<svg[^>]*?)\s+width='[^']*'\s+height='[^']*'", r"\1", svg, count=1)
    return svg


def main():
    data = json.loads(DATA.read_text(encoding="utf-8"))
    OUT.mkdir(exist_ok=True)
    n = 0
    for key, m in data.items():
        if not m.get("smiles"):
            continue
        (OUT / f"{key}.svg").write_text(render(m["smiles"]), encoding="utf-8")
        n += 1
    print(f"{n} structures -> {OUT}")


if __name__ == "__main__":
    main()
