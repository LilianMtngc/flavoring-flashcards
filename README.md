# Flavoring flashcards

Flash cards pour apprendre à reconnaître les molécules aromatiques à partir de leur structure chimique,
avec un focus sur les **notes brunes** (caramel, cacao, café, grillé, noisette, érable, vanille, fumé…).

**▶ Ouvrir les cartes : https://lilianmtngc.github.io/flavoring-flashcards/**

- **Recto** : la structure chimique (ou, pour un produit naturel, sa ou ses molécules marqueurs).
- **Verso** : le nom, la famille chimique, la formule, le numéro CAS, les descripteurs, le dosage, les usages et les sources.

446 cartes : 278 molécules, 32 produits naturels, 74 molécules théoriques et 62 comparaisons de « sosies », dont 73 notes brunes (hors sosies).

Quatre catégories :
- **Molécules** et **produits naturels** du jeu de départ, avec leurs descripteurs et leur dosage d'évaluation (ppm) ;
- **Théoriques** : molécules de référence ajoutées parce qu'il serait gênant de les ignorer quand on apprend l'aromatique
  (2-acétyl-1-pyrroline, thiols du pamplemousse et du cassis, défauts comme le TCA, la géosmine ou le (E)-2-nonénal,
  énantiomères de la carvone, terpènes et esters fruités…). Leurs descripteurs viennent de leur page Wikipédia, citation à l'appui
  ([`data/theoriques_wiki.csv`](data/theoriques_wiki.csv)), ou à défaut du profil FEMA traduit terme à terme
  ([`data/glossaire_fema.csv`](data/glossaire_fema.csv)). Elles n'ont pas de dosage ;
- **Sosies** : cartes de comparaison de 2 à 4 molécules aux structures très proches (furanéol / sotolon,
  eugénol / isoeugénol, ionones / damascones, géraniol / nérol, carvones…), différence surlignée.
  Chaque molécule concernée affiche aussi un encart « À ne pas confondre avec ».

## Utilisation
- **Navigateur** : [la page en ligne](https://lilianmtngc.github.io/flavoring-flashcards/), ou `index.html`
  téléchargé et ouvert localement. Le fichier est autonome et fonctionne hors connexion.
  - Clic ou `Espace` : retourner la carte ; `←` / `→` : naviguer ; `1` : à revoir ; `2` : je sais.
  - Filtres : complexité, molécule / produit naturel / théorique / sosies, notes brunes, à revoir, je sais, famille chimique, recherche.
  - Les cartes sues sont teintées en vert, celles à revoir en rouge. Bouton clair / sombre en haut à droite.
  - La progression est enregistrée dans le navigateur, sur cet appareil uniquement.
- **Anki** : [`flavoring-flashcards.apkg`](flavoring-flashcards.apkg), à importer avec Fichier → Importer.
  Les cartes sont taguées `note_brune`, `complexite_*`, `famille_*` et `molecule` / `naturel` / `theorique`.

## D'où viennent les informations

| Information | Source |
|---|---|
| Structure (SMILES, InChIKey, formule, masse) | [PubChem](https://pubchem.ncbi.nlm.nih.gov/), recoupée avec [NCI CACTUS](https://cactus.nci.nih.gov/) et [OPSIN](https://opsin.ch.cam.ac.uk/) |
| N° FEMA, profil aromatique, n° JECFA | [FEMA Flavor Library](https://www.femaflavor.org/flavor-library), lien sur chaque carte |
| Descriptions d'odeur et de goût | PubChem, en majorité Hazardous Substances Data Bank (HSDB), lien sur chaque carte |
| Marqueurs des produits naturels | Wikipédia (anglais), lien sur chaque carte |
| Descripteurs et dosages (ppm) | Évaluation sensorielle interne ; pour les molécules théoriques, Wikipédia ou profil FEMA |
| Complexité | Calculée à partir de la structure (voir plus bas) |
| Usages et origine des notes | Synthèse rédigée avec assistance IA, non sourcée ligne à ligne |

### Vérification des structures
Chaque molécule est cherchée dans 5 sources : PubChem par le nom, PubChem par le n° CAS, NCI CACTUS par le nom
et par le CAS, et OPSIN. On compare la connectivité (1er bloc de l'InChIKey).

La structure retenue est celle de la majorité. Elle doit être confirmée par au moins 3 sources, dont une hors PubChem.

Des contrôles chimiques automatiques s'ajoutent : cycle des γ- et δ-lactones, groupes fonctionnels attendus d'après le nom,
présence de N et S, stéréo E/Z ou énantiomère quand le nom l'exige.

Les rares désaccords ont été tranchés à la main, avec la raison tracée. Le détail figure dans
[`verification.html`](https://lilianmtngc.github.io/flavoring-flashcards/verification.html).

### Marqueurs des produits naturels
Un couple « produit → molécule » n'est retenu que dans deux cas :
- la page Wikipédia du produit cite la molécule ;
- la page de la molécule cite le produit dans une phrase qui indique sa présence (« found in », « component of », « responsible for »…).

Une molécule qui « sent comme » le produit ne suffit pas.

Les liens retenus ont été relus. Les produits sans marqueur confirmé n'ont pas de carte.

### Complexité (Facile / Moyen / Difficile)
Elle est calculée à partir de l'indice de complexité topologique de Bertz
(S. H. Bertz, *J. Am. Chem. Soc.* 1981, 103, 3599), via RDKit (`BertzCT`). Les trois niveaux correspondent
aux terciles de cet indice sur l'ensemble des molécules :
- Facile : indice < 116 ;
- Moyen : de 116 à 196 ;
- Difficile : ≥ 196.

Pour une carte de sosies, on retient la molécule la plus complexe.

Pour un produit naturel, on retient le marqueur le plus complexe.

### Une carte par structure
Deux cartes qui auraient exactement le même recto (même structure, stéréochimie comprise, comparée par l'InChIKey complète)
sont fusionnées. La carte porte le nom de la molécule ; les produits naturels dont elle est l'unique marqueur sont cités
(« Marqueur de : … »), avec leurs propres descripteurs et dosages. Exemple : filbertone et noisette.

### Sosies
Calculés avec RDKit ([`scripts/sosies.py`](scripts/sosies.py)). Deux molécules sont sosies si :
- elles ont la même formule brute, ou ne diffèrent que d'une double liaison C=C en gardant les mêmes fonctions ;
- leur plus grande sous-structure commune (ordres de liaison ignorés) couvre tous leurs atomes sauf un
  (sauf deux pour des isomères à cycles non aromatiques, comme furanéol / sotolon) ;
- elles sont toutes deux des acides carboxyliques, ou aucune des deux.

Les groupes comptent au plus 4 molécules, toutes sosies deux à deux ; les stéréo-isomères forment un groupe à part.
Les dessins sont alignés sur le squelette commun. Sont surlignés : les atomes placés ailleurs, les liaisons dont l'ordre
change, ou les centres stéréogènes et doubles liaisons E / Z. La nature de la différence est décrite automatiquement.

### Notes brunes
Une carte est classée « note brune » dans deux cas :
- un descripteur brun figure parmi ses deux premiers descripteurs ;
- elle compte au moins deux descripteurs bruns.

Descripteurs bruns : caramel, cacao, café, grillé, noisette, nutty, praliné, érable, vanillé, fumé, malt, biscuité, pain, pop-corn,
riz basmati, brûlé, sucre brun, tabac, noix, fenugrec…

### Ce qui n'est pas publié
- **Correspondances fondées sur un nom commercial** que PubChem ou CACTUS ne confirment pas.
- **Mélanges de composition non communiquée** et produits sans marqueur sourcé.

## Reproduire
```bash
pip install -r requirements.txt
python scripts/02_resolve.py      # PubChem, CACTUS, OPSIN (cache : data/resolve_cache.json)
python scripts/03_verify.py       # vote entre sources et contrôles chimiques -> verification.html
python scripts/04_render.py       # dessins 2D (RDKit) -> structures/*.svg
python scripts/06_sources.py      # FEMA et PubChem (cache : data/sources_cache.json)
python scripts/07_markers.py      # sources Wikipédia des marqueurs -> data/marqueurs_sources.json
python scripts/05_build_cards.py  # index.html, flavoring-flashcards.apkg, data/cartes.json
```

| Fichier | Contenu |
|---|---|
| `data/descripteurs.csv` | descripteurs et dosage (ppm) par molécule ou produit, et par forme (naturelle / non précisée) |
| `data/glossaire_fema.csv` | traduction des termes du profil FEMA (molécules théoriques) |
| `data/molecules_ref.csv` | nom affiché, nom de requête et n° CAS de chaque molécule (`note=theorique` : molécule théorique) |
| `data/naturels_ref.csv`, `data/naturels_wiki.csv` | produits naturels, marqueurs, pages Wikipédia consultées |
| `data/revue_manuelle.csv`, `data/overrides.csv` | décisions manuelles de vérification, avec leur justification |
| `data/usages.csv` | usages et origine des notes |

## Credit

- Lilian Montagnac (Bioinformatics Engineer)
- https://github.com/LilianMtngc
- AI Assisted (Opus 5.5)
