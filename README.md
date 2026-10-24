<h1><img src="icon-192.png" alt="" width="48" height="48" align="absmiddle">&nbsp;AromaCards</h1>

Flash cards pour apprendre à reconnaître les molécules aromatiques à partir de leur structure chimique,
avec un focus sur les **notes brunes** (caramel, cacao, café, grillé, noisette, érable, vanille, fumé…).

**▶ Ouvrir les cartes : https://lilianmtngc.github.io/aromacards/**

- **Recto** : la structure chimique (ou, pour un produit naturel, sa ou ses molécules marqueurs).
- **Verso** : le nom, la famille chimique, la formule, le numéro CAS, les descripteurs, le dosage, les usages et les sources.

459 cartes : 278 molécules, 32 produits naturels, 86 molécules théoriques et 63 comparaisons d'analogues, dont 73 notes brunes (hors analogues), 31 molécules de défaut et 8 cartes avec une mention CMR / CIRC.

Quatre catégories :
- **Molécules** et **produits naturels** du jeu de départ, avec leurs descripteurs et leur dosage d'évaluation (ppm) ;
- **Théoriques** : molécules de référence ajoutées parce qu'il serait gênant de les ignorer quand on apprend l'aromatique
  (2-acétyl-1-pyrroline, thiols du pamplemousse et du cassis, défauts comme le TCA, la géosmine ou le (E)-2-nonénal,
  énantiomères de la carvone et du limonène, terpènes, esters fruités, agents rafraîchissants, chauffants et picotants…). Leurs descripteurs viennent de leur page Wikipédia, citation à l'appui
  ([`data/theoriques_wiki.csv`](data/theoriques_wiki.csv)), ou à défaut du profil FEMA traduit terme à terme
  ([`data/glossaire_fema.csv`](data/glossaire_fema.csv)). Elles n'ont pas de dosage ;
- **Analogues** (onglet à part, à côté d'Étude et Catalogue) : comparaisons de 2 à 4 molécules aux structures très proches
  (furanéol / sotolon, eugénol / isoeugénol, ionones / damascones, géraniol / nérol, carvones…), différence surlignée,
  en fiches ou en liste. Chaque molécule concernée affiche aussi un encart « À ne pas confondre avec ».

## Utilisation
- **Navigateur** : [la page en ligne](https://lilianmtngc.github.io/aromacards/), ou `index.html`
  téléchargé et ouvert localement. Le fichier est autonome et fonctionne hors connexion.
  - Clic ou `Espace` : retourner la carte ; `←` / `→` : naviguer ; `1` : à revoir ; `2` : je sais.
  - En Étude, l'ordre est aléatoire : un tirage sans répétition parmi les cartes filtrées, puis un nouveau tirage ;
    les cartes vues récemment passent en fin de tirage. « Précédent » revient sur les cartes déjà vues.
    Le Catalogue garde son tri alphabétique ou par famille.
  - Sur ordinateur, les filtres sont dans la colonne de gauche et le verso s'affiche sur deux colonnes ;
    sur mobile, ils s'ouvrent avec le bouton « Filtres ».
  - Filtres : type de carte, complexité, progression (pas encore vues, à revoir, je sais, à réviser aujourd'hui),
    défauts, CMR / CIRC, famille chimique, familles olfactives (dont notes brunes et sensations en bouche), étiquettes.
    Les familles olfactives et les étiquettes se plient et se déplient. Plusieurs familles se combinent en **OU**
    (au moins une) ou en **ET** (toutes à la fois). « Sensations en bouche » se sélectionne en entier (toutes les sensations)
    ou se déplie pour choisir une sensation précise.
    Les filtres actifs s'affichent en pastilles, effaçables d'un clic.
  - **Recto** au choix (en haut à droite de la carte) : la structure (par défaut), le nom (retrouver la structure parmi 4)
    ou les descripteurs (retrouver la molécule parmi 4 ; les leurres sont les analogues et les molécules aux descripteurs proches).
  - **Répétition espacée** : « Je sais » repousse la carte (1, 3, 7, 14, 30 puis 60 jours), « À revoir » la garde pour
    aujourd'hui. Le filtre « À réviser aujourd'hui » (ou le lien de l'onglet Statistiques) ouvre les cartes du jour.
  - **Statistiques** : onglet sur le bord droit (progression par type, complexité et famille olfactive, cartes les plus ratées).
  - **Catalogue** : tri alphabétique (par défaut) ou par famille chimique, bouton en haut à droite.
  - **Vocabulaire** : onglet sur le bord gauche, par ordre alphabétique. Noms usuels et pièges (valéraldéhyde = pentanal, caproïque / caprylique /
    caprique, « aldéhydes C-n »…), stéréochimie, préfixes, fonctions, hétérocycles, terpènes, chimie des arômes, identifiants.
    Chaque entrée renvoie vers les cartes concernées ([`data/glossaire.csv`](data/glossaire.csv), rédigé avec assistance IA).
  - **Raccourcis** : `/` rechercher, `F` filtres, `E` étiquette, `M` changer de recto, `S` statistiques, `V` vocabulaire, `D` descripteur personnel,
    `?` liste des raccourcis (aussi sous la carte : « Autres raccourcis »).
  - L'adresse de la page suit la carte affichée et les filtres : on peut la copier pour partager une carte.
  - **Installable** sur mobile (« Ajouter à l'écran d'accueil ») et utilisable hors connexion.
  - **Étiquettes personnelles** : « + Étiquette » au verso d'une carte. Nom et couleur au choix ; on peut ensuite filtrer,
    renommer, recolorer ou supprimer une étiquette depuis le panneau des filtres.
  - **Descripteurs personnels** : bouton « + » après les descripteurs d'une carte. Ils s'affichent en couleur sous les
    descripteurs communs (rubrique « Mes descripteurs »), se modifient ou se suppriment d'un clic, et la recherche les trouve.
    Ils ne changent ni les descripteurs communs ni leurs sources, et restent distincts des étiquettes.
  - Boutons ↻ « À revoir » et ✓ « Je sais » sous la carte : le contour des cartes sues devient vert, celui des cartes
    à revoir rouge. Bouton clair / sombre en haut à droite.
  - La progression, les révisions, les étiquettes et les descripteurs personnels sont enregistrés dans le navigateur, sur cet appareil uniquement.
    « Exporter » (onglet Statistiques, rubrique « Mes données ») les sauvegarde dans un fichier JSON, à réimporter ailleurs.
- **Anki** : [`aromacards.apkg`](aromacards.apkg), à importer avec Fichier → Importer.
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
| Classification CMR et CIRC | PubChem : classification harmonisée européenne (règlement CLP, annexe VI) et groupes du CIRC |
| Défauts | Liste éditoriale ([`data/defauts.csv`](data/defauts.csv)), rédigée avec assistance IA |
| Familles olfactives | Calculées à partir des descripteurs ([`data/familles_olfactives.csv`](data/familles_olfactives.csv)) |

### Vérification des structures
Chaque molécule est cherchée dans 5 sources : PubChem par le nom, PubChem par le n° CAS, NCI CACTUS par le nom
et par le CAS, et OPSIN. On compare la connectivité (1er bloc de l'InChIKey).

La structure retenue est celle de la majorité. Elle doit être confirmée par au moins 3 sources, dont une hors PubChem.

Des contrôles chimiques automatiques s'ajoutent : cycle des γ- et δ-lactones, groupes fonctionnels attendus d'après le nom,
présence de N et S, stéréo E/Z ou énantiomère quand le nom l'exige.

Les rares désaccords ont été tranchés à la main, avec la raison tracée. Le détail figure dans
[`verification.html`](https://lilianmtngc.github.io/aromacards/verification.html).

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
- Facile : indice < 119 ;
- Moyen : de 119 à 200 ;
- Difficile : ≥ 200.

Pour une comparaison d'analogues, on retient la molécule la plus complexe.

Pour un produit naturel, on retient le marqueur le plus complexe.

### Une carte par structure
Deux cartes qui auraient exactement le même recto (même structure, stéréochimie comprise, comparée par l'InChIKey complète)
sont fusionnées. La carte porte le nom de la molécule ; les produits naturels dont elle est l'unique marqueur sont cités
(« Marqueur de : … »), avec leurs propres descripteurs et dosages. Exemple : filbertone et noisette.

### Analogues
Calculés avec RDKit ([`scripts/sosies.py`](scripts/sosies.py)). Deux molécules sont analogues si :
- elles ont la même formule brute, ou ne diffèrent que d'une double liaison C=C en gardant les mêmes fonctions ;
- leur plus grande sous-structure commune (ordres de liaison ignorés) couvre tous leurs atomes sauf un
  (sauf deux pour des isomères à cycles non aromatiques, comme furanéol / sotolon) ;
- elles sont toutes deux des acides carboxyliques, ou aucune des deux.

Les groupes comptent au plus 4 molécules, toutes analogues deux à deux ; les stéréo-isomères forment un groupe à part.
Les dessins sont alignés sur le squelette commun. Un halo orange semi-transparent souligne les liaisons vers les atomes
placés ailleurs, les liaisons dont l'ordre change, ou les coins de stéréochimie et doubles liaisons E / Z. La nature de la différence est décrite automatiquement.

### Familles olfactives
Chaque descripteur est rattaché à une famille d'odeur (hespéridé, fruité, floral, vert, lacté, épicé, amandé, boisé, camphré, soufré,
terreux, phénolique…) ou à une sensation en bouche (sucré, salé, amer, umami, acide, piquant, froid, autres sensations).
« Froid » regroupe la fraîcheur mentholée ; le camphré est une famille d'odeur (voir l'entrée du vocabulaire).
Une famille à laquelle aucune carte n'appartient n'apparaît pas dans les filtres.
Une carte appartient à une famille d'odeur selon la même règle que les notes brunes : un descripteur de la famille parmi ses deux premiers,
ou au moins deux descripteurs de la famille. Pour une sensation en bouche, il suffit que la carte porte le descripteur, quelle que soit sa position.

### Mention CMR / CIRC
Une carte porte la mention dans deux cas :
- la substance a une **classification harmonisée européenne** (règlement CLP, annexe VI) comme cancérogène (H350, H351),
  mutagène (H340, H341) ou toxique pour la reproduction (H360, H361) ;
- le **CIRC** la classe dans le groupe 1, 2A ou 2B.

Ces données viennent de PubChem ([`scripts/08_securite.py`](scripts/08_securite.py)). Il s'agit du danger de la substance pure,
pas d'une évaluation du risque aux doses d'usage alimentaire. Les autoclassifications d'autres pays ne sont pas retenues.

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
python scripts/08_securite.py     # classifications CLP harmonisées et CIRC (cache : data/securite_cache.json)
python scripts/05_build_cards.py  # index.html, aromacards.apkg, data/cartes.json
```

| Fichier | Contenu |
|---|---|
| `data/descripteurs.csv` | descripteurs et dosage (ppm) par molécule ou produit, et par forme (naturelle / non précisée) |
| `data/glossaire_fema.csv` | traduction des termes du profil FEMA (molécules théoriques) |
| `data/molecules_ref.csv` | nom affiché, nom de requête et n° CAS de chaque molécule (`note=theorique` : molécule théorique) |
| `data/naturels_ref.csv`, `data/naturels_wiki.csv` | produits naturels, marqueurs, pages Wikipédia consultées |
| `data/revue_manuelle.csv`, `data/overrides.csv` | décisions manuelles de vérification, avec leur justification |
| `data/usages.csv` | usages et origine des notes |
| `data/identifiants.csv` | identifiants des molécules et produits évalués, par forme (et version quand il y en a plusieurs) |
| `fonts/` | police du titre : [Young Serif](https://github.com/noirblancrouge/YoungSerif), réduite aux lettres de « AromaCards » (licence SIL OFL 1.1) |
