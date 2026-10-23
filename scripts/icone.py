"""Icône du site : deux flash cards, un p-cymène vertical sur celle du dessus.

- version « application » (fond en dégradé, foncé en bas à gauche, plus clair en haut à droite) : icon-192.png, icon-512.png ;
- version « onglet » sans fond : favicon.svg (net à toutes les tailles) et favicon-32.png / favicon-64.png.
Le dessin est décrit une seule fois (coordonnées sur 512) et rendu en SVG et en PNG (Pillow, suréchantillonné).
"""
import io
import math

from PIL import Image, ImageDraw, ImageFilter

DARK, LIGHT = (22, 11, 5), (150, 92, 45)          # dégradé du fond : bas gauche → haut droit
BACK, FRONT, INK = (234, 165, 90), (246, 232, 207), (58, 34, 16)


def _rot(cx, cy, a, pts):
    """Rotation (degrés, sens horaire à l'écran) autour de (cx, cy)."""
    r = math.radians(a)
    return [(cx + x * math.cos(r) - y * math.sin(r), cy + x * math.sin(r) + y * math.cos(r)) for x, y in pts]


def geometry(fill=1.0):
    """Cartes et traits de la molécule ; `fill` agrandit l'ensemble (version sans fond)."""
    s = fill
    c0 = (256 + 30 * s, 256 - 6 * s)                 # carte de derrière
    c1 = (256 - 18 * s, 262 + 4 * s)                 # carte du dessus
    back = dict(c=c0, w=250 * s, h=330 * s, r=34 * s, a=-12, color=BACK)
    front = dict(c=c1, w=262 * s, h=344 * s, r=34 * s, a=6, color=FRONT)
    # p-cymène vertical, coordonnées relatives au centre de la carte du dessus
    R, b = 62 * s, 54 * s                            # rayon du cycle, longueur de liaison
    cy = -22 * s
    ring = [(R * math.cos(math.radians(-90 + 60 * k)), cy + R * math.sin(math.radians(-90 + 60 * k))) for k in range(6)]
    top, bottom = ring[0], ring[3]
    methyl = [top, (top[0], top[1] - b)]
    ch = (bottom[0], bottom[1] + b)
    iso = [[bottom, ch], [ch, (ch[0] - b * math.sin(math.radians(60)), ch[1] + b * math.cos(math.radians(60)))],
           [ch, (ch[0] + b * math.sin(math.radians(60)), ch[1] + b * math.cos(math.radians(60)))]]
    doubles = []
    for k in (1, 3, 5):                              # doubles liaisons alternées, côté intérieur
        (x1, y1), (x2, y2) = ring[k], ring[(k + 1) % 6]
        f, t = 0.36, 0.2
        ix1, iy1 = x1 + (0 - x1) * f, y1 + (cy - y1) * f
        ix2, iy2 = x2 + (0 - x2) * f, y2 + (cy - y2) * f
        doubles.append([(ix1 + (ix2 - ix1) * t, iy1 + (iy2 - iy1) * t), (ix2 + (ix1 - ix2) * t, iy2 + (iy1 - iy2) * t)])
    w = 13 * s
    lines = [(ring + [ring[0]], w), (methyl, w)] + [(l, w) for l in iso] + [(l, w * 0.75) for l in doubles]
    lines = [(_rot(c1[0], c1[1], front["a"], pts), width) for pts, width in lines]
    return back, front, lines


# ---------------------------------------------------------------- SVG

def svg(background=True):
    back, front, lines = geometry(1.0 if background else 1.22)
    out = ["<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 512 512'>"]
    if background:
        out.append("<defs><linearGradient id='bg' x1='0' y1='1' x2='1' y2='0'>"
                   f"<stop offset='0' stop-color='rgb{DARK}'/><stop offset='1' stop-color='rgb{LIGHT}'/></linearGradient>"
                   "<filter id='sh' x='-20%' y='-20%' width='140%' height='140%'><feDropShadow dx='6' dy='12' stdDeviation='12' "
                   "flood-color='#000' flood-opacity='.35'/></filter></defs>"
                   "<rect width='512' height='512' rx='112' fill='url(#bg)'/>")
    for card in (back, front):
        (cx, cy), w, h = card["c"], card["w"], card["h"]
        out.append(f"<rect x='{cx - w / 2:.1f}' y='{cy - h / 2:.1f}' width='{w:.1f}' height='{h:.1f}' rx='{card['r']:.1f}' "
                   f"fill='rgb{card['color']}' transform='rotate({card['a']} {cx:.1f} {cy:.1f})'"
                   + (" filter='url(#sh)'" if background else "") + "/>")
    for pts, width in lines:
        d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        out.append(f"<path d='{d}' fill='none' stroke='rgb{INK}' stroke-width='{width:.1f}' stroke-linecap='round' stroke-linejoin='round'/>")
    out.append("</svg>")
    return "".join(out)


# ---------------------------------------------------------------- PNG

def png(size, background=True):
    k = 4 * size / 512                               # suréchantillonnage x4 puis réduction
    S = int(512 * k)
    back, front, lines = geometry(1.0 if background else 1.22)
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    if background:
        grad = Image.new("RGBA", (S, S))
        px = grad.load()
        for y in range(S):
            for x in range(S):
                t = (x + (S - 1 - y)) / (2 * (S - 1))
                px[x, y] = tuple(int(DARK[i] + (LIGHT[i] - DARK[i]) * t) for i in range(3)) + (255,)
        mask = Image.new("L", (S, S), 0)
        ImageDraw.Draw(mask).rounded_rectangle([0, 0, S - 1, S - 1], radius=int(112 * k), fill=255)
        img.paste(grad, (0, 0), mask)

    def card_layer(card, color, offset=(0, 0)):
        layer = Image.new("RGBA", (S, S), (0, 0, 0, 0))
        (cx, cy), w, h = card["c"], card["w"], card["h"]
        cx, cy = (cx + offset[0]) * k, (cy + offset[1]) * k
        ImageDraw.Draw(layer).rounded_rectangle([cx - w * k / 2, cy - h * k / 2, cx + w * k / 2, cy + h * k / 2],
                                                radius=int(card["r"] * k), fill=color)
        return layer.rotate(-card["a"], resample=Image.BICUBIC, center=(cx, cy))

    for card in (back, front):
        if background:
            sh = card_layer(card, (0, 0, 0, 95), (6, 12)).filter(ImageFilter.GaussianBlur(12 * k))
            img = Image.alpha_composite(img, sh)
        img = Image.alpha_composite(img, card_layer(card, card["color"] + (255,)))
    d = ImageDraw.Draw(img)
    for pts, width in lines:
        p = [(x * k, y * k) for x, y in pts]
        d.line(p, fill=INK + (255,), width=int(width * k), joint="curve")
        for x, y in (p[0], p[-1]):
            r = width * k / 2
            d.ellipse([x - r, y - r, x + r, y + r], fill=INK + (255,))
    return img.resize((size, size), Image.LANCZOS)


def write_all(root):
    """Écrit les icônes du site dans `root`."""
    for n in (192, 512):
        png(n, True).save(root / f"icon-{n}.png")
    png(180, True).save(root / "apple-touch-icon.png")
    for n in (32, 64):
        png(n, False).save(root / f"favicon-{n}.png")
    (root / "favicon.svg").write_text(svg(False), encoding="utf-8")
