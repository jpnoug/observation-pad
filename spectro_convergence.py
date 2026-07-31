#!/usr/bin/env python3
"""
Spectro Convergence — visualisation rapide de l'empilement d'une série de poses.

Compagnon de l'Observation Pad. On sélectionne à la main les FITS d'une cible,
le module détecte tout seul l'inclinaison du spectre et les zones d'extraction,
puis affiche :

  1. le profil spatial, avec les zones détectées (cible + fond) surlignées
  2. la convergence du signal : spectre moyen et dispersion ±2 sigma,
     avec les deux estimations de SNR reportées sur le graphe

Les poses longues sur objet faible sont truffées de pixels chauds et de
cosmiques, qui peuvent dépasser cent fois le signal utile. Deux protections :

  - le fond de ciel est estimé par la MÉDIANE des lignes des bandes de fond.
    Avec une moyenne, un seul pixel chaud dans le fond était multiplié par le
    nombre de lignes de la cible lors de la soustraction et creusait un cratère
    de plusieurs milliers d'ADU dans le spectre.
  - dans la zone cible, les aberrants sont rejetés sur critère de profil
    spatial (principe de Horne 1986) : le signal réel vérifie f(x)·P(y), donc
    une raie d'émission — qui monte sur toute la hauteur — est intégralement
    conservée, alors qu'un pixel chaud, limité à une ligne, est écarté.

Les deux estimateurs de SNR sont EMPIRIQUES : ils se mesurent directement sur
les données, sans gain ni bruit de lecture. Aucun réglage de caméra n'est donc
nécessaire, le module marche indifféremment avec l'Atik ou la ZWO.

  DER_SNR  — Stoehr et al. 2008, ST-ECF Newsletter 42. Robuste aux raies et aux
             pixels chauds (il n'utilise que des médianes).
  Buil     — soustraction de deux poses consécutives : le bruit réel apparaît
             dans la différence (variance comptée double, d'où le /racine(2)).

Dépendances : numpy, astropy, matplotlib, tkinter. Scipy n'est pas utilisé.

Usage :
    python spectro_convergence.py                 # ouvre le sélecteur de fichiers
    python spectro_convergence.py f1.fits f2.fits # charge directement

Depuis l'Observation Pad (bouton SNR), la fenêtre s'ouvre en Toplevel fille et
le bouton « Log » écrit la synthèse directement dans le journal d'observation.
"""

import os
import sys
import datetime
import tkinter as tk
from tkinter import filedialog, messagebox

import numpy as np

try:
    from astropy.io import fits
except ImportError:
    print("ERREUR : astropy est requis.  pip install astropy")
    sys.exit(1)

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# ─── Thème sombre/rouge, identique à l'Observation Pad ────────────────────────
# Ces valeurs doivent rester synchronisées avec observation_pad.py :
# même rouge vif sur noir, mêmes polices, mêmes graisses.
BG      = "#000000"
FG      = "#ff3333"   # rouge vif, bon contraste sur fond noir
FG_DIM  = "#c0392b"   # rouge moyen (labels, axes, statut)
BTN_BG  = "#1a0a08"
FONT_UI  = ("DejaVu Sans", 9)
FONT_TT  = ("DejaVu Sans Mono", 10)

# Couleurs matplotlib
C_LINE  = "#e05545"   # courbes principales
C_BAND  = "#c0392b"   # bandeau de dispersion
C_TGT   = "#f0a000"   # zone cible
C_BG    = "#4060a0"   # zones de fond
C_AXIS  = "#c0392b"   # axes/cadres : même rouge moyen que FG_DIM du pad

# ─── Géométrie automatique des zones, en multiples de la FWHM ────────────────
K_TGT_DEFAUT  = 1.25   # demi-hauteur de la zone cible      (±1.25 FWHM ≈ 99.7 %)
K_GAP_DEFAUT  = 1.00   # garde entre le bord cible et le fond
K_BG_WIDTH    = 1.50   # épaisseur de chaque bande de fond

FITS_TYPES = [("FITS files", "*.fits *.fit *.fts *.FITS *.FIT *.FTS"),
              ("All files", "*.*")]


# ══════════════════════════════════════════════════════════════════════════════
#  Calcul
# ══════════════════════════════════════════════════════════════════════════════

def lire_fits(chemin):
    """Retourne (data float32, en-tête utile) ou None si illisible."""
    try:
        with fits.open(chemin, memmap=False) as hdul:
            hdu = next((h for h in hdul if h.data is not None), None)
            if hdu is None:
                return None
            data = np.asarray(hdu.data, dtype=np.float32)
            if data.ndim == 3:            # cube couleur ou monoplan : on aplatit
                data = data.mean(axis=0)
            if data.ndim != 2:
                return None
            h = hdu.header
            return {
                "data":     data,
                "fichier":  os.path.basename(chemin),
                "chemin":   chemin,
                "objet":    str(h.get("OBJECT", "")).strip(),
                "date_obs": str(h.get("DATE-OBS", "")),
                "exptime":  h.get("EXPTIME", h.get("EXPOSURE", None)),
            }
    except Exception:
        return None


def detecter_inclinaison(data, n_blocs=16):
    """
    Mesure l'inclinaison du spectre par régression du centroïde spatial
    calculé sur des tranches successives de colonnes.

    Retourne la pente en pixels de décalage vertical par pixel horizontal.
    """
    h, w = data.shape
    # Centre approximatif : profil des colonnes centrales, peu affectées par l'inclinaison
    c0, c1 = int(w * 0.40), int(w * 0.60)
    prof = data[:, c0:c1].sum(axis=1)
    prof = prof - np.median(prof)
    if prof.max() <= 0:
        return 0.0
    y_approx = int(np.argmax(prof))

    # Fenêtre de recherche généreuse autour du centre approximatif
    demi = max(20, h // 20)
    y0, y1 = max(0, y_approx - demi), min(h, y_approx + demi + 1)
    lignes = np.arange(y0, y1)

    xs, ys = [], []
    bord = int(w * 0.05)                       # on écarte les bords, souvent bruités
    bloc = max(1, (w - 2 * bord) // n_blocs)
    for x in range(bord, w - bord - bloc, bloc):
        tranche = data[y0:y1, x:x + bloc].sum(axis=1)
        tranche = tranche - np.median(tranche)
        tranche[tranche < 0] = 0.0
        s = tranche.sum()
        if s <= 0:
            continue
        xs.append(x + bloc / 2.0)
        ys.append(float((lignes * tranche).sum() / s))

    if len(xs) < 4:
        return 0.0
    xs, ys = np.asarray(xs), np.asarray(ys)
    # Ajustement robuste : on rejette les blocs à plus de 2 sigma du premier fit
    pente, ordonnee = np.polyfit(xs, ys, 1)
    residus = ys - (pente * xs + ordonnee)
    bons = np.abs(residus) < 2.0 * max(residus.std(), 1e-6)
    if bons.sum() >= 4:
        pente, _ = np.polyfit(xs[bons], ys[bons], 1)
    return float(pente)


def redresser(data, pente):
    """
    Redresse l'image en décalant chaque colonne d'un nombre entier de lignes,
    de façon à rendre le spectre horizontal. Pas d'interpolation, donc pas de
    lissage parasite du bruit — ce qui compte pour une mesure de SNR.
    """
    if abs(pente) < 1e-4:
        return data
    h, w = data.shape
    xs = np.arange(w)
    decal = np.round(pente * (xs - w / 2.0)).astype(int)
    yy = np.arange(h)[:, None] + decal[None, :]
    np.clip(yy, 0, h - 1, out=yy)
    return data[yy, xs[None, :]]


def profil_spatial(data):
    """Somme de l'image le long de l'axe de dispersion : une valeur par ligne."""
    return data.sum(axis=1).astype(np.float64)


def mesurer_source(prof):
    """
    Sur le profil spatial, retourne (centre, fwhm) de la source.
    Le centre est un centroïde pondéré, la FWHM une mesure directe à mi-hauteur.
    """
    base = np.median(prof)
    p = prof - base
    pic = float(p.max())
    if pic <= 0:
        return len(prof) / 2.0, 10.0
    i_pic = int(np.argmax(p))

    # Largeur à mi-hauteur : on descend de part et d'autre du pic, puis on
    # interpole linéairement le croisement pour une précision sous-pixel
    moitie = pic / 2.0
    g = i_pic
    while g > 0 and p[g] > moitie:
        g -= 1
    d = i_pic
    while d < len(p) - 1 and p[d] > moitie:
        d += 1
    x_g = float(g)
    if p[g + 1] > p[g]:
        x_g = g + (moitie - p[g]) / (p[g + 1] - p[g])
    x_d = float(d)
    if p[d - 1] > p[d]:
        x_d = d - (moitie - p[d]) / (p[d - 1] - p[d])
    fwhm = max(x_d - x_g, 2.0)

    # Centroïde sur ±3 FWHM autour du pic, plus stable que le simple argmax
    a = max(0, int(i_pic - 3 * fwhm))
    b = min(len(p), int(i_pic + 3 * fwhm) + 1)
    fen = np.clip(p[a:b], 0, None)
    s = fen.sum()
    centre = float((np.arange(a, b) * fen).sum() / s) if s > 0 else float(i_pic)
    return centre, fwhm


def zones_auto(prof, k_tgt=K_TGT_DEFAUT, k_gap=K_GAP_DEFAUT):
    """
    Place automatiquement la zone cible et les deux bandes de fond à partir du
    profil spatial. Retourne un dictionnaire de bornes entières + un diagnostic.
    """
    h = len(prof)
    centre, fwhm = mesurer_source(prof)

    demi_t = k_tgt * fwhm
    t0, t1 = int(round(centre - demi_t)), int(round(centre + demi_t))

    depart = demi_t + k_gap * fwhm
    largeur = K_BG_WIDTH * fwhm
    ha0, ha1 = int(round(centre - depart - largeur)), int(round(centre - depart))
    hb0, hb1 = int(round(centre + depart)), int(round(centre + depart + largeur))

    avert = []
    # Contraste de la source : hauteur du pic rapportée à la dispersion robuste
    # du profil loin du pic. En dessous de 10, il n'y a probablement pas de
    # spectre exploitable (pose sur le mauvais champ, étoile hors de la fente,
    # ou sélection de darks/bias par erreur).
    p = prof - np.median(prof)
    loin = np.ones(h, dtype=bool)
    a = max(0, int(centre - 4 * fwhm)); b = min(h, int(centre + 4 * fwhm) + 1)
    loin[a:b] = False
    if loin.sum() > 10:
        mad = np.median(np.abs(p[loin] - np.median(p[loin])))
        bruit = 1.4826 * mad
        if bruit > 0 and p.max() / bruit < 10.0:
            avert.append("no clear source detected — check the selection")

    if t0 < 0 or t1 > h:
        avert.append("target region truncated by image edge")
    if ha0 < 0 or hb1 > h:
        avert.append("background bands truncated by image edge")

    z = {
        "centre": centre, "fwhm": fwhm,
        "tgt": (max(0, t0), min(h, t1)),
        "bga": (max(0, ha0), max(0, ha1)),
        "bgb": (min(h, hb0), min(h, hb1)),
        "avert": avert,
    }
    if z["bga"][1] - z["bga"][0] < 3 and z["bgb"][1] - z["bgb"][0] < 3:
        z["avert"].append("no usable background band")
    return z


def rejet_aberrants(S, profil, k=6.0, plancher=5.0):
    """
    Rejette les pixels aberrants d'une zone déjà soustraite du fond de ciel.

    Le signal réel obéit au modèle  pixel(y, x) = f(x) · P(y)  : le profil
    spatial P est le même à toutes les longueurs d'onde, seule l'intensité f(x)
    change. Une raie d'émission fait monter f(x) et respecte donc le modèle sur
    toute la hauteur — elle est intégralement conservée. Un pixel chaud ou un
    cosmique ne touche qu'une ligne et s'écarte violemment du modèle.

    C'est le principe du rejet de Horne (1986). Un simple filtre médian 2D, lui,
    raboterait le sommet des raies fines.

    Retourne (zone corrigée, nombre de pixels rejetés).
    """
    P = np.clip(np.asarray(profil, dtype=np.float64), 1e-6, None)[:, None]
    utiles = P[:, 0] > 0.2 * P.max()
    if utiles.sum() < 3:
        return S, 0
    # Intensité par colonne, estimée robustement sur les lignes où la source domine
    f = np.median(S[utiles] / P[utiles], axis=0)
    modele = f[None, :] * P
    res = S - modele
    mad = np.median(np.abs(res - np.median(res, axis=0)), axis=0)
    bruit = np.maximum(1.4826 * mad, plancher)
    masque = res > k * bruit[None, :]
    return np.where(masque, modele, S), int(masque.sum())


def extraire(data, zones, profil_norm=None, rejet=True):
    """
    Extraction sur image déjà redressée.

    Le fond est estimé par la MÉDIANE des lignes des bandes de fond, et non par
    leur moyenne : un seul pixel chaud dans le fond, multiplié par le nombre de
    lignes de la cible au moment de la soustraction, creusait sinon un cratère
    de plusieurs milliers d'ADU dans le spectre.

    Retourne (spectre soustrait du fond, fond par ligne, nb de lignes, nb rejets).
    """
    (t0, t1) = zones["tgt"]
    if t1 - t0 < 1:
        return None, None, 0, 0
    n_lignes = t1 - t0

    bandes = []
    for cle in ("bga", "bgb"):
        a, b = zones[cle]
        if b - a >= 3:
            bandes.append(data[a:b, :])
    if bandes:
        fond_par_ligne = np.median(np.concatenate(bandes, axis=0),
                                   axis=0).astype(np.float64)
    else:
        fond_par_ligne = np.zeros(data.shape[1], dtype=np.float64)

    S = data[t0:t1, :].astype(np.float64) - fond_par_ligne[None, :]

    n_rej = 0
    if rejet and profil_norm is not None:
        S, n_rej = rejet_aberrants(S, profil_norm)

    return S.sum(axis=0), fond_par_ligne, n_lignes, n_rej


def masque_utile(spec, seuil=0.15):
    """
    Colonnes exploitables pour les mesures de SNR : celles où le signal dépasse
    une fraction du maximum. Écarte les extrémités mortes du spectre, où le
    rapport signal/bruit n'a pas de sens.
    """
    if spec is None or spec.size == 0:
        return None
    pic = float(np.nanmax(spec))
    if pic <= 0:
        return None
    return spec > seuil * pic


def der_snr(flux, masque=None):
    """
    DER_SNR (Stoehr et al. 2008) : estimation empirique du SNR à partir du
    spectre lui-même.

        signal = médiane(flux)
        bruit  = 1.482602 / racine(6) * médiane(|2·f_i - f_{i-2} - f_{i+2}|)

    Le terme en différences secondes s'annule sur un continuum lisse et ne
    retient que le bruit haute fréquence : l'estimateur est donc peu sensible
    à la présence de raies et aux pixels chauds isolés.
    """
    f = flux if masque is None else flux[masque]
    f = f[np.isfinite(f)]
    if f.size < 5:
        return None
    signal = float(np.median(f))
    if signal <= 0:
        return None
    diff = np.abs(2.0 * f[2:-2] - f[:-4] - f[4:])
    bruit = 1.482602 / np.sqrt(6.0) * float(np.median(diff))
    return signal / bruit if bruit > 0 else None


def snr_buil(a, b, masque=None):
    """
    SNR par soustraction de deux poses de la même cible (méthode C. Buil / ISIS).
    La différence ne contient que du bruit ; sa variance compte deux fois celle
    d'une pose, d'où la division par racine(2).
    """
    if a is None or b is None or a.size != b.size:
        return None
    a, b = a.astype(float), b.astype(float)
    if masque is not None:
        a, b = a[masque], b[masque]
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if a.size < 5:
        return None
    signal = float(np.mean((a + b) / 2.0))
    bruit = float(np.std(a - b)) / np.sqrt(2.0)
    if bruit <= 0 or signal <= 0:
        return None
    return signal / bruit


def analyser(chemins, k_tgt=K_TGT_DEFAUT, k_gap=K_GAP_DEFAUT):
    """
    Charge, redresse, extrait et mesure. Retourne un dictionnaire de résultats
    ou lève ValueError avec un message lisible.
    """
    poses = [p for p in (lire_fits(c) for c in chemins) if p is not None]
    if not poses:
        raise ValueError("No readable 2D FITS in the selection.")

    poses.sort(key=lambda p: (p["date_obs"], p["fichier"]))

    formes = {p["data"].shape for p in poses}
    if len(formes) > 1:
        raise ValueError("Frames don't all have the same dimensions.\n"
                         "Select a single series with a single binning.")

    objets = {p["objet"] for p in poses if p["objet"]}
    objet = objets.pop() if len(objets) == 1 else (
        " / ".join(sorted(objets)) if objets else "(sans OBJECT)")
    melange = len(objets) > 0     # il en restait après le pop : plusieurs cibles

    # Géométrie mesurée sur l'empilement médian : plus robuste qu'une seule pose
    pile = np.median(np.stack([p["data"] for p in poses], axis=0), axis=0)
    pente = detecter_inclinaison(pile)
    pile_droite = redresser(pile, pente)
    prof = profil_spatial(pile_droite)
    zones = zones_auto(prof, k_tgt=k_tgt, k_gap=k_gap)

    # Profil spatial normalisé de la source, mesuré sur l'empilement médian et
    # soustrait du niveau de ciel. Empirique plutôt que gaussien : une galaxie
    # étendue comme une étoile ponctuelle sont traitées correctement.
    t0, t1 = zones["tgt"]
    niveaux_fond = []
    for cle in ("bga", "bgb"):
        a, b = zones[cle]
        if b - a >= 3:
            niveaux_fond.append(prof[a:b])
    ciel = float(np.median(np.concatenate(niveaux_fond))) if niveaux_fond else 0.0
    profil_norm = np.clip(prof[t0:t1] - ciel, 0.0, None)
    if profil_norm.max() > 0:
        profil_norm = profil_norm / profil_norm.max()
    else:
        profil_norm = None

    spectres, n_rejets = [], 0
    for p in poses:
        spec, _fond, _n, nr = extraire(redresser(p["data"], pente), zones,
                                       profil_norm=profil_norm)
        if spec is not None:
            spectres.append(spec)
            n_rejets += nr
    if not spectres:
        raise ValueError("Extraction failed: empty target region.")

    S = np.vstack(spectres)
    moyenne = S.mean(axis=0)
    sigma = S.std(axis=0, ddof=1) if S.shape[0] > 1 else np.zeros_like(moyenne)
    masque = masque_utile(moyenne)

    # DER cumulé : DER de la pose moyenne des n premières, croissance en racine(n)
    der_cumul = []
    for n in range(1, S.shape[0] + 1):
        d = der_snr(S[:n].mean(axis=0), masque)
        der_cumul.append(d if d is not None else np.nan)

    der_pose = [der_snr(s, masque) for s in S]
    der_pose_moy = float(np.nanmean([d for d in der_pose if d is not None])) \
        if any(d is not None for d in der_pose) else None

    paires = [snr_buil(S[i], S[i + 1], masque) for i in range(S.shape[0] - 1)]
    paires = [v for v in paires if v is not None]
    buil_pose = float(np.mean(paires)) if paires else None
    buil_pile = buil_pose * np.sqrt(S.shape[0]) if buil_pose is not None else None

    exps = [p["exptime"] for p in poses if p["exptime"] is not None]
    exp = float(exps[0]) if exps and len(set(map(float, exps))) == 1 else None

    return {
        "poses": poses, "objet": objet, "melange": melange,
        "pente": pente, "prof": prof, "zones": zones,
        "moyenne": moyenne, "sigma": sigma, "masque": masque,
        "n": S.shape[0], "exptime": exp, "n_rejets": n_rejets,
        "der_cumul": der_cumul, "der_pose": der_pose_moy,
        "buil_pose": buil_pose, "buil_pile": buil_pile,
    }


def ligne_resume(r):
    """Ligne de synthèse, copiable dans le journal d'observation."""
    bouts = [r["objet"], f"{r['n']} frames"]
    if r["exptime"]:
        bouts.append(f"{r['exptime']:.0f} s")
    if r["der_cumul"] and np.isfinite(r["der_cumul"][-1]):
        bouts.append(f"DER {r['der_cumul'][-1]:.0f}")
    if r["buil_pile"]:
        bouts.append(f"Buil {r['buil_pile']:.0f}")
    bouts.append(f"FWHM {r['zones']['fwhm']:.1f} px")
    return "  —  ".join(bouts)


# ══════════════════════════════════════════════════════════════════════════════
#  Interface
# ══════════════════════════════════════════════════════════════════════════════

class Application(tk.Toplevel):
    """Fenêtre d'analyse. Toplevel (et non Tk) pour pouvoir s'ouvrir comme
    fille de l'Observation Pad : le bouton « Log » écrit alors directement
    dans le journal via le rappel `on_log`. En usage autonome, main() crée
    une racine masquée et le rappel retombe sur le presse-papiers."""

    def __init__(self, master=None, chemins_initiaux=None, on_log=None):
        super().__init__(master)
        self.title("Spectro Convergence")
        self.configure(bg=BG)
        self.geometry("620x680")
        self.minsize(560, 460)

        self.resultat = None
        self.on_log = on_log
        self.dernier_dossier = os.path.expanduser("~")

        self._barre()
        self._figure()
        self._statut()

        if chemins_initiaux:
            self._charger(list(chemins_initiaux))

    # ── construction ─────────────────────────────────────────────────────────
    def _bouton(self, parent, texte, commande):
        b = tk.Button(parent, text=texte, command=commande,
                      bg=BTN_BG, fg=FG, activebackground=FG_DIM,
                      activeforeground=BG, font=FONT_UI,
                      relief="flat", bd=1, highlightthickness=0,
                      padx=10, pady=3)
        b.pack(side="left", padx=4)
        return b

    def _barre(self):
        barre = tk.Frame(self, bg=BG)
        barre.pack(side="top", fill="x", padx=8, pady=(8, 4))

        self._bouton(barre, "Files…", self._choisir)
        self._bouton(barre, "Log", self._logger)

        tk.Frame(barre, bg=FG_DIM, width=1).pack(side="left", fill="y", padx=10)

        tk.Label(barre, text="Target ±", bg=BG, fg=FG_DIM,
                 font=FONT_UI).pack(side="left")
        self.var_tgt = tk.DoubleVar(value=K_TGT_DEFAUT)
        self._spin(barre, self.var_tgt)
        tk.Label(barre, text="FWHM   Gap", bg=BG, fg=FG_DIM,
                 font=FONT_UI).pack(side="left", padx=(4, 2))
        self.var_gap = tk.DoubleVar(value=K_GAP_DEFAUT)
        self._spin(barre, self.var_gap)
        tk.Label(barre, text="FWHM", bg=BG, fg=FG_DIM,
                 font=FONT_UI).pack(side="left", padx=(4, 0))

    def _spin(self, parent, var):
        s = tk.Spinbox(parent, from_=0.25, to=6.0, increment=0.25,
                       textvariable=var, width=5, font=FONT_UI,
                       bg=BTN_BG, fg=FG, buttonbackground=BTN_BG,
                       insertbackground=FG, relief="flat",
                       highlightthickness=1, highlightbackground=FG_DIM,
                       command=self._recalculer)
        s.pack(side="left")
        s.bind("<Return>", lambda _e: self._recalculer())
        return s

    def _figure(self):
        # Figure resserrée : la fenêtre doit rester étroite à côté de CCDciel.
        # top/hspace élargis car les titres tiennent désormais sur deux lignes.
        self.fig = Figure(figsize=(5.8, 5.3), dpi=100, facecolor=BG)
        gs = self.fig.add_gridspec(2, 1, height_ratios=[1.0, 2.3],
                                   hspace=0.55,
                                   left=0.13, right=0.97, top=0.89, bottom=0.10)
        self.ax_prof = self.fig.add_subplot(gs[0])
        self.ax_conv = self.fig.add_subplot(gs[1])
        for ax in (self.ax_prof, self.ax_conv):
            self._styler(ax)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().configure(bg=BG, highlightthickness=0)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True,
                                         padx=8, pady=4)
        self._message("Select the frames of a target with « Files… »")
        self.canvas.draw()

    def _styler(self, ax):
        ax.set_facecolor(BG)
        for sp in ax.spines.values():
            sp.set_edgecolor(C_AXIS)
        ax.tick_params(colors=FG, labelsize=7.5, length=3)

    def _statut(self):
        self.var_statut = tk.StringVar(value="")
        tk.Label(self, textvariable=self.var_statut, bg=BG, fg=FG,
                 font=FONT_TT, anchor="w", justify="left"
                 ).pack(side="bottom", fill="x", padx=10, pady=(0, 8))

    def _message(self, texte):
        for ax in (self.ax_prof, self.ax_conv):
            ax.clear()
            self._styler(ax)
            ax.set_xticks([]); ax.set_yticks([])
        self.ax_conv.text(0.5, 0.5, texte, color=FG_DIM, fontsize=10,
                          ha="center", va="center",
                          transform=self.ax_conv.transAxes)

    # ── actions ──────────────────────────────────────────────────────────────
    def _choisir(self):
        chemins = filedialog.askopenfilenames(
            title="Target frames", initialdir=self.dernier_dossier,
            filetypes=FITS_TYPES)
        if chemins:
            self.dernier_dossier = os.path.dirname(chemins[0])
            self._charger(list(chemins))

    def _charger(self, chemins):
        self.chemins = chemins
        self.var_statut.set(f"Reading {len(chemins)} file(s)…")
        self.update_idletasks()
        self._recalculer()

    def _recalculer(self):
        if not getattr(self, "chemins", None):
            return
        try:
            k_t, k_g = float(self.var_tgt.get()), float(self.var_gap.get())
            self.resultat = analyser(self.chemins, k_tgt=k_t, k_gap=k_g)
        except ValueError as e:
            self.resultat = None
            self._message(str(e))
            self.var_statut.set("")
            self.canvas.draw()
            return
        except Exception as e:
            self.resultat = None
            messagebox.showerror("Error", f"{type(e).__name__}: {e}")
            return
        self._tracer()

    # ── tracé ────────────────────────────────────────────────────────────────
    def _tracer(self):
        r = self.resultat
        for ax in (self.ax_prof, self.ax_conv):
            ax.clear()
            self._styler(ax)

        # 1 — profil spatial et zones détectées
        prof, z = r["prof"], r["zones"]
        ys = np.arange(len(prof))
        self.ax_prof.plot(ys, prof, color=C_LINE, lw=0.9)
        self.ax_prof.axvspan(*z["tgt"], color=C_TGT, alpha=0.22, lw=0,
                             label="target")
        premier = True
        for cle in ("bga", "bgb"):
            a, b = z[cle]
            if b > a:
                self.ax_prof.axvspan(a, b, color=C_BG, alpha=0.28, lw=0,
                                     label="background" if premier else None)
                premier = False
        # Cadrage sur l'emprise réelle des zones, pour qu'elles restent visibles
        # quel que soit le réglage de largeur
        bords = [z["tgt"][0], z["tgt"][1]]
        for cle in ("bga", "bgb"):
            a, b = z[cle]
            if b > a:
                bords += [a, b]
        marge = max(15.0, 1.5 * z["fwhm"])
        self.ax_prof.set_xlim(max(0, min(bords) - marge),
                              min(len(prof), max(bords) + marge))
        leg = self.ax_prof.legend(loc="upper right", fontsize=7,
                                  facecolor=BTN_BG, edgecolor=C_AXIS)
        for t in leg.get_texts():
            t.set_color(FG)
        t0, t1 = z["tgt"]
        # Titre sur deux lignes : réduit la largeur minimale de la fenêtre
        self.ax_prof.set_title(
            "Spatial profile\n"
            f"centre {z['centre']:.1f} px · FWHM {z['fwhm']:.1f} px · "
            f"tilt {np.degrees(np.arctan(r['pente'])):.2f}° · "
            f"width {t1 - t0} px",
            color=FG, fontsize=9, pad=4)
        self.ax_prof.set_xlabel("Y pixel", color=FG, fontsize=8)
        self.ax_prof.set_yticks([])

        # 2 — convergence : moyenne et dispersion ±2 sigma
        moy, sig = r["moyenne"], r["sigma"]
        xs = np.arange(moy.size)
        haut, bas = moy + 2 * sig, moy - 2 * sig
        self.ax_conv.fill_between(xs, bas, haut, color=C_BAND, alpha=0.30, lw=0,
                                  label="±2σ")
        # bords fins : le bandeau reste repérable même quand il est étroit
        self.ax_conv.plot(xs, haut, color=C_BAND, lw=0.5, alpha=0.9)
        self.ax_conv.plot(xs, bas, color=C_BAND, lw=0.5, alpha=0.9)
        self.ax_conv.plot(xs, moy, color=C_LINE, lw=0.9, label="mean")

        # Échelle robuste : un artefact résiduel ne doit pas écraser le spectre
        vals = np.concatenate([bas, haut])
        vals = vals[np.isfinite(vals)]
        if vals.size:
            lo, hi = np.percentile(vals, 0.2), np.percentile(vals, 99.8)
            marge = 0.08 * max(hi - lo, 1.0)
            self.ax_conv.set_ylim(lo - marge, hi + marge)

        der = r["der_cumul"][-1] if r["der_cumul"] else np.nan
        der_txt = f"DER_SNR {der:.0f}" if np.isfinite(der) else "DER_SNR —"
        buil_txt = f"Buil {r['buil_pile']:.0f}" if r["buil_pile"] else "Buil —"
        self.ax_conv.set_title(
            "Signal convergence\n"
            f"{der_txt} · {buil_txt}",
            color=FG, fontsize=9, pad=4)
        self.ax_conv.set_xlabel("X pixel", color=FG, fontsize=8)
        self.ax_conv.set_ylabel("Summed ADU", color=FG, fontsize=8)
        self.ax_conv.set_xlim(0, moy.size)
        leg = self.ax_conv.legend(loc="upper right", fontsize=7.5,
                                  facecolor=BTN_BG, edgecolor=C_AXIS)
        for t in leg.get_texts():
            t.set_color(FG)

        self.canvas.draw()

        # bandeau de statut
        lignes = [ligne_resume(r)]
        detail = []
        if r["der_pose"]:
            detail.append(f"DER per frame {r['der_pose']:.0f}")
        if r["buil_pose"]:
            detail.append(f"Buil per frame {r['buil_pose']:.0f}")
        if r["n_rejets"]:
            detail.append(f"{r['n_rejets']} outlier pixels rejected "
                          f"(hot / cosmic rays)")
        if detail:
            lignes.append("   ".join(detail))
        av = list(r["zones"]["avert"])
        if r["melange"]:
            av.append("multiple OBJECT values in the selection")
        if av:
            lignes.append("⚠ " + " ; ".join(av))
        self.var_statut.set("\n".join(lignes))

    def _logger(self):
        """Envoie une ligne de synthèse au journal de l'Observation Pad.
        Sans rappel (usage autonome), retombe sur le presse-papiers."""
        if self.resultat is None:
            messagebox.showinfo("Log", "Nothing to log yet — load frames first.")
            return
        r = self.resultat
        z = r["zones"]
        tu = datetime.datetime.now(datetime.timezone.utc).strftime("%H:%M")
        seg = [f"{r['objet'] or '?'} ({r['n']} frames)",
               f"centre {z['centre']:.1f} px",
               f"FWHM {z['fwhm']:.1f} px",
               f"tilt {np.degrees(np.arctan(r['pente'])):.2f}°",
               f"width {z['tgt'][1] - z['tgt'][0]} px"]
        der = r["der_cumul"][-1] if r["der_cumul"] else np.nan
        seg.append(f"DER_SNR {der:.0f}" if np.isfinite(der) else "DER_SNR —")
        seg.append(f"Buil {r['buil_pile']:.0f}" if r["buil_pile"] else "Buil —")
        ligne = f"[{tu} TU] SNR check : " + " / ".join(seg)

        if callable(self.on_log):
            self.on_log(ligne)
            self.var_statut.set("Written to the observation log.")
        else:
            try:
                self.clipboard_clear()
                self.clipboard_append(ligne)
                self.var_statut.set("No log available — copied to clipboard.")
            except tk.TclError:
                self.var_statut.set(ligne)


def main():
    # Usage autonome : racine masquée, la fenêtre visible reste le Toplevel.
    args = [a for a in sys.argv[1:] if os.path.isfile(a)]
    root = tk.Tk()
    root.withdraw()
    app = Application(root, args or None)
    app.protocol("WM_DELETE_WINDOW", root.destroy)
    root.mainloop()


if __name__ == "__main__":
    main()
