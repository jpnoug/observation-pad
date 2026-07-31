#!/usr/bin/env python3
"""
Observation Pad — notepad for spectroscopy observing sessions.

Free-form notepad: the main text area is written however you like. Buttons
INSERT text at the cursor to save typing in the dark.

Buttons:
  Files           — choose the acquisition folder (where the log is written)
  New obs         — recap of FITS acquired since the last New obs (count +
                    names), then separator and new "new obs" / target block
  Insert TU       — inserts the UTC timestamp: [HH:MM TU]
  Check last fits — reads the last FITS: target, type/exp, RA/DEC, altitude,
                    airmass, CCD temperature, pier side, header time
  MTO             — Moon/Sun altitude + illumination and current Open-Meteo
                    weather (T, T-Td, gusts, cirrus). Live and fast.
  AOD             — aerosol optical depth from two sources: CAMS model
                    (always available) and AERONET ground measurement from
                    the nearest photometer. The AERONET search is anchored
                    on the last frame's DATE-OBS, so it still works when the
                    pad is relaunched after the run. Slower — click once or
                    twice a night.
  Calib           — one line: bias (0 sec) dark (600 sec) neon (13 sec) ...
  Airmass         — popup: airmass summary (Light frames) grouped by target,
                    with per-target average — check reference star placement
  SNR             — launches spectro_convergence.py (separate window) on the
                    current target's frames: stacking convergence and SNR
  Save            — saves the log
  End             — (confirmation) recaps all FITS of the night, grouped by
                    target with type and exposure time

Plain-text log obs_YYYY-MM-DD.txt, written to the acquisition folder, using
an "observing night" rule (before local noon = previous day).

EVENT-DRIVEN saving: written on every button click + on close.
No timer, no periodic CPU wake-up.

To be installed in the star_finder_gui venv (astropy is already present there).
"""

import os
import sys
import re
import math
import datetime
import subprocess
import configparser
import json
import urllib.request
import urllib.error
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# --- Lecture FITS / calculs astro (astropy présent dans le venv star_finder) --
try:
    from astropy.io import fits
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle, get_body
    from astropy.time import Time
    import astropy.units as u
    HAVE_ASTROPY = True
except ImportError:
    HAVE_ASTROPY = False

# --- Réglages par défaut (éditables ici, et à la main dans l'en-tête) ---------
SITE_DEFAUT = "PBO"
SETUP_DEFAUT = "Newton 200 / AM5 / Alpy 600 / Atik 414ex / Guide 174MM"
TEMP_CCD_DEFAUT = "-10"

# Coordonnées du site par défaut (PBO : Mairie de Pibrac, Haute-Garonne).
#   43°37'06.73"N  1°17'01.17"E
# Elles ne servent que de repli : si un FITS est disponible, ses SITELAT/SITELONG
# font foi (utile en nomade, où le site réel diffère du site par défaut).
SITE_LAT = 43.618536   # degrés, + Nord
SITE_LON = 1.283658    # degrés, + Est
SITE_ELEV = 190        # mètres (approx. ; sans effet notable sur Moon/Sol)

# Thème sombre/rouge pour préserver la vision nocturne
BG = "#000000"
FG = "#ff3333"          # rouge vif, bon contraste sur fond noir
FG_DIM = "#c0392b"      # rouge moyen (labels, statut, chemin par défaut)
BTN_BG = "#1a0a08"
FONT_TXT = ("DejaVu Sans Mono", 10)
FONT_BTN = ("DejaVu Sans", 9)

FITS_EXTS = (".fits", ".fit", ".fts")


# --- Utilitaires date/heure ---------------------------------------------------
def now_utc_hm():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%H:%M")


def observing_date():
    """Date de la nuit d'observation : avant midi local -> jour précédent."""
    now = datetime.datetime.now()
    if now.hour < 12:
        now = now - datetime.timedelta(days=1)
    return now.date()


def open_in_file_manager(path):
    """Ouvre un dossier avec le gestionnaire de fichiers système (non bloquant)."""
    try:
        if sys.platform.startswith("linux"):
            subprocess.Popen(["xdg-open", path])
        elif sys.platform == "darwin":
            subprocess.Popen(["open", path])
        elif os.name == "nt":
            os.startfile(path)  # noqa: S606
        else:
            return False
        return True
    except Exception:
        return False


# --- Utilitaires FITS ---------------------------------------------------------
def _activity_time(path):
    # Activité la plus récente du fichier. Sous Windows, une copie conserve le
    # mtime de la source mais fixe le ctime (création) à l'instant de la copie ;
    # prendre le max des deux garantit que le dernier fichier ajouté au dossier
    # est bien détecté comme le plus récent, quelle que soit sa date d'origine.
    try:
        return max(os.path.getmtime(path), os.path.getctime(path))
    except OSError:
        return 0.0


def list_fits(acq_dir):
    out = [os.path.join(acq_dir, n) for n in os.listdir(acq_dir)
           if n.lower().endswith(FITS_EXTS)]
    return sorted(out, key=_activity_time)


def _recency_key(path):
    # Pour départager un lot de fichiers arrivés ensemble : DATE-OBS (heure
    # d'observation, immuable dans le header) en priorité, activité fichier sinon.
    h = read_header(path)
    if h is not None and h.get("DATE-OBS"):
        try:
            return Time(str(h["DATE-OBS"]), format="isot", scale="utc").unix
        except Exception:
            pass
    return _activity_time(path)


def read_header(path):
    if not HAVE_ASTROPY:
        return None
    try:
        with fits.open(path, memmap=False) as hdul:
            return hdul[0].header
    except Exception:
        return None


def default_browse_dir():
    """Dossier proposé à l'ouverture du sélecteur : Documents s'il existe
    (Linux avec xdg-user-dirs, et Windows), sinon le dossier personnel."""
    home = os.path.expanduser("~")
    for name in ("Documents", "Documenten", "Dokumente"):
        d = os.path.join(home, name)
        if os.path.isdir(d):
            return d
    return home


def airmass_ky(alt_deg):
    """Airmass de Kasten-Young (1989) : robuste jusqu'aux basses hauteurs."""
    if alt_deg is None or alt_deg <= 0:
        return None
    a = math.radians(alt_deg)
    return 1.0 / (math.sin(a) + 0.50572 * (alt_deg + 6.07995) ** -1.6364)


def default_location():
    """Position du site par défaut (constantes en haut du script)."""
    if not HAVE_ASTROPY or SITE_LAT is None or SITE_LON is None:
        return None
    try:
        return EarthLocation(lat=SITE_LAT * u.deg, lon=SITE_LON * u.deg,
                             height=(SITE_ELEV or 0) * u.m)
    except Exception:
        return None


# Nom du site lu dans observer.ini, pour l'écrire dans l'en-tête / le statut.
# Rempli par location_from_ini(), None si le fichier est absent ou muet.
INI_SITE_NAME = None


def _find_observer_ini():
    """Cherche observer.ini à côté du script (dossier SpectroStarFinder)."""
    here = os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(here, "observer.ini")
    return p if os.path.isfile(p) else None


def location_from_ini():
    """Position lue dans observer.ini (section [observer]), source partagée
    avec SpectroStarFinder. Renvoie None si absent/illisible. Met à jour
    INI_SITE_NAME au passage."""
    global INI_SITE_NAME
    if not HAVE_ASTROPY:
        return None
    path = _find_observer_ini()
    if path is None:
        return None
    try:
        cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
        cfg.read(path, encoding="utf-8")
        sec = cfg["observer"]
        lat = float(sec["latitude"])
        lon = float(sec["longitude"])
        elev = float(sec.get("elevation", "0") or 0)
        loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=elev * u.m)
        INI_SITE_NAME = (sec.get("name", "") or "").strip() or None
        return loc
    except Exception:
        return None


def location_from_header(h):
    """Position écrite dans le header (INDI : SITELAT/SITELONG), sinon None."""
    if not HAVE_ASTROPY or h is None:
        return None
    if "SITELAT" in h and "SITELONG" in h:
        try:
            return EarthLocation(
                lat=Angle(str(h["SITELAT"]), u.deg),
                lon=Angle(str(h["SITELONG"]), u.deg),
                height=float(h.get("SITEELEV", 0)) * u.m,
            )
        except Exception:
            return None
    return None


def moon_illumination():
    """Fraction éclairée du disque lunaire (0–1) à l'instant courant.
    Indépendante du site (variation topocentrique négligeable). Formule de
    Meeus : fraction = (1 + cos i)/2, i = angle de phase."""
    if not HAVE_ASTROPY:
        return None
    try:
        t = Time.now()
        sun = get_body("sun", t)
        moon = get_body("moon", t)
        elong = sun.separation(moon)          # élongation géocentrique
        # angle de phase i : atan2 tenant compte des distances Terre-Soleil/Lune
        i = math.atan2(
            sun.distance.to("km").value * math.sin(elong.radian),
            moon.distance.to("km").value
            - sun.distance.to("km").value * math.cos(elong.radian),
        )
        return (1 + math.cos(i)) / 2
    except Exception:
        return None


def fetch_weather(loc, timeout=3):
    """Interroge Open-Meteo pour le site `loc`. Renvoie un dict
    {temperature, t_td, wind_gust, cloud_high} ou None (réseau absent, etc.).
    Appel synchrone bref : la fenêtre peut se figer jusqu'à `timeout` s."""
    if loc is None:
        return None
    lat = loc.lat.deg
    lon = loc.lon.deg
    url = (
        "https://api.open-meteo.com/v1/forecast"
        f"?latitude={lat:.4f}&longitude={lon:.4f}"
        "&current=temperature_2m,dew_point_2m,wind_gusts_10m,cloud_cover_high"
        "&timezone=UTC"
    )
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        cur = data.get("current", {})
        temp = cur.get("temperature_2m")
        dew = cur.get("dew_point_2m")
        return {
            "temperature": temp,
            "t_td": (temp - dew) if (temp is not None and dew is not None) else None,
            "wind_gust": cur.get("wind_gusts_10m"),
            "cloud_high": cur.get("cloud_cover_high"),
        }
    except (urllib.error.URLError, OSError, ValueError, KeyError):
        return None


def fetch_aod(loc, when=None, timeout=3):
    """Épaisseur optique des aérosols à 550 nm (AOD), via l'API Air Quality
    d'Open-Meteo, alimentée par CAMS (11 km sur l'Europe, 45 km ailleurs).
    Utile pour l'extinction atmosphérique lors du traitement des spectres.

    L'AOD n'est exposée qu'en série HORAIRE. On retient l'échéance la plus
    proche de `when` — la DATE-OBS des poses en dépouillement a posteriori,
    l'instant courant sinon. `past_days` permet de remonter dans le passé,
    sans quoi seule la journée en cours serait couverte.

    Renvoie (valeur, horodatage) ou (None, None).
    """
    if loc is None:
        return None, None
    lat = loc.lat.deg
    lon = loc.lon.deg
    ref = when or datetime.datetime.now(datetime.timezone.utc)
    today = datetime.datetime.now(datetime.timezone.utc).date()
    # Nombre de jours à remonter pour couvrir `ref` (l'API accepte 92 max)
    past = min(92, max(1, (today - ref.date()).days + 1))
    url = (
        "https://air-quality-api.open-meteo.com/v1/air-quality"
        f"?latitude={lat:.4f}&longitude={lon:.4f}"
        "&hourly=aerosol_optical_depth"
        f"&past_days={past}&forecast_days=1&timezone=UTC"
    )
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        hourly = data.get("hourly", {})
        times = hourly.get("time") or []
        vals = hourly.get("aerosol_optical_depth") or []
        if not times or not vals:
            return None, None
        # Échéance horaire la plus proche de la référence
        best_v, best_t, best_gap = None, None, None
        for ts, v in zip(times, vals):
            if v is None:
                continue
            try:
                t = datetime.datetime.fromisoformat(ts).replace(
                    tzinfo=datetime.timezone.utc)
            except ValueError:
                continue
            gap = abs((t - ref).total_seconds())
            if best_gap is None or gap < best_gap:
                best_v, best_t, best_gap = float(v), t, gap
        return best_v, best_t
    except (urllib.error.URLError, OSError, ValueError, KeyError, TypeError):
        return None, None


def _haversine_km(lat1, lon1, lat2, lon2):
    """Distance orthodromique approchée, en km."""
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = (math.sin(dp / 2) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2)
    return 2 * r * math.asin(math.sqrt(a))


def fetch_aeronet(loc, when=None, timeout=6, box_deg=1.0, days=7):
    """AOD MESURÉE au sol par le réseau AERONET (photomètres CIMEL).

    Interroge le service web NASA sur une boîte englobante autour du site,
    retient la station la plus proche puis sa mesure la plus récente.
    `lunar_merge=1` inclut la photométrie LUNAIRE : des mesures de nuit sont
    possibles quand la Lune est assez brillante (pas près de la nouvelle Lune).

    Fenêtre de `days` jours : le niveau 1.5 a une latence de publication, et
    une nuit couverte ne produit aucune donnée. On préfère une mesure un peu
    ancienne — dont l'horodatage est reporté dans le log, à toi de juger — à
    aucune mesure du tout.

    L'AOD est ramenée à 550 nm via le coefficient d'Ångström, pour être
    directement comparable à la valeur CAMS :  AOD(550) = AOD(l) (550/l)^-alpha

    Renvoie (record | None, url). L'URL est toujours renvoyée pour pouvoir
    l'écrire dans le journal quand aucune donnée n'est trouvée.
    """
    if loc is None:
        return None, None
    lat = loc.lat.deg
    lon = loc.lon.deg
    today = datetime.datetime.now(datetime.timezone.utc).date()
    if when is None:
        # Session en direct : on remonte `days` jours depuis aujourd'hui.
        end = today
        start = today - datetime.timedelta(days=days)
    else:
        # Dépouillement a posteriori : fenêtre CENTRÉE sur la date des poses,
        # pour trouver la mesure contemporaine de l'observation.
        half = max(1, days // 2)
        start = when.date() - datetime.timedelta(days=half)
        end = min(today, when.date() + datetime.timedelta(days=half))
    url = (
        "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3"
        f"?lat1={lat - box_deg:.4f}&lon1={lon - box_deg:.4f}"
        f"&lat2={lat + box_deg:.4f}&lon2={lon + box_deg:.4f}"
        f"&year={start.year}&month={start.month}&day={start.day}"
        f"&year2={end.year}&month2={end.month}&day2={end.day}"
        "&AOD15=1&AVG=10&lunar_merge=1&if_no_html=1"
    )
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            text = resp.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, OSError):
        return None, url

    lines = [l for l in text.splitlines() if l.strip()]
    # En-tête = la ligne qui nomme les colonnes. On repère n'importe quelle
    # colonne "AOD_<nnn>nm" : selon la nuit et l'instrument, toutes les
    # longueurs d'onde ne sont pas restituées (beaucoup de -999).
    rx_aod = re.compile(r"^AOD_(\d+)nm$")
    hdr_i = next((i for i, l in enumerate(lines)
                  if any(rx_aod.match(c.strip().replace(" ", ""))
                         for c in l.split(","))), None)
    if hdr_i is None:
        return None, url
    cols = [c.strip().replace(" ", "") for c in lines[hdr_i].split(",")]

    def col(pred):
        return next((i for i, c in enumerate(cols) if pred(c)), None)

    # toutes les colonnes AOD disponibles, avec leur longueur d'onde
    aod_cols = []
    for i, c in enumerate(cols):
        m = rx_aod.match(c)
        if m:
            aod_cols.append((i, float(m.group(1))))
    i_500 = col(lambda c: c == "AOD_500nm")
    i_ang = col(lambda c: "Angstrom_Exponent" in c and "440-870" in c)
    i_lat = col(lambda c: c.startswith("Site_Latitude"))
    i_lon = col(lambda c: c.startswith("Site_Longitude"))
    i_elev = col(lambda c: c.startswith("Site_Elevation"))
    i_site = col(lambda c: "AERONET_Site" in c or c == "Site")
    i_date = col(lambda c: c.startswith("Date("))
    i_time = col(lambda c: c.startswith("Time("))
    if not aod_cols or i_lat is None or i_lon is None:
        return None, url

    def num(row, i):
        if i is None or i >= len(row):
            return None
        try:
            v = float(row[i])
        except ValueError:
            return None
        return None if v <= -900 else v   # -999 = valeur absente

    ref = when or datetime.datetime.now(datetime.timezone.utc)
    best = None
    for line in lines[hdr_i + 1:]:
        row = [c.strip() for c in line.split(",")]
        slat, slon = num(row, i_lat), num(row, i_lon)
        if slat is None or slon is None:
            continue
        # canal utilisable le plus proche de 550 nm
        avail = [(abs(w - 550.0), w, num(row, i)) for i, w in aod_cols]
        avail = [(d, w, v) for d, w, v in avail if v is not None]
        if not avail:
            continue
        _, lam, aod = min(avail)
        alpha = num(row, i_ang)
        # extrapolation d'Ångström vers 550 nm (identité si lam == 550)
        aod550 = aod * (550.0 / lam) ** (-alpha) if alpha is not None else None
        d = _haversine_km(lat, lon, slat, slon)
        when = None
        if i_date is not None and i_time is not None and i_time < len(row):
            try:
                dd, mm, yy = row[i_date].split(":")
                when = datetime.datetime(
                    int(yy), int(mm), int(dd),
                    *[int(x) for x in row[i_time].split(":")],
                    tzinfo=datetime.timezone.utc)
            except (ValueError, IndexError):
                when = None
        kind = next((v for v in row if v.lower() in ("lunar", "solar")), None)
        rec = {
            "site": row[i_site] if (i_site is not None and i_site < len(row))
                    else "?",
            "dist_km": d, "aod": aod, "lam": lam, "aod550": aod550,
            "aod500": num(row, i_500), "elev": num(row, i_elev),
            "angstrom": alpha, "kind": kind, "when": when,
        }
        # Station la plus proche ; à station égale, la mesure la plus proche
        # TEMPORELLEMENT de `ref` — soit l'instant des poses (dépouillement
        # a posteriori), soit maintenant (session en direct).
        def _gap(r):
            return abs((r["when"] - ref).total_seconds()) if r["when"] else 9e18
        key_new = (rec["dist_km"], _gap(rec))
        key_old = (best["dist_km"], _gap(best)) if best else None
        if best is None or key_new < key_old:
            best = rec
    return best, url


def moon_sun_alt(loc=None):
    """Hauteurs Lune et Soleil à l'heure TU courante, depuis `loc`.
    Si `loc` est None, on retombe sur le site par défaut."""
    if not HAVE_ASTROPY:
        return None, None
    if loc is None:
        loc = default_location()
    if loc is None:
        return None, None
    try:
        t = Time.now()
        altaz = AltAz(obstime=t, location=loc)
        moon_alt = float(get_body("moon", t, loc).transform_to(altaz).alt.deg)
        sun_alt = float(get_body("sun", t, loc).transform_to(altaz).alt.deg)
        return moon_alt, sun_alt
    except Exception:
        return None, None


def _sky_from_header(h):
    if "OBJCTRA" in h and "OBJCTDEC" in h:
        try:
            return SkyCoord(str(h["OBJCTRA"]), str(h["OBJCTDEC"]),
                            unit=(u.hourangle, u.deg))
        except Exception:
            pass
    if "RA" in h and "DEC" in h:
        try:
            return SkyCoord(float(h["RA"]) * u.deg, float(h["DEC"]) * u.deg)
        except Exception:
            pass
    return None


def _location(h):
    # Priorité : header FITS (site réel écrit par la monture) > observer.ini
    # (réglé pour SpectroStarFinder) > constantes en dur (PBO). Le header
    # domine car il ne peut pas être "oublié" ; l'ini couvre le début de nuit
    # avant toute pose ; les constantes ne servent que si l'ini est absent.
    return location_from_header(h) or location_from_ini() or default_location()


def target_info(path):
    """Infos du dernier FITS : nom, RA/DEC, hauteur, airmass (header ou calcul)."""
    info = {"file": os.path.basename(path), "name": None, "dateobs": None,
            "ra": None, "dec": None, "alt": None, "airmass": None,
            "alt_src": None, "airmass_src": None,
            "imgtyp": None, "exptime": None, "ccdtemp": None, "pier": None}
    h = read_header(path)
    if h is None:
        info["name"] = os.path.splitext(info["file"])[0]
        return info

    info["name"] = str(h.get("OBJECT", "")).strip() or \
        os.path.splitext(info["file"])[0]

    # Type de pose (Light/Dark/Flat/Calib) et temps de pose
    info["imgtyp"] = str(h.get("IMAGETYP", h.get("FRAME", ""))).strip() or None
    try:
        info["exptime"] = float(h.get("EXPTIME", h.get("EXPOSURE")))
    except (TypeError, ValueError):
        pass
    # Température CCD réelle : contrôle que le refroidissement tient la consigne
    try:
        info["ccdtemp"] = float(h["CCD-TEMP"])
    except (KeyError, TypeError, ValueError):
        pass
    # Côté du pied : signale un retournement au méridien
    info["pier"] = str(h.get("PIERSIDE", "")).strip() or None

    # DATE-OBS : on garde la valeur BRUTE du header (début de pose, précision
    # milliseconde) — pas de troncature, l'étiquette correspond au contenu.
    dateobs = h.get("DATE-OBS")
    if dateobs:
        info["dateobs"] = str(dateobs).strip()

    sky = _sky_from_header(h)
    if sky is not None:
        # precision=0 : OBJCTRA/OBJCTDEC ne portent que des secondes entières
        # ('10 38 58'), inutile d'ajouter une décimale qui n'existe pas.
        info["ra"] = sky.ra.to_string(unit=u.hour, sep=":", precision=0, pad=True)
        info["dec"] = sky.dec.to_string(unit=u.deg, sep=":", precision=0,
                                        alwayssign=True, pad=True)

    # Hauteur / airmass : priorité au header, sinon calcul. On mémorise
    # l'origine pour ne pas étiqueter OBJCTALT/AIRMASS une valeur calculée.
    if "OBJCTALT" in h:
        try:
            info["alt"] = float(h["OBJCTALT"])
            info["alt_src"] = "header"
        except Exception:
            pass
    if "AIRMASS" in h:
        try:
            info["airmass"] = float(h["AIRMASS"])
            info["airmass_src"] = "header"
        except Exception:
            pass

    if info["alt"] is None and sky is not None:
        loc = _location(h)
        dateobs = h.get("DATE-OBS")
        if loc is not None and dateobs:
            try:
                t = Time(str(dateobs), format="isot", scale="utc")
                altaz = sky.transform_to(AltAz(obstime=t, location=loc))
                info["alt"] = float(altaz.alt.deg)
                info["alt_src"] = "calc"
            except Exception:
                pass

    if info["airmass"] is None and info["alt"] is not None:
        info["airmass"] = airmass_ky(info["alt"])
        info["airmass_src"] = "calc"

    return info


# --- Application --------------------------------------------------------------
class ObservationPad:
    def __init__(self, root):
        self.root = root
        self.acq_dir = None
        self.save_path = None
        self.recapped = set()   # fichiers FITS déjà inclus dans un récap (New)
        self.seen = set()       # fichiers FITS déjà observés par le programme
        self.latest_fits = None  # dernier FITS arrivé, suivi par observation

        root.title("Observation Pad")
        root.configure(bg=BG)

        self._init_style()
        self._build_toolbar()
        self._build_pathbar()
        self._build_status()
        self._build_text()

        self._insert_header()
        self._fit_to_buttons()
        root.protocol("WM_DELETE_WINDOW", self._on_close)

    def _fit_to_buttons(self):
        # Largeur par défaut = juste de quoi afficher la plus large des deux
        # rangées de boutons. Fenêtre volontairement compacte : elle cohabite
        # avec CCDciel et star_finder sur le PC d'acquisition.
        self.root.update_idletasks()
        w = max(r.winfo_reqwidth() for r in self.toolbar_rows) + 16
        self.root.geometry(f"{w}x680")
        self.root.minsize(w, 320)

    # -- Interface -------------------------------------------------------------
    def _init_style(self):
        # Thème 'clam' : entièrement redessinable -> scrollbar sombre même
        # sous Windows (le tk.Scrollbar natif ignorait les couleurs et
        # apparaissait blanc/éblouissant).
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        style.configure(
            "Dark.Vertical.TScrollbar",
            background=BTN_BG, troughcolor=BG, bordercolor=BG,
            arrowcolor=FG_DIM, darkcolor=BTN_BG, lightcolor=BTN_BG,
            relief="flat",
        )
        style.map(
            "Dark.Vertical.TScrollbar",
            background=[("active", FG_DIM)],
            arrowcolor=[("active", FG)],
        )

    def _build_toolbar(self):
        # Deux rangées : écriture dans le journal en haut, contrôles et
        # session en bas. Garde la fenêtre étroite sur le PC d'acquisition.
        holder = tk.Frame(self.root, bg=BG)
        holder.pack(side=tk.TOP, fill=tk.X, padx=6, pady=(6, 4))
        row1 = tk.Frame(holder, bg=BG)
        row1.pack(side=tk.TOP, fill=tk.X)
        row2 = tk.Frame(holder, bg=BG)
        row2.pack(side=tk.TOP, fill=tk.X, pady=(4, 0))
        self.toolbar_rows = (row1, row2)

        rows = [
            # rangée 1 — ce qui écrit dans le journal
            (row1, [
                ("Files",     self.choose_folder),
                ("New obs",   self.btn_block),
                ("Insert TU", self.btn_time),
                ("MTO",       self.btn_mto),
                ("AOD",       self.btn_aod),
                ("Calib",     self.btn_calib),
            ]),
            # rangée 2 — contrôles et fin de session
            (row2, [
                ("Check last fits", self.btn_fits),
                ("Airmass",         self.btn_airmass),
                ("SNR",             self.btn_convergence),
                ("Save",            self.save),
                ("End",             self.btn_end),
            ]),
        ]
        for parent, buttons in rows:
            for label, cmd in buttons:
                b = tk.Button(
                    parent, text=label, command=self._wrap(cmd),
                    bg=BTN_BG, fg=FG, activebackground=FG_DIM,
                    activeforeground=BG,
                    font=FONT_BTN, relief=tk.FLAT, padx=10, pady=4,
                    highlightbackground=FG_DIM, highlightthickness=1,
                )
                b.pack(side=tk.LEFT, padx=3)

    def _wrap(self, cmd):
        # Efface le statut avant chaque clic : un message ("FITS lu : ...",
        # "Enregistré : ...") ne reste affiché que jusqu'au bouton suivant,
        # au lieu de traîner indéfiniment en bas de fenêtre.
        def handler():
            self._set_status("")
            cmd()
        return handler

    def _build_text(self):
        frame = tk.Frame(self.root, bg=BG)
        frame.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 6))
        self.text = tk.Text(
            frame, bg=BG, fg=FG, insertbackground=FG,
            font=FONT_TXT, wrap=tk.WORD, undo=False,
            selectbackground=FG_DIM, selectforeground=BG,
            padx=8, pady=8, relief=tk.FLAT,
            # Liseré rouge tout autour, pour détacher la zone de texte du fond
            # noir de la fenêtre. Même couleur avec/sans focus (pas de saut
            # visuel quand on clique dedans).
            highlightthickness=2, highlightbackground=FG_DIM, highlightcolor=FG,
            bd=0,
        )
        # undo=False désactive déjà Ctrl-Z, mais certains bindings système
        # (<<Undo>>) peuvent persister selon la plateforme : on les neutralise
        # explicitement pour ne jamais risquer d'effacer une partie du journal
        # par un Ctrl-Z accidentel en pleine nuit.
        for seq in ("<Control-z>", "<Control-Z>", "<<Undo>>",
                    "<Control-y>", "<Control-Y>", "<<Redo>>"):
            self.text.bind(seq, lambda e: "break")
        scroll = ttk.Scrollbar(frame, orient="vertical",
                               command=self.text.yview,
                               style="Dark.Vertical.TScrollbar")
        self.text.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.text.focus_set()

    def _build_status(self):
        self.status = tk.Label(
            self.root, text="No acquisition folder selected.",
            bg=BG, fg=FG_DIM, font=FONT_BTN, anchor="w", padx=8,
        )
        self.status.pack(side=tk.BOTTOM, fill=tk.X, pady=(2, 4))

    def _set_status(self, msg):
        self.status.configure(text=msg)

    def _build_pathbar(self):
        # Bottom bar: small open button + acquisition folder path
        bar = tk.Frame(self.root, bg=BG)
        bar.pack(side=tk.BOTTOM, fill=tk.X, padx=6, pady=(0, 6))
        self.open_btn = tk.Button(
            bar, text="Open", command=self._wrap(self.open_folder),
            bg=BTN_BG, fg=FG, activebackground=FG_DIM, activeforeground=BG,
            font=FONT_BTN, relief=tk.FLAT, padx=8, pady=2,
            highlightbackground=FG_DIM, highlightthickness=1,
        )
        self.open_btn.pack(side=tk.LEFT, padx=(0, 6))
        self.path_label = tk.Label(
            bar, text="(no folder selected)",
            bg=BG, fg=FG_DIM, font=FONT_BTN, anchor="w",
        )
        self.path_label.pack(side=tk.LEFT, fill=tk.X, expand=True)

    def open_folder(self):
        if not self.acq_dir:
            self._set_status("No folder to open (click Files first).")
            return
        if not open_in_file_manager(self.acq_dir):
            self._set_status("Could not open the folder on this system.")

    # -- Pre-filled header ------------------------------------------------------
    def _insert_header(self):
        d = observing_date().isoformat()
        header = (
            f"Date : {d}\n"
            f"Site : {SITE_DEFAUT}\n"
            f"Setup : {SETUP_DEFAUT}\n"
            f"CCD Temperature : {TEMP_CCD_DEFAUT}\n"
            f"Polar Align : \n"
            f"{'=' * 23}\n"
        )
        self.text.insert("1.0", header)
        # Cursor right after "Polar Align : " (line 5): field to fill in
        # right away when the pad opens, before the first New obs.
        self.text.mark_set(tk.INSERT, "5.end")

    # -- Insertion au curseur --------------------------------------------------
    def _insert(self, s):
        self.text.insert(tk.INSERT, s)
        self.text.see(tk.INSERT)
        self.text.focus_set()
        self._quiet_save()

    def _insert_line(self, s):
        line_start = self.text.index("insert linestart")
        prefix = "" if self.text.get(line_start, "insert") == "" else "\n"
        self._insert(prefix + s)

    def _require_dir(self):
        if self.acq_dir:
            return True
        self._set_status("Choose the acquisition folder first (Files).")
        return False

    # -- Boutons simples -------------------------------------------------------
    def btn_time(self):
        self._insert_line(f"[{now_utc_hm()} TU] ")

    def _site_location(self):
        """Position de calcul et libellé de sa source. Priorité : header du
        dernier FITS > observer.ini > constantes en dur (PBO)."""
        if self.acq_dir:
            self._observe()
            if self.latest_fits:
                loc = location_from_header(read_header(self.latest_fits))
                if loc is not None:
                    return loc, "FITS"
        loc = location_from_ini()
        if loc is not None:
            return loc, (INI_SITE_NAME or "ini")
        return default_location(), SITE_DEFAUT

    def btn_block(self):
        # New obs closes the previous target with a recap of the FITS
        # acquired since the last New obs, then opens a new target block.
        self._observe()   # keeps the seen-files tracking up to date (for Fits)
        body = []
        tu = now_utc_hm()
        sep = "-" * 40
        if self.acq_dir:
            new = [f for f in list_fits(self.acq_dir) if f not in self.recapped]
            if new:
                names = [os.path.basename(f) for f in new]
                # "new" : these are files that appeared since the last New obs,
                # not the whole content of the folder.
                mot = "new files" if len(names) > 1 else "new file"
                # Aggregate a few header values over the new frames, to close
                # the target with its exposure budget and mean airmass.
                exps, ccds, ams = [], [], []
                for f in new:
                    h = read_header(f)
                    if h is None:
                        continue
                    try:
                        exps.append(float(h.get("EXPTIME", h.get("EXPOSURE"))))
                    except (TypeError, ValueError):
                        pass
                    try:
                        ccds.append(float(h["CCD-TEMP"]))
                    except (KeyError, TypeError, ValueError):
                        pass
                    try:
                        ams.append(float(h["AIRMASS"]))
                    except (KeyError, TypeError, ValueError):
                        pass
                extra = []
                if exps:
                    unit = exps[0] if len(set(round(e, 1) for e in exps)) == 1 \
                        else None
                    total = sum(exps)
                    if unit is not None:
                        extra.append(f"exp {int(round(unit))}s ×{len(exps)} "
                                     f"= {total / 60:.1f} min")
                    else:
                        extra.append(f"exp total {total / 60:.1f} min "
                                     f"({len(exps)} frames)")
                if ams:
                    extra.append(f"airmass {sum(ams) / len(ams):.2f} mean")
                if ccds:
                    extra.append(f"CCD {sum(ccds) / len(ccds):.1f}°C")
                tail = ("  |  " + " / ".join(extra)) if extra else ""
                body.append(f"[{tu} TU] Recap : {len(names)} {mot} ; "
                            + ", ".join(names) + tail)
                self.recapped.update(new)
                recap_msg = (f"Recap : {len(names)} {mot} "
                             f"({len(self.recapped)} total tonight)")
            else:
                recap_msg = ("No new FITS since the last New obs "
                             f"— {len(self.recapped)} already recapped")
        else:
            recap_msg = "No folder selected (Files): recap unavailable"
        body.append(sep)
        body.append("")
        body.append(f"[{tu} TU] new obs")
        self._set_status(recap_msg)
        body.append("*** Target: ")
        # leading "\n" = one blank line before the block (before the recap)
        self._insert_line("\n" + "\n".join(body))

    def btn_mto(self):
        # Moon/Sun altitude + illumination, and current Open-Meteo weather.
        # Deliberately LIVE and lightweight: a single HTTP request, so the
        # window is never frozen for long. Aerosols moved to the AOD button.
        tu = now_utc_hm()
        loc, src = self._site_location()
        lines = []
        moon_alt, sun_alt = moon_sun_alt(loc)
        illum = moon_illumination()
        if moon_alt is None:
            lines.append(f"[{tu} TU] Moon illum. ?  h ?  Sun h ?")
        else:
            illum_txt = f"{illum * 100:.0f}%" if illum is not None else "?"
            lines.append(f"[{tu} TU] Moon illum. {illum_txt} h {moon_alt:.0f}° "
                         f"Sun h {sun_alt:.0f}°")
        wx = fetch_weather(loc)
        if wx is not None:
            seg = []
            if wx["temperature"] is not None:
                seg.append(f"T {wx['temperature']:.0f}°C")
            if wx["t_td"] is not None:
                seg.append(f"T-Td {wx['t_td']:.1f}°C")
            if wx["wind_gust"] is not None:
                seg.append(f"gusts {wx['wind_gust']:.0f} km/h")
            if wx["cloud_high"] is not None:
                seg.append(f"cirrus {wx['cloud_high']:.0f}%")
            if seg:
                lines.append(f"[{tu} TU] open_meteo : " + " / ".join(seg))
        self._insert_line("\n".join(lines))
        if moon_alt is not None:
            self._set_status(f"MTO inserted  |  site {src}")

    def btn_aod(self):
        """Aerosol optical depth, from two independent sources.

        CAMS (Open-Meteo Air Quality) is a MODEL, always available.
        AERONET is a ground MEASUREMENT from the nearest CIMEL photometer,
        far more accurate but only when the station published data.

        Unlike MTO, the AERONET search is anchored on the DATE-OBS of the
        last frame when there is one — so relaunching the pad the day after
        an observing run still retrieves the measurement contemporaneous
        with the exposures.

        Two HTTP requests: slower than MTO, meant to be clicked once or
        twice per night rather than continuously.
        """
        tu = now_utc_hm()
        loc, _src = self._site_location()
        lines = []

        # Référence temporelle commune aux DEUX sources : la DATE-OBS du
        # dernier FITS s'il y en a un, sinon maintenant. Sans cela, CAMS
        # décrirait l'instant du clic et AERONET la nuit d'observation.
        ref_when = None
        if self.acq_dir and self.latest_fits:
            hh = read_header(self.latest_fits)
            if hh is not None and hh.get("DATE-OBS"):
                try:
                    ref_when = datetime.datetime.fromisoformat(
                        str(hh["DATE-OBS"]).strip()).replace(
                            tzinfo=datetime.timezone.utc)
                except ValueError:
                    ref_when = None
        base = ref_when or datetime.datetime.now(datetime.timezone.utc)
        tag = "vs frame" if ref_when else "vs now"

        # 1) CAMS model — always something, low accuracy (MAE ~0.08)
        aod, aod_t = fetch_aod(loc, when=ref_when)
        if aod is not None:
            note = "CAMS model"
            if aod_t is not None:
                gap = (aod_t - base).total_seconds() / 3600.0
                note += (f", {aod_t:%Y-%m-%d %H:%M} TU "
                         f"({gap:+.1f} h {tag})")
            lines.append(f"[{tu} TU] AOD 550nm : {aod:.3f} ({note})")

        # 2) AERONET ground measurement, same time reference
        ae, ae_url = fetch_aeronet(loc, when=ref_when)
        if ae is not None:
            head = f"{ae['site']} ({ae['dist_km']:.0f} km"
            if ae["elev"] is not None:
                head += f", {ae['elev']:.0f} m"
            head += ")"
            seg = [head]
            if ae["aod550"] is not None:
                seg.append(f"AOD550 {ae['aod550']:.3f} "
                           f"(from {ae['lam']:.0f}nm)")
            else:
                seg.append(f"AOD{ae['lam']:.0f} {ae['aod']:.3f}")
            if ae["aod500"] is not None:
                seg.append(f"AOD500 {ae['aod500']:.3f}")
            if ae["angstrom"] is not None:
                seg.append(f"Angstrom {ae['angstrom']:.2f}")
            if ae["kind"]:
                seg.append(ae["kind"])
            if ae["when"] is not None:
                gap = (ae["when"] - base).total_seconds() / 3600.0
                seg.append(ae["when"].strftime("%Y-%m-%d %H:%M TU")
                           + f" ({gap:+.1f} h {tag})")
            # A station much higher up does not "see" the aerosol layer below
            # it: its AOD then underestimates extinction at the observing
            # site. Flag it rather than let it pass unnoticed.
            if ae["elev"] is not None and loc is not None:
                dz = ae["elev"] - loc.height.to_value(u.m)
                if abs(dz) > 500:
                    seg.append(f"!! alt {dz:+.0f} m vs site")
            lines.append(f"[{tu} TU] AERONET : " + " / ".join(seg))
        elif ae_url:
            # No usable measurement in the window: log the query link so you
            # can check for yourself rather than being left with nothing.
            lines.append(f"[{tu} TU] AERONET : no data — {ae_url}")

        if lines:
            self._insert_line("\n".join(lines))
            self._set_status("AOD inserted"
                             + (" (AERONET anchored on frame time)"
                                if ref_when and ae is not None else ""))
        else:
            self._set_status("No AOD available (no network?).")

    def btn_calib(self):
        # A single line, with the usual exposure times in parentheses
        # (edit by hand if this session needs different values).
        tu = now_utc_hm()
        self._insert_line(
            f"[{tu} TU] bias (0 sec) dark (600 sec) neon (13 sec) flat (6 sec)"
        )

    def btn_airmass(self):
        # Airmass summary, in a POPUP (a quick check, not logged).
        # Light frames only (calibrations have no meaningful airmass),
        # grouped by target in acquisition order.
        if not self._require_dir():
            return
        files = list_fits(self.acq_dir)   # acquisition order (activity)
        order, groups = [], {}
        for f in files:
            h = read_header(f)
            if h is None:
                continue
            typ = str(h.get("IMAGETYP", h.get("FRAME", ""))).strip().lower()
            if typ and typ != "light":
                continue   # excludes bias / dark / flat / calibration
            am = h.get("AIRMASS")
            try:
                am = float(am)
            except (TypeError, ValueError):
                am = None
            obj = str(h.get("OBJECT", "")).strip() or \
                os.path.splitext(os.path.basename(f))[0]
            if obj not in groups:
                groups[obj] = []
                order.append(obj)
            groups[obj].append((os.path.basename(f), am))

        self._show_airmass_popup(order, groups)

    def _show_airmass_popup(self, order, groups):
        win = tk.Toplevel(self.root)
        win.title("Airmass Summary")
        win.configure(bg=BG)
        win.geometry("560x520")

        txt = tk.Text(win, bg=BG, fg=FG, insertbackground=FG,
                      font=FONT_TXT, wrap=tk.NONE, relief=tk.FLAT,
                      padx=10, pady=10,
                      highlightthickness=2, highlightbackground=FG_DIM,
                      highlightcolor=FG_DIM, bd=0)
        scroll = ttk.Scrollbar(win, orient="vertical", command=txt.yview,
                               style="Dark.Vertical.TScrollbar")
        txt.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        txt.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        if not order:
            txt.insert("1.0", "No Light frames with airmass in this folder.")
        else:
            out = []
            for obj in order:
                rows = groups[obj]
                vals = [am for _, am in rows if am is not None]
                avg = sum(vals) / len(vals) if vals else None
                head = f"{obj}"
                if avg is not None:
                    head += f"   (average {avg:.2f})"
                out.append(head)
                for name, am in rows:
                    am_txt = f"{am:.2f}" if am is not None else "  —"
                    out.append(f"    {name:<34s} {am_txt}")
                out.append("")
            txt.insert("1.0", "\n".join(out))
        txt.configure(state=tk.DISABLED)

        btn = tk.Button(win, text="Close", command=win.destroy,
                        bg=BTN_BG, fg=FG, activebackground=FG_DIM,
                        activeforeground=BG, font=FONT_BTN, relief=tk.FLAT,
                        padx=12, pady=4, highlightbackground=FG_DIM,
                        highlightthickness=1)
        btn.pack(side=tk.BOTTOM, pady=6)
        win.transient(self.root)

    def btn_convergence(self):
        # Spectro Convergence s'ouvre en fenêtre FILLE (Toplevel), ce qui
        # permet à son bouton « Log » d'écrire directement dans ce journal
        # via le rappel ci-dessous. L'import est PARESSEUX : numpy et
        # matplotlib ne sont chargés qu'au premier clic, le démarrage du pad
        # reste léger sur le PC de terrain.
        here = os.path.dirname(os.path.abspath(__file__))
        if here not in sys.path:
            sys.path.insert(0, here)
        try:
            from spectro_convergence import Application as ConvergenceWindow
        except ImportError as e:
            self._set_status(f"spectro_convergence not available: {e}")
            return
        # Série de la cible en cours : tous les Light partageant l'OBJECT du
        # dernier FITS. Le module garde son bouton « Files… » pour changer.
        frames = []
        if self.acq_dir:
            self._observe()
            if self.latest_fits:
                h = read_header(self.latest_fits)
                obj = str(h.get("OBJECT", "")).strip() if h is not None else ""
                if obj:
                    for f in list_fits(self.acq_dir):
                        hh = read_header(f)
                        if hh is None:
                            continue
                        typ = str(hh.get("IMAGETYP",
                                         hh.get("FRAME", ""))).strip().lower()
                        if typ and typ != "light":
                            continue
                        if str(hh.get("OBJECT", "")).strip() == obj:
                            frames.append(f)
        try:
            win = ConvergenceWindow(self.root, frames or None,
                                    on_log=self._log_from_tool)
            win.transient(self.root)
        except Exception as e:
            self._set_status(f"Could not open Spectro Convergence: {e}")
            return
        self._set_status("Spectro Convergence opened"
                         + (f" — {len(frames)} frames preloaded" if frames else
                            " — select frames in its window"))

    def _log_from_tool(self, line):
        """Rappel offert aux fenêtres outils : écrit une ligne dans le journal
        et la sauvegarde, comme si un bouton du pad l'avait insérée."""
        self._insert_line(line)
        self._set_status("Line added from Spectro Convergence.")

    def btn_end(self):
        # Asks for confirmation, then recaps ALL FITS of the night
        # (whole folder content, not just the new ones).
        if not messagebox.askyesno(
            "End session",
            "End the observing session?\n\n"
            "A recap of all FITS from the night will be inserted "
            "into the log.",
        ):
            return
        tu = now_utc_hm()
        lines = [f"[{tu} TU] End session"]
        if self.acq_dir:
            files = list_fits(self.acq_dir)
            if files:
                order, groups = [], {}
                for f in files:
                    h = read_header(f)
                    if h is not None:
                        obj = str(h.get("OBJECT", "")).strip() or \
                            os.path.splitext(os.path.basename(f))[0]
                        typ = str(h.get("IMAGETYP", h.get("FRAME", ""))).strip() \
                            or "?"
                        try:
                            e = float(h.get("EXPTIME", h.get("EXPOSURE")))
                            expo = f"{int(round(e))}s" if abs(e - round(e)) < 0.05 \
                                else f"{e:.1f}s"
                        except (TypeError, ValueError):
                            expo = "?s"
                    else:
                        obj = os.path.splitext(os.path.basename(f))[0]
                        typ, expo = "?", "?s"
                    key = (typ, expo)
                    if obj not in groups:
                        groups[obj] = {}
                        order.append(obj)
                    groups[obj][key] = groups[obj].get(key, 0) + 1
                lines.append(f"Night summary: {len(files)} files")
                for obj in order:
                    parts = []
                    for (typ, expo), n in groups[obj].items():
                        tag = f" ({typ})" if typ != "?" else ""
                        parts.append(f"{n}×{expo}{tag}")
                    lines.append(f"  {obj} — " + " + ".join(parts))
            else:
                lines.append("Night summary: no FITS file found")
        else:
            lines.append("Night summary: no folder selected (Files)")
        self._insert_line("\n" + "\n".join(lines))
        self._set_status("Session ended — recap inserted.")

    # -- Suivi du dernier FITS arrivé (indépendant des horodatages) -----------
    def _observe(self):
        """Compare le contenu du dossier à ce que le programme a déjà vu, et
        détermine le dernier FITS arrivé. Fiable même si la copie des fichiers
        préserve des dates trompeuses : c'est l'ordre d'apparition constaté par
        le programme (à chaque clic) qui fait foi, pas le mtime/ctime."""
        if not self.acq_dir:
            return
        current = list_fits(self.acq_dir)
        cur_set = set(current)
        new = [f for f in current if f not in self.seen]
        if new:
            # dernier arrivé du lot : DATE-OBS la plus récente, sinon activité
            self.latest_fits = max(new, key=_recency_key)
            self.seen.update(new)
        self.seen &= cur_set   # purge des fichiers supprimés du dossier
        if self.latest_fits not in cur_set:
            self.latest_fits = max(current, key=_recency_key) if current else None

    # -- Fits button: snapshot of the last FITS to arrive ----------------------
    def btn_fits(self):
        if not self._require_dir():
            return
        self._observe()
        if not self.latest_fits:
            self._set_status("No FITS file in the folder.")
            return
        info = target_info(self.latest_fits)
        # Assumed mix: DATE-OBS is reproduced RAW under its keyword name (the
        # label exactly matches the header content); other values are
        # formatted/rounded for readability, so labelled in plain words
        # rather than a keyword name that would be misleading.
        parts = []
        if info["name"]:
            parts.append(f"OBJECT : {info['name']}")
        if info["imgtyp"]:
            parts.append(f"type : {info['imgtyp']}")
        if info["exptime"] is not None:
            e = info["exptime"]
            val = f"{int(round(e))}s" if abs(e - round(e)) < 0.05 else f"{e:.1f}s"
            parts.append(f"exp : {val}")
        if info["ra"]:
            parts.append(f"RA {info['ra']} DEC {info['dec']}")
        if info["alt"] is not None:
            suf = "" if info["alt_src"] == "header" else " (calc)"
            parts.append(f"alt : {info['alt']:.1f}°{suf}")
        if info["airmass"] is not None:
            suf = "" if info["airmass_src"] == "header" else " (calc)"
            parts.append(f"airmass : {info['airmass']:.2f}{suf}")
        if info["ccdtemp"] is not None:
            parts.append(f"CCD : {info['ccdtemp']:.1f}°C")
        if info["pier"]:
            parts.append(f"pier : {info['pier']}")
        if info["dateobs"]:
            parts.append(f"DATE-OBS : {info['dateobs']}")
        line = (f"[{now_utc_hm()} TU] Check FITS : {info['file']} ; "
                + " / ".join(parts))
        self._insert_line(line + "\n")
        self._set_status(f"FITS read: {info['file']}")

    # -- Folder / saving --------------------------------------------------------
    def choose_folder(self):
        d = filedialog.askdirectory(
            title="Acquisition Folder",
            initialdir=self.acq_dir or default_browse_dir(),
        )
        if not d:
            return
        self.acq_dir = d
        self.save_path = os.path.join(d, f"obs_{observing_date().isoformat()}.txt")
        self.recapped = set()   # new folder: first recap covers everything there
        # Files already present form the baseline: any file added afterwards
        # is detected as new, whatever its original date.
        existing = list_fits(d)
        self.seen = set(existing)
        self.latest_fits = max(existing, key=_recency_key) if existing else None
        if os.path.exists(self.save_path):
            if messagebox.askyesno(
                "Existing Log",
                "A log already exists for tonight in this folder.\n"
                "Reload it (replaces the current content)?",
            ):
                with open(self.save_path, "r", encoding="utf-8") as fh:
                    self.text.delete("1.0", tk.END)
                    self.text.insert("1.0", fh.read())
                self.text.mark_set(tk.INSERT, tk.END)
                # Log resumed: files already there were most likely logged already
                self.recapped = set(list_fits(d))
        self.path_label.configure(text=self.acq_dir, fg=FG)
        self._set_status(f"Log: {os.path.basename(self.save_path)}")
        self._insert_line(f"[{now_utc_hm()} TU] Files : {self.acq_dir}")

    def save(self):
        if not self.save_path:
            p = filedialog.asksaveasfilename(
                title="Save Log",
                initialdir=self.acq_dir or default_browse_dir(),
                defaultextension=".txt",
                initialfile=f"obs_{observing_date().isoformat()}.txt",
            )
            if not p:
                return
            self.save_path = p
        try:
            with open(self.save_path, "w", encoding="utf-8") as fh:
                fh.write(self.text.get("1.0", "end-1c"))
            self._set_status(f"Saved: {self.save_path}  ({now_utc_hm()} TU)")
        except Exception as e:
            self._set_status(f"Save failed: {e}")

    def _quiet_save(self):
        # Sauvegarde silencieuse, uniquement lors d'une action (clic bouton).
        # Aucun timer, aucun réveil CPU périodique.
        if not self.save_path:
            return
        try:
            with open(self.save_path, "w", encoding="utf-8") as fh:
                fh.write(self.text.get("1.0", "end-1c"))
        except Exception:
            pass

    def _on_close(self):
        if self.save_path:
            self.save()
        self.root.destroy()


def main():
    root = tk.Tk()
    ObservationPad(root)
    root.mainloop()


if __name__ == "__main__":
    main()
