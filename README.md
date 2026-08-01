# Observation Pad

English version below.

Bloc-notes en texte brut pour les sessions de spectroscopie, accompagné d'un outil pour vérifier la convergence d'un empilement de spectres bruts. Conçu pour une utilisation sur le terrain. C'est essentiellement une zone de texte libre. Les boutons ne font qu'insérer du texte pré-formaté (horodatage, nom de fichier, etc.) à l'emplacement du curseur.
Le journal est écrit sous la forme `obs_YYYY-MM-DD.txt` dans le répertoire d'acquisition. La sauvegarde est événementielle, à chaque clic de bouton et à la fermeture.


---

## Fichiers

| Fichier | Rôle |
|---|---|
| `observation_pad.py` | le bloc-notes lui-même |
| `spectro_convergence.py` | convergence d'empilement / fenêtre SNR |
| `observer.ini` | coordonnées du site d'observation |

---

## Prérequis

Python 3.9+ avec Tkinter (`sudo apt install python3-tk` sous Debian/Ubuntu).

```
astropy
numpy
matplotlib
```

```bash
pip install -r requirements.txt
python observation_pad.py
```

---

## Le site d'observation

Trois sources, dans cet ordre de priorité :

1. `SITELAT` / `SITELONG` dans l'en-tête FITS.
2. `observer.ini`, à côté du script.
3. Constantes en tête de `observation_pad.py`.

---

## Boutons

**Files** : choisir le répertoire d'acquisition. C'est là que le journal est écrit, et où les fichiers FITS sont recherchés.

**New obs** : clôt la cible en cours et passe à la suivante. Insère un récapitulatif des poses acquises depuis le dernier New obs, puis un séparateur et une nouvelle ligne `*** Target:`.

**Insert TU** : insère un horodatage sur une nouvelle ligne. Tapez votre note à la suite.

**MTO** : hauteur de la Lune et du Soleil, plus météo courante. Voir ci-dessous.

**AOD** : épaisseur optique des aérosols, à partir de deux sources indépendantes. Voir ci-dessous.

**Calib** : insère une ligne avec les poses de calibration habituelles.

**Check last fits** : lit le dernier FITS écrit et inscrit son en-tête dans le journal. `DATE-OBS` est reproduit tel quel sous son nom de mot-clé. Si l'altitude ou la masse d'air ont dû être calculées parce qu'elles manquaient dans l'en-tête, elles sont marquées `(calc)`.

**Airmass** : tableau popup de la masse d'air de chaque pose Light, pour vérification rapide.

**SNR** : ouvre `spectro_convergence.py` sur les poses de la cible en cours (toutes les Lights partageant le même `OBJECT` que la dernière). Le bouton **Log** écrit une ligne de résumé dans ce journal : centre, FWHM, inclinaison, largeur d'extraction, DER_SNR et Buil SNR.

**Save** : sauvegarde le journal, rarement nécessaire puisque chaque bouton sauvegarde déjà.

**End** : après confirmation, ajoute un résumé de la nuit complète.

---

## MTO : ciel et météo

Insère deux lignes décrivant le ciel courant :

```
[21:30 TU] Moon illum. 99% h 11° Sun h -12°
[21:30 TU] open_meteo : T 27°C / T-Td 14.8°C / gusts 21 km/h / cirrus 0%
```

**Éclairement et altitude de la Lune** calculés localement avec astropy.

**Altitude du Soleil** −12° correspond au crépuscule nautique, −18° à la nuit astronomique complète.

**T − Td**, l'écart entre la température de l'air et le point de rosée, = la marge, en degrés C, avant condensation.

**Gusts** les fichues rafales qui décrochent le guidage.

**Cirrus** correspond à `cloud_cover_high`, les nuages hauts et fins, tueur silencieux de la spectrophotométrie. Ils atténuent sans avoir l'air nuageux, on les remarque rarement la nuit, et ils compliquent la courbe de réponse.

La météo vient d'[Open-Meteo](https://open-meteo.com/) (pas de clé API). Au-dessus de la France, cela correspond au modèle Météo-France AROME à 1–2 km. C'est un **modèle météo**, pas une mesure sur votre site.

---

## AOD : épaisseur optique des aérosols

Les aérosols gouvernent l'extinction atmosphérique, l'AOD est donc important quand on corrige les spectres. Deux sources indépendantes sont consignées :

```
[21:31 TU] AOD 550nm : 0.137 (CAMS model, 2026-07-28 22:00 TU (+0.4 h vs frame))
[21:31 TU] AERONET : Toulouse_MF (21 km, 160 m) / AOD550 0.122 (from 500nm) / AOD500 0.144 / Angstrom 1.75 / lunar / 2026-07-28 21:05 TU (-0.6 h vs frame)
```

**CAMS** (via air-quality d'Open-Meteo) est un **modèle**. Il est toujours disponible, quel que soit le temps, mais sa précision est modeste (erreur absolue moyenne d'environ 0.08).

**AERONET** est une **mesure au sol** du photomètre CIMEL soleil/lune le plus proche, précise à environ ±0.03. La station la plus proche est trouvée automatiquement à partir des coordonnées du site. `lunar_merge` est activé : avec une Lune suffisamment brillante, ce sont de vraies mesures *nocturnes* potentiellement contemporaines de vos poses.

Trois éléments sont rapportés pour juger si la valeur est exploitable :

- **Distance et altitude de la station.** Un photomètre 2700 m plus haut que vous ne voit pas la couche d'aérosols en dessous et sous-estimera votre extinction. Au-delà de 500 m de différence, la ligne est marquée `!! alt … vs site`.
- **Heure de la mesure et son décalage** par rapport à vos poses.
- **Quelle longueur d'onde a réellement été mesurée.** AERONET n'a pas de canal à 550 nm, et en mode lunaire de nombreux canaux ne renvoient rien, donc le canal utilisable le plus proche est converti à 550 nm via l'exposant d'Ångström, ce qui le rend directement comparable au CAMS. L'`AOD500` brut est affiché à côté quand il est disponible.

Si rien n'est publié pour votre fenêtre, le journal reçoit l'URL de la requête pour que vous puissiez vérifier vous-même.

**Les deux sources sont à la date de l'observation** : le `DATE-OBS` du dernier FITS quand il y en a un, sinon l'heure courante. Cela signifie que vous pouvez relancer le pad le lendemain d'une session et retrouver les valeurs d'aérosols contemporaines de vos poses. Contrairement à MTO, AOD fait deux requêtes réseau et peut être lent, il est prévu pour être cliqué une ou deux fois par nuit, ce qui est acceptable car la charge en aérosols évolue lentement.

---

## Convergence Spectro

Sélectionnez une série de spectres bruts 2D du même objet. La fenêtre montre le profil spatial de la fente avec les zones d'extraction et de ciel, et la convergence du signal empilé avec sa bande de dispersion.

Deux estimateurs de SNR sont donnés : **DER_SNR**, calculé à partir du spectre lui-même, et l'estimateur **Buil** basé sur les statistiques de l'empilement.

Le bouton **Log** écrit le résumé dans le journal d'observation.

---

## Sources de données

- Météo et aérosols CAMS : [Open-Meteo](https://open-meteo.com/) — CC BY 4.0
- Mesures d'aérosols : [AERONET](https://aeronet.gsfc.nasa.gov/) (NASA GSFC). [Citer les données](https://aeronet.gsfc.nasa.gov/new_web/data_usage.html).

## Licence
MIT

---
---

# Observation Pad

A plain-text observing notepad for spectroscopy sessions, plus a companion tool to check the convergence of a stack of raw spectra. Written for field use. Basically a free-form text area. The buttons only insert pre-formatted text (timestamp filename etc) at the cursor. 
The log is written as `obs_YYYY-MM-DD.txt` inside the acquisition folder. Saving is event-driven, on every button click and on close.

## Files

| File | Role |
|---|---|
| `observation_pad.py` | the notepad itself |
| `spectro_convergence.py` | stacking convergence / SNR window |
| `observer.ini` | observing site coordinates |

---

## Requirements

Python 3.9+ with Tkinter (`sudo apt install python3-tk` on Debian/Ubuntu).

```
astropy
numpy
matplotlib
```

```bash
pip install -r requirements.txt
python observation_pad.py
```

---

## The observing site

Three sources, in this order of priority:

1. `SITELAT` / `SITELONG` in the FITS header.
2. `observer.ini`, next to the script.
3. Constants at the top of `observation_pad.py`.

---

## Buttons


**Files** : choose the acquisition folder. This is where the log is written,
and where FITS files are looked for.

**New obs** : closes the current target and opens the next one. Inserts a
recap of the frames acquired since the last New obs, then a separator and a fresh
`*** Target:` line.

**Insert TU** : inserts a timestamp on a new line. Type your note after it.

**MTO** : Moon and Sun hight, plus current weather. See below.

**AOD** : aerosol optical depth, from two independent sources. See below.

**Calib** : inserts one line with the usual calibration exposures.

**Check last fits** : reads the most recently written FITS and logs its header. `DATE-OBS` is reproduced verbatim under its keyword name. If altitude or airmass had to be computed because the header lacked them, they are tagged `(calc)`.

**Airmass** : popup table of the airmass of every Light frame, for quick sanity check.

**SNR** : opens `spectro_convergence.py` on the frames of the current target (all Lights sharing the last frame's `OBJECT`). Its **Log** button writes a summary line back into this log: centre, FWHM, tilt, extraction width,
DER_SNR and Buil SNR.

**Save** : save the log, rarely needed, since every button already saves.

**End** : after confirmation, appends a summary of the whole night.

---

## MTO — sky and weather

Inserts two lines describing the current sky:

```
[21:30 TU] Moon illum. 99% h 11° Sun h -12°
[21:30 TU] open_meteo : T 27°C / T-Td 14.8°C / gusts 21 km/h / cirrus 0%
```

**Moon illumination and altitude** computed locally with astropy.

**Sun altitude**  −12° is nautical twilight, −18° is full astronomical darkness.

**T − Td**, the spread between air temperature and dew point, = the margin, in C degrees, before condensation.

**Gusts** damn gusts that throw the guiding off. 

**Cirrus** is `cloud_cover_high`, high thin cloud silent killer for spectrophotometry. Attenuates without looking cloudy, you rarely notice it at night, and complicates your response curve.

Weather comes from [Open-Meteo](https://open-meteo.com/) (no API key). Over France that means Météo-France AROME at 1–2 km. It is a **model on a grid cell**, not a measurement at your site.

---

## AOD — aerosol optical depth

Aerosols drive atmospheric extinction, so AOD matters when you correct spectra. Two independent sources are logged side by side:

```
[21:31 TU] AOD 550nm : 0.137 (CAMS model, 2026-07-28 22:00 TU (+0.4 h vs frame))
[21:31 TU] AERONET : Toulouse_MF (21 km, 160 m) / AOD550 0.122 (from 500nm) / AOD500 0.144 / Angstrom 1.75 / lunar / 2026-07-28 21:05 TU (-0.6 h vs frame)
```

**CAMS** (via Open-Meteo's air-quality endpoint) is a **model**. It is always available, whatever the weather, but its accuracy is modest (mean absolute error around 0.08).

**AERONET** is a **ground measurement** from the nearest CIMEL sun/moon photometer, accurate to about ±0.03 at level 1.5. The nearest station is found automatically from the site coordinates. `lunar_merge` is enabled: with a bright enough Moon, these are real *night-time* measurements possibly contemporaneous with your exposures.

Three things are reported so you can judge whether the value is usable:

- **Distance and elevation of the station.** A photometer 2700 m higher than you does not see the aerosol layer below it and will underestimate your extinction. Beyond 500 m of difference the line is flagged `!! alt … vs site`.
- **Measurement time and its offset** from your frames.
- **Which wavelength was actually measured.** AERONET has no 550 nm channel, and in lunar mode many channels return no value, so the nearest usable one is converted to 550 nm through the Ångström exponent, making it directly comparable to CAMS. The raw `AOD500` is shown alongside when available.

If nothing is published for your window, the log gets the query URL so you can look yourself.

**Both sources are anchored on the same reference time**: the `DATE-OBS` of the last FITS when there is one, otherwise the current time. This means you can relaunch the pad the day *after* a run and still retrieve the aerosol values contemporaneous with your exposures. Unlike MTO, AOD makes two network requests and can be slow — it is meant to be clicked once or twice a night, ok because aerosol load changes slowly.

---

## Spectro Convergence

Select a series of raw 2D spectra of the same target. The window shows the slit spatial profile with the extraction and sky zones, and the convergence of the stacked signal with its dispersion band.

Two SNR estimates are given: **DER_SNR**, computed from the spectrum itself, and the **Buil** estimator based on the stack statistics.

The **Log** button writes the summary into the observation log.

---

## Data sources

- Weather and CAMS aerosols: [Open-Meteo](https://open-meteo.com/) — CC BY 4.0
- Aerosol measurements: [AERONET](https://aeronet.gsfc.nasa.gov/) (NASA GSFC). [Cite the data](https://aeronet.gsfc.nasa.gov/new_web/data_usage.html).

## License
MIT
