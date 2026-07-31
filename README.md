# Observation Pad

A plain-text observing notepad for spectroscopy sessions, plus a companion
tool to check the convergence of a stack of raw spectra. Written for field use. Observation Pad is a **free-form text area**: you type whatever you want. The buttons only *insert* pre-formatted text at the cursor, so you never have to type a
timestamp or copy a filename by hand.

The log is written as `obs_YYYY-MM-DD.txt` **inside the acquisition folder**.
Saving is **event-driven** — on every button click and on close.
---

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

1. **`SITELAT` / `SITELONG` in the FITS header** — the real site of the
   night, as recorded by the mount.
2. **`observer.ini`**, next to the script. Covers the start of the night, before the first frame.
3. **Constants at the top of `observation_pad.py`** — last resort.

Edit `observer.ini` when you travel; the calculations follow automatically.

---

## Buttons


**Files** — choose the acquisition folder. This is where the log is written,
and where FITS files are looked for. Logs the path.

**New obs** — closes the current target and opens the next one. Inserts a
recap of the frames acquired *since the last New obs* — count, names,
exposure budget, mean airmass, CCD temperature — then a separator and a fresh
`*** Target:` line, with the cursor ready for the target name.

**Insert TU** — inserts `[HH:MM TU]` on a new line. Type your note after it.
This is the workhorse: one button, then free text.

**MTO** — Moon and Sun geometry, plus current weather. See below.

**AOD** — aerosol optical depth, from two independent sources. See below.

**Calib** — inserts one line with the usual calibration exposures:
`bias (0 sec) dark (600 sec) neon (13 sec) flat (6 sec)`. Edit by hand if the
night differs.

**Check last fits** — reads the most recently written FITS and logs its
header

*A note on labels:* `DATE-OBS` is reproduced **verbatim** under its keyword
name. If altitude or airmass had to be computed because
the header lacked them, they are tagged `(calc)`.

**Airmass** — popup table of the airmass of every Light frame, grouped by
target with a per-target average. Calibration frames are excluded. Its point
is a quick sanity check: compare the average airmass of your reference star
to that of the science target, and see at a glance whether the reference was
badly placed in the sky.

**SNR** — opens `spectro_convergence.py` on the frames of the current target
(all Lights sharing the last frame's `OBJECT`). Its **Log** button writes a
summary line back into this log: centre, FWHM, tilt, extraction width,
DER_SNR and Buil SNR.

**Save** — explicit save. Rarely needed, since every button already saves.

**End** — after confirmation, appends a summary of the whole night: every
FITS in the folder, grouped by target, with type and exposure time.

---

## MTO — sky and weather

Inserts two lines describing the current sky:

```
[21:30 TU] Moon illum. 99% h 11° Sun h -12°
[21:30 TU] open_meteo : T 27°C / T-Td 14.8°C / gusts 21 km/h / cirrus 0%
```

**Moon illumination and altitude** omputed locally with astropy.

**Sun altitude**  −12° is nautical twilight, −18° is full astronomical darkness.

**T − Td**, the spread between air temperature and dew point, the margin, in degrees, before condensation.
**Gusts** rather than mean wind, because gusts are what throw the guiding
off. **Cirrus** is `cloud_cover_high` specifically: high thin cloud is the
silent killer for spectrophotometry — it attenuates without looking cloudy,
you rarely notice it at night, and weeks later it explains a response curve
that will not fit.

Weather comes from [Open-Meteo](https://open-meteo.com/) (no API key). Over
France that means Météo-France AROME at 1–2 km. It is a **model on a grid
cell**, not a measurement at your site.

---

## AOD — aerosol optical depth

Aerosols drive atmospheric extinction, so AOD matters when you correct
spectra. Two independent sources are logged side by side:

```
[21:31 TU] AOD 550nm : 0.137 (CAMS model, 2026-07-28 22:00 TU (+0.4 h vs frame))
[21:31 TU] AERONET : Toulouse_MF (21 km, 160 m) / AOD550 0.122 (from 500nm) / AOD500 0.144 / Angstrom 1.75 / lunar / 2026-07-28 21:05 TU (-0.6 h vs frame)
```

**CAMS** (via Open-Meteo's air-quality endpoint) is a **model**. It is always
available, whatever the weather, but its accuracy is modest — mean absolute
error around 0.08, which is the same order as a clear-night AOD. Treat it as
a qualitative indication.

**AERONET** is a **ground measurement** from the nearest CIMEL sun/moon
photometer, accurate to about ±0.03 at level 1.5. The nearest station is found
automatically from the site coordinates, so it follows you when you travel.
Crucially, `lunar_merge` is enabled: with a bright enough Moon, these are real
*night-time* measurements, contemporaneous with your exposures.

Three things are reported so you can judge whether the value is usable:

- **Distance and elevation of the station.** A photometer 2700 m higher than
  you does not see the aerosol layer *below* it and will underestimate your
  extinction. Beyond 500 m of difference the line is flagged `!! alt … vs site`.
- **Measurement time and its offset** from your frames.
- **Which wavelength was actually measured.** AERONET has no 550 nm channel,
  and in lunar mode many channels return no value, so the nearest usable one
  is converted to 550 nm through the Ångström exponent — making it directly
  comparable to CAMS. The raw `AOD500` is shown alongside when available.

If nothing is published for your window, the log gets the query URL so you can
go and look yourself.

**Both sources are anchored on the same reference time**: the `DATE-OBS` of
the last FITS when there is one, otherwise the current time. This means you
can relaunch the pad the day *after* a run and still retrieve the aerosol
values contemporaneous with your exposures. Unlike MTO, AOD makes two network
requests and can be slow — it is meant to be clicked once or twice a night,
which is plenty since aerosol load changes slowly.

---

## Spectro Convergence

Select a series of raw 2D spectra of the same target. The window shows the
spatial profile with the extraction and sky zones, and the convergence of the
stacked signal with its dispersion band — so you can see whether adding more
frames is still buying you signal.

Two SNR estimates are given: **DER_SNR**, computed from the spectrum itself,
and the **Buil** estimator based on the stack statistics.

The **Log** button writes the summary into the observation log. Run standalone
(`python spectro_convergence.py f1.fits f2.fits …`) it copies the line to the
clipboard instead.

---

## Data sources

- Weather and CAMS aerosols: [Open-Meteo](https://open-meteo.com/) — CC BY 4.0
- Aerosol measurements: [AERONET](https://aeronet.gsfc.nasa.gov/) (NASA GSFC).
  Please [cite the data](https://aeronet.gsfc.nasa.gov/new_web/data_usage.html)
  and contact the site PIs if you publish anything based on them.

## License

MIT
