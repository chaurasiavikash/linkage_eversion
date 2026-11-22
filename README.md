# Linkage eversion visualizer

Explore cyclic linkages with three links each of types A, B, and C. Choose a
preset or enter a custom sequence, generate its motion, and use the playback
controls to inspect it. The bottom legend identifies links, hinges, hinge axes,
and hinge trajectories.

## Guided tutorial

Click **Tutorial** at the top of the right sidebar to watch an automated guide.
The on-screen cursor first moves down the left controls, demonstrating the four
base sequences, a custom sequence, adjacent swaps, playback, and speed. It then
follows the right sidebar from top to bottom: fullscreen, swept volume, hinge
trajectories, export, chirality, camera views, and background color. Export, background color, and fullscreen controls
are explained without opening native dialogs.

Steps advance automatically. Use **Pause/Resume**, **Next**, or **Stop** in the
tutorial card. **C** also pauses/resumes outside text inputs; **Esc** closes the
tour. Completing or stopping restores your starting motion, custom input, phase,
chirality, overlays, playback speed, playback state, and camera view. The tutorial
can be run again at any time after it closes.

## Export and plot the motion

Click **Export data** under **Data** in the right sidebar. Playback pauses and
the browser saves a ZIP containing the current configuration and full motion.
Extract the ZIP and keep its seven files together.

From a terminal in the extracted folder:

```bash
python3 -m pip install matplotlib
python3 plot_motion.py
```

The plot combines the exported snapshot and all nine hinge-center trajectories
in one 3D view, with the same A/B/C palette as the visualizer.

## Animate the full cycle

```bash
python3 animate_motion.py
python3 animate_motion.py --no-paths
python3 animate_motion.py --fps 15
```

Animation starts at the exported frame, visits every saved frame, and repeats.
Axes limits stay fixed. FPS changes playback speed; the exported phase is a
fraction of the saved cycle, not elapsed physical time.

Save an animation instead of opening a window:

```bash
python3 animate_motion.py --save motion.gif --fps 30
python3 animate_motion.py --save motion.mp4 --fps 30
```

MP4 requires `ffmpeg` installed separately and available on `PATH`. GIF saving
uses Matplotlib's Pillow writer and can consume substantial memory for large
animations. Saved paths are relative to the terminal's working directory;
existing outputs with the same name may be overwritten. Use
`python3 animate_motion.py --help` for options. Interactive playback requires a
desktop GUI backend; use `--save` on a headless machine.

## Files inside each export

- `<sequence>.js`: full-precision hinge-direction arrays for the entire motion.
- `configuration.csv`: exported frame's nine hinge centers and endpoint pairs.
- `eversion_paths.csv`: nine hinge-center trajectories over every saved frame.
- `metadata.json`: sequence, chirality, frame and phase information, hinge length,
  and coordinate convention.
- `plot_motion.py`: combined plot and safe JavaScript-assignment reader.
- `animate_motion.py`: full-cycle playback and GIF/MP4 saving.
- `README.txt`: detailed commands, array shapes, indexing, and coordinate mapping.

Python uses the exported scene coordinates. Camera transforms are excluded,
so its default view differs from the browser camera. These scripts draw line
geometry rather than shaded cylinders.
