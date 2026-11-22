// Local ZIP export. fflate 0.8.2 is vendored with its MIT license.
function normalizedFrame(frame, count) {
  return ((frame % count) + count) % count;
}

function motionExportFilename(sequence, frame, count) {
  const fraction = Number((normalizedFrame(frame, count) / count).toFixed(6));
  return sequence + '_' + fraction + 'T.zip';
}

const motionPlotScript = String.raw`"""Plot the exported snapshot and full hinge trajectories on one 3D axes.

Run: python3 plot_motion.py. Install: python3 -m pip install matplotlib
All input files are resolved relative to this script, not the working directory.
"""
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
LINK_COLORS = {"A": "#E25C5C", "B": "#4A6B82", "C": "#D4A337"}
PATH_COLOR = "#B4A000"  # Dark yellow stays legible on Matplotlib's white background.


def read_motion(path):
    """Read exported numData assignment as JSON, without executing JavaScript."""
    text = Path(path).read_text(encoding="ascii")
    prefix = "var numData = "
    start = text.index(prefix) + len(prefix)
    data, end = json.JSONDecoder().raw_decode(text[start:])
    if text[start + end:].strip() != ";":
        raise ValueError("Unexpected trailing JavaScript")
    return data


def read_table(name):
    with (ROOT / name).open(newline="", encoding="ascii") as stream:
        return [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(stream)]


def load_bundle():
    meta = json.loads((ROOT / "metadata.json").read_text(encoding="ascii"))
    sequence = meta["sequence"]
    motion = read_motion(ROOT / (sequence + ".js"))[sequence]
    count, hinges = meta["frame_count"], meta["hinge_count"]
    if count < 1 or hinges != len(sequence) or set(sequence)-set(LINK_COLORS):
        raise ValueError("Invalid frame count, hinge count, or sequence")
    for key in ("bx_data", "by_data", "bz_data"):
        if len(motion[key]) != hinges+1 or any(len(row) != count for row in motion[key]):
            raise ValueError("Unexpected direction array shape: " + key)
        if not all(math.isfinite(value) for row in motion[key] for value in row):
            raise ValueError("Non-finite direction coordinates")
    configuration = sorted(read_table("configuration.csv"), key=lambda row: row["hinge_index"])
    if [row["hinge_index"] for row in configuration] != list(range(1, hinges+1)):
        raise ValueError("Expected one configuration row per hinge")
    frames = [[None]*hinges for _ in range(count)]
    for row in read_table("eversion_paths.csv"):
        frame, hinge = int(row["frame_index"]), int(row["hinge_index"])-1
        if not 0 <= frame < count or not 0 <= hinge < hinges or frames[frame][hinge] is not None:
            raise ValueError("Invalid or duplicate trajectory sample")
        point = tuple(row[axis] for axis in "xyz")
        if not all(math.isfinite(value) for value in point):
            raise ValueError("Non-finite trajectory coordinates")
        frames[frame][hinge] = point
    if any(point is None for frame in frames for point in frame):
        raise ValueError("Missing trajectory samples")
    return meta, motion, configuration, frames


def frame_geometry(meta, motion, frames, frame):
    """Return centers and hinge endpoints, in the exported scene coordinates."""
    frame %= meta["frame_count"]
    centers = frames[frame]
    half_length = meta["hinge_length_parameter"]*1.7/2
    endpoints = []
    for hinge, center in enumerate(centers, start=1):
        # Browser scene mapping: (model x, model z, -model y).
        direction = (motion["bx_data"][hinge][frame], motion["bz_data"][hinge][frame],
                     -motion["by_data"][hinge][frame])
        a = tuple(c-half_length*d for c, d in zip(center, direction))
        b = tuple(c+half_length*d for c, d in zip(center, direction))
        endpoints.append((a, b))
    return centers, endpoints


def make_axes(meta, motion, frames, show_paths=True):
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    if show_paths:
        for hinge in range(meta["hinge_count"]):
            path = [frame[hinge] for frame in frames]
            ax.plot(*zip(*path), color=PATH_COLOR, linewidth=1, alpha=0.55)
    # One fixed cubic envelope includes every frame and both hinge endpoints.
    points = []
    for frame in range(meta["frame_count"]):
        _, endpoints = frame_geometry(meta, motion, frames, frame)
        points.extend(point for pair in endpoints for point in pair)
    low = [min(point[axis] for point in points) for axis in range(3)]
    high = [max(point[axis] for point in points) for axis in range(3)]
    half_span = max(max(b-a for a, b in zip(low, high))*0.55, 1e-6)
    for setter, a, b in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), low, high):
        midpoint = (a+b)/2
        setter(midpoint-half_span, midpoint+half_span)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.view_init(elev=25, azim=-55)
    keys = [Line2D([], [], color=color, linewidth=3, label=letter)
            for letter, color in LINK_COLORS.items()]
    keys.append(Line2D([], [], color="black", linewidth=3, label="Hinge"))
    if show_paths:
        keys.append(Line2D([], [], color=PATH_COLOR, linewidth=1, label="Hinge trajectory"))
    ax.legend(handles=keys, loc="upper left", fontsize=9)
    fig.subplots_adjust(left=0.02, right=0.96, bottom=0.08, top=0.9)
    return fig, ax


def plot():
    meta, motion, configuration, frames = load_bundle()
    fig, ax = make_axes(meta, motion, frames)
    centers = [tuple(row["center_"+axis] for axis in "xyz") for row in configuration]
    for row in configuration:
        a, b = (tuple(row[prefix+axis] for axis in "xyz") for prefix in ("a_", "b_"))
        ax.plot(*zip(a, b), color="black", linewidth=3)
    for i, center in enumerate(centers):
        following = centers[(i+1) % len(centers)]
        ax.plot(*zip(center, following), color=LINK_COLORS[meta["sequence"][(i+1) % len(centers)]], linewidth=3)
    ax.set_title("{} at {:.6f}T ({})\nSnapshot and full-cycle hinge trajectories".format(
        meta["sequence"], meta["phase_fraction"], meta["chirality"]))
    return fig


if __name__ == "__main__":
    plot()
    plt.show()
`;

const motionAnimationScript = String.raw`"""Animate every saved frame, with optional full-cycle hinge trajectories.

python3 animate_motion.py
python3 animate_motion.py --save motion.gif --fps 30
python3 animate_motion.py --save motion.mp4 --fps 30  # Requires ffmpeg.
python3 animate_motion.py --no-paths
Requires plot_motion.py and the data files beside this script.
FPS controls playback speed only; exported T is not physical time in seconds.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter

from plot_motion import LINK_COLORS, frame_geometry, load_bundle, make_axes


def animate(fps=30, show_paths=True):
    if fps <= 0:
        raise ValueError("FPS must be positive")
    meta, motion, _, frames = load_bundle()
    fig, ax = make_axes(meta, motion, frames, show_paths=show_paths)
    hinges = [ax.plot([], [], [], color="black", linewidth=3)[0]
              for _ in range(meta["hinge_count"])]
    links = [ax.plot([], [], [], color=LINK_COLORS[meta["sequence"][(i+1) % meta["hinge_count"]]], linewidth=3)[0]
             for i in range(meta["hinge_count"])]
    title = ax.set_title("")

    def update(frame):
        centers, endpoints = frame_geometry(meta, motion, frames, frame)
        for i, (a, b) in enumerate(endpoints):
            hinges[i].set_data_3d(*zip(a, b))
            links[i].set_data_3d(*zip(centers[i], centers[(i+1) % len(centers)]))
        title.set_text("{} at {:.6f}T ({})".format(
            meta["sequence"], frame/meta["frame_count"], meta["chirality"]))
        return hinges+links+[title]

    start = meta["displayed_frame_index"]
    order = [(start+i) % meta["frame_count"] for i in range(meta["frame_count"])]
    animation = FuncAnimation(fig, update, frames=order, init_func=lambda: update(start),
                              interval=1000/fps, repeat=True, blit=False, cache_frame_data=False)
    return fig, animation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fps", type=float, default=30, help="Playback frames per second (default: 30)")
    parser.add_argument("--no-paths", action="store_true", help="Hide full-cycle trajectories")
    parser.add_argument("--save", type=Path, help="Write a .gif or .mp4 instead of opening a window")
    args = parser.parse_args()
    if not 0 < args.fps < float("inf"):
        parser.error("--fps must be finite and positive")
    writer = None
    if args.save:
        if args.save.suffix.lower() == ".gif":
            writer = PillowWriter(fps=args.fps)
        elif args.save.suffix.lower() == ".mp4":
            if not FFMpegWriter.isAvailable():
                parser.error("MP4 export requires ffmpeg; use .gif or install ffmpeg")
            writer = FFMpegWriter(fps=args.fps)
        else:
            parser.error("--save filename must end with .gif or .mp4")
    fig, animation = animate(fps=args.fps, show_paths=not args.no_paths)
    if args.save:
        animation.save(str(args.save), writer=writer, dpi=100)
        plt.close(fig)
        print("Saved:", args.save.resolve())
    else:
        plt.show()


if __name__ == "__main__":
    main()
`;

const motionExportReadme = String.raw`LINKAGE EVERSION EXPORT

QUICK START
1. Extract the ZIP before running Python. Keep all seven files together.
2. Open a terminal in the extracted folder.
3. Install the plotting dependency:
   python3 -m pip install matplotlib
4. Display the exported configuration and all hinge trajectories in one 3D plot:
   python3 plot_motion.py
5. Play the full saved motion:
   python3 animate_motion.py

Close the plot window to exit. The scripts find input files beside themselves,
so you can also run them by their full paths from another working directory.
Use python instead of python3 if that is your Python 3 command.

ANIMATION OPTIONS
python3 animate_motion.py --fps 15
  Slower playback. Default is 30 frames per second; every saved frame is used.
python3 animate_motion.py --no-paths
  Show the moving linkage without the full-cycle trajectory overlay.
python3 animate_motion.py --save motion.gif --fps 30
  Save a GIF instead of opening a plot window.
python3 animate_motion.py --save motion.mp4 --fps 30
  Save an MP4. Requires ffmpeg installed separately and available on PATH.
python3 animate_motion.py --help
  Show command-line options. Options can be combined.

Saved animation paths are relative to your terminal's working directory unless
you provide an absolute path. Existing output files with the same name may be
overwritten. GIF saving can use substantial memory for long, large animations.

WHAT THE PLOTS SHOW
The static plot overlays the frozen exported configuration and nine full-cycle
hinge-center trajectories. The animation starts at that exported frame, visits
all saved frames, then repeats. Axis limits and scale stay fixed throughout.
A/B/C links use the visualizer palette; hinge segments are black. Trajectories
are dark yellow for readability on Matplotlib's white background.

FPS controls playback speed only, not physical time. T means one full saved
cycle. Phase is frame_index/frame_count. Each cycle contains frames 0 through
frame_count - 1, without an extra copy of frame zero at the end. Frame zero is
not necessarily the mechanism's special geometric reference configuration.

PACKAGE CONTENTS
<sequence>.js
  Full-precision direction arrays in the original format:
  var numData = {"<sequence>": {"bx_data": ..., "by_data": ..., "bz_data": ...}};
  Each array is [hinge_row][frame_index], shape 10 x frame_count.
  Rows are 0..9; terminal row 9 is negative row 0. read_motion() in
  plot_motion.py parses this assignment as JSON without executing JavaScript.
configuration.csv
  The frozen displayed configuration: nine hinges, one per row. Includes
  hinge_index, center_x/y/z, and the two hinge endpoints a_x/y/z and b_x/y/z.
eversion_paths.csv
  All nine hinge-center trajectories across every saved frame. Columns:
  frame_index,phase_fraction,hinge_index,x,y,z. Frames are zero-based;
  hinge indices are 1..9. These are full-cycle paths, unlike the browser overlay
  that currently samples about half of the saved cycle.
metadata.json
  Sequence, chirality, hinge count, frame count, exported frame, exact sampled
  phase fraction, displayed clock phase, hinge length, and coordinate convention.
plot_motion.py
  Combined static plot and shared data-reading helpers.
animate_motion.py
  Full-cycle Matplotlib playback and GIF/MP4 saving.
README.txt
  These instructions.

COORDINATES AND READING THE MOTION
Coordinates match exported scene geometry, including chirality, centering,
length normalization to 6, and the existing model z offset -0.25. The mapping
from model to scene coordinates is (model x, model z, -model y). Camera
transforms are excluded; Python does not reproduce the browser's current view.

Animation reads exact centers from eversion_paths.csv. Hinge directions come
from rows 1..9 of the .js arrays, mapped as (bx, bz, -by). Endpoints are center
plus/minus direction times 1.7*hinge_length_parameter/2. Consecutive centers are
connected cyclically, including last to first. Snapshot plotting reads the
frozen endpoints directly from configuration.csv.

This is a numerical line plot, not the browser or Blender's cylinder rendering.
A visible interactive animation requires a desktop Matplotlib GUI backend;
use --save on a headless machine. MP4 saving reports an error if ffmpeg is
unavailable; GIF saving uses Matplotlib's Pillow writer.
`;

function buildMotionExport(snapshot) {
  const {sequence, motion, chirality, hingeLength, displayedPhase} = snapshot;
  const count = motion.bx_data[0].length;
  const frame = normalizedFrame(snapshot.frame, count);
  const geometry = evaluateFrameGeometry(motion, frame, chirality, hingeLength);
  const configuration = ['hinge_index,center_x,center_y,center_z,a_x,a_y,a_z,b_x,b_y,b_z'];
  geometry.forEach((hinge, index) => {
    configuration.push([index + 1, ...hinge.center, ...hinge.a, ...hinge.b].join(','));
  });
  const paths = ['frame_index,phase_fraction,hinge_index,x,y,z'];
  for (let index = 0; index < count; index++) {
    evaluateFrameGeometry(motion, index, chirality, hingeLength).forEach((hinge, i) => {
      paths.push([index, index/count, i + 1, ...hinge.center].join(','));
    });
  }
  if ([...configuration, ...paths].some(line => /NaN|Infinity/.test(line))) {
    throw new Error('Motion contains non-finite geometry.');
  }
  const metadata = {
    format_version: 1, sequence,
    chirality: chirality === -1 ? 'left-handed' : 'right-handed',
    chirality_parameter: chirality, hinge_count: geometry.length,
    frame_count: count, displayed_frame_index: frame, phase_fraction: frame/count,
    phase_fraction_numerator: frame, phase_fraction_denominator: count,
    displayed_clock_phase_radians: displayedPhase,
    displayed_clock_phase_fraction: ((displayedPhase / (2*Math.PI)) % 1 + 1) % 1,
    hinge_length_parameter: hingeLength,
    coordinate_convention: 'Scene coordinates (x,z,-y) after model centering, length normalization to 6, and model z offset -0.25. Camera transform excluded.',
    period_definition: 'T is all saved frames, one playback-clock revolution; not physical seconds.',
    indexing: 'Frames zero-based; hinge indices 1 through 9 match displayed geometry. Direction array row 9 is terminal, negative row 0.',
  };
  const comments = '// Full saved cycle for this exact sequence, including manual switches.\n' +
    '// Each coordinate array is [hinge_row][frame_index], shape 10 x frame_count.\n' +
    '// Read safely in Python using read_motion() in plot_motion.py.\n' +
    '// Chirality and frozen frame are in metadata.json, not encoded in direction arrays.\n';
  const files = {
    [sequence + '.js']: comments + 'var numData = ' + JSON.stringify({[sequence]: motion}) + ';\n',
    'configuration.csv': configuration.join('\n') + '\n',
    'eversion_paths.csv': paths.join('\n') + '\n',
    'metadata.json': JSON.stringify(metadata, null, 2) + '\n',
    'plot_motion.py': motionPlotScript,
    'animate_motion.py': motionAnimationScript,
    'README.txt': motionExportReadme,
  };
  return {filename: motionExportFilename(sequence, frame, count), files};
}

async function saveMotionExport(snapshot) {
  const filename = motionExportFilename(snapshot.sequence, snapshot.frame, snapshot.motion.bx_data[0].length);
  // Picker must be invoked before any deferred work to retain user activation.
  let handle = null;
  if (typeof window.showSaveFilePicker === 'function' && window.isSecureContext) {
    try {
      handle = await window.showSaveFilePicker({suggestedName: filename,
        types: [{description: 'Linkage motion ZIP', accept: {'application/zip': ['.zip']}}]});
    } catch (error) {
      if (error.name === 'AbortError') return false;
      if (error.name !== 'SecurityError' && error.name !== 'NotAllowedError') throw error;
      // Embedded/restricted contexts can expose an unusable picker. Download instead.
    }
  }
  const bundle = buildMotionExport(snapshot);
  const entries = {};
  for (const [name, text] of Object.entries(bundle.files)) entries[name] = fflate.strToU8(text);
  const zip = fflate.zipSync(entries, {level: 6});
  if (handle) {
    const writable = await handle.createWritable();
    try {
      await writable.write(zip);
      await writable.close();
    } catch (error) {
      try { await writable.abort(); } catch (_) { /* Original failure is reported. */ }
      throw error;
    }
  } else {
    const url = URL.createObjectURL(new Blob([zip], {type: 'application/zip'}));
    const anchor = document.createElement('a');
    anchor.href = url;
    anchor.download = filename;
    document.body.appendChild(anchor);
    try { anchor.click(); } finally {
      anchor.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    }
  }
  return true;
}
