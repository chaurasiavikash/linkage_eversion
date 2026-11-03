// Local ZIP export. fflate 0.8.2 is vendored with its MIT license.
function normalizedFrame(frame, count) {
  return ((frame % count) + count) % count;
}

function motionExportFilename(sequence, frame, count) {
  const fraction = Number((normalizedFrame(frame, count) / count).toFixed(6));
  return sequence + '_' + fraction + 'T.zip';
}

const motionPlotScript = String.raw`"""Run: python3 plot_motion.py. Requires matplotlib (pip install matplotlib)."""
import csv
import json
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


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


def plot():
    meta = json.loads((ROOT / "metadata.json").read_text(encoding="ascii"))
    motion = read_motion(ROOT / (meta["sequence"] + ".js"))[meta["sequence"]]
    assert len(motion["bx_data"][0]) == meta["frame_count"]
    configuration = read_table("configuration.csv")
    paths = read_table("eversion_paths.csv")
    fig = plt.figure(figsize=(12, 7))
    current = fig.add_subplot(121, projection="3d")
    trajectories = fig.add_subplot(122, projection="3d")
    centers = [[row[axis] for axis in ("center_x", "center_y", "center_z")]
               for row in configuration]
    for row in configuration:
        a = [row["a_" + axis] for axis in "xyz"]
        b = [row["b_" + axis] for axis in "xyz"]
        current.plot(*zip(a, b), color="black")
    current.plot(*zip(*(centers + centers[:1])), color="tab:blue")
    for hinge in range(1, meta["hinge_count"] + 1):
        rows = [row for row in paths if int(row["hinge_index"]) == hinge]
        trajectories.plot(*[[row[axis] for row in rows] for axis in "xyz"],
                          label=str(hinge))
    current.set_title("{} at {:.6f}T ({})".format(
        meta["sequence"], meta["phase_fraction"], meta["chirality"]))
    trajectories.set_title("Hinge-center paths over one saved cycle")
    for ax in (current, trajectories):
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_aspect("equal")
    fig.subplots_adjust(left=0.04, right=0.96, bottom=0.12, top=0.84, wspace=0.18)
    return fig


if __name__ == "__main__":
    plot()
    plt.show()
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
    'README.txt': 'LINKAGE EVERSION EXPORT\n\n' +
      'Install plotting dependency: python3 -m pip install matplotlib\n' +
      'Run: python3 plot_motion.py (files resolve relative to the script).\n\n' +
      'The sequence .js file preserves the numData assignment and full precision.\n' +
      'Arrays bx_data/by_data/bz_data are direction components, [row][frame].\n' +
      'There are 10 direction rows (0..9); terminal row 9 = -row 0.\n' +
      'configuration.csv contains the frozen displayed hinge centers and endpoints.\n' +
      'eversion_paths.csv contains nine center paths over every saved frame.\n' +
      'Frames are zero-based; geometry hinge indices are 1..9.\n' +
      'T means all saved frames, not seconds. Phase = frame_index/frame_count.\n' +
      'Frame zero is not duplicated at the end. Paths sample the full cycle,\n' +
      'unlike the current on-screen path overlay, which samples about half.\n' +
      'metadata.json records exact frame fraction, clock phase, chirality, and\n' +
      'hinge length. Coordinates match scene geometry, before camera transforms.\n' +
      'Scene mapping is (model x, model z, -model y), after existing centering,\n' +
      'length normalization to 6 and model z offset -0.25.\n' +
      'The Python figure connects centers to show linkage connectivity.\n',
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
