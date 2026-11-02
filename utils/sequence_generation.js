// Sequence distance uses linear adjacent swaps, matching the switch button.
function validateSequence(input) {
  const sequence = input.trim().toUpperCase();
  if (!/^[ABC]{9}$/.test(sequence) ||
      ['A', 'B', 'C'].some(letter => sequence.split(letter).length - 1 !== 3)) {
    throw new Error('Enter 9 letters with exactly 3 A, 3 B, and 3 C.');
  }
  return sequence;
}

function sequenceSwitchPath(source, target) {
  const letters = source.split('');
  const switches = [];
  for (let position = 0; position < target.length; position++) {
    const found = letters.indexOf(target[position], position);
    if (found < 0) throw new Error('Source and target must have the same letter counts.');
    for (let index = found - 1; index >= position; index--) {
      [letters[index], letters[index + 1]] = [letters[index + 1], letters[index]];
      switches.push(index);
    }
  }
  return switches;
}

function generateSequenceMotion(input, database) {
  const sequence = validateSequence(input);
  let base;
  let switches;
  for (const candidate of Object.keys(database)) {
    const path = sequenceSwitchPath(candidate, sequence);
    if (!switches || path.length < switches.length) {
      base = candidate;
      switches = path;
    }
  }
  if (!base) throw new Error('No saved motion data is available.');
  const motion = {};
  const keys = ['bx_data', 'by_data', 'bz_data'];
  for (const key of keys) motion[key] = database[base][key].map(row => row.slice());
  const arrays = keys.map(key => motion[key]);
  for (const index of switches) {
    for (let frame = 0; frame < arrays[0][0].length; frame++) {
      const vectors = arrays.map(array => array.map(row => row[frame]));
      sequence_switch(...vectors, index + 1);
      arrays.forEach((array, axis) => { array[index + 1][frame] = vectors[axis][index + 1]; });
    }
  }
  return { sequence, base, switches, motion };
}
