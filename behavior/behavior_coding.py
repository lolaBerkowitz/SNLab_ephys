import csv
from functools import partial
from pathlib import Path

import napari
from napari.utils.notifications import show_info

# this is where data will be written
CSV_OUT = Path('~/Desktop/data.csv').expanduser()
if not CSV_OUT.exists():
    CSV_OUT.write_text("file,frame,action\n")

# adjust keybindings to your liking
KEYMAP = {
    'd': 'object1_start',
    'f': 'object1_stop',
    'j': 'object2_start',
    'k': 'object2_stop',
}

viewer = napari.Viewer()

# this writes the frame, layer source, and action each time you press a key
def on_keypress(key, viewer):
    action = KEYMAP[key]
    frame = viewer.dims.current_step[0]
    layer = viewer.layers.selection.active or viewer.layers[-1]

    show_info(action)  # if you want some visual feedback
    with open(CSV_OUT, 'a') as f:
        csv.writer(f).writerow([layer.source.path, frame, action])

for key in KEYMAP:
    viewer.bind_key(key, partial(on_keypress, key))

napari.run()