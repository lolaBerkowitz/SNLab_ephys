import csv
from functools import partial
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
import os

import napari
from napari.utils.notifications import show_info

# this is where data will be written
currdir = os.getcwd()
tempdir = filedialog.askdirectory(parent= tk.Tk(), initialdir=currdir, title='Please select a directory')
basename = os.path.basename(tempdir)
save_path = os.path.join(tempdir, basename +'_ol_scoring.csv')

CSV_OUT = Path(save_path)
if not CSV_OUT.exists():
    CSV_OUT.write_text("file,frame,object,action\n")

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
    if (key == 'd') | (key == 'f'):
        obj_id  = '1'
    elif (key == 'j') | (key == 'k'):
        obj_id  = '2'
    if (key == 'd') | (key == 'j'):
        action = 'start'
    elif (key == 'f') | (key == 'k'):
        action = 'stop'
    frame = viewer.dims.current_step[0]
    layer = viewer.layers.selection.active or viewer.layers[-1]

    show_info('object :'+obj_id+' '+action+' Frame: '+str(frame))  # if you want some visual feedback
    with open(CSV_OUT, 'a') as f:
        source_path = getattr(layer.source, 'path', 'unknown')
        csv.writer(f).writerow([source_path, frame, obj_id, action])

def make_callback(k):
    return lambda: on_keypress(k, viewer)

for key in KEYMAP:
    viewer.bind_key(key, lambda viewer, key=key: on_keypress(key, viewer))
napari.run()