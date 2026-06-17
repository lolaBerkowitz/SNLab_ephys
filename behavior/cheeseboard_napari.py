import glob
import os
from itertools import compress
from pathlib import Path
from sys import path

import napari
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import cv2
from napari_video.napari_video import VideoReaderNP


def save_layers_to_video_directory(viewer, video_path, point_layers):
    video_path_obj = Path(video_path)
    video_dir = video_path_obj.parent
    video_name = video_path_obj.stem  # Extracts the filename without the extension

    for layer in point_layers:
        points = layer.data  # Extract point coordinates (N, 2)
        frames = layer.metadata.get("axis-0", np.zeros(len(points), dtype=int))
        visible = layer.metadata.get("visible", np.ones(len(points), dtype=bool))

        if len(points) == 0:
            print(f"No data to save for layer '{layer.name}'.")
            continue  # Skip empty layers

        # Create DataFrame with required columns
        df = pd.DataFrame(
            {
                "index": np.arange(len(points)),  # Index of each point
                "axis-0": frames,  # Frame number
                "axis-1": points[:, 0],  # X-coordinates
                "axis-2": points[:, 1],  # Y-coordinates
            }
        )

        # Updated naming convention: videoname_layername.csv
        save_path = video_dir / f"{video_name}_{layer.name}.csv"
        df.to_csv(save_path, index=False)
        print(f"Saved {layer.name} points to {save_path}")


def update_point_visibility(layer, viewer):
    if "frames" not in layer.metadata:
        return

    current_frame = viewer.dims.current_step[0]
    frames = layer.metadata["frames"]

    # Update visibility based on current frame
    visible = frames <= current_frame
    layer.metadata["visible"] = visible

    # Update the displayed points
    if len(layer.data) > 0:
        layer.shown = visible


def store_frame_metadata(layer, viewer):
    current_frame = viewer.dims.current_step[0]  # Get current frame number
    num_points = len(layer.data)  # Total points in the layer

    if "frames" not in layer.metadata:
        layer.metadata["frames"] = np.zeros(num_points, dtype=int)  # Initialize
        layer.metadata["visible"] = np.ones(
            num_points, dtype=bool
        )  # Initialize visibility

    existing_frames = layer.metadata["frames"]  # Get stored frame numbers

    # If new points were added, store their frame numbers
    if num_points > len(existing_frames):
        new_frames = np.full(
            num_points - len(existing_frames), current_frame
        )  # Assign current frame to new points
        layer.metadata["frames"] = np.concatenate(
            [existing_frames, new_frames]
        )  # Update frames
        # New points should be visible if current frame >= their frame
        new_visible = np.ones(num_points - len(existing_frames), dtype=bool)
        layer.metadata["visible"] = np.concatenate(
            [layer.metadata["visible"], new_visible]
        )

    # Update visibility for all points
    update_point_visibility(layer, viewer)


def annotate_video(video_path):
    vr = VideoReaderNP(video_path)
    viewer = napari.view_image(vr, name=video_path)

    # Add point layers
    rewards_layer = viewer.add_points(name="rewards", face_color="red", size=10)
    start_layer = viewer.add_points(name="start", face_color="green", size=10)
    trials_layer = viewer.add_points(name="trials", face_color="blue", size=10)
    diameter_layer = viewer.add_points(name="diameter", face_color="yellow", size=10)

    point_layers = [rewards_layer, start_layer, trials_layer, diameter_layer]

    # Attach event listeners
    for layer in point_layers:
        # Update frame metadata when points are added
        layer.events.data.connect(
            lambda event, l=layer: store_frame_metadata(l, viewer)
        )

        # Update visibility when frame changes
        viewer.dims.events.current_step.connect(
            lambda event: update_point_visibility(layer, viewer)
        )

    # Bind save function to 'S' key
    @viewer.bind_key("Shift-S")
    def save_on_keypress(viewer):
        save_layers_to_video_directory(viewer, video_path, point_layers)

    napari.run()

def verify_manual_annotation(video_path, fs=None, trial_wiggle_room=10):
    if not is_annotated(video_path):
        print(f"{video_path} --- not annotated")
        return
    
    if fs is None:
        # load video and get fps 
        cap = cv2.VideoCapture(video_path)
        fs = cap.get(cv2.CAP_PROP_FPS)
        cap.release()

    video_path_obj = Path(video_path)
    video_dir = video_path_obj.parent
    video_name = video_path_obj.stem  # Extracts the filename without extension

    # Read the updated CSV naming schemes
    trials = pd.read_csv(video_dir / f"{video_name}_trials.csv")
    rewards = pd.read_csv(video_dir / f"{video_name}_rewards.csv")
    start = pd.read_csv(video_dir / f"{video_name}_start.csv")
    diameter = pd.read_csv(video_dir / f"{video_name}_diameter.csv")

    # check if any trial durations
    trials = trials.sort_values("axis-0").reset_index(drop=True)
    start_ind = trials[trials.index % 2 == 0]["axis-0"]
    stop_ind = trials[trials.index % 2 == 1]["axis-0"]
    starts = start_ind / fs
    stops = stop_ind / fs
    
    # check if same number of starts and stops
    assert len(starts) == len(
        stops
    ), f"Found {len(starts)} starts and {len(stops)} stops"
    
    # get durations per trial
    durations = stops.values - starts.values

    # check if any trial is longer than x seconds
    assert durations.max() < 60 * 6, f"Found {durations.max()} seconds"
    if durations.max() > 120:
        print(f"{video_path} --- Found {durations.max()} seconds")

    # check if any trial is shorter than 5 seconds
    if durations.min() < 5:
        test = 1
    assert durations.min() > 2, f"Found {durations.min()} seconds"

    # check if more than 30 trials
    assert durations.shape[0] <= 35, f"Found {durations.shape[0]} trials"

    # check if at least 2 rewards
    assert rewards.shape[0] >= 2, f"Found {rewards.shape[0]} rewards"

    # check if 1 start
    assert start.shape[0] == 1, f"Found {start.shape[0]} start points"

    # check if 1 diameter
    assert diameter.shape[0] == 2, f"Found {diameter.shape[0]} diameter points"

    print(f"{video_path} --- verified")


def is_annotated(video_path):
    video_path_obj = Path(video_path)
    video_dir = video_path_obj.parent
    video_name = video_path_obj.stem

    return (
        os.path.exists(video_dir / f"{video_name}_rewards.csv")
        and os.path.exists(video_dir / f"{video_name}_start.csv")
        and os.path.exists(video_dir / f"{video_name}_trials.csv")
        and os.path.exists(video_dir / f"{video_name}_diameter.csv")
    )


def which_is_missing(video_path):
    video_path_obj = Path(video_path)
    video_dir = video_path_obj.parent
    video_name = video_path_obj.stem
    
    missing = []
    
    if not os.path.exists(video_dir / f"{video_name}_rewards.csv"):
        missing.append(f"{video_name}_rewards.csv")
    if not os.path.exists(video_dir / f"{video_name}_start.csv"):
        missing.append(f"{video_name}_start.csv")
    if not os.path.exists(video_dir / f"{video_name}_trials.csv"):
        missing.append(f"{video_name}_trials.csv")
    if not os.path.exists(video_dir / f"{video_name}_diameter.csv"):
        missing.append(f"{video_name}_diameter.csv")
        
    return missing

def select_csv_file():
    """Opens a file dialog for the user to select a CSV file."""
    # Hide the main tkinter root window
    root = tk.Tk()
    root.withdraw()
    
    # Make the popup appear on top of other windows
    root.attributes('-topmost', True)

    # Open the file dialog
    file_path = filedialog.askopenfilename(
        title="Select the CSV file containing basepaths",
        filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
    )
    
    return file_path

if __name__ == "__main__":
    # Prompt user to select the CSV file
    print("Please select your CSV file in the popup window...")
    csv_path = select_csv_file()

    # Check if the user actually selected a file or just closed the window
    if not csv_path:
        print("No file selected. Exiting script.")
        exit()

    print(f"Loading data from: {csv_path}")
    print("discovering files...")

    # Load the unique basepaths from the user-selected CSV
    basepaths = pd.read_csv(csv_path).basepath.unique()
    
    files = []

    for folder in basepaths:
        files.extend(glob.glob(os.path.join(folder, "**", "*.avi"), recursive=True))

    # print("filter to just probe sessions to start")
    files_series = pd.Series(files)
    idx = ~files_series.str.contains("backup") & ~files_series.str.contains("test")
    files = list(compress(files, idx))

failed = []

for video_path_i, video_path in enumerate(files):
    print(f"{video_path} --- {video_path_i} of {len(files)} files")

    video_dir = Path(video_path).parent

    if is_annotated(video_path):
        try:
            verify_manual_annotation(video_path)
        except Exception as e:
            print(f"  !! ERROR ({type(e).__name__}): {e}")
            failed.append(video_path)
        continue

    missing = which_is_missing(video_path)
    if missing:
        for m in missing:
            print(f"Missing {m}")
    annotate_video(video_path)
    try:
        verify_manual_annotation(video_path)
    except Exception as e:
        print(f"  !! ERROR ({type(e).__name__}): {e}")
        failed.append(video_path)

if failed:
    print(f"\n{len(failed)} videos failed — fix them then re-verifying...")
    input("Press Enter when ready to re-verify failed videos...")
    for video_path in failed:
        try:
            verify_manual_annotation(video_path)
            print(f"  ✓ fixed: {video_path}")
        except Exception as e:
            print(f"  !! still failing ({type(e).__name__}): {e}")
