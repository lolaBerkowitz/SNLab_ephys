import glob
import os
from itertools import compress
from pathlib import Path

import napari
import numpy as np
import pandas as pd
from napari_video.napari_video import VideoReaderNP


# Goes through all the layers of the video, extracting the coordinates,
# frames, and the visible data and creates a DataFrame with the according
# coordinates for each frame.
#
def save_layers_to_video_directory(viewer, video_path, point_layers):
    video_dir = Path(video_path).parent

    #A layer is a level at which an aspect of the video can be manipulated (i.e.,
    #can contain videos, audio, images, text, or effects). Layers are used to
    #superimpose on video clip over another or to add special efects.
    for layer in point_layers:
        # Extract point coordinates (N, 2)
        points = layer.data
        # Loads the frame data and, if not found, replaces with 0s
        frames = layer.metadata.get("frames", np.zeros(len(points), dtype=int))
        # Loads the visible data and, if not found, replaces with 1s (always True)
        visible = layer.metadata.get("visible", np.ones(len(points), dtype=bool))

        if len(points) == 0:
            print(f"No data to save for layer '{layer.name}'.")
            continue  # Skip empty layers

        # Create DataFrame with required columns
        df = pd.DataFrame(
            {
                "index": np.arange(len(points)),  # Index of each point
                "axis-0": frames,  # Frame number
                "axis-1": points[:, 0],  # X-coordinates -> first column of array
                "axis-2": points[:, 1],  # Y-coordinates -> second column of array
            }
        )

        save_path = video_dir / f"{layer.name}.csv"
        df.to_csv(save_path, index=False)
        print(f"Saved {layer.name} points to {save_path}")

#
def update_point_visibility(layer, viewer):
    if "frames" not in layer.metadata:
        return

    # Gets the current frame from the video, starting from where the viewer is
    current_frame = viewer.dims.current_step[0]
    # Gets the metadata associated with frames
    frames = layer.metadata["frames"]

    # Update visibility based on current -> so visibility is true for the frames
    # before and at the frame the viewer is on
    visible = frames <= current_frame
    layer.metadata["visible"] = visible

    # Update the displayed points -> points that are displayed are the ones
    # that are listed as visible
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


# Goes through each layer and modifies data based on certain events, saving the
# data to the layer when press 'S' key
def annotate_video(video_path):
    vr = VideoReaderNP(video_path)
    viewer = napari.view_image(vr, name=video_path)

    # Add point layers
    rewards_layer = viewer.add_points(name="rewards", face_color="red", size=10)
    start_layer = viewer.add_points(name="start", face_color="green", size=10)
    trials_layer = viewer.add_points(name="trials", face_color="blue", size=10)

    point_layers = [rewards_layer, start_layer, trials_layer]

    # Attach event listeners -> so whenever the data changes (i.e., points
    # are removed, added, or modified) or the viewer changes a frame,
    # updates metadata and changes visuals (i.e., visible points), respectively
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

# Checks trial durations, number of trials, number of rewards, and each trial
# has only 1 start time
def verify_manual_annotation(video_path, fs=40, trial_wiggle_room=10):
    if not is_annotated(video_path):
        print(f"{video_path} --- not annotated")
        return

    video_dir = Path(video_path).parent

    trials = pd.read_csv(os.path.join(video_dir, "trials.csv"))
    rewards = pd.read_csv(os.path.join(video_dir, "rewards.csv"))
    start = pd.read_csv(os.path.join(video_dir, "start.csv"))

    # check if any trial durations
    trials = trials.sort_values("axis-0").reset_index(drop=True)
    # fills the even indexes of trials with start times
    start_ind = trials[trials.index % 2 == 0]["axis-0"]
    # fills the odd indexes of trials with end times
    stop_ind = trials[trials.index % 2 == 1]["axis-0"]
    # converts the indices into seconds by dividing by sampling frequency
    starts = start_ind / fs
    stops = stop_ind / fs
    # check if same number of starts and stops -> sanity check ig?
    assert len(starts) == len(stops), (
        f"Found {len(starts)} starts and {len(stops)} stops"
    )
    # get durations per trial
    durations = stops.values - starts.values

    # check if any trial is longer than 360 seconds
    assert durations.max() < 60 * 6, f"Found {durations.max()} seconds"
    if durations.max() > 90:
        print(f"{video_path} --- Found {durations.max()} seconds")

    # check if any trial is shorter than 90 seconds
    if durations.min() < 90:
        test = 1
    assert durations.min() > 5, f"Found {durations.min()} seconds"

    # check if more than 25 trials
    assert durations.shape[0] <= 30, f"Found {durations.shape[0]} trials"

    # check if at least 2 rewards
    assert rewards.shape[0] >= 2, f"Found {rewards.shape[0]} rewards"

    # check if 1 start
    assert start.shape[0] == 1, f"Found {start.shape[0]} start points"

    print(f"{video_path} --- verified")


def is_annotated(video_path):
    video_dir = Path(video_path).parent
    print("hey")
    if os.path.exists(os.path.join(video_dir, "rewards.csv")) & os.path.exists(os.path.join(video_dir, "start.csv")) & os.path.exists(os.path.join(video_dir, "trials.csv")):
        print("treldkjf")
    else:
        print("suckjerja;ldkjf;alskjdfa;lsdkjf")
    return (
        os.path.exists(os.path.join(video_dir, "rewards.csv"))
        & os.path.exists(os.path.join(video_dir, "start.csv"))
        & os.path.exists(os.path.join(video_dir, "trials.csv"))
    )


if __name__ == "__main__":
    print("discovering files...")

    # files = glob.glob(r"U:\data\hpc_ctx_project\**\*.avi", recursive=True)

    # files_series = pd.Series(files)
    # idx = (
    #     files_series.str.contains("cheeseboard|cheesboard|open_field|acquisition|probe")
    #     & ~files_series.str.contains("backup")
    #     & ~files_series.str.contains("test")
    # )

    # files = list(compress(files, idx))
    basepaths = pd.read_csv("/Users/shuwend/Desktop/SN_Research/sample_data.csv").basepath.unique()
    files = []

    for folder in basepaths:
        print('hei')
        files.extend(glob.glob(os.path.join(folder, "**", "*.avi"), recursive=True))

    print(f"Initial file count: {len(files)}")
    print(files)

    print("filter to just probe sessions to start")
    files_series = pd.Series(files)
    idx = ~files_series.str.contains("backup") & ~files_series.str.contains("test")
    files = list(compress(files, idx))
    print(files)
    for video_path_i, video_path in enumerate(files):
        print(f"{video_path} --- {video_path_i} of {len(files)} files")

        video_dir = Path(video_path).parent

        if is_annotated(video_path):
            verify_manual_annotation(video_path)
            continue

        annotate_video(video_path)
        verify_manual_annotation(video_path)
