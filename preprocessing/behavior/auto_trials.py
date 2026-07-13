import glob
import nelpy as nel
import neuro_py as npy
import numpy as np
import pandas as pd
import os
import scipy.io as sio
from matplotlib import patches as patches
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema
from neuro_py.process.intervals import find_interval  # type: ignore
from neuro_py.io import loading  # type: ignore


def get_trial_intervals(
    positions: pd.DataFrame,
    dist_from_start: int = 30,
    min_trial_duration: int = 2,
    inter_trial_gap: int = 15,
):
    """
    Identify trial intervals based on the animal's distance from the start position
    Source: Ryan's hpc_ctx\behavior\cheeseboard.py

    ----------
    Inputs:
        positions : pd.DataFrame
            DataFrame containing the x,y positions and time of the rat
        dist_from_start : int, optional
            Distance from the start position to consider the rat in the cheeseboard in cm, by default 30
        min_trial_duration : int, optional
            Minimum duration of a trial in seconds, by default 2
        inter_trial_gap : int, optional
            Maximum gap between trials in seconds (how fast can you reset maze?), by default 15
    ----------
    Outputs:
        trial_intervals : nel.EpochArray
            Intervals when the rat is in the cheeseboard (raw estimate of trial intervals)
    """
    # find duration that rat is away from start
    intervals = np.array(
        find_interval(positions.distance_to_start_smooth.values > dist_from_start)
    )
    if intervals.size == 0:
        return nel.EpochArray(np.empty((0, 2)))
    trial_intervals = positions.time.values[intervals.astype(int)]

    # put trial intervals into nelpy format and remove trials that are too short
    trial_intervals = nel.EpochArray(trial_intervals)
    trial_intervals = nel.EpochArray(
        trial_intervals.data[trial_intervals.lengths > min_trial_duration, :]
    )

    # merge trials that are less than 15 seconds apart
    trial_intervals = trial_intervals.merge(gap=inter_trial_gap)

    return trial_intervals


def reward_visit_time(
    reward: np.ndarray,
    animal_positions: pd.DataFrame,
    distance_cap: float = 8.0,
    min_duration: float = 0.01,
    show_plot: bool = False,
):
    """
    Find the time interval of the animal visits a reward.
    The animal is considered visiting the reward when it is spending the longest
    duration within a certain distance to the reward in a trial.
    ----------
    Inputs:
        reward: np.ndarray
            coordinate of a given reward location
        animal_positions: pd.DataFrame
            DataFrame containing the x,y positions and time of the rat
        distance_cap: float, optional
            limit of how close the animal should be to the reward for it to be considered as a valid visit in cm, by default 8.0
        min_duration: float, optional
            minimal duration of the animal needs to stay around the reward in seconds, by default 0.01
        show_plot: bool, optional
            show animal distance to reward (raw & smoothed) vs. time, by default False
    ----------
    Outputs:
        visit_time: nel.EpochArray
            time interval of animal visiting the reward
    """
    # compute distance locally without mutating input DataFrame
    dist_raw = np.linalg.norm(
        animal_positions[["x", "y"]].values - reward, axis=1
    ).astype(float)

    dist_series = pd.Series(dist_raw).interpolate(method="linear")
    dist_smooth = dist_series.rolling(window=40, min_periods=1).mean().values

    if show_plot:
        plt.figure(figsize=(12, 6))
        plt.plot(animal_positions["time"].values, dist_raw, alpha=0.4, label="Raw")
        plt.plot(animal_positions["time"].values, dist_smooth, label="Smoothed")
        plt.axhline(y=distance_cap, color="red", linestyle="--")
        plt.xlabel("Time (s)")
        plt.ylabel("Distance to reward (cm)")
        plt.legend()
        plt.show()

    minima = argrelextrema(dist_smooth, np.less, order=10)[0]

    if len(minima) == 0:
        return nel.EpochArray(np.empty((0, 2)))

    visit_time = None
    max_duration = 0
    for idx in minima:
        if dist_smooth[idx] > distance_cap:
            continue

        left, right = idx, idx
        while left > 0 and dist_smooth[left] <= distance_cap:
            left -= 1
        while right < len(dist_smooth) - 1 and dist_smooth[right] <= distance_cap:
            right += 1

        duration = (
            animal_positions["time"].iloc[right] - animal_positions["time"].iloc[left]
        )
        if duration > max_duration:
            max_duration = duration
            visit_time = (left, right)

    if visit_time is None or max_duration < min_duration:
        return nel.EpochArray(np.empty((0, 2)))

    start = animal_positions["time"].iloc[visit_time[0]]
    stop = animal_positions["time"].iloc[visit_time[1]]
    return nel.EpochArray(np.array([[start, stop]]))


def find_active_trials_cheeseboard(
    current_positions: pd.DataFrame,
    rewards: np.ndarray,
    start: np.ndarray,
    dist_from_start: int = 30,
    min_trial_duration: int = 5,
    inter_trial_gap: int = 15,
    max_trial_duration: int = 120,
):
    """
    Find active trial intervals for cheeseboard task.

    ----------
    Inputs:
        current_positions: pd.DataFrame
            DataFrame containing the x,y positions and time of the rat
        rewards: np.ndarray
            Array of reward coordinates, shape (n_rewards, 2)
        start: np.ndarray
            Start position coordinate, shape (2,)
        dist_from_start: int, optional
            Distance from start to consider animal on cheeseboard in cm, by default 30
        min_trial_duration: int, optional
            Minimum trial duration in seconds, by default 5
        inter_trial_gap: int, optional
            Maximum gap to merge adjacent trials in seconds, by default 15
        max_trial_duration: int, optional
            Maximum allowed trial duration in seconds, by default 120
    ----------
    Outputs:
        trials_df: pd.DataFrame
            DataFrame with columns:
                trial_start     : float, trial start time (s)
                trial_stop      : float, trial stop time (s) — time animal reached last reward,
                                  or capped at trial_start + max_trial_duration
                duration        : float, trial duration (s)
                n_rewards_visited: int, number of rewards visited
                reward_visits   : np.ndarray of shape (n_visited, 2), [start, stop] per reward visit
    """
    current_positions = current_positions.copy()
    rewards = np.atleast_2d(rewards)
    start = np.array(start).flatten()

    current_positions["distance_to_start"] = np.linalg.norm(
        current_positions[["x", "y"]].values - start, axis=1
    )
    current_positions["distance_to_start"] = current_positions[
        "distance_to_start"
    ].interpolate(method="linear")
    current_positions["distance_to_start_smooth"] = (
        current_positions["distance_to_start"]
        .rolling(window=200, center=True)
        .mean()
        .copy()
    )

    if current_positions["distance_to_start_smooth"].isna().all():
        return pd.DataFrame()

    trial_intervals = get_trial_intervals(
        current_positions, dist_from_start, min_trial_duration, inter_trial_gap
    )

    records = []
    for i, trial in enumerate(trial_intervals):
        trial_idx = current_positions.time.between(trial.starts[0], trial.stops[0])
        trial_positions = current_positions.loc[trial_idx].copy()

        # cap trial duration
        capped = trial.duration >= max_trial_duration

        # find reward visit time for each reward
        reward_visit_epochs = []
        not_visited = 0
        rewards_intervals = nel.EpochArray(np.empty((0, 2)))
        for reward in rewards:
            visit_time = reward_visit_time(reward, trial_positions)
            if visit_time.isempty:
                not_visited += 1
                continue
            rewards_intervals += visit_time
            reward_visit_epochs.append(
                [float(visit_time.starts[0]), float(visit_time.stops[0])]
            )

        # skip trials where no rewards were visited
        if not_visited == len(rewards):
            continue

        if capped:
            stop_time = trial.starts[0] + max_trial_duration
        else:
            stop_time = float(rewards_intervals.starts.max())

        records.append(
            {
                "trial_start": float(trial.starts[0]),
                "trial_stop": stop_time,
                "duration": stop_time - float(trial.starts[0]),
                "n_rewards_visited": len(rewards) - not_visited,
                "reward_visits": np.array(
                    reward_visit_epochs
                ),  # TODO: consider converting to average time spent at reward site
            }
        )

    return pd.DataFrame(records)


def run(
    basepath: str,
    dist_from_start: int = 30,
    min_trial_duration: int = 5,
    inter_trial_gap: int = 15,
    max_trial_duration: int = 120,
):
    """
    Run auto trial detection for all cheeseboard epochs in a session.

    ----------
    Inputs:
        basepath: str
            Path to the session directory
        dist_from_start: int, optional
            Distance from start to consider animal on cheeseboard in cm, by default 30
        min_trial_duration: int, optional
            Minimum trial duration in seconds, by default 5
        inter_trial_gap: int, optional
            Maximum gap to merge adjacent trials in seconds, by default 15
        max_trial_duration: int, optional
            Maximum allowed trial duration in seconds, by default 120
    ----------
    Outputs:
        df: pd.DataFrame
            DataFrame with one row per trial, columns:
                trial_start, trial_stop, duration, n_rewards_visited,
                reward_visits, basepath, epoch_name, task_manipulation, task_notes
    """

    # load session epochs from session file
    epoch_df = loading.load_epoch(basepath)
    # make string columns string
    epoch_df["environment"] = epoch_df["environment"].astype(str)
    epoch_df["manipulation"] = epoch_df["manipulation"].astype(str)
    epoch_df["notes"] = epoch_df["notes"].astype(str)

    # load behavior file
    filename = os.path.join(
        basepath, os.path.basename(basepath) + ".session.mat"
    )
    data = sio.loadmat(filename, simplify_cells=True)

    # load position from behavior file
    positions = loading.load_animal_behavior(basepath)

    # loading epochs from behavior file
    beh_epoch = data["behavior"]["epochs"]
    if isinstance(beh_epoch, dict):
        beh_epoch = [beh_epoch]
    beh_epoch_df = pd.DataFrame(beh_epoch)

    if "environment" not in beh_epoch_df.columns:
        beh_epoch_df["environment"] = epoch_df["environment"].values

    if beh_epoch_df.shape[0] != epoch_df.shape[0]:
        raise ValueError("Epochs do not match")

    df = pd.DataFrame()
    for cheeseboard_ind in epoch_df.query("environment=='circular_field'").index:
        current_idx = positions.time.between(
            epoch_df.loc[cheeseboard_ind].startTime,
            epoch_df.loc[cheeseboard_ind].stopTime,
        )
        current_positions = positions.loc[current_idx]

        # TODO: Check and see how to get 'rewards in the behavior file epochs
        if isinstance(data['session']['behavioralTracking'], dict):
            video_name = data['session']['behavioralTracking']['notes'].split('.')[0]
            assert data['session']['behavioralTracking']['epoch'] - 1 == cheeseboard_ind
            # TODO load the rewards csv that also contains video_name
            # TODO load the start 
        elif isinstance(data['session']['behavioralTracking'], list):
            assert data['session']['behavioralTracking'][cheeseboard_ind]['epoch'] - 1 == cheeseboard_ind
            video_name = data['session']['behavioralTracking'][cheeseboard_ind]['notes'].split('.')[0]
            # TODO load the rewards csv that also contains video_name
            # TODO load the start 

            # skip if there is no reward csv & print a warning that no rewards / start exist

            # format rewards (n_rewards,2) and start as ndarray (2,) # xy coordinates in cm


            # # MUST HAVE REWARDS FOR CODE TO COMPLETE
            # # skip epochs without rewards
            # if (
            #     "rewards" not in beh_epoch_df.columns
            #     or np.isnan(beh_epoch_df.iloc[cheeseboard_ind]["rewards"]).all()
            # ):
            #     continue

            # # format rewards
            # rewards = np.atleast_2d(beh_epoch_df.iloc[cheeseboard_ind]["rewards"])
            # start = np.array(beh_epoch_df.iloc[cheeseboard_ind]["start"]).reshape(-1, 2)[0]

        # TEMP TO GET RID OF WARNINGS
        rewards = np.empty((0,2))
        start = np.empty((0,2))

        trials_df = find_active_trials_cheeseboard(
            current_positions,
            rewards,
            start,
            dist_from_start=dist_from_start,
            min_trial_duration=min_trial_duration,
            inter_trial_gap=inter_trial_gap,
            max_trial_duration=max_trial_duration,
        )

        if trials_df.empty:
            continue

        # add metadata to dataframe
        trials_df["basepath"] = basepath
        trials_df["epoch_name"] = epoch_df.iloc[cheeseboard_ind]["name"]
        trials_df["task_manipulation"] = epoch_df.iloc[cheeseboard_ind].manipulation
        trials_df["task_notes"] = str(epoch_df.iloc[cheeseboard_ind].notes)

        df = pd.concat([df, trials_df], axis=0, ignore_index=True)

    return df


def load_results(save_path: str, verbose: bool = False) -> pd.DataFrame:
    """
    Load and concatenate auto trial results across all sessions in save_path.

    ----------
    Inputs:
        save_path: str
            Directory containing per-session pickle files produced by run()
        verbose: bool, optional
            Print each loaded session path, by default False
    ----------
    Outputs:
        results: pd.DataFrame
            Concatenated trial DataFrame across all sessions
    """
    sessions = glob.glob(save_path + os.sep + "*.pkl")
    if not sessions:
        raise FileNotFoundError(f"No pickle files found in {save_path}")

    dfs = []
    for session_file in sessions:
        try:
            df = npy.process.batch_analysis.load_specific_data(session_file)
            if verbose:
                print(f"Loaded: {session_file}")
            dfs.append(df)
        except Exception as e:
            print(f"Failed to load {session_file}: {e}")

    return pd.concat(dfs, axis=0, ignore_index=True)


def inspect_trials(
    basepath: str,
    ncols: int = 5,
):
    """
    Plot animal trajectory per trial for all cheeseboard epochs in a session.

    ----------
    Inputs:
        basepath: str
            Path to the session directory
        ncols: int, optional
            Number of columns in the subplot grid, by default 5
    ----------
    Outputs:
        figs: list of matplotlib.figure.Figure
            One figure per cheeseboard epoch
    """
    # TODO: Update the reward loading like you do in the run function

    # load data
    filename = os.path.join(
        basepath, os.path.basename(basepath) + ".animal.behavior.mat"
    )

    data = sio.loadmat(filename, simplify_cells=True)
    epoch_df = loading.load_epoch(basepath)
    positions = loading.load_animal_behavior(basepath)

    beh_epoch = data["behavior"]["epochs"]
    if isinstance(beh_epoch, dict):
        beh_epoch = [beh_epoch]
    beh_epoch_df = pd.DataFrame(beh_epoch)

    if "environment" not in beh_epoch_df.columns:
        beh_epoch_df["environment"] = epoch_df["environment"].values

    figs = []
    for cheeseboard_ind in epoch_df.query("environment=='circular_field'").index:
        current_idx = positions.time.between(
            epoch_df.loc[cheeseboard_ind].startTime,
            epoch_df.loc[cheeseboard_ind].stopTime,
        )
        current_positions = positions.loc[current_idx]
        rewards = np.atleast_2d(beh_epoch_df.iloc[cheeseboard_ind]["rewards"])
        start = np.array(beh_epoch_df.iloc[cheeseboard_ind]["start"]).reshape(-1, 2)[0]

        trials_df = find_active_trials_cheeseboard(current_positions, rewards, start)

        if trials_df.empty:
            continue

        n = len(trials_df)
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3, nrows * 3))
        fig.suptitle(epoch_df.loc[cheeseboard_ind]["name"], fontsize=11)
        axes = axes.flatten()

        for i, (_, trial) in enumerate(trials_df.iterrows()):
            ax = axes[i]
            mask = current_positions.time.between(trial.trial_start, trial.trial_stop)
            traj = current_positions.loc[mask]

            ax.plot(traj.x, traj.y, linewidth=1, color="steelblue", alpha=0.8)
            ax.scatter(rewards[:, 0], rewards[:, 1], color="orange", zorder=10, s=60)
            ax.scatter(
                start[0], start[1], color="green", s=80, zorder=11, edgecolors="black"
            )
            ax.set_title(f"Trial {i} | {trial.duration:.1f}s", fontsize=9)
            ax.set_aspect("equal")
            ax.axis("off")

        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        plt.tight_layout()
        figs.append(fig)

    return figs
