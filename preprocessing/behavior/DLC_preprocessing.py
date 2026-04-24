import os
from typing import Optional, Tuple, Union

import glob  # for file names
import pickle
import re
import shutil

import cv2
import numpy as np
import pandas as pd  # to create data frames
import ruamel.yaml
import seaborn as sns
from matplotlib import pyplot as plt  # visualization
from matplotlib.figure import Figure

def strip_dlc_suffix(key: str) -> str:
    """Remove trailing DLC suffix from a key string.

    Parameters
    ----------
    key : str
        Input key containing a ``DLC_`` token.

    Returns
    -------
    str
        Key truncated before ``DLC_`` with trailing underscores removed.
    """
    return key.split("DLC_")[0].rstrip("_")


def update_config_cropping(
    config_path: str,
    xmin: Union[int, float],
    xmax: Union[int, float],
    ymin: Union[int, float],
    ymax: Union[int, float],
) -> Optional[dict]:
    """Update DeepLabCut config cropping parameters on disk.

    Parameters
    ----------
    config_path : str
        Path to the DLC ``config.yaml`` file.
    xmin : int or float
        Left crop bound (x1).
    xmax : int or float
        Right crop bound (x2).
    ymin : int or float
        Top crop bound (y1).
    ymax : int or float
        Bottom crop bound (y2).

    Returns
    -------
    dict or None
        Existing config mapping when no update is needed, otherwise ``None``
        after writing updated values.

    Raises
    ------
    ValueError
        If the YAML file is empty or malformed.
    """

    yaml = ruamel.yaml.YAML()
    yaml.preserve_quotes = True  # keeps formatting stable if desired

    # ---- Read existing config ----
    with open(config_path, "r") as fp:
        data = yaml.load(fp)

    if data is None:
        raise ValueError("YAML file is empty or malformed.")

    # ---- Early exit if already set ----
    if (
        data.get("x1") == xmin
        and data.get("x2") == xmax
        and data.get("y1") == ymin
        and data.get("y2") == ymax
        and data.get("cropping") is True
    ):
        print(
            f"Cropping parameters already set to ({xmin}, {xmax}, {ymin}, {ymax}). No update needed."
        )
        return data

    # ---- Update cropping parameters ----
    data["cropping"] = True
    data["x1"] = int(xmin)
    data["x2"] = int(xmax)
    data["y1"] = int(ymin)
    data["y2"] = int(ymax)

    # ---- Write safely ----
    with open(config_path, "w") as f:
        yaml.dump(data, f)
        f.flush()
        os.fsync(f.fileno())  # ensure write completes to disk


def percent_good(df: pd.DataFrame, cutoff: float = 0.80) -> float:
    """Compute proportion of frames above likelihood cutoff.

    Parameters
    ----------
    df : pandas.DataFrame
        DLC body-part dataframe containing a ``likelihood`` column.
    cutoff : float, optional
        Likelihood threshold, by default 0.80.

    Returns
    -------
    float
        Fraction of likelihood values greater than ``cutoff``.
    """
    lh = df["likelihood"]
    return np.mean(lh > cutoff)


def plot_dlc(
    df: pd.DataFrame, varname: str, cutoff: float = 0.90, fs: float = 30
) -> Figure:
    """Plot coordinate quality diagnostics for a single body part.

    Parameters
    ----------
    df : pandas.DataFrame
        DLC dataframe for one body part with ``x``, ``y``, and ``likelihood``.
    varname : str
        Body-part name used in figure title.
    cutoff : float, optional
        Likelihood threshold for "good" points, by default 0.90.
    fs : float, optional
        Sampling frequency in Hz, by default 30.

    Returns
    -------
    matplotlib.figure.Figure
        Figure with trajectory and time-series quality panels.
    """
    # Index coordinate and likelihood values.
    x = df["x"]
    y = df["y"]
    lh = df["likelihood"]

    # Create timestamps to align coordinates over time.
    ts = np.linspace(0, len(x) / fs, len(x))

    # Easy index for all "good" coordinates.
    x_g = x.loc[lh > cutoff]
    y_g = y.loc[lh > cutoff]
    ts_g = ts[lh > cutoff]

    # Set up figure and subplots.
    fig, ax = plt.subplots(4, 1, figsize=(15.0, 15.0))

    # First plot, shows all tracked coordinates
    ax[0].scatter(x, y, s=0.05)
    ax[0].set_aspect("equal", "box")
    ax[0].set_title(f"{cutoff * 100} th percentile likelihood for {varname}")

    # Second plot, these are only coordinates that pass cutoff
    ax[1].scatter(x_g, y_g, s=0.05)
    ax[1].set_aspect("equal", "box")
    ax[1].set_title(f"{cutoff * 100} th percentile likelihood for {varname}")

    # Third plot, shows x coords over time
    ax[2].plot(ts, x, color="r", label="bad")
    ax[2].plot(ts_g, x_g, label="good")
    ax[2].set_title("All X coords over time")
    ax[2].legend()

    # Fourth plot, shows y coords over time
    ax[3].plot(ts, y, color="r", label="bad")
    ax[3].plot(ts_g, y_g, label="good")
    ax[3].set_title("All Y coords over time")
    ax[3].legend()

    plt.tight_layout()

    return fig


def plot_keypoints(
    file_name: str,
    cutoff: float = 0.80,
    savefig: bool = False,
    save_path: Optional[str] = None,
) -> Figure:
    """Plot all keypoint trajectories for a DLC CSV file.

    Parameters
    ----------
    file_name : str
        Path to DLC CSV file.
    cutoff : float, optional
        Likelihood threshold used for restricted trajectories, by default 0.80.
    savefig : bool, optional
        Whether to save the generated figure, by default False.
    save_path : str or None, optional
        Destination directory for saving figure. If ``None``, uses CSV parent
        directory.

    Returns
    -------
    matplotlib.figure.Figure
        Figure comparing likelihood-restricted and raw trajectories.
    """

    df = pd.read_csv(file_name, header=[1, 2])
    header = df.columns.get_level_values(0).drop_duplicates().to_numpy()
    vidname = os.path.basename(file_name)
    body_parts = header[1:]

    fig, ax = plt.subplots(
        2,
        len(body_parts),
        figsize=set_size("paper", 2, subplots=(5, len(body_parts))),
        sharex=True,
        sharey=True,
        dpi=150,
    )
    ax = ax.ravel()

    for i, var in enumerate(body_parts):
        x = df[var]["x"]
        y = df[var]["y"]
        lh = df[var]["likelihood"]

        # Easy index for all "good" coordinates.
        x_g = x.loc[lh > cutoff]
        y_g = y.loc[lh > cutoff]

        ax[i].plot(
            x_g,
            y_g,
            label="likelihood restricted",
            color="red",
            alpha=0.5,
            linewidth=0.5,
        )
        ax[i + len(body_parts)].plot(
            x,
            y,
            label="raw",
            color="k",
            alpha=0.5,
            linewidth=0.5,
        )
        percent = percent_good(df[var], cutoff)
        ax[i].set_title(f"{var} ({percent:.2f})", fontsize=8)
        ax[i].set_aspect("equal")
        sns.despine()

    fig.suptitle(vidname, fontsize=8)
    ax[len(body_parts)].legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    ax[len(body_parts) + len(body_parts) - 1].legend(
        bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8
    )
    if save_path is None:
        save_path = os.path.dirname(file_name)

    if savefig:
        plt.savefig(
            os.path.join(
                save_path,
                f"{vidname}.png",
            ),
            bbox_inches="tight",
        )

    return fig


def plot_all_checks(file_name: str, cutoff: float = 0.80) -> pd.DataFrame:
    """Generate per-body-part diagnostic plots and summary table.

    Parameters
    ----------
    file_name : str
        Path to DLC CSV file.
    cutoff : float, optional
        Likelihood threshold used for summary metrics, by default 0.80.

    Returns
    -------
    pandas.DataFrame
        Summary dataframe with proportion of good samples per body part.
    """

    df = pd.read_csv(file_name, header=[1, 2])
    header = df.columns.get_level_values(0).drop_duplicates().to_numpy()
    vidname = os.path.basename(file_name)
    body_parts = header[1:]
    good_prop = []

    for var in body_parts:
        fig = plot_dlc(df[var], var, cutoff=cutoff)
        plt.title(vidname)
        # fig.savefig(file_name.replace('.csv','')+'_'+var+'.png')
        plt.show()
        good_prop.append(percent_good(df[var], cutoff=cutoff))

    summary = pd.DataFrame(columns=["proportion_good", "body_part"])
    summary["proportion_good"] = good_prop
    summary["body_part"] = body_parts

    return summary


# Backward-compatible alias.
plot_DLC = plot_dlc


def set_size(
    width: Union[float, str] = "double_col",
    fraction: float = 1,
    subplots: Tuple[int, int] = (1, 1),
    ratio: Optional[float] = None,  # override golden ratio
) -> Tuple[float, float]:
    """
    Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width : float or str
        Document width in points (float) or predefined document type (str).
        Supported types: see WIDTHS dictionary for all presets (nature_single, nature_double,
        science_single, science_double, science_triple, cell_single, cell_1p5, cell_double,
        single_col, double_col, beamer, thesis, textwidth, paper).
    fraction : float, optional
        Fraction of the width which you wish the figure to occupy, by default 1.
    subplots : tuple of int, optional
        Number of rows and columns of subplots, by default (1, 1). Used to adjust height accordingly.
    ratio : float, optional
        The aspect ratio of the figure (height/width), by default None (uses golden ratio).

    Returns
    -------
    tuple of float
        Dimensions of the figure in inches (width, height).
    """
    WIDTHS = {
        # Nature
        "nature_single": 255,  # 90mm
        "nature_double": 510,  # 180mm
        # Science
        "science_single": 162,  # 5.7cm
        "science_double": 343,  # 12.1cm
        "science_triple": 521,  # 18.4cm
        # Cell/Neuron
        "cell_single": 241,  # 8.5cm
        "cell_1p5": 323,  # 11.4cm
        "cell_double": 493,  # 17.4cm
        # Generic aliases
        "single_col": 255,
        "double_col": 510,
        "beamer": 307.28987,
        "thesis": 426.79135,
        "textwidth": 418,
        "paper": 595.276,
    }

    if isinstance(width, str):
        if width not in WIDTHS:
            raise ValueError(
                f"Unknown width preset '{width}'. Choose from: {list(WIDTHS)}"
            )

        width_pt = WIDTHS[width]
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    aspect = ratio if ratio is not None else (5**0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * aspect * (subplots[0] / subplots[1])

    return fig_width_in, fig_height_in


def extract_video_name(file_path: str) -> str:
    """Extract canonical video stem used to match related DLC files.

    Parameters
    ----------
    file_path : str
        Path to CSV, H5, pickle, or video file.

    Returns
    -------
    str
        Normalized video stem used for cross-file matching.
    """
    name = os.path.basename(file_path)
    name = os.path.splitext(name)[0]
    name = re.split(r"DLC", name, flags=re.IGNORECASE)[0]
    name = re.sub(r"_meta$", "", name, flags=re.IGNORECASE)
    return name.rstrip("_")


def add_original_suffix(file_path: str) -> str:
    """Build backup filename using ``_original`` suffix.

    Parameters
    ----------
    file_path : str
        Source file path.

    Returns
    -------
    str
        Basename with ``_original`` inserted before extension.
    """
    stem, ext = os.path.splitext(os.path.basename(file_path))
    return f"{stem}_original{ext}"


def backup_original_file(file_path: str, backup_dir: str) -> str:
    """Copy a source file to backup directory with original suffix.

    Parameters
    ----------
    file_path : str
        Source file path to back up.
    backup_dir : str
        Destination backup directory.

    Returns
    -------
    str
        Full path to created backup file.
    """
    backup_name = add_original_suffix(file_path)
    backup_path = os.path.join(backup_dir, backup_name)
    shutil.copy2(file_path, backup_path)
    return backup_path


def get_xy_columns(df: pd.DataFrame) -> Tuple[pd.Index, pd.Index]:
    """Return x and y coordinate columns from DLC dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        DLC dataframe expected to have a 3-level MultiIndex header.

    Returns
    -------
    tuple of pandas.Index
        Tuple ``(x_cols, y_cols)`` containing matching coordinate columns.

    Raises
    ------
    ValueError
        If dataframe does not have expected MultiIndex structure or no x/y
        columns are found.
    """
    if not isinstance(df.columns, pd.MultiIndex) or df.columns.nlevels < 3:
        raise ValueError("Expected DLC dataframe with 3-level MultiIndex columns")

    coord_level = df.columns.get_level_values(2).astype(str).str.lower()
    x_cols = df.columns[coord_level == "x"]
    y_cols = df.columns[coord_level == "y"]

    if len(x_cols) == 0 or len(y_cols) == 0:
        raise ValueError("No x/y columns found in DLC dataframe")
    return x_cols, y_cols


def apply_xy_offset_to_dataframe(
    df: pd.DataFrame, x_min: Union[int, float], y_min: Union[int, float]
) -> pd.DataFrame:
    """Apply x/y offsets to all coordinate columns in a DLC dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        DLC dataframe with coordinate columns.
    x_min : float
        X offset to add.
    y_min : float
        Y offset to add.

    Returns
    -------
    pandas.DataFrame
        Dataframe with shifted x/y coordinate columns.
    """
    x_cols, y_cols = get_xy_columns(df)
    df[x_cols] = df[x_cols].astype(float) + float(x_min)
    df[y_cols] = df[y_cols].astype(float) + float(y_min)
    return df


def inspect_h5_keys(h5_file: str) -> list[str]:
    """List dataset keys stored in an H5 file.

    Parameters
    ----------
    h5_file : str
        Path to H5 file.

    Returns
    -------
    list of str
        Dataset keys present in the HDF store.
    """
    with pd.HDFStore(h5_file, mode="r") as store:
        return list(store.keys())


def correct_h5_in_place(
    h5_file: str, x_min: Union[int, float], y_min: Union[int, float]
) -> None:
    """Shift H5-stored DLC coordinates while preserving key and storage format.

    Parameters
    ----------
    h5_file : str
        Path to DLC H5 file.
    x_min : float
        X offset to add.
    y_min : float
        Y offset to add.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If H5 file does not contain exactly one dataset key.
    """
    keys = inspect_h5_keys(h5_file)
    if len(keys) != 1:
        raise ValueError(f"Expected exactly one dataset key in {h5_file}, found {keys}")

    key = keys[0]
    df = pd.read_hdf(h5_file, key=key)
    df = apply_xy_offset_to_dataframe(df, x_min, y_min)

    with pd.HDFStore(h5_file, mode="r+") as store:
        storer = store.get_storer(key)
        format_type = getattr(storer, "format_type", "fixed") or "fixed"
        data_columns = getattr(storer, "data_columns", None)
        store.remove(key)
        if format_type == "table":
            put_kwargs = {"format": "table"}
            if data_columns is not None:
                put_kwargs["data_columns"] = list(data_columns)
            store.put(key, df, **put_kwargs)
        else:
            store.put(key, df, format="fixed")


def find_matching_video_file(basepath: str, video_name: str) -> Optional[str]:
    """Find a video file by stem in a base directory.

    Parameters
    ----------
    basepath : str
        Directory containing video files.
    video_name : str
        Video stem name to match.

    Returns
    -------
    str or None
        Matching video path if found, otherwise ``None``.
    """
    preferred_ext = [".avi", ".mp4", ".mov", ".mkv"]
    for ext in preferred_ext:
        candidate = os.path.join(basepath, f"{video_name}{ext}")
        if os.path.exists(candidate):
            return candidate

    wildcard = glob.glob(os.path.join(basepath, f"{video_name}.*"))
    for candidate in wildcard:
        if os.path.splitext(candidate)[1].lower() in preferred_ext:
            return candidate
    return None


def get_video_dimensions(
    basepath: str, video_name: str
) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    """Get video path and frame dimensions for a matched stem.

    Parameters
    ----------
    basepath : str
        Directory containing video files.
    video_name : str
        Video stem to resolve.

    Returns
    -------
    tuple
        ``(video_file, width, height)`` where width/height may be ``None`` if
        video cannot be read.
    """
    video_file = find_matching_video_file(basepath, video_name)
    if video_file is None:
        return None, None, None

    cap = cv2.VideoCapture(video_file)
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    cap.release()

    if width <= 0 or height <= 0:
        return video_file, None, None

    return video_file, width, height


def update_pickle_with_video_dimensions(
    pickle_file: str, basepath: str, video_name: str
) -> Optional[str]:
    """Update pickle metadata from matching video dimensions.

    Parameters
    ----------
    pickle_file : str
        Path to pickle metadata file.
    basepath : str
        Directory containing related files.
    video_name : str
        Video stem used to locate matching video.

    Returns
    -------
    str or None
        Matched video path when update succeeds, otherwise ``None``.
    """
    video_file, width, height = get_video_dimensions(basepath, video_name)
    if video_file is None:
        print(f"Warning: no video file found for {video_name}; skipping pickle update")
        return None

    if width is None or height is None:
        print(f"Warning: failed to read video dimensions from {video_file}")
        return None

    with open(pickle_file, "rb") as f:
        data = pickle.load(f)

    full_bounds = [0, width, 0, height]
    data.setdefault("data", {})["video_dimensions"] = full_bounds
    data.setdefault("data", {})["cropping_parameters"] = full_bounds

    with open(pickle_file, "wb") as f:
        pickle.dump(data, f)

    return video_file


def apply_crop_parameters_to_dlc(
    basepath: str, backup_folder: str = "dlc_originals"
) -> list[dict[str, object]]:
    """
    Apply crop offsets from pickle files to DLC csv/h5 files in-place.

    This function creates one-time backups with '_original' suffix in
    `basepath/backup_folder`, then overwrites the original csv/h5 files.

    Parameters
    ----------
    basepath : str
        Path containing DLC csv, DLC h5, pickle, and video files.
    backup_folder : str, optional
        Backup folder inside basepath. Default is 'dlc_originals'.

    Returns
    -------
    list of dict
        Per-video summary of corrected files.
    """
    backup_dir = os.path.join(basepath, backup_folder)
    backup_complete = os.path.join(backup_dir, "_backup_complete")

    if os.path.exists(backup_complete):
        msg = (
            "Crop correction already applied for this basepath. "
            "Backups exist in 'dlc_originals'. Re-running would double-correct coordinates."
        )
        raise RuntimeError(msg)

    dlc_csv_files = glob.glob(os.path.join(basepath, "*DLC*.csv"))
    dlc_h5_files = glob.glob(os.path.join(basepath, "*DLC*.h5"))
    pickle_files = glob.glob(os.path.join(basepath, "*.pickle"))

    csv_map = {extract_video_name(p): p for p in dlc_csv_files}
    h5_map = {extract_video_name(p): p for p in dlc_h5_files}
    pickle_map = {extract_video_name(p): p for p in pickle_files}

    common_video_names = sorted(set(csv_map) & set(h5_map) & set(pickle_map))

    if not common_video_names:
        raise FileNotFoundError(
            "No matching csv/h5/pickle sets were found in basepath. "
            "Ensure files share the same video stem."
        )

    corrected_files = []

    os.makedirs(backup_dir, exist_ok=False)

    for video_name in common_video_names:
        csv_file = csv_map[video_name]
        h5_file = h5_map[video_name]
        pickle_file = pickle_map[video_name]

        backup_original_file(csv_file, backup_dir)
        backup_original_file(h5_file, backup_dir)
        backup_original_file(pickle_file, backup_dir)

        pkl = pd.read_pickle(pickle_file)
        cropping = pkl["data"]["cropping_parameters"]
        x_min = cropping[0]
        y_min = cropping[2]

        h5_keys = inspect_h5_keys(h5_file)
        if len(h5_keys) != 1:
            raise ValueError(f"Expected exactly one dataset key in {h5_file}, found {h5_keys}")

        csv_df = pd.read_csv(csv_file, header=[0, 1, 2], low_memory=False)
        csv_df = apply_xy_offset_to_dataframe(csv_df, x_min, y_min)
        csv_df.to_csv(csv_file, index=False)
        correct_h5_in_place(h5_file, x_min, y_min)

        video_file = update_pickle_with_video_dimensions(
            pickle_file, basepath, video_name
        )

        corrected_files.append(
            {
                "video_name": video_name,
                "csv": csv_file,
                "h5": h5_file,
                "pickle": pickle_file,
                "video": video_file,
                "x_min": x_min,
                "y_min": y_min,
                "h5_key": h5_keys[0],
            }
        )

        print(f"Saved in-place: {csv_file}")
        print(f"Saved in-place: {h5_file}")

    with open(backup_complete, "w", encoding="utf-8") as f:
        f.write("backup complete\n")

    return corrected_files


# Backward-compatible alias.
apply_crop_parameters_to_DLC = apply_crop_parameters_to_dlc


def infer_applied_offsets_from_csv_pair(
    original_csv: str, corrected_csv: str
) -> Tuple[float, float]:
    """Infer applied coordinate offsets from original/corrected CSV pair.

    Parameters
    ----------
    original_csv : str
        Path to original DLC CSV.
    corrected_csv : str
        Path to corrected DLC CSV.

    Returns
    -------
    tuple of float
        Inferred ``(x_offset, y_offset)`` from median coordinate shifts.
    """
    original_df = pd.read_csv(original_csv, header=[0, 1, 2], low_memory=False)
    corrected_df = pd.read_csv(corrected_csv, header=[0, 1, 2], low_memory=False)

    x_cols_orig, y_cols_orig = get_xy_columns(original_df)
    x_cols_corr, y_cols_corr = get_xy_columns(corrected_df)

    x_offset = float(
        np.nanmedian(corrected_df[x_cols_corr[0]].astype(float).values)
        - np.nanmedian(original_df[x_cols_orig[0]].astype(float).values)
    )
    y_offset = float(
        np.nanmedian(corrected_df[y_cols_corr[0]].astype(float).values)
        - np.nanmedian(original_df[y_cols_orig[0]].astype(float).values)
    )

    return x_offset, y_offset


def robust_fix_pickle_cropping_parameters(
    basepath: str, backup_folder: str = "dlc_originals", tolerance: float = 0.5
) -> list[dict[str, object]]:
    """
    Repair stale pickle cropping parameters by inspecting original vs corrected CSV files.

    This function does NOT modify CSV/H5 coordinates. It only reconciles pickle metadata.

    Behavior:
    - Sets pickle cropping parameters to full video bounds [0, width, 0, height]
      when corrections appear already applied or pickle metadata is stale.
    - Updates video_dimensions alongside cropping_parameters for consistency.
    - Reports unresolved mismatches for manual review.

    Parameters
    ----------
    basepath : str
        Path containing corrected DLC files and pickle files.
    backup_folder : str, optional
        Folder containing '_original' backups. Default is 'dlc_originals'.
    tolerance : float, optional
        Absolute tolerance (pixels) for offset comparisons. Default is 0.5.

    Returns
    -------
    list of dict
        Per-video status with inferred offsets and update action.
    """
    backup_dir = os.path.join(basepath, backup_folder)
    if not os.path.isdir(backup_dir):
        raise FileNotFoundError(f"Backup folder not found: {backup_dir}")

    pickle_files = glob.glob(os.path.join(basepath, "*.pickle"))
    if len(pickle_files) == 0:
        raise FileNotFoundError(f"No pickle files found in {basepath}")

    pickle_map = {extract_video_name(p): p for p in pickle_files}
    results = []

    for video_name, pickle_file in sorted(pickle_map.items()):
        corrected_csv = glob.glob(os.path.join(basepath, f"{video_name}*DLC*.csv"))
        original_csv = glob.glob(os.path.join(backup_dir, f"{video_name}*DLC*_original.csv"))

        if len(corrected_csv) != 1 or len(original_csv) != 1:
            results.append(
                {
                    "video_name": video_name,
                    "status": "skipped_missing_or_ambiguous_csv",
                    "corrected_csv_count": len(corrected_csv),
                    "original_csv_count": len(original_csv),
                }
            )
            continue

        corrected_csv = corrected_csv[0]
        original_csv = original_csv[0]

        with open(pickle_file, "rb") as f:
            pkl_data = pickle.load(f)

        crop = pkl_data.get("data", {}).get("cropping_parameters", [0, 0, 0, 0])
        x_pickle = float(crop[0])
        y_pickle = float(crop[2])

        x_inferred, y_inferred = infer_applied_offsets_from_csv_pair(original_csv, corrected_csv)

        offsets_match_pickle = (
            abs(x_inferred - x_pickle) <= tolerance and abs(y_inferred - y_pickle) <= tolerance
        )
        offsets_are_zero = abs(x_inferred) <= tolerance and abs(y_inferred) <= tolerance
        pickle_is_zero = abs(x_pickle) <= tolerance and abs(y_pickle) <= tolerance

        video_file, width, height = get_video_dimensions(basepath, video_name)
        if video_file is None or width is None or height is None:
            results.append(
                {
                    "video_name": video_name,
                    "pickle": pickle_file,
                    "original_csv": original_csv,
                    "corrected_csv": corrected_csv,
                    "x_pickle": x_pickle,
                    "y_pickle": y_pickle,
                    "x_inferred": x_inferred,
                    "y_inferred": y_inferred,
                    "action": "skipped_video_missing_or_invalid",
                }
            )
            continue

        full_bounds = [0, width, 0, height]

        if offsets_match_pickle and not pickle_is_zero:
            pkl_data.setdefault("data", {})["cropping_parameters"] = full_bounds
            pkl_data.setdefault("data", {})["video_dimensions"] = full_bounds
            with open(pickle_file, "wb") as f:
                pickle.dump(pkl_data, f)
            action = "updated_pickle_to_full_video_bounds"
        elif pickle_is_zero and not offsets_are_zero:
            pkl_data.setdefault("data", {})["cropping_parameters"] = full_bounds
            pkl_data.setdefault("data", {})["video_dimensions"] = full_bounds
            with open(pickle_file, "wb") as f:
                pickle.dump(pkl_data, f)
            action = "fixed_zero_pickle_to_full_video_bounds"
        elif offsets_are_zero and not pickle_is_zero:
            action = "mismatch_offsets_zero_but_pickle_nonzero"
        elif pickle_is_zero and offsets_are_zero:
            pkl_data.setdefault("data", {})["cropping_parameters"] = full_bounds
            pkl_data.setdefault("data", {})["video_dimensions"] = full_bounds
            with open(pickle_file, "wb") as f:
                pickle.dump(pkl_data, f)
            action = "normalized_zero_to_full_video_bounds"
        else:
            action = "mismatch_manual_review"

        results.append(
            {
                "video_name": video_name,
                "pickle": pickle_file,
                "original_csv": original_csv,
                "corrected_csv": corrected_csv,
                "x_pickle": x_pickle,
                "y_pickle": y_pickle,
                "x_inferred": x_inferred,
                "y_inferred": y_inferred,
                "video_bounds": full_bounds,
                "action": action,
            }
        )

    print("\nRobust pickle crop repair summary")
    for r in results:
        label = r.get("action", r.get("status"))
        print(f"{r['video_name']}: {label}")

    return results