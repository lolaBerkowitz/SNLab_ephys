import cv2
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import os
import numpy as np
import pandas as pd
import sys

sys.path.append(r"C:\Users\schafferlab\github\SNLab_ephys\preprocessing\behavior")
from DLC_preprocessing import update_config_cropping


def _raise_window_to_front(fig=None, window_name=None):
    """Best-effort raise/focus for matplotlib or OpenCV windows."""
    if fig is not None:
        try:
            manager = fig.canvas.manager
            if hasattr(manager, "window"):
                window = manager.window
                # Tk backend
                if hasattr(window, "wm_attributes"):
                    window.wm_attributes("-topmost", 1)
                    window.wm_attributes("-topmost", 0)
                    if hasattr(window, "focus_force"):
                        window.focus_force()
                elif hasattr(window, "attributes"):
                    window.attributes("-topmost", True)
                    window.attributes("-topmost", False)
                    if hasattr(window, "focus_force"):
                        window.focus_force()
                # Qt backend
                if hasattr(window, "raise_"):
                    window.raise_()
                if hasattr(window, "activateWindow"):
                    window.activateWindow()
        except Exception:
            pass

    if window_name is not None:
        try:
            cv2.setWindowProperty(window_name, cv2.WND_PROP_TOPMOST, 1)
        except Exception:
            pass

def _normalize_path(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return None
    return os.path.normpath(str(p))


def _video_stem_from_dlc_file(dlc_file):
    """Infer video stem from dlc_file name using existing strip_dlc_suffix helper."""
    if dlc_file is None or (isinstance(dlc_file, float) and np.isnan(dlc_file)):
        return None
    dlc_name = os.path.basename(str(dlc_file))
    return dlc_name.split("DLC_")[0].rstrip("_")


def load_preview_image(video_file, mode="middle", n_samples=100):
    """Load preview image from video as middle frame, mean frame, or both."""
    video_file = _normalize_path(video_file)
    if video_file is None or not os.path.exists(video_file):
        raise FileNotFoundError(f"Video not found: {video_file}")

    cap = cv2.VideoCapture(video_file)
    if not cap.isOpened():
        raise RuntimeError(f"Could not open video: {video_file}")

    frame_count = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    if frame_count <= 0:
        cap.release()
        raise RuntimeError(f"Video has no frames: {video_file}")

    middle_idx = max(frame_count // 2, 0)
    cap.set(cv2.CAP_PROP_POS_FRAMES, middle_idx)
    ok, mid_bgr = cap.read()
    if not ok:
        cap.release()
        raise RuntimeError(f"Could not read middle frame from: {video_file}")
    mid_rgb = cv2.cvtColor(mid_bgr, cv2.COLOR_BGR2RGB)

    if mode == "middle":
        cap.release()
        return {"middle": mid_rgb, "height": height, "width": width, "frame_count": frame_count}

    sample_idx = np.linspace(0, frame_count - 1, num=min(n_samples, frame_count), dtype=int)
    acc = None
    valid = 0
    for idx in sample_idx:
        cap.set(cv2.CAP_PROP_POS_FRAMES, int(idx))
        ok, bgr = cap.read()
        if not ok:
            continue
        rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB).astype(np.float32)
        acc = rgb if acc is None else acc + rgb
        valid += 1

    cap.release()
    if valid == 0:
        raise RuntimeError(f"Could not read sampled frames from: {video_file}")

    mean_rgb = (acc / valid).astype(np.uint8)
    if mode == "mean":
        return {"mean": mean_rgb, "height": height, "width": width, "frame_count": frame_count}
    if mode == "both":
        return {
            "middle": mid_rgb,
            "mean": mean_rgb,
            "height": height,
            "width": width,
            "frame_count": frame_count,
        }
    raise ValueError("mode must be one of: 'middle', 'mean', 'both'")


def show_preview(preview_dict, title_prefix="Preview"):
    """Display one or two preview images."""
    has_middle = "middle" in preview_dict
    has_mean = "mean" in preview_dict

    if has_middle and has_mean:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        axes[0].imshow(preview_dict["middle"])
        axes[0].set_title(f"{title_prefix}: middle frame")
        axes[0].axis("off")
        axes[1].imshow(preview_dict["mean"])
        axes[1].set_title(f"{title_prefix}: mean sampled frame")
        axes[1].axis("off")
        plt.tight_layout()
        plt.show()
    elif has_middle:
        plt.figure(figsize=(7, 6))
        plt.imshow(preview_dict["middle"])
        plt.title(f"{title_prefix}: middle frame")
        plt.axis("off")
        plt.show()
    elif has_mean:
        plt.figure(figsize=(7, 6))
        plt.imshow(preview_dict["mean"])
        plt.title(f"{title_prefix}: mean sampled frame")
        plt.axis("off")
        plt.show()


def _pick_crop_bounds_cv2(image):
    """Collect four clicks in an OpenCV window and return bounding crop limits."""
    points = []
    canvas = cv2.cvtColor(image.copy(), cv2.COLOR_RGB2BGR)
    window_name = "Select crop (click 4 points, ESC to cancel)"

    def _on_mouse(event, x, y, flags, param):
        if event == cv2.EVENT_LBUTTONDOWN and len(points) < 4:
            points.append((x, y))
            cv2.circle(canvas, (x, y), 5, (0, 0, 255), -1)
            cv2.imshow(window_name, canvas)

    cv2.namedWindow(window_name, cv2.WINDOW_NORMAL)
    cv2.imshow(window_name, canvas)
    _raise_window_to_front(window_name=window_name)
    cv2.setMouseCallback(window_name, _on_mouse)

    while len(points) < 4:
        key = cv2.waitKey(20) & 0xFF
        if key == 27:
            cv2.destroyWindow(window_name)
            raise RuntimeError("Crop selection cancelled by user (ESC)")

    cv2.destroyWindow(window_name)

    xs = [int(round(p[0])) for p in points]
    ys = [int(round(p[1])) for p in points]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    return x_min, x_max, y_min, y_max


def pick_crop_bounds(image, title="Click four points around the crop region"):
    """Pick four points with backend-aware fallback and return bounding limits."""
    backend = matplotlib.get_backend().lower()
    is_inline_backend = "inline" in backend

    if not is_inline_backend:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.imshow(image)
        ax.set_title(title)
        ax.axis("on")
        _raise_window_to_front(fig=fig)

        canvas_name = type(fig.canvas).__name__.lower()
        is_agg_canvas = "agg" in canvas_name

        if not is_agg_canvas:
            pts = plt.ginput(4, timeout=-1)
            plt.close(fig)
            if len(pts) == 4:
                xs = [int(round(p[0])) for p in pts]
                ys = [int(round(p[1])) for p in pts]
                x_min, x_max = min(xs), max(xs)
                y_min, y_max = min(ys), max(ys)
                return x_min, x_max, y_min, y_max
        else:
            plt.close(fig)

    try:
        return _pick_crop_bounds_cv2(image)
    except Exception as exc:
        print(f"OpenCV click fallback failed: {exc}")
        print("Enter coordinates manually as x1,y1,x2,y2,x3,y3,x4,y4")
        raw = input("coords: ").strip()
        vals = [int(v) for v in raw.split(",")]
        if len(vals) != 8:
            raise ValueError("Expected 8 integers: x1,y1,x2,y2,x3,y3,x4,y4")
        xs = vals[0::2]
        ys = vals[1::2]
        x_min, x_max = min(xs), max(xs)
        y_min, y_max = min(ys), max(ys)
        return x_min, x_max, y_min, y_max


def validate_and_clip_bounds(x_min, x_max, y_min, y_max, width, height):
    """Clip crop bounds to image dimensions and enforce non-zero region."""
    x_min = int(np.clip(x_min, 0, width - 1))
    x_max = int(np.clip(x_max, 0, width - 1))
    y_min = int(np.clip(y_min, 0, height - 1))
    y_max = int(np.clip(y_max, 0, height - 1))

    if x_min >= x_max or y_min >= y_max:
        raise ValueError(
            f"Invalid crop after clipping: x=({x_min},{x_max}), y=({y_min},{y_max}), frame=({width},{height})"
        )
    return x_min, x_max, y_min, y_max


def _resolve_target_row(df, basepath=None, video_file=None):
    """Resolve a unique row in check_dlc by exact video_file, then basepath+stem fallback."""
    basepath_norm = _normalize_path(basepath)
    video_file_norm = _normalize_path(video_file)

    if video_file_norm is not None and "video_file" in df.columns:
        row_mask = df["video_file"].astype(str).map(_normalize_path) == video_file_norm
        if row_mask.sum() == 1:
            return df.index[row_mask][0]
        if row_mask.sum() > 1:
            raise ValueError(f"Multiple rows matched video_file={video_file_norm}")

    if video_file_norm is None and basepath_norm is None:
        raise ValueError("Provide either video_file or basepath")

    if video_file_norm is None:
        avi_files = sorted(glob.glob(os.path.join(basepath_norm, "*.avi")))
        if len(avi_files) == 0:
            raise FileNotFoundError(f"No .avi files found in basepath={basepath_norm}")
        video_file_norm = _normalize_path(avi_files[0])

    target_stem = Path(video_file_norm).stem

    if "basepath" in df.columns and "dlc_file" in df.columns:
        base_mask = df["basepath"].astype(str).map(_normalize_path) == basepath_norm if basepath_norm else True
        stem_mask = df["dlc_file"].map(_video_stem_from_dlc_file) == target_stem
        row_mask = base_mask & stem_mask
        if row_mask.sum() == 1:
            return df.index[row_mask][0]
        if row_mask.sum() > 1:
            raise ValueError(
                f"Multiple rows matched basepath+stem for basepath={basepath_norm}, stem={target_stem}"
            )

    raise ValueError("Could not uniquely match a row in check_dlc.csv")


def select_crop_and_update_csv(
    check_csv_path,
    basepath=None,
    video_file=None,
    preview_mode="both",
    n_samples=100,
    save=True,
    backup=True,
    show_overlay=True,
    update_dlc_config=False,
    config_col="config_path",
):
    """Interactive crop picker + CSV updater for check_dlc rows."""
    check_csv_path = _normalize_path(check_csv_path)
    df = pd.read_csv(check_csv_path)

    row_idx = _resolve_target_row(df, basepath=basepath, video_file=video_file)

    if video_file is None and "video_file" in df.columns and pd.notna(df.loc[row_idx, "video_file"]):
        video_file = df.loc[row_idx, "video_file"]
    elif video_file is None:
        basepath_val = _normalize_path(df.loc[row_idx, "basepath"])
        dlc_stem = _video_stem_from_dlc_file(df.loc[row_idx, "dlc_file"])
        video_file = os.path.join(basepath_val, f"{dlc_stem}.avi")

    video_file = _normalize_path(video_file)
    preview = load_preview_image(video_file, mode=preview_mode, n_samples=n_samples)
    show_preview(preview, title_prefix=os.path.basename(video_file))

    click_image = preview["mean"] if preview_mode == "mean" else preview["middle"]
    x_min, x_max, y_min, y_max = pick_crop_bounds(click_image)
    x_min, x_max, y_min, y_max = validate_and_clip_bounds(
        x_min, x_max, y_min, y_max, width=preview["width"], height=preview["height"]
    )

    if show_overlay:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.imshow(click_image)
        rect = plt.Rectangle(
            (x_min, y_min), x_max - x_min, y_max - y_min,
            fill=False, linewidth=2, edgecolor="red"
        )
        ax.add_patch(rect)
        ax.set_title(
            f"Selected crop: x_min={x_min}, x_max={x_max}, y_min={y_min}, y_max={y_max}"
        )
        ax.axis("on")
        plt.show()

    df.loc[row_idx, "x_min"] = x_min
    df.loc[row_idx, "x_max"] = x_max
    df.loc[row_idx, "y_min"] = y_min
    df.loc[row_idx, "y_max"] = y_max
    if "video_file" in df.columns:
        df.loc[row_idx, "video_file"] = video_file
    else:
        df["video_file"] = np.nan
        df.loc[row_idx, "video_file"] = video_file

    if save:
        if backup:
            backup_path = check_csv_path.replace(".csv", ".bak.csv")
            if not os.path.exists(backup_path):
                df_backup = pd.read_csv(check_csv_path)
                df_backup.to_csv(backup_path, index=False)
        df.to_csv(check_csv_path, index=False)

    if update_dlc_config and config_col in df.columns and pd.notna(df.loc[row_idx, config_col]):
        update_config_cropping(
            _normalize_path(df.loc[row_idx, config_col]),
            int(x_min),
            int(x_max),
            int(y_min),
            int(y_max),
        )

    return {
        "row_index": int(row_idx),
        "video_file": video_file,
        "x_min": int(x_min),
        "x_max": int(x_max),
        "y_min": int(y_min),
        "y_max": int(y_max),
        "saved": bool(save),
    }

def get_video_bounds(video_path):
    """Return full-frame crop bounds for a video.

    Parameters
    ----------
    video_path : str
        Path to video file.

    Returns
    -------
    tuple[int, int, int, int]
        (y_min, y_max, x_min, x_max) in pixel units.
    """
    import os
    import cv2

    if not os.path.exists(video_path):
        raise FileNotFoundError(f"Video not found: {video_path}")

    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        raise RuntimeError(f"Could not open video: {video_path}")

    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    cap.release()

    if width <= 0 or height <= 0:
        raise RuntimeError(f"Invalid video dimensions for: {video_path}")

    y_min, y_max = 0, height
    x_min, x_max = 0, width
    return y_min, y_max, x_min, x_max