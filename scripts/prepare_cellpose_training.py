"""Prepare Cellpose training data from Cell Tracking Challenge silver truth.

Copies the raw frames and silver-truth masks of one sequence into the folder layout expected by
`cellpose.io.load_train_test_data`: `<name>.tif` + `<name>_masks.tif`.

Every k-th frame goes to a separate `test/` folder. It is only used to monitor training;
the actual evaluation is done on the gold truth of the *other* sequence (see notebook 04).

Usage (from the repository root):
    python scripts/prepare_cellpose_training.py
    python scripts/prepare_cellpose_training.py --seq 02 --test-every 10
"""
import argparse
import shutil
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="data/raw/Fluo-N2DL-HeLa", help="CTC dataset folder")
    ap.add_argument("--seq", default="02", help="sequence used for training")
    ap.add_argument("--out", default=None, help="output folder (default: data/training/cellpose_seq<SEQ>)")
    ap.add_argument("--test-every", type=int, default=10, help="every k-th frame goes to test/")
    args = ap.parse_args()

    data = Path(args.data)
    out = Path(args.out or f"data/training/cellpose_seq{args.seq}")
    raw_dir = data / args.seq
    st_dir = data / f"{args.seq}_ST" / "SEG"
    if not raw_dir.exists() or not st_dir.exists():
        raise FileNotFoundError(f"expected {raw_dir} and {st_dir} – run from the repository root")

    train_dir, test_dir = out / "train", out / "test"
    train_dir.mkdir(parents=True, exist_ok=True)
    test_dir.mkdir(parents=True, exist_ok=True)

    n = {"train": 0, "test": 0, "skipped": 0}
    for img in sorted(raw_dir.glob("t*.tif")):
        t = int(img.stem[1:])                          # t012.tif -> 12
        mask = st_dir / f"man_seg{t:03d}.tif"
        if not mask.exists():
            print(f"no silver-truth mask for frame {t}, skipped")
            n["skipped"] += 1
            continue
        split = "test" if t % args.test_every == 0 else "train"
        dest = test_dir if split == "test" else train_dir
        name = f"seq{args.seq}_t{t:03d}"
        shutil.copy2(img, dest / f"{name}.tif")
        shutil.copy2(mask, dest / f"{name}_masks.tif")
        n[split] += 1

    print(f"train: {n['train']} images | test: {n['test']} images | skipped: {n['skipped']}")
    print(f"output: {out.resolve()}")


if __name__ == "__main__":
    main()
