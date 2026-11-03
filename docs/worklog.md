# Work log

Chronological record of what was done outside the notebooks (setup, terminal commands, data,
decisions). Analysis steps themselves are documented inside the notebooks.

Format per entry: date – what – why – commands (if any) – result.

---

## 2026-10-29 – Python environment

**What:** Removed an old, broken Anaconda installation (conda 23.7.4) and installed Miniforge.
**Why:** The old install had plugin errors; Miniforge uses only the community channel
`conda-forge` and avoids conflicts between two conda installations.

```bash
# remove old Anaconda (shell config first, then folders)
conda init --reverse --all
rm -rf ~/anaconda3 ~/.conda ~/.condarc ~/.continuum

# install Miniforge (Apple Silicon)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
bash Miniforge3-MacOSX-arm64.sh

# do not activate the base environment automatically
conda config --set auto_activate false
```

**Result:** conda 26.9.4, channel `conda-forge` only.

## 2026-10-29 – Project environment `bia`

**Why a separate environment:** keeps project packages isolated from the system and from other projects.
StarDist (TensorFlow) is deliberately **not** installed here to avoid conflicts with PyTorch
(Cellpose); StarDist will run via the Fiji plugin or Colab.

```bash
conda create -y -n bia python=3.12
conda activate bia
python -m pip install "napari[all]" cellpose tifffile scikit-image pandas matplotlib jupyterlab

# checks
python -c "import torch; print('MPS:', torch.backends.mps.is_available())"
python -c "import cellpose, napari; print(cellpose.version, napari.__version__)"
```

**Result:** Python 3.12.16, cellpose 4.2.1.1, napari 0.9.2, PyTorch MPS (Apple GPU) available.

**Reminder:** every new terminal starts without an environment → `conda activate bia` first.

## 2026-10-29 – Repository

```bash
cd ~/Projects/ctc-hela-dl-tracking
git remote add origin https://github.com/Tritium13/ctc-hela-dl-tracking.git
git push -u origin main
```

Private for now; will be made public once the README and first results are presentable.

## 2026-10-29 / 30 – Dataset download (pending)

**Dataset:** Fluo-N2DL-HeLa (training), chosen for many nuclei per frame, visible mitoses,
and available silver-truth annotations.

```bash
cd ~/Projects/ctc-hela-dl-tracking/data/raw
curl -L -O https://data.celltrackingchallenge.net/training-datasets/Fluo-N2DL-HeLa.zip
unzip -q Fluo-N2DL-HeLa.zip
```

**Result:** server not reachable (timeout via https and http, also from a mobile network →
not a local network problem). Retry later.

## 2026-10-30 – First notebook on example data

`notebooks/01_cellpose_sample.ipynb` – Cellpose-SAM on `skimage.data.human_mitosis()`
while the CTC server is down. 333 nuclei; some small fragments and one elongated debris object.

Starting JupyterLab:

```bash
conda activate bia
cd ~/Projects/ctc-hela-dl-tracking
jupyter lab
```

## 2026-10-30 – Data

Fluo-N2DL-HeLa training data obtained from a Kaggle mirror (xunbch/hela-train) while the CTC server was down;
treated as the original. Location: `data/raw/Fluo-N2DL-HeLa/`. Verify by checksum against the official zip when available.

## 2026-10-30 – Fine-tuning on Kaggle

**Why Kaggle:** free GPU, dataset already available there (no upload), no Google account needed.
Requires phone verification for GPU and internet access.

**Setup:** Kaggle notebook imported from `notebooks/04_finetune_cellpose_kaggle.ipynb`;
input dataset `xunbch/hela-train` (mounted at `/kaggle/input/datasets/xunbch/hela-train/`);
accelerator 2 × Tesla T4 (15 GB each, Cellpose uses one); cellpose 4.2.1.1, Python 3.13, torch 2.11 (CUDA 12.8).

**Run:** version "v1 – seq02 ST, 100 epochs, lr 1e-5", started via *Save Version → Save & Run All (Commit)*.
Training data: `scripts/prepare_cellpose_training.py` logic inline (seq 02 silver truth, 82 train / 10 test frames).
Outputs (Kaggle notebook output tab): `models/cpsam_hela_seq02`, `learning_curves.png/.csv`,
`metrics_pretrained_vs_finetuned_seq01.csv`.

## 2026-11-03 – CTC server still unreachable

data.celltrackingchallenge.net still times out (checked 2026-10-29, 10-30, 11-03). Checksum comparison
with the Kaggle copy postponed.
