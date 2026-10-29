# ctc-hela-dl-tracking

Deep-learning segmentation and tracking of HeLa nuclei from the
[Cell Tracking Challenge](https://celltrackingchallenge.net), followed by a
statistical analysis of nuclear motion in R.

> **Status:** work in progress. Results and figures will be added as each stage is completed.

## Question (draft)

1. How well do pretrained deep-learning models (Cellpose-SAM, StarDist `2D_versatile_fluo`)
   segment H2B-GFP-labelled HeLa nuclei compared with the challenge reference annotations,
   and does training a model on this dataset improve the result?
2. Does nuclear migration speed change around mitosis?

## Data

**Fluo-N2DL-HeLa** (training set) from the Cell Tracking Challenge:
HeLa cells stably expressing H2B-GFP, Mitocheck Consortium.
Olympus IX81, 10x/0.4, pixel size 0.645 × 0.645 µm, time step 30 min.

The data are **not** included in this repository. Download them from the
[2D+Time datasets page](https://celltrackingchallenge.net/2d-datasets/) and unpack into `data/raw/`.

## Workflow

| Stage | Tools | Status |
|---|---|---|
| 1. Orientation: inspect data | Python, napari, Fiji | in progress |
| 2. Segmentation: pretrained models, evaluation against reference, own training | Cellpose, StarDist, ZeroCostDL4Mic (Colab) | planned |
| 3. Tracking and statistics | TrackMate, R (lme4) | planned |
| 4. Documentation | — | planned |

## Repository layout

```
data/        raw data (not tracked by git)
notebooks/   Jupyter notebooks, one per analysis step
scripts/     reusable Python scripts
R/           statistical analysis in R
results/     tables and figures
```

## Environment

Local analysis runs in a conda environment (Miniforge, Python 3.12) on macOS (Apple Silicon):

```bash
conda env create -f environment.yml
conda activate bia
```

Model training runs in Google Colab.

## License and acknowledgements

Code in this repository: MIT License (see `LICENSE`).
The dataset belongs to its providers; see the Cell Tracking Challenge terms of use.
Pretrained Cellpose models are trained on data licensed CC-BY-NC and are not redistributed here.
