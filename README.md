# multievolve-baselines

This repository contains scripts and notebooks that compare MULTI-evolve's
neural network predictions against additive (linear) baselines for the APEX,
dCasRx, and HuABC2 engineering tasks and reanalyze the DMS benchmarking results.

## Setup

### 1. Clone with submodule

```bash
git clone --recurse-submodules https://github.com/dewitt-lab/multievolve-baselines.git
cd multievolve-baselines
```

### 2. Apply patch

Our analysis extends MULTI-evolve's `Fcn` class to support zero hidden layers
(a pure linear map, tantamount to linear regression). Apply the patch:

```bash
git -C MULTI-evolve apply ../patches/0001-support-zero-hidden-layers.patch
```

### 3. Create the conda environment

Use the environment file from MULTI-evolve:

```bash
conda env create -f MULTI-evolve/env.yml
conda activate multievolve
pip install -e MULTI-evolve
pip install openpyxl
```

### 4. Download data from Zenodo

The analysis requires supplementary data tables from the original MULTI-evolve
publication. These are accessible from [Zenodo: 10.5281/zenodo.18240028](https://doi.org/10.5281/zenodo.18240028):

Download and move the following files into the `data/` directory:

- `Table S14 APEX variants.xlsx` – APEX engineering campaign data
- `Table S15 dCasRx variants.xlsx` — dCasRx engineering campaign data
- `Table S16 HuABC2 variants.xlsx` — HuABC2 engineering campaign data
- `Table S3 *.xlsx` — DMS benchmarking results

DMS benchmarking reanalysis requires data from [Zenodo: 10.5281/zenodo.17620759](https://doi.org/10.5281/zenodo.17620759). Download these files and move to `data/benchmark/`.

## Reproducing the analysis

All commands are run from the repository root with the `multievolve` conda
environment activated. The pipeline has three stages: data preparation,
model training/scoring, and figure generation.

### Stage 1: Prepare training data

Extract training CSVs from the supplementary tables:

```bash
python analysis/prep_training_data.py all
```

This creates `output/{apex,dcasrx,huabc2}/` with training CSVs
and mutation lists (the `output/` directory is gitignored).

### Stage 2: Train models and score combinatorial variants

Train the MULTI-evolve FCNN ensemble (10-fold CV, matching the original
hyperparameters) and score all combinatorial variants for each protein. Also
train with zero hidden layers (linear model) for comparison:

```bash
# Standard MULTI-evolve FCNN
python analysis/run_multievolve.py apex
python analysis/run_multievolve.py dcasrx
python analysis/run_multievolve.py huabc2 --measurement expr
python analysis/run_multievolve.py huabc2 --measurement bind

# Zero hidden layers (linear baseline)
python analysis/run_multievolve.py apex --zerolayer
python analysis/run_multievolve.py dcasrx --zerolayer
python analysis/run_multievolve.py huabc2 --measurement expr --zerolayer
python analysis/run_multievolve.py huabc2 --measurement bind --zerolayer
```

Output: `output/{protein}/{protein}_multievolve_proposals.csv` and
`output/{protein}/{protein}_0layer_multievolve_proposals.csv`.

Run the hyperparameter sweep (extended grid including zero hidden layers and
ridge regression):

```bash
python analysis/multievolve_hyperparameter_sweep.py apex
python analysis/multievolve_hyperparameter_sweep.py dcasrx
python analysis/multievolve_hyperparameter_sweep.py huabc2 --measurement expr
python analysis/multievolve_hyperparameter_sweep.py huabc2 --measurement bind
```

Output: `output/{protein}/{protein}_log2_hparams_results.jsonl`.

### Stage 3: Generate figures

Open and run the Jupyter notebooks in `analysis/`. Each notebook reads the
intermediate CSVs produced in Stage 2 and writes figures to `figures/`.

| Notebook | Figures produced |
|----------|-----------------|
| `additive-comparison-APEX.ipynb` | APEX additive vs. MULTI-evolve correlation |
| `additive-comparison-dCasRx.ipynb` | dCasRx additive vs. MULTI-evolve correlation |
| `additive-comparison-HuABC2.ipynb` | HuABC2 additive vs. MULTI-evolve correlation, Pareto frontiers |
| `epistasis.ipynb` | Epistasis in training data vs. NN predictions (set `PROTEIN` and `MEASUREMENT` in the config cell) |
| `additive-proteingym_fig2ef.ipynb` | DMS benchmark reanalysis |
| `analyze_hparam_sweep.ipynb` | Hyperparameter sweep results |


## MULTI-evolve dependency

This repository uses [MULTI-evolve](https://github.com/ArcInstitute/MULTI-evolve)
(Tran et al., *Science*, 2026) as a git submodule, pinned at commit `5cdb074`.
The submodule provides the neural network training code, featurizers, and
proposer infrastructure. Our analysis scripts call into this codebase to train
models.

A small patch (`patches/0001-support-zero-hidden-layers.patch`) modifies the
`Fcn` class to accept an empty architecture (zero hidden layers), producing a
pure linear map from one-hot features to the output. This is used to test
whether MULTI-evolve's nonlinear layers contribute beyond a linear baseline.

## License

MIT
