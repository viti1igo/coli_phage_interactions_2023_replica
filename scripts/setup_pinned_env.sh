#!/bin/bash
# Env for the pinned rerun (T81). The authors state scikit-learn 1.1.2 and give no Python version for the model step (3.12.2 is only the kernel metadata of the plotting notebook visualize_predictions.ipynb).
# scikit-learn 1.1.2 has no Python 3.12 build (conda-forge osx-arm64 has py310 only; the release predates 3.12), so Python is 3.10.
# scipy 1.11.2 is the version the paper states (scipy.stats). numpy, pandas and joblib versions are not stated anywhere: set to releases contemporary with scikit-learn 1.1.2.
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba-phage}
micromamba create -y -n ml_pinned --override-channels -c conda-forge python=3.10 scikit-learn=1.1.2 numpy=1.23.5 pandas=1.5.3 scipy=1.11.2 joblib=1.2.0
