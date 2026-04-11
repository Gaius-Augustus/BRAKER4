#!/bin/bash
# Per-scenario overrides for scenario_08_es.
# Sourced by run_test.sh after compute_profile.sh, so anything exported here
# wins over the shared HPC defaults.

# This is the optional-QC test scenario — the fastest HPC scenario (ES mode,
# pre-masked O. tauri), used to exercise the optional protein-level annotation
# steps that depend on the BRAKER protein output.

# OMArk proteome quality assessment.
export BRAKER4_RUN_OMARK=1

# FANTASIA-Lite functional annotation (GPU-only). The shared HPC config.ini
# defines the SIF path, HuggingFace cache, and SLURM GPU partition; this flag
# only flips the toggle for this one scenario. See README "run_fantasia" for
# the fragility caveats.
export BRAKER4_RUN_FANTASIA=1
