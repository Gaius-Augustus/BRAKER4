#!/bin/bash
# Per-scenario overrides for scenario_10_et_varus.

# This is the VARUS auto-download test.
export BRAKER4_USE_VARUS=1

# VARUS download + masking + alignment takes hours; bump runtime.
export BRAKER4_MAX_RUNTIME=4320
