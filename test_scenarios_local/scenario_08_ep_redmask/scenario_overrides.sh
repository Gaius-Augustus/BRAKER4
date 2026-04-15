#!/bin/bash
# Per-scenario overrides for scenario_08_ep_redmask (local).

# Use Red instead of RepeatModeler+RepeatMasker for repeat masking.
export BRAKER4_MASKING_TOOL=red

# Tiny test data — BUSCO is not informative; skip it for speed.
export BRAKER4_SKIP_BUSCO=1

# Keep intermediate files for inspection.
export BRAKER4_NO_CLEANUP=1
