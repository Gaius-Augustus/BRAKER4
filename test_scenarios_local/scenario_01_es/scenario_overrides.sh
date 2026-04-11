#!/bin/bash
# Per-scenario overrides for scenario_01_es (local).

# Tiny test data — BUSCO is not informative; skip it for speed.
export BRAKER4_SKIP_BUSCO=1

# Keep intermediate files so the developer can inspect a failed local run.
export BRAKER4_NO_CLEANUP=1
