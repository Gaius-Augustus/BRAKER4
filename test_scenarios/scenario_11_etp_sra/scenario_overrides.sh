#!/bin/bash
# Per-scenario overrides for scenario_11_etp_sra.

# This scenario explicitly exercises the AUGUSTUS optimisation step
# (skip_optimize_augustus = 0), which is skipped in the other test scenarios
# for speed.
export BRAKER4_SKIP_OPTIMIZE_AUGUSTUS=0
