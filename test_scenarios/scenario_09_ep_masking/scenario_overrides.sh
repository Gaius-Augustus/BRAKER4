#!/bin/bash
# Per-scenario overrides for scenario_09_ep_masking.

# Repeat masking + full ProtHint takes hours; bump the per-job runtime well
# above the shared 120 min default.
export BRAKER4_MAX_RUNTIME=4320
