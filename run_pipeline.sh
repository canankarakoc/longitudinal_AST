#!/usr/bin/env bash
# Run the AMR pipeline, filtering rEDM's C++ NaN chatter out of stderr.
#
# rEDM prints "WARNING: ComputeError(): nan found..." once per embedding vector
# it refuses to build across a patient boundary. That is the NA-separator design
# working as intended - suppressWarnings() and capture.output() cannot catch it,
# because it is written from C++ straight to the process stderr. Real errors
# still come through.
set -euo pipefail
cd "$(dirname "$0")"
STEP="${1:-all}"
run() { echo "=== $1 ==="; Rscript "$1" 2> >(grep -v 'ComputeError(): nan found' >&2); }
case "$STEP" in
  prep)     run analytical_data_prep.R ;;
  summary)  run data_inspection_summary.R ;;
  pipeline) run AMR_pipeline.R ;;
  all)      run analytical_data_prep.R; run data_inspection_summary.R; run AMR_pipeline.R ;;
  *) echo "usage: $0 [prep|summary|pipeline|all]"; exit 1 ;;
esac
