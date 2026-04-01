#!/bin/bash

COHORT="CA"   # (CA/CT)

ROOT=$(cd "$(dirname "$0")/.." && pwd)

IMAGE="chgyi/binbash-r-base:4.2.2-twosamplemr"

docker run --rm \
  -v "$ROOT":/work \
  -w /work \
  "$IMAGE" \
  Rscript R/run_mr.R $COHORT