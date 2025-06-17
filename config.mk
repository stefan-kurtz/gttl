GTTL ?= $(shell git rev-parse --show-toplevel 2>/dev/null || pwd | awk -F'gttl' '/gttl/ {print $$1 "gttl"}')
export PYTHONPATH=${GTTL}/scripts

include ${GTTL}/compile_flags.mk
include ${GTTL}/data_files.mk
include ${GTTL}/default_targets.mk
