# If no path is defined in the environment, we attempt to use the git root, or otherwise
# the folder named "gttl" in pwd.
GTTL ?= $(shell git rev-parse --show-toplevel 2>/dev/null || pwd | awk -F'gttl' '/gttl/ {print $$1 "gttl"}')
# Then, no matter what the variable previously was, we forcibly convert it
# to a relative path. This fixes issues where file-paths are diff'd in testing scripts.
GTTL := $(shell python3 -c "import os; print(os.path.relpath('$(GTTL)', os.getcwd()).replace('\\\\\\\\', '/'))")

export PYTHONPATH=${GTTL}/scripts

include ${GTTL}/compile_flags.mk
include ${GTTL}/data_files.mk
include ${GTTL}/default_targets.mk
