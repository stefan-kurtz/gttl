GTTL ?= $(shell sh -c 'git rev-parse --show-toplevel 2>/dev/null || pwd | sed -n "s#^\(.*gttl\)/.*#\1#p"')
export PYTHONPATH=${GTTL}/scripts

include ${GTTL}/compile_flags.mk
include ${GTTL}/data_files.mk
include ${GTTL}/default_targets.mk
