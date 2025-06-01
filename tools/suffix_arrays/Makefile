CC?=gcc
CXX?=g++
CXX_PURE=$(shell echo ${CXX} | sed -e 's/^ccache //')
CXX_VERSION=$(shell basename ${CXX_PURE})
# PROFFLAGS= -DTEXT_OUTPUT -DSAINSHOWSTATEBU -DSAINSHOWSTATE
# PROFFLAGS= -DTEXT_OUTPUT
# PROFFLAGS= -DSAINSHOWSTATE
PROFFLAGS=
SYSTEM?=$(shell uname -s)
MACHINE?=$(shell uname -m)

CFLAGS=-Wall -Werror -Wextra -pedantic -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -DLARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DCXXOPTS_NO_REGEX -funroll-loops -O3 -g
LDFLAGS=-g ${PROFFLAGS} -ldl
CPPFLAGS=-I ${GTTL}/src
CXXFLAGS=-std=c++20
LDLIBS=-lz -lstdc++

ifeq ($(debug),yes)
  ifneq (${sanitize},no)
    ifneq (${MACHINE},arm64)
      ifneq (,$(findstring clang++,$(CXX)),clang++)
        CFLAGS += -fsanitize=address -fno-omit-frame-pointer
        LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
      endif
    endif
  endif
else
  CFLAGS += -DNDEBUG
  ifeq ($(SYSTEM),Darwin)
  else
    LDLIBS+=-ldl
  endif
endif

ifeq (${divsufsort},yes)
  DIVSUFSORTPATH=${PACKAGES}/libdivsufsort
  CFLAGS += -DWITHDIVSUFSORT -I${DIVSUFSORTPATH}/include
  LDFLAGS += -L${DIVSUFSORTPATH}/lib -ldivsufsort
endif

SA_INDUCED_OBJ=sa_induced.o\
               sa_induced_options.o

INPUTFILES=${shell ./enum_gttl_fasta_files.py ../../testdata}

all:sa_induced.x
ifndef GTTL
	$(error GTTL is undefined. Please define the Path to the gttl sources)
endif

# generic compilation rule which creates dependency file on the fly
%.o: %.cpp
	$(CXX) -DCXX_VERSION=\"${CXX_VERSION}\" -DCXX_FLAGS=\"'${CXXFLAGS} ${CFLAGS}'\" -c $< -o $@ $(CXXFLAGS) $(CFLAGS) ${CPPFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d)

sa_induced.x:${SA_INDUCED_OBJ}
	${CXX} -o $@ ${LDFLAGS} ${SA_INDUCED_OBJ} ${LDLIBS}

sa_reader.x:sa_reader.o
	${CXX} -o $@ ${LDFLAGS} $< ${LDLIBS}

sa_exact_pattern_search.x:sa_exact_pattern_search.o
	${CXX} -o $@ ${LDFLAGS} $< ${LDLIBS}

.PHONY:install
install:sa_induced.x
	mv $< ${HOME}/bin-ops/${CONFIGGUESS}/.

.PHONY:code_check
code_check:
	code_check.py -wt `find . -name '*.[ch]pp'`
	${GTTL}/scripts/check_ifndef.py `find . -name '*.hpp'`	

.PHONY:tags
tags:
	@ctags -w --c++-kinds=+p --fields=+iaSKlm --extra=+q \
	`find ${GTTL}/src -name '*.hpp'` *.[hc]pp > $@

.PHONY:text
test:test_error \
     test_plain_mixed_gzip \
     test_skylines \
     test_check_guess_and_multi \
     test_saincheck \
     test_plcp \
     test_sfx_compare \
     test_big
	@echo "$@ passed"

.PHONY:text_error
test_error:sa_induced.x
	@${GTTL}/testsuite/check_err.py -p ${GTTL} ./sa_induced.x
	@./check_sain_err.py ${GTTL}
	@echo "$@ passed"

.PHONY:test_plain_mixed_gzip
test_plain_mixed_gzip:sa_induced.x
	$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@gzip -9 -c ${GTTL}/testdata/at1MB.fna > ${TMPFILE}
	@./sa_induced.x --plain_input_format --indexname at1MB2Xpp -t ${GTTL}/testdata/at1MB.fna  ${GTTL}/testdata/at1MB.fna
	@./sa_induced.x --plain_input_format --indexname at1MB2Xpz -t ${GTTL}/testdata/at1MB.fna  ${TMPFILE}
	@./sa_induced.x --plain_input_format --indexname at1MB2Xzp -t ${TMPFILE} ${GTTL}/testdata/at1MB.fna
	@./sa_induced.x --plain_input_format --indexname at1MB2Xzz -t ${TMPFILE} ${TMPFILE}
	@cmp -s at1MB2Xpp.tis at1MB2Xpz.tis
	@cmp -s at1MB2Xpp.tis at1MB2Xzp.tis
	@cmp -s at1MB2Xpp.tis at1MB2Xzz.tis
	@${RM} ${TMPFILE} at1MB2X[pz][pz].tis at1MB2X[pz][pz].prj
	@echo "$@ passed"

.PHONY:text_skylines
test_skylines:sa_induced.x
	@./runskylines.sh
	@echo "$@ passed"

.PHONY:test_check_guess_and_multi
test_check_guess_and_multi:sa_induced.x
	@./check_guess_and_multi.sh
	@echo "$@ passed"

.PHONY:test_saincheck
test_saincheck:sa_induced.x
	@./saincheck.py --absolute_suftab
	@./saincheck.py --relative_suftab
	@./saincheck.py --lcptab kasai13n
	@./saincheck.py --absolute_suftab --lcptab kasai9n
	@./saincheck.py --absolute_suftab --lcptab plcp5n
	@./saincheck.py --reverse_complement --lcptab kasai13n
	@./saincheck.py --absolute_suftab --reverse_complement --lcptab kasai9n
	@./saincheck.py --absolute_suftab --reverse_complement --lcptab plcp5n
	@echo "$@ passed"

test_plcp:sa_induced.x
	@for filename in ${INPUTFILES} ; do \
	   echo "test_plcp.sh $${filename}"; \
	   ./test_plcp.sh $${filename}; \
	done
	@for filename in `enum_big_files.py ECOLI HIQ SOIL_UNIQUE`; do\
	  echo "test_plcp.sh $${filename}"; \
	  ./test_plcp.sh $${filename}; \
	done
	@echo "$@ passed"

test_sfx_compare:sa_induced.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@for filename in ${INPUTFILES}; do \
	   cp $$filename ${TMPFILE}; \
	   echo "sfx-compare.sh $$filename"; \
	   ./sfx-compare.sh ${TMPFILE}; \
	done
	@${RM} ${TMPFILE}
	@echo "$@ passed"

test_big:sa_induced.x
	@for filename in `enum_big_files.py ECOLI HIQ SOIL_UNIQUE`; do\
	  ./sa_induced.x --relative_suftab --verbose --check_suftab $${filename}; \
	  ${RM} `basename $${filename}`.bsf `basename $${filename}`.prj; \
	done
	@echo "$@ passed"

.PHONY:run_big
run_big:sa_induced.x
	@./run_big.py Protocols/*.tsv

# read dependencies
-include $(wildcard *.d)

.PHONY:clean
clean:
	@${RM} *.[odx] *.tis *.suf *.prj *.bsf .*.prj .*.bsf tmp.*
	@${RM} -r sa_induced.x.dSYM __pycache__
	@../../scripts/cleanpp.sh
