GTTL ?= ../..
include ${GTTL}/config.mk

# PROFFLAGS= -DTEXT_OUTPUT -DSAINSHOWSTATEBU -DSAINSHOWSTATE
# PROFFLAGS= -DTEXT_OUTPUT
# PROFFLAGS= -DSAINSHOWSTATE
PROFFLAGS=

LDFLAGS += ${PROFFLAGS}

ifeq (${divsufsort},yes)
  DIVSUFSORTPATH=${PACKAGES}/libdivsufsort
  CFLAGS += -DWITHDIVSUFSORT -I${DIVSUFSORTPATH}/include
  LDFLAGS += -L${DIVSUFSORTPATH}/lib -ldivsufsort
endif

SA_INDUCED_OBJ=sa_induced.o\
               sa_induced_options.o

.PHONY:all
all:sa_induced.x lcp_checker.x
ifndef GTTL
	$(error GTTL is undefined. Please define the Path to the gttl sources)
endif

sa_induced.x:${SA_INDUCED_OBJ}
	${CXX} -o $@ ${LDFLAGS} ${SA_INDUCED_OBJ} ${LDLIBS}

lcp_checker.x:lcp_checker.o
	${CXX} -o $@ ${LDFLAGS} $< ${LDLIBS}

sa_exact_pattern_search.x:sa_exact_pattern_search.o
	${CXX} -o $@ ${LDFLAGS} $< ${LDLIBS}

.PHONY:install
install:sa_induced.x
	mv $< ${HOME}/bin-ops/${CONFIGGUESS}/.

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
	@FASTA_FILES="$(GTTL_FASTA_FILES)" ./saincheck.py --absolute_suftab
	@FASTA_FILES="$(GTTL_FASTA_FILES)" ./saincheck.py --relative_suftab
	@FASTA_FILES="$(GTTL_FASTA_FILES)" ./saincheck.py --lcptab kasai13n
	@FASTA_FILES="$(GTTL_FASTA_FILES)" ./saincheck.py --absolute_suftab --lcptab kasai9n
	@FASTA_FILES="$(GTTL_FASTA_FILES)" ./saincheck.py --absolute_suftab --lcptab plcp5n
	@FASTA_FILES="$(GTTL_FASTA_FILES)" ./saincheck.py --reverse_complement --lcptab kasai13n
	@FASTA_FILES="$(GTTL_FASTA_FILES)" ./saincheck.py --absolute_suftab --reverse_complement --lcptab kasai9n
	@FASTA_FILES="$(GTTL_FASTA_FILES)" ./saincheck.py --absolute_suftab --reverse_complement --lcptab plcp5n
	@echo "$@ passed"

test_plcp:sa_induced.x
	@for filename in ${GTTL_FASTA_FILES} ; do \
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
	@for filename in ${GTTL_FASTA_FILES}; do \
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
