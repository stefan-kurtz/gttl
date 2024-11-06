CXX?=g++
LD=${CXX}
CFLAGS=-m64 -Wall -Werror -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -funroll-loops -O3 -pedantic
GTTL?=../..
CPPFLAGS=-I ${GTTL}/src
CXXFLAGS=-std=c++17
LDFLAGS=-g -m64
SYSTEM?=$(shell uname -s)
MACHINE?=$(shell uname -m)
OBJ=unwords_mn.o unwords_opt.o
AT=${GTTL}/testdata/at1MB.fna
VAC=${GTTL}/testdata/vaccg.fna
Y3=${GTTL}/testdata/ychrIII.fna
PROTEIN_FILE=${GTTL}/testdata/protein.fsa
SW_175=${GTTL}/testdata/sw175.fna
IRC_OPTION?=--ignore_rc
PYTHONPATH=../../scripts

ifeq ($(OS),Windows_NT)
	CFLAGS += -D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_WARNINGS
	LDLIBS+=-L C:\\vcpkg\\packages\\zlib_x64-windows\\lib -lzlib
else
	LDLIBS+=-lm -lz -lpthread -lstdc++
  SHELL=/bin/bash
endif

export PYTHONPATH=../../scripts

ifeq (,$(findstring g++-,$(CXX)))
  CFLAG += -march=native
  TIME_OPTIONS=-time
endif

ifeq ($(debug),yes)
  CFLAGS += -g
  LDFLAGS += -g
  ifneq (${sanitize},no)
    ifneq (${MACHINE},arm64)
      CFLAGS += -fsanitize=address -fno-omit-frame-pointer
      LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
    endif
  endif
else
  CFLAGS += -DNDEBUG
  ifeq ($(SYSTEM),Darwin)
  else
    ifeq ($(OS),Windows_NT)
    else
      LDLIBS+=-ldl
    endif
  endif
endif

ifeq ($(prof),yes)
  CFLAGS += -pg
  LDFLAGS += -pg
endif

SRC=${wildcard *.cpp}

unwords_mn.x:${OBJ}
	$(LD) ${LDFLAGS} ${OBJ} -o $@ ${LDLIBS}

%.o:%.cpp
	$(CXX) ${TIME_OPTION} -o $@ -c $< ${CXXFLAGS} ${CFLAGS} ${CPPFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d)

.PHONY:test
test: check_unwords_sw175 test_unwords_length check_unwords_p check_unwords_vac check_unwords_y3 check_unwords_at \
      compare_unwords_vac compare_unwords_y3 compare_unwords_at check_deBruijn
	@echo "$@ passed"

.PHONY:check_deBruijn
check_deBruijn: check_deBruijn1 check_deBruijn3 check_deBruijn10
	@echo "Congratulations. $@ passed!"

.PHONY:check_deBruijn3
check_deBruijn3:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@echo "> de Bruijn sequence for q=3" > ${TMPFILE}
	@${GTTL}/scripts/de_Bruijn.py -q 3 --dna >> ${TMPFILE}
	@echo "AA" >> ${TMPFILE}
	@./test_de_bruijn.py -q 3 ${TMPFILE}
	@${RM} ${TMPFILE}

.PHONY:check_deBruijn1
check_deBruijn1:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@echo "> de Bruijn sequence for q=1" > ${TMPFILE}
	@${GTTL}/scripts/de_Bruijn.py -q 1 --dna >> ${TMPFILE}
	@echo "AAA" >> ${TMPFILE}
	@./test_de_bruijn.py -q 1 ${TMPFILE}
	@${RM} ${TMPFILE}

.PHONY:check_deBruijn10
check_deBruijn10:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@echo "> de Bruijn sequence for q=10" > ${TMPFILE}
	@${GTTL}/scripts/de_Bruijn.py -q 10 --dna >> ${TMPFILE}
	@echo "AAAAAAAAA" >> ${TMPFILE}
	@./test_de_bruijn.py -q 10 ${TMPFILE}
	@${RM} ${TMPFILE}

.PHONY:check_unwords_vac
check_unwords_vac:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@./unwords_mn.x ${IRC_OPTION} ${VAC} > ${TMPFILE}
	@./check_unwords.py ${TMPFILE} ${VAC}
	@${RM} ${TMPFILE}

.PHONY:check_unwords_y3
check_unwords_y3:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@./unwords_mn.x ${IRC_OPTION} ${Y3} > ${TMPFILE}
	@./check_unwords.py ${TMPFILE} ${Y3}
	@${RM} ${TMPFILE}

.PHONY:check_unwords_at
check_unwords_at:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@./unwords_mn.x ${IRC_OPTION} ${AT} > ${TMPFILE}
	@./check_unwords.py ${TMPFILE} ${AT}
	@${RM} ${TMPFILE}

.PHONY:check_unwords_p
check_unwords_p:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@./unwords_mn.x ${PROTEIN_FILE} > ${TMPFILE}
	@./check_unwords.py ${TMPFILE} ${PROTEIN_FILE}
	@${RM} ${TMPFILE}

.PHONY:check_unwords_sw175
check_unwords_sw175:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@./unwords_mn.x ${SW_175} > ${TMPFILE}
	@./check_unwords.py ${TMPFILE} ${SW_175}
	@${RM} ${TMPFILE}

.PHONY:compare_unwords_vac
compare_unwords_vac:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE1 := $(shell mktemp --tmpdir=.))
	@./unwords_mn.x ${IRC_OPTION} ${VAC} > ${TMPFILE1}
	@$(eval TMPFILE2 := $(shell mktemp --tmpdir=.))
	@diff --strip-trailing-cr -I '^[#>]' -i ${TMPFILE1} ${GTTL}/testdata/unwords_vac.txt
	@${RM} ${TMPFILE1}
	@echo "Congratulations. $@ passed"

.PHONY:compare_unwords_y3
compare_unwords_y3:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE1 := $(shell mktemp --tmpdir=.))
	@./unwords_mn.x ${IRC_OPTION} ${Y3} > ${TMPFILE1}
	@$(eval TMPFILE2 := $(shell mktemp --tmpdir=.))
	@diff --strip-trailing-cr -I '^[#>]' -i ${TMPFILE1} ${GTTL}/testdata/unwords_y3.txt
	@${RM} ${TMPFILE1}
	@echo "Congratulations. $@ passed"

.PHONY:compare_unwords_at
compare_unwords_at:unwords_mn.x
	@echo "$@"
	@$(eval TMPFILE1 := $(shell mktemp --tmpdir=.))
	@./unwords_mn.x ${IRC_OPTION} ${AT} > ${TMPFILE1}
	@$(eval TMPFILE2 := $(shell mktemp --tmpdir=.))
	@diff --strip-trailing-cr -I '^[#>]' -i ${TMPFILE1} ${GTTL}/testdata/unwords_at.txt
	@${RM} ${TMPFILE1}
	@echo "Congratulations. $@ passed"

.PHONY:test_runtime
test_runtime: unwords_mn.x
	@echo "$@"
	@./test_runtime.py --computer_data --unwords_mn

.PHONY:test_runtime_protein
test_runtime_protein: unwords_mn.x
	@echo "$@"
	@./test_runtime.py --computer_data --unwords_mn --protein

.PHONY:compare_runtime
compare_runtime: unwords_mn.x
	@echo "$@"
	@./test_runtime.py --computer_data --unwords_mn --unwords_bits

.PHONY:test_unwords_length
test_unwords_length: unwords_mn.x
	@./test_unwords_length.py ${GTTL}/testdata

.PHONY:tags
tags:
	@ctags -w --c++-kinds=+p --fields=+iaSKlm --extra=+q `ls *.[hc]pp` \
        `find ${GTTL}/src -name '*.hpp'` > $@

.PHONY:code_check
code_check:
	@${GTTL}/scripts/code_check.py -wt `find . -name '*.[hc]pp'`
	@${GTTL}/scripts/code_check.py -wt `find . -name '*.py'`

.PHONY:clean
clean:
	${RM} -r *.[oxd] *.x.dSYM/ tmp.*

-include ${wildcard *.d}
