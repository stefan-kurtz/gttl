CXX?=g++
CFLAGS=-g -m64 -Wall -Werror -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -funroll-loops -O3 -pedantic
CPPFLAGS=-I../src/
CXXFLAGS=-std=c++17
LDFLAGS=-g -m64
LDLIBS=-lm -lz -lpthread
SYSTEM?=$(shell uname -s)

ifeq ($(debug),yes)
  CFLAGS += -g
else
  CFLAGS += -DNDEBUG -march=native
endif

ifeq ($(SYSTEM),Darwin)
  ifeq ($(debug),yes)
    CFLAGS += -fsanitize=address -fno-omit-frame-pointer
    LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
  endif
endif

SRC=${wildcard  *.cpp}
EXECS=${subst .cpp,.x,${SRC}}

all:${EXECS}

test:${EXECS}
	@./alphabet.x | diff - ../testdata/DNA_AA.tsv
	@./check_err.py --fail_for_protein_sequence enum_nthash.x
	@./check_err.py --fail_for_protein_sequence count_non_wc_ranges.x
	@./check_err.py multiseq_mn.x
	@./multiseq_mn.x -w 0 ../testdata/small.fna | diff -I '^#' - ../testdata/small.fna
	@./multiseq_mn.x -w 60 ../testdata/sw175.fna | diff -I '^#' - ../testdata/sw175.fna
	@./multiseq_mn.x ../testdata/sw175.fna | grep -v '^# TIME' | diff - ../testdata/sw175_stat.tsv
	@./multiseq_mn.x -w 70 ../testdata/at1MB.fna | diff -I '^#' - ../testdata/at1MB.fna
	@./multiseq_mn.x ../testdata/at1MB.fna | grep -v '^# TIME' | diff - ../testdata/at1MB_stat.tsv
	@./enum_nthash.x ../testdata/at1MB.fna | diff - ../testdata/at1MB_nthash.tsv
	@./count_non_wc_ranges.x ../testdata/at1MB.fna | diff - ../testdata/at1MB_non_wc.tsv 
	@echo "Congratulations. Tests passed."

%.o:%.cpp
	$(CXX) $(CXXFLAGS) ${CFLAGS} ${CPPFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d) $< -c -o $@

%.x:%.o
	$(CXX) ${LDFLAGS} $< -o $@ ${LDLIBS}

.PHONY:code_check
code_check:
	code_check.py -wt `find . -name '*.cpp'`

.PHONY:clean
clean:
	${RM} -r *.[odx] *.x.dSYM/

-include ${wildcard *.d}
