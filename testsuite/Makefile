CXX?=g++
CFLAGS=-Wall -Werror -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -DLARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -funroll-loops -O3
CPPFLAGS=-I ../src
LD=${CXX}
LDLIBS=-lm -lz -lpthread
CXXFLAGS=-std=c++17
SYSTEM?=$(shell uname -s)
MACHINE?=$(shell uname -m)

ifeq ($(debug),yes)
  CFLAGS += -g
  LDFLAGS += -g
  ifneq (${sanitize},no)
    ifneq (${MACHINE},arm64)
      ifneq ($(findstring clang++,$(CXX)),clang++)
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

ifeq ($(findstring clang++,$(CXX)),clang++)
  CFLAGS += -Wno-psabi
else
  CFLAGS += -time -Wno-psabi
endif

SRC=${wildcard  *.cpp}
EXECS=${subst .cpp,.x,${SRC}}

all:${EXECS}
	echo ${compiler_type}

test:${EXECS}
	@./alphabet.x | diff - ../testdata/DNA_AA.tsv
	@./check_err.py --fail_for_protein_sequence ./enum_nthash.x
	@./check_err.py --fail_for_protein_sequence ./count_non_wc_ranges.x
	@./check_err.py ./multiseq_mn.x
	@./multiseq_mn.x -w 0 ../testdata/small.fna | diff -I '^#' - ../testdata/small.fna
	@./multiseq_mn.x -w 60 ../testdata/sw175.fna | diff -I '^#' - ../testdata/sw175.fna
	@./multiseq_mn.x -p ../testdata/sw175.fna | grep -v '^# TIME' | diff - ../testdata/sw175_stat.tsv
	@./multiseq_mn.x -w 70 ../testdata/at1MB.fna | diff -I '^#' - ../testdata/at1MB.fna
	@./multiseq_mn.x ../testdata/at1MB.fna | grep -v '^# TIME' | diff - ../testdata/at1MB_stat.tsv
	@./enum_nthash.x ../testdata/at1MB.fna | diff - ../testdata/at1MB_nthash.tsv
	@./count_non_wc_ranges.x ../testdata/at1MB.fna | diff - ../testdata/at1MB_non_wc.tsv
	@./bitpacker_main.x 100
	@./endian.x
	@./char_range_compare.sh ../testdata/small.fna
	@./char_range_compare.sh ../testdata/simple.fna
	@./char_range_compare.sh ../testdata/at1MB.fna
	@echo "Congratulations. $@ passed."

%.o:%.cpp
	$(CXX) $(CXXFLAGS) ${CFLAGS} ${CPPFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d) $< -c -o $@

%.x:%.o
	$(CXX) ${LDFLAGS} $< -o $@ ${LDLIBS}

.PHONY:code_check
code_check:
	code_check.py -wt `find . -name '*.cpp'`

.PHONY:clean
clean:
	${RM} -r *.[odx] *.x.dSYM/ __pycache__

-include ${wildcard *.d}
