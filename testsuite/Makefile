CXX?=g++
CFLAGS=-Wall -Werror -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -DLARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -funroll-loops -O3
CPPFLAGS=-I ../src
LD=${CXX}
LDLIBS=-lm -lz -lpthread
CXXFLAGS=-std=c++17
AT1MB=../testdata/at1MB.fna
SYSTEM?=$(shell uname -s)
MACHINE?=$(shell uname -m)

ifeq ($(debug),yes)
  CFLAGS += -g
  LDFLAGS += -g
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

ifeq (,$(findstring g++,$(CXX)))
  CFLAGS += -time -Wno-psabi
endif

SRC=${wildcard  *.cpp}
EXECS=${subst .cpp,.x,${SRC}}

all:${EXECS}

test:${EXECS} test_matrix_partition test_line_iterator
	@./alphabet.x | diff - ../testdata/DNA_AA.tsv
	@./check_err.py --fail_for_protein_sequence ./enum_nthash.x
	@./check_err.py --fail_for_protein_sequence ./count_non_wc_ranges.x
	@./check_err.py ./multiseq_mn.x
	@./multiseq_mn.x -w 0 ../testdata/small.fna | diff -I '^#' - ../testdata/small.fna
	@./multiseq_mn.x -w 60 ../testdata/sw175.fna | diff -I '^#' - ../testdata/sw175.fna
	@./multiseq_mn.x -r -p ../testdata/sw175.fna | grep -v '^# TIME' | diff - ../testdata/sw175_stat.tsv
	@./multiseq_mn.x -w 70 ${AT1MB} | diff -I '^#' - ${AT1MB}
	@./multiseq_mn.x -r ${AT1MB} | grep -v '^# TIME' | diff - ../testdata/at1MB_stat.tsv
	@./enum_nthash.x ${AT1MB} | diff - ../testdata/at1MB_nthash.tsv
	@./count_non_wc_ranges.x ${AT1MB} | diff - ../testdata/at1MB_non_wc.tsv
	@./bitpacker_main.x 100
	@./endian.x
	@./char_range_compare.sh ../testdata/small.fna
	@./char_range_compare.sh ../testdata/simple.fna
	@./char_range_compare.sh ${AT1MB}
	@echo "Congratulations. $@ passed."

test_line_iterator:line_iterator.x
	@./line_iterator.py `find .. -name '*.[ch]pp'` | sh -s
	@./line_iterator.x ${AT1MB} | diff - ${AT1MB}
	@$(eval TMPFILE := $(shell mktemp))
	@cat ${AT1MB} | tr -d '\n' > ${TMPFILE}
	@echo "" >> ${TMPFILE}
	@./line_iterator.x ${TMPFILE} | diff - ${TMPFILE}
	@cat ${AT1MB} ${AT1MB} > ${TMPFILE}
	@./line_iterator.x --all ${AT1MB} ${AT1MB} | diff - ${TMPFILE}
	@./line_iterator.x ${AT1MB} ${AT1MB} | diff - ${TMPFILE}
	@cat ${AT1MB} ${AT1MB} ${AT1MB} > ${TMPFILE}
	@./line_iterator.x --all ${AT1MB} ${AT1MB} ${AT1MB} | diff - ${TMPFILE}
	@./line_iterator.x ${AT1MB} ${AT1MB} ${AT1MB} | diff - ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed."

test_matrix_partition:matrix_partition.py matrix_partition.x
	./matrix_partition_cmp.sh 10 100 99	
	./matrix_partition_cmp.sh 17 511 137
	./matrix_partition_cmp.sh 97 10000 16000
	./matrix_partition_cmp.sh 250 50000 43001
	./matrix_partition_cmp.sh 10 100 0
	./matrix_partition_cmp.sh 17 400 0
	./matrix_partition_cmp.sh 19 998 0
	./matrix_partition_cmp.sh 500 13000 0
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
