CXX?=g++
CFLAGS=-Wall -Werror -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -DLARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -funroll-loops -O3
CPPFLAGS=-I ../src
LD=${CXX}
LDLIBS=-lm -lz -lpthread -lstdc++
CXXFLAGS=-std=c++17
AT1MB=../testdata/at1MB.fna
SYSTEM?=$(shell uname -s)
MACHINE?=$(shell uname -m)

ifeq ($(prof),yes)
  CFLAGS += -pg
  LDFLAGS += -pg
endif

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
  CFLAGS += -Wno-psabi
endif

MAIN_SRC=${shell grep -wl main *.cpp}
EXECS=${subst .cpp,.x,${MAIN_SRC}}

all:${EXECS}

test:test_alphabet \
     test_cycle_of_numbers \
     test_endian \
     test_nthash \
     test_bitpacker \
     test_char_range \
     test_matrix_partition \
     test_non_wc_ranges \
     test_line_iterator \
     test_fastq_iterator \
     test_multiseq \
     test_thread_pool \
     test_sort_kvt \
     test_sort_kvp \
     test_sort_i \
     test_eoplist \
     test_invint \
     test_line_iterator_mapped \
     test_nt_minimizer_count \
     test_guess_if_protein_seq \
     test_fs_prio_store
	@echo "Congratulations. $@ passed."

test_endian:./endian.x
	@./endian.x
	@echo "Congratulations. $@ passed."

test_bitpacker:./bitpacker_main.x
	@./bitpacker_main.x 100
	@echo "Congratulations. $@ passed."

test_non_wc_ranges:./count_non_wc_ranges.x
	@./check_err.py --fail_for_protein_sequence ./count_non_wc_ranges.x
	@./count_non_wc_ranges.x ${AT1MB} | diff - ../testdata/at1MB_non_wc.tsv
	@echo "Congratulations. $@ passed."

test_alphabet:./alphabet.x
	@./alphabet.x | diff - ../testdata/DNA_AA.tsv
	@echo "Congratulations. $@ passed."

test_cycle_of_numbers:./cycle_of_numbers_mn.x
	@$(eval TMPFILE := $(shell mktemp))
	@./cycle_of_numbers.py > ${TMPFILE}
	@./cycle_of_numbers_mn.x | sort -nu | diff - ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed."

test_char_range:char_range_mn.x
	@./char_range_compare.sh ../testdata/small.fna
	@./char_range_compare.sh ../testdata/simple.fna
	@./char_range_compare.sh ${AT1MB}
	@echo "Congratulations. $@ passed."

test_nthash:enum_nthash.x test_nthash_wc_palindrome
	@./check_err.py --fail_for_protein_sequence ./enum_nthash.x
	@./enum_nthash.x --bytes_unit ${AT1MB} | \
           grep -v '^# TIME' | diff - ../testdata/at1MB_nthash.tsv
	@./enum_nthash.x --with_rc --bytes_unit ${AT1MB} | \
           grep -v '^# TIME' | grep -v '_rc' | \
           diff - ../testdata/at1MB_nthash.tsv
	@./enum_nthash.x --with_rc --bytes_unit ${AT1MB} | \
           grep '_rc' | diff - ../testdata/at1MB_nthash_rc.tsv
	@./enum_nthash.x --bytes_unit ../testdata/vaccg.fna | \
           grep -v '^# TIME' | diff - ../testdata/vaccg_nthash.tsv
	@./enum_nthash.x --with_rc --bytes_unit ../testdata/vaccg.fna | \
           grep -v '^# TIME' | grep -v '_rc' | \
           diff - ../testdata/vaccg_nthash.tsv
	@./enum_nthash.x --bytes_unit ../testdata/ychrIII.fna | \
           grep -v '^# TIME' | diff - ../testdata/ychrIII_nthash.tsv
	@./enum_nthash.x --with_rc --bytes_unit ../testdata/ychrIII.fna | \
           grep -v '^# TIME' | grep -v '_rc' | \
           diff - ../testdata/ychrIII_nthash.tsv
	@echo "Congratulations. $@ passed."

.PHONY:test_nthash_wc_palindrome
test_nthash_wc_palindrome:enum_nthash.x
	@$(eval TMPFILE := $(shell mktemp))
	@./enum_nthash.x -k 8 -s --with_rc ../testdata/wc_palindromes8.fna | grep -v '^#' | cut -f 1 > ${TMPFILE}
	@./enum_nthash.x -k 8 -s --with_rc ../testdata/wc_palindromes8.fna | grep -v '^#' | cut -f 2 | diff - ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed."

.PHONY:test_multiseq
test_multiseq:./multiseq_mn.x
	@./check_err.py ./multiseq_mn.x
	@./multiseq_mn.x -w 0 ../testdata/small.fna | diff -I '^#' - ../testdata/small.fna
	@./multiseq_mn.x -w 60 ../testdata/sw175.fna | diff -I '^#' - ../testdata/sw175.fna
	@./multiseq_mn.x -r -p ../testdata/sw175.fna | grep -v '^# TIME' | diff - ../testdata/sw175_stat.tsv
	@./multiseq_mn.x -w 70 ${AT1MB} | diff -I '^#' - ${AT1MB}
	@./multiseq_mn.x -r ${AT1MB} | grep -v '^# TIME' | diff - ../testdata/at1MB_stat.tsv
	@./multiseq_mn.x -z -w 0 ../testdata/varlen_paired_1.fastq ../testdata/varlen_paired_2.fastq | grep -v '^#' | diff - ../testdata/varlen_paired_both.fasta
	@./multiseq_mn.x -w 70 -s ${AT1MB} | grep '^>' | diff - ../testdata/at1MB_short_header.txt
	@echo "Congratulations. $@ passed."

.PHONY:test_line_iterator
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

.PHONY:test_matrix_partition
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

.PHONY:test_fastq_iterator
test_fastq_iterator:./fastq_mn.x
	@for num_threads in 1 2 3 4 5; do \
	  ./fastq_mn.x --threads $${num_threads} ../testdata/70x_161nt_phred64.fastq | grep '^#' | diff - ../testdata/70x_161nt_phred64_dist.tsv;\
	done
	@./test_fastq_iterator.sh
	@./check_split.sh ../testdata/70x_161nt_phred64.fastq
	@echo "Congratulations. $@ passed."

.PHONY:test_thread_pool
test_thread_pool:thread_pool_mn.x
	@./thread_pool_mn.x 4 40
	@echo "Congratulations. $@ passed."

.PHONY:test_sort_kvt
test_sort_kvt:sort_key_value_pairs.x
	@for num in 66 666 66666 666666 6666666 66666666; do \
	   ./sort_key_value_pairs.x t rad $$num;\
	   ./sort_key_value_pairs.x t std $$num;\
	done
	@echo "Congratulations. $@ passed."

.PHONY:test_sort_kvp
test_sort_kvp:sort_key_value_pairs.x
	@for num in 100 1000 10000 100000 1000000 10000000 100000000; do \
	   ./sort_key_value_pairs.x p rad $$num;\
	   ./sort_key_value_pairs.x p std $$num;\
	done
	@echo "Congratulations. $@ passed."

.PHONY:test_sort_i
test_sort_i:sort_key_value_pairs.x
	@for num in 200 2000 20000 200000 2000000 20000000 200000000; do \
	   ./sort_key_value_pairs.x i rad $$num;\
	   ./sort_key_value_pairs.x i std $$num;\
	done
	@echo "Congratulations. $@ passed."

.PHONY:test_eoplist
test_eoplist:eoplist_mn.x
	@./eoplist_mn.x silent 100000
	@./eoplist_mn.x display 1000
	@echo "Congratulations. $@ passed."

.PHONY:test_invint
test_invint:enum_invint.x
	@./enum_invint.x 15 ../testdata/vaccg.fna | diff - ../testdata/vaccg_invint15.tsv
	@./enum_invint.x 15 ../testdata/ychrIII.fna | diff - ../testdata/ychrIII_invint15.tsv
	@./enum_invint.x 14 ../testdata/at1MB.fna | diff - ../testdata/at1MB_invint14.tsv
	@echo "Congratulations. $@ passed."

.PHONY:test_line_iterator_mapped
test_line_iterator_mapped:line_iterator_mapped.x
	@./line_iterator_mapped.x `ls ../testdata/*.fna | egrep -v 'empty.fna|\.gz'`
	@./line_iterator_mapped.x *.[hc]pp
	@echo "Congratulations. $@ passed."

.PHONY:test_nt_minimizer_count
test_nt_minimizer_count:nt_minimizer.x
	@$(eval TMPFILE := $(shell mktemp))
	@echo "# number of hashed kmers	37465" > ${TMPFILE}
	@./nt_minimizer.x -l 47 -k 18 -s ${AT} | grep 'number of hashed kmers' | diff - ${TMPFILE}
	@./nt_minimizer.x -l 47 -k 18 -t 2 -s ${AT} | grep 'number of hashed kmers' | diff - ${TMPFILE}
	@echo "# number of hashed kmers	37577" > ${TMPFILE}
	@./nt_minimizer.x -l 47 -k 18 -c -s ${AT} | grep 'number of hashed kmers' | diff - ${TMPFILE}
	@./nt_minimizer.x -l 47 -k 18 -c -t 2 -s ${AT} | grep 'number of hashed kmers' | diff - ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed"

.PHONY:test_guess_if_protein_seq
test_guess_if_protein_seq:./guess_if_protein_seq.x
	./test_guess_if_protein_seq.sh
	@echo "Congratulations. $@ passed"

.PHONY:test_fs_prio_store
test_fs_prio_store:./fs_prio_store_mn.x
	@./fs_prio_store_mn.x 1000 1000 1000 100
	@./fs_prio_store_mn.x 1000 1000 1000 10000
	@./fs_prio_store_mn.x 1000 2000 10000000 1000
	@echo "Congratulations. $@ passed"

fastq_mn.x:fastq_mn.o fastq_opt.o
	$(CXX) ${LDFLAGS} fastq_mn.o fastq_opt.o -o $@ ${LDLIBS}

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
