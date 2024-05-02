CXX?=g++
CFLAGS=-Wall -Werror -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -DLARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -funroll-loops -O3
CPPFLAGS=-I ../src
LD=${CXX}
LDLIBS=-lm -lz -lpthread -lstdc++
CXXFLAGS=-std=c++17
AT1MB=../testdata/at1MB.fna
SW175=../testdata/sw175.fna
Y3=../testdata/ychrIII.fna
VAC=../testdata/vaccg.fna
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
     test_fasta_iterator \
     test_multiseq \
     test_thread_pool \
     test_sort_kvt \
     test_sort_kvp \
     test_sort_i \
     test_eoplist \
     test_invint \
     test_line_iterator_mapped \
     test_minimizer_count \
     test_guess_if_protein_seq \
     test_fs_prio_store \
     test_split_string
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
	@./multiseq_mn.x --width 0 ../testdata/small.fna | diff -I '^#' - ../testdata/small.fna
	@./multiseq_mn.x --width 60 ${SW175} | diff -I '^#' - ${SW175}
	@./multiseq_mn.x --statistics --rankdist --protein ${SW175} | grep -v '^# TIME' | diff - ../testdata/sw175_stat.tsv
	@./multiseq_mn.x --width 70 ${AT1MB} | diff -I '^#' - ${AT1MB}
	@./multiseq_mn.x --statistics --rankdist ${AT1MB} | grep -v '^# TIME' | diff - ../testdata/at1MB_stat.tsv
	@./multiseq_mn.x --zipped --width 0 ../testdata/varlen_paired_1.fastq ../testdata/varlen_paired_2.fastq | grep -v '^#' | diff - ../testdata/varlen_paired_both.fasta
	@./multiseq_mn.x --width 70 --short_header ${AT1MB} | grep '^>' | diff - ../testdata/at1MB_short_header.txt
	@$(eval TMPFILE := $(shell mktemp))
	@echo "42" > ${TMPFILE}
	@./multiseq_mn.x --sample 42 --seed 34824347 --short_header --width 0 ${AT1MB} | grep -c '^>' | diff - ${TMPFILE}
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

.PHONY:test_fasta_iterator
test_fasta_iterator:./fasta_mn.x
	@for input_mode in "" --mapped; do \
	  ./fasta_mn.x $${input_mode} --width 70 ${AT1MB} | diff - ${AT1MB};\
	  ./fasta_mn.x $${input_mode} --width 60 ${SW175} | diff - ${SW175};\
	  ./fasta_mn.x $${input_mode} --width 61 ${Y3} | diff - ${Y3};\
	  ./fasta_mn.x $${input_mode} --width 61 ${VAC} | diff - ${VAC};\
	done
	@echo "Congratulations. $@ passed."

.PHONY:test_thread_pool
test_thread_pool:thread_pool_mn.x
	@./thread_pool_mn.x 4 40
	@echo "Congratulations. $@ passed."

.PHONY:test_sort_kvt
test_sort_kvt:sort_key_value_pairs.x
	@for num in 66 666 66666 666666 6666666 66666666; do \
	   ./sort_key_value_pairs.x -d t -m lsb-radix $$num;\
	   ./sort_key_value_pairs.x -d t -m mergesort $$num;\
	   ./sort_key_value_pairs.x -d t -m stdsort $$num;\
	done
	@echo "Congratulations. $@ passed."

.PHONY:test_sort_kvp
test_sort_kvp:sort_key_value_pairs.x
	@for num in 100 1000 10000 100000 1000000 10000000 100000000; do \
	   ./sort_key_value_pairs.x -d p -m lsb-radix $$num;\
	   ./sort_key_value_pairs.x -d p -m mergesort $$num;\
	   ./sort_key_value_pairs.x -d p -m stdsort $$num;\
	done
	@echo "Congratulations. $@ passed."

.PHONY:test_sort_i
test_sort_i:sort_key_value_pairs.x
	@for num in 200 2000 20000 200000 2000000 20000000 200000000; do \
	   ./sort_key_value_pairs.x -d i -m lsb-radix $$num;\
	   ./sort_key_value_pairs.x -d i -m mergesort $$num;\
	   ./sort_key_value_pairs.x -d i -m stdsort $$num;\
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

.PHONY:test_minimizer_count
test_minimizer_count:minimizer.x
	@./check_err.py ./minimizer.x
	@$(eval TMPFILE := $(shell mktemp))
	@echo "# number of hashed kmers	37465" > ${TMPFILE}
	@./minimizer.x -w 30 -k 18 -s ${AT} | grep 'number of hashed kmers' | diff - ${TMPFILE}
	@./minimizer.x -w 30 -k 18 -t 2 -s ${AT} | grep 'number of hashed kmers' | diff - ${TMPFILE}
	@echo "# number of hashed kmers	37577" > ${TMPFILE}
	@./minimizer.x -w 30 -k 18 -c -s ${AT} | grep 'number of hashed kmers' | diff - ${TMPFILE}
	@./minimizer.x -w 30 -k 18 -c -t 2 -s ${AT} | grep 'number of hashed kmers' | diff - ${TMPFILE}
	@./minimizer.x -w 30 -k 18 --show_mode 1 ${AT} > ${TMPFILE}
	@./minimizer.x -w 30 -k 18 --show_mode 2 ${AT} | diff -I '^#' - ${TMPFILE}
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

.PHONY:test_split_string
test_split_string:split_string_mn.x
	@$(eval TMPFILE := $(shell mktemp))
	@echo "A" > ${TMPFILE}
	@./split_string_mn.x '	' ${TMPFILE}
	@echo "A\tB" > ${TMPFILE}
	@./split_string_mn.x '	' ${TMPFILE}
	@echo "A\tB\nC" > ${TMPFILE}
	@./split_string_mn.x '	' ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed"

fastq_mn.x:fastq_mn.o seq_reader_options.o
	$(CXX) ${LDFLAGS} fastq_mn.o seq_reader_options.o -o $@ ${LDLIBS}

fasta_mn.x:fasta_mn.o seq_reader_options.o
	$(CXX) ${LDFLAGS} fasta_mn.o seq_reader_options.o -o $@ ${LDLIBS}

minimizer.x:minimizer_opt.o minimizer_mn.o
	$(CXX) ${LDFLAGS} minimizer_opt.o minimizer_mn.o -o $@ ${LDLIBS}

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
