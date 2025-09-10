GTTL := ..
include ${GTTL}/config.mk

MAIN_SRC=${shell grep -wl main *.cpp}
EXECS=${subst .cpp,.x,${MAIN_SRC}}

all:${EXECS}

ifneq ($(OS),Windows_NT)
test: test_untar
endif

test:test_alphabet \
     test_multibitvector \
     test_cycle_of_numbers \
     test_endian \
     test_nthash \
     test_plane_sweep_filter \
     test_bitpacker \
     test_char_range \
     test_matrix_partition \
     test_non_wc_ranges \
     test_line_generator \
     test_fastq_generator \
     test_fasta_generator \
     test_multiseq \
     test_thread_pool \
     test_sort_kvt \
     test_sort_kvp \
     test_sort_i \
     test_eoplist \
     test_invint \
     test_line_generator_mapped \
     test_minimizer_count \
     test_guess_if_protein_seq \
     test_fs_prio_store \
     test_rdbuf \
     test_split_string \
     test_split_files \
     test_multiseq_factory \
     test_read_vector \
     test_string_of_digits \
     test_binary_iterator \
     test_multiple_options
	@echo "Congratulations. $@ passed."
ifeq ($(valgrind),yes)
	@echo "Valgrind could find no issues in the executables it tested"
endif

test_endian:./endian.x
	@${VALGRIND} ./endian.x
	@echo "Congratulations. $@ passed."

test_bitpacker:./bitpacker_main.x
	@./bitpacker_main.x 100
	@echo "Congratulations. $@ passed."

test_non_wc_ranges:./count_non_wc_ranges.x
	@./check_err.py --fail_for_protein_sequence ./count_non_wc_ranges.x
	@${VALGRIND} ./count_non_wc_ranges.x ${AT1MB} | diff --strip-trailing-cr - ../testdata/at1MB_non_wc.tsv
	@echo "Congratulations. $@ passed."

test_alphabet:./alphabet.x
	@${VALGRIND} ./alphabet.x | diff --strip-trailing-cr - ../testdata/DNA_AA.tsv
	@echo "Congratulations. $@ passed."

test_plane_sweep_filter:./plane_sweep_filter.x
	@${VALGRIND} ./plane_sweep_filter.x | diff --strip-trailing-cr - ../testdata/non_redundant_matches.txt
	@echo "Congratulations. $@ passed."

test_cycle_of_numbers:./cycle_of_numbers_mn.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@./cycle_of_numbers.py > ${TMPFILE}
	@${VALGRIND} ./cycle_of_numbers_mn.x | sort -nu | diff --strip-trailing-cr - ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed."

test_char_range:char_range_mn.x
	@./char_range_compare.sh ../testdata/small.fna
	@./char_range_compare.sh ../testdata/simple.fna
	@./char_range_compare.sh ${AT1MB}
	@echo "Congratulations. $@ passed."

test_nthash:enum_nthash.x test_nthash_wc_palindrome
	@./check_err.py ./enum_nthash.x
	@${VALGRIND} ./enum_nthash.x --bytes_unit ${AT1MB} | \
           grep -v '^# TIME' | diff --strip-trailing-cr - ../testdata/at1MB_nthash.tsv
	@${VALGRIND} ./enum_nthash.x --with_rc --bytes_unit ${AT1MB} | \
           grep -v '^# TIME' | grep -v '_rc' | \
           diff --strip-trailing-cr - ../testdata/at1MB_nthash.tsv
	@${VALGRIND} ./enum_nthash.x --with_rc --bytes_unit ${AT1MB} | \
           grep '_rc' | diff --strip-trailing-cr - ../testdata/at1MB_nthash_rc.tsv
	@${VALGRIND} ./enum_nthash.x --bytes_unit ../testdata/vaccg.fna | \
           grep -v '^# TIME' | diff --strip-trailing-cr - ../testdata/vaccg_nthash.tsv
	@${VALGRIND} ./enum_nthash.x --with_rc --bytes_unit ../testdata/vaccg.fna | \
           grep -v '^# TIME' | grep -v '_rc' | \
           diff --strip-trailing-cr - ../testdata/vaccg_nthash.tsv
	@${VALGRIND} ./enum_nthash.x --bytes_unit ../testdata/ychrIII.fna | \
           grep -v '^# TIME' | diff --strip-trailing-cr - ../testdata/ychrIII_nthash.tsv
	@${VALGRIND} ./enum_nthash.x --with_rc --bytes_unit ../testdata/ychrIII.fna | \
           grep -v '^# TIME' | grep -v '_rc' | \
           diff --strip-trailing-cr - ../testdata/ychrIII_nthash.tsv
	@${VALGRIND} ./enum_nthash.x --bytes_unit ../testdata/protein.fsa | \
                     grep -v '^# TIME' | \
                     diff --strip-trailing-cr - ../testdata/protein_nthash.tsv
	@echo "Congratulations. $@ passed."

.PHONY:test_nthash_wc_palindrome
test_nthash_wc_palindrome:enum_nthash.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@${VALGRIND} ./enum_nthash.x -k 8 -s --with_rc ../testdata/wc_palindromes8.fna | grep -v '^#' | cut -f 1 > ${TMPFILE}
	@${VALGRIND} ./enum_nthash.x -k 8 -s --with_rc ../testdata/wc_palindromes8.fna | grep -v '^#' | cut -f 2 | diff --strip-trailing-cr - ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed."

.PHONY:test_multiseq
test_multiseq:./multiseq_mn.x test_sorted_by_header_multiseq
	@./check_err.py ./multiseq_mn.x
	@${VALGRIND} ./multiseq_mn.x --width 0 ../testdata/small.fna | diff --strip-trailing-cr -I '^#' - ../testdata/small.fna
	@${VALGRIND} ./multiseq_mn.x --width 60 ${SW175} | diff --strip-trailing-cr -I '^#' - ${SW175}
	@${VALGRIND} ./multiseq_mn.x --statistics --rankdist --protein ${SW175} | grep -v '^# TIME' | diff --strip-trailing-cr - ../testdata/sw175_stat.tsv
	@${VALGRIND} ./multiseq_mn.x --width 70 ${AT1MB} | diff --strip-trailing-cr -I '^#' - ${AT1MB}
	@${VALGRIND} ./multiseq_mn.x --statistics --rankdist ${AT1MB} | grep -v '^# TIME' | diff --strip-trailing-cr - ../testdata/at1MB_stat.tsv
	@${VALGRIND} ./multiseq_mn.x --zipped --width 0 ../testdata/varlen_paired_1.fastq ../testdata/varlen_paired_2.fastq | grep -v '^#' | diff --strip-trailing-cr - ../testdata/varlen_paired_both.fasta
	@${VALGRIND} ./multiseq_mn.x --width 70 --short_header ${AT1MB} | grep '^>' | diff --strip-trailing-cr - ../testdata/at1MB_short_header.txt
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@echo "42" > ${TMPFILE}
	@${VALGRIND} ./multiseq_mn.x --sample 42 --seed 34824347 --short_header --width 0 ${AT1MB} | grep -c '^>' | diff --strip-trailing-cr - ${TMPFILE}
	@${VALGRIND} ./multiseq_mn.x --length_dist 1 ../testdata/at1MB.fna | diff --strip-trailing-cr -I '^#' - ../testdata/at1MB_length_dist_1.tsv
	@${VALGRIND} ./multiseq_mn.x --length_dist 10 ../testdata/at1MB.fna | diff --strip-trailing-cr -I '^#' - ../testdata/at1MB_length_dist_10.tsv
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed."

.PHONY:test_sorted_by_header_multiseq
test_sorted_by_header_multiseq:./multiseq_mn.x
	@./cmp_sorted_by_header.sh ../testdata/small.fna
	@./cmp_sorted_by_header.sh ../testdata/varlen_paired_1.fasta
	@./cmp_sorted_by_header.sh ${VAC}
	@./cmp_sorted_by_header.sh ${Y3}
	@./cmp_sorted_by_header.sh ${SW175}
	@./cmp_sorted_by_header.sh ${AT1MB}
	@echo "Congratulations. $@ passed."

.PHONY:test_line_generator
test_line_generator:line_generator.x
	@./line_generator.py `find .. -name '*.[ch]pp'` | sh -s
	@${VALGRIND} ./line_generator.x ${AT1MB} | diff --strip-trailing-cr - ${AT1MB}
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@cat ${AT1MB} | tr -d '\n' > ${TMPFILE}
	@echo "" >> ${TMPFILE}
	@${VALGRIND} ./line_generator.x ${TMPFILE} | diff --strip-trailing-cr - ${TMPFILE}
	@cat ${AT1MB} ${AT1MB} > ${TMPFILE}
	@${VALGRIND} ./line_generator.x ${AT1MB} ${AT1MB} | diff --strip-trailing-cr - ${TMPFILE}
	@cat ${AT1MB} ${AT1MB} ${AT1MB} > ${TMPFILE}
	@${VALGRIND} ./line_generator.x ${AT1MB} ${AT1MB} ${AT1MB} | diff --strip-trailing-cr - ${TMPFILE}
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

.PHONY:test_fastq_generator
test_fastq_generator:./fastq_mn.x
	@for num_threads in 1 2 3 4 5; do \
	  ${VALGRIND} ./fastq_mn.x --threads $${num_threads} ../testdata/70x_161nt_phred64.fastq | grep -E '# char|# total_count_entries' | diff --strip-trailing-cr - ../testdata/70x_161nt_phred64_dist.tsv || exit 1;\
	done
	@sh ./test_fastq_generator.sh
	@./check_split.sh ../testdata/70x_161nt_phred64.fastq
	@echo "Congratulations. $@ passed."

.PHONY:test_fasta_generator
test_fasta_generator:./fasta_mn.x
	@for input_mode in "" --mapped; do \
	  ${VALGRIND} ./fasta_mn.x $${input_mode} --width 70 ${AT1MB} | diff --strip-trailing-cr - ${AT1MB} || exit 1;\
	  ${VALGRIND} ./fasta_mn.x $${input_mode} --width 60 ${SW175} | diff --strip-trailing-cr - ${SW175} || exit 1;\
	  ${VALGRIND} ./fasta_mn.x $${input_mode} --width 61 ${Y3} | diff --strip-trailing-cr - ${Y3} || exit 1;\
	  ${VALGRIND} ./fasta_mn.x $${input_mode} --width 61 ${VAC} | diff --strip-trailing-cr - ${VAC} || exit 1;\
	done
	@echo "Congratulations. $@ passed."

.PHONY:test_thread_pool
test_thread_pool:thread_pool_mn.x
	@./thread_pool_mn.x 4 40
	@echo "Congratulations. $@ passed."

.PHONY:test_sort_kvt
test_sort_kvt:sort_key_value_pairs.x
	@for num in 66 666 66666 666666 6666666; do \
	   ./sort_key_value_pairs.x -d t -m lsb-radix $$num || exit 1;\
	   ./sort_key_value_pairs.x -d t -m mergesort $$num || exit 1;\
	   ./sort_key_value_pairs.x -d t -m stdsort $$num || exit 1;\
	done
	@echo "Congratulations. $@ passed."

.PHONY:test_sort_kvp
test_sort_kvp:sort_key_value_pairs.x
	@for num in 100 1000 10000 100000 1000000 10000000; do \
	   ./sort_key_value_pairs.x -d p -m lsb-radix $$num || exit 1;\
	   ./sort_key_value_pairs.x -d p -m mergesort $$num || exit 1;\
	   ./sort_key_value_pairs.x -d p -m stdsort $$num || exit 1;\
	done
	@echo "Congratulations. $@ passed."

.PHONY:test_sort_i
test_sort_i:sort_key_value_pairs.x
	@for num in 200 2000 20000 200000 2000000 20000000; do \
	   ./sort_key_value_pairs.x -t 1 -d i -m lsb-radix $$num || exit 1;\
	   ./sort_key_value_pairs.x -t 2 -d i -m lsb-radix $$num || exit 1;\
	   ./sort_key_value_pairs.x -t 3 -d i -m lsb-radix $$num || exit 1;\
	   ./sort_key_value_pairs.x -d i -m mergesort $$num || exit 1;\
	   ./sort_key_value_pairs.x -d i -m stdsort $$num || exit 1;\
	done
	@echo "Congratulations. $@ passed."

.PHONY:test_eoplist
test_eoplist:eoplist_mn.x
	@${VALGRIND} ./eoplist_mn.x silent 100000
	@${VALGRIND} ./eoplist_mn.x display 10
	@echo "Congratulations. $@ passed."

.PHONY:test_invint
test_invint:enum_invint.x
	@${VALGRIND} ./enum_invint.x 15 ../testdata/vaccg.fna | diff --strip-trailing-cr - ../testdata/vaccg_invint15.tsv
	@${VALGRIND} ./enum_invint.x 15 ../testdata/ychrIII.fna | diff --strip-trailing-cr - ../testdata/ychrIII_invint15.tsv
	@${VALGRIND} ./enum_invint.x 14 ../testdata/at1MB.fna | diff --strip-trailing-cr - ../testdata/at1MB_invint14.tsv
	@echo "Congratulations. $@ passed."

.PHONY:test_line_generator_mapped
test_line_generator_mapped:line_generator_mapped.x
	@${VALGRIND} ./line_generator_mapped.x `ls ../testdata/*.fna | grep -E -v 'empty.fna|\.gz'`
	@${VALGRIND} ./line_generator_mapped.x *.[hc]pp
	@echo "Congratulations. $@ passed."

.PHONY:test_minimizer_count
test_minimizer_count:minimizer_mn.x
#	@./check_err.py ./minimizer_mn.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@echo "# number of hashed kmers	37465" > ${TMPFILE}
	@${VALGRIND} ./minimizer_mn.x -w 30 -k 18 -s ${AT1MB} | grep 'number of hashed kmers' | diff --strip-trailing-cr - ${TMPFILE}
	@${VALGRIND} ./minimizer_mn.x -w 30 -k 18 -t 2 -s ${AT1MB} | grep 'number of hashed kmers' | diff --strip-trailing-cr - ${TMPFILE}
	@echo "# number of hashed kmers	37577" > ${TMPFILE}
	@${VALGRIND} ./minimizer_mn.x -w 30 -k 18 -c -s ${AT1MB} | grep 'number of hashed kmers' | diff --strip-trailing-cr - ${TMPFILE}
	@${VALGRIND} ./minimizer_mn.x -w 30 -k 18 -c -t 2 -s ${AT1MB} | grep 'number of hashed kmers' | diff --strip-trailing-cr - ${TMPFILE}
	@${VALGRIND} ./minimizer_mn.x -w 30 -k 18 --show_mode 1 ${AT1MB} > ${TMPFILE}
	@${VALGRIND} ./minimizer_mn.x -w 30 -k 18 --show_mode 2 ${AT1MB} | diff --strip-trailing-cr -I '^#' - ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed"

.PHONY:test_guess_if_protein_seq
test_guess_if_protein_seq:./guess_if_protein_seq.x
	./test_guess_if_protein_seq.sh
	@echo "Congratulations. $@ passed"

.PHONY:test_fs_prio_store
test_fs_prio_store:./fs_prio_store_mn.x
	@${VALGRIND} ./fs_prio_store_mn.x 1000 1000 1000 100
	@${VALGRIND} ./fs_prio_store_mn.x 1000 1000 1000 10000
	@${VALGRIND} ./fs_prio_store_mn.x 1000 2000 10000000 1000
	@echo "Congratulations. $@ passed"

.PHONY:test_rdbuf
test_rdbuf:./rdbuf_mn.x test_rdbuf_zlib
	@${VALGRIND} ./rdbuf_mn.x ${AT1MB} | diff --strip-trailing-cr - ${AT1MB}
	@${VALGRIND} ./rdbuf_mn.x ${SW175} | diff --strip-trailing-cr - ${SW175}
	@${VALGRIND} ./rdbuf_mn.x ${Y3} | diff --strip-trailing-cr - ${Y3}
	@${VALGRIND} ./rdbuf_mn.x ${VAC} | diff --strip-trailing-cr - ${VAC}

.PHONY:test_rdbuf_zlib
test_rdbuf_zlib:./rdbuf_mn.x
ifneq ($(without_zlib),yes)
	@${VALGRIND} ./rdbuf_mn.x ../testdata/SRR19536726_1_1000.fastq.gz | zdiff --strip-trailing-cr - ../testdata/SRR19536726_1_1000.fastq.gz
	@${VALGRIND} ./rdbuf_mn.x ../testdata/varlen_paired_both.fasta.gz | zdiff --strip-trailing-cr - ../testdata/varlen_paired_both.fasta.gz
	@echo "Congratulations. $@ passed"
endif

.PHONY:test_split_string
test_split_string:split_string_mn.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@echo "A" > ${TMPFILE}
	@${VALGRIND} ./split_string_mn.x '	' ${TMPFILE}
	@echo "A\tB" > ${TMPFILE}
	@${VALGRIND} ./split_string_mn.x '	' ${TMPFILE}
	@echo "A\tB\nC" > ${TMPFILE}
	@${VALGRIND} ./split_string_mn.x '	' ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed"

.PHONY:test_split_files
test_split_files:test_split_fastq_files test_split_fasta_files test_split_files_multithread
		@echo "Congratulations. $@ passed"

.PHONY:test_split_files_multithread
test_split_files_multithread:split_files_mn.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@${RM} ${TMPFILE}
	@${VALGRIND} ./split_files_mn.x -t 4 -p 5 -o ${TMPFILE} -c 0 ../testdata/varlen_paired_both.fasta.gz
	@cat ${TMPFILE}*.fasta | diff --strip-trailing-cr - <(gunzip -c ../testdata/varlen_paired_both.fasta.gz)
	@${RM} ${TMPFILE}*.fasta
	@${VALGRIND} ./split_files_mn.x -t 4 -p 5 -o ${TMPFILE} -c 7 ../testdata/varlen_paired_both.fasta.gz
	@gunzip -c ${TMPFILE}*.fasta.gz | diff --strip-trailing-cr - <(gunzip -c ../testdata/varlen_paired_both.fasta.gz)
	@${RM} ${TMPFILE}*.fasta.gz
	@${VALGRIND} ./split_files_mn.x -t 8 -l 20000 -o ${TMPFILE} ../testdata/SRR19536726_1_1000.fastq.gz
	@gunzip -c ${TMPFILE}*.fastq.gz | diff --strip-trailing-cr -I "^\+" - <(gunzip -c ../testdata/SRR19536726_1_1000.fastq.gz)
	@${RM} ${TMPFILE}*.fastq.gz
	@echo "Congratulations. $@ passed"

# We check uncompressed, custom higher compression level and all 3 types of splitting
.PHONY:test_split_fasta_files
test_split_fasta_files:split_files_mn.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@${RM} ${TMPFILE}
	@${VALGRIND} ./split_files_mn.x -p 5 -o ${TMPFILE} -c 0 ../testdata/varlen_paired_both.fasta.gz
	@cat ${TMPFILE}*.fasta | diff --strip-trailing-cr - <(gunzip -c ../testdata/varlen_paired_both.fasta.gz)
	@${RM} ${TMPFILE}*.fasta
	@${VALGRIND} ./split_files_mn.x -p 5 -o ${TMPFILE} -c 7 ../testdata/varlen_paired_both.fasta.gz
	@gunzip -c ${TMPFILE}*.fasta.gz | diff --strip-trailing-cr - <(gunzip -c ../testdata/varlen_paired_both.fasta.gz)
	@${RM} ${TMPFILE}*.fasta.gz
	@${VALGRIND} ./split_files_mn.x -l 120 -o ${TMPFILE} ../testdata/varlen_paired_both.fasta.gz
	@gunzip -c ${TMPFILE}*.fasta.gz | diff --strip-trailing-cr - <(gunzip -c ../testdata/varlen_paired_both.fasta.gz)
	@${RM} ${TMPFILE}*.fasta.gz
	@${VALGRIND} ./split_files_mn.x -n 2 -o ${TMPFILE} ../testdata/varlen_paired_both.fasta.gz
	@gunzip -c ${TMPFILE}*.fasta.gz | diff --strip-trailing-cr - <(gunzip -c ../testdata/varlen_paired_both.fasta.gz)
	@${RM} ${TMPFILE}*.fasta.gz
	@echo "Congratulations. $@ passed"


.PHONY:test_split_fastq_files
test_split_fastq_files:split_files_mn.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@${RM} ${TMPFILE}
	@${VALGRIND} ./split_files_mn.x -p 5 -o ${TMPFILE} -c 0 ../testdata/SRR19536726_1_1000.fastq.gz
	@cat ${TMPFILE}*.fastq | diff --strip-trailing-cr -I "^\+" - <(gunzip -c ../testdata/SRR19536726_1_1000.fastq.gz)
	@${RM} ${TMPFILE}*.fastq
	@${VALGRIND} ./split_files_mn.x -p 5 -o ${TMPFILE} -c 7 ../testdata/SRR19536726_1_1000.fastq.gz
	@gunzip -c ${TMPFILE}*.fastq.gz | diff --strip-trailing-cr -I "^\+" - <(gunzip -c ../testdata/SRR19536726_1_1000.fastq.gz)
	@${RM} ${TMPFILE}*.fastq.gz
	@${VALGRIND} ./split_files_mn.x -l 20000 -o ${TMPFILE} ../testdata/SRR19536726_1_1000.fastq.gz
	@gunzip -c ${TMPFILE}*.fastq.gz | diff --strip-trailing-cr -I "^\+" - <(gunzip -c ../testdata/SRR19536726_1_1000.fastq.gz)
	@${RM} ${TMPFILE}*.fastq.gz
	@${VALGRIND} ./split_files_mn.x -n 200 -o ${TMPFILE} ../testdata/SRR19536726_1_1000.fastq.gz
	@gunzip -c ${TMPFILE}*.fastq.gz | diff --strip-trailing-cr -I "^\+" - <(gunzip -c ../testdata/SRR19536726_1_1000.fastq.gz)
	@${RM} ${TMPFILE}*.fastq.gz
	@echo "Congratulations. $@ passed"

.PHONY:test_multiseq_factory
test_multiseq_factory:multiseq_factory_mn.x fastq_mn.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@echo "# number of parts	20" > ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x -p 20 ../testdata/at1MB.fna | diff --strip-trailing-cr - ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x --width 70 -p 20 ../testdata/at1MB.fna | diff -I '^#' --strip-trailing-cr - ../testdata/at1MB.fna
	@${VALGRIND} ./multiseq_factory_mn.x -l 39000 ../testdata/at1MB.fna | diff --strip-trailing-cr - ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x -n 100 ../testdata/at1MB.fna | diff --strip-trailing-cr - ${TMPFILE}
	@echo "# number of parts	9" > ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x -p 9 ../testdata/sw175.fna | diff - --strip-trailing-cr ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x --width 60 -p 9 ../testdata/sw175.fna | diff -I '^#' - --strip-trailing-cr ../testdata/sw175.fna
	@${VALGRIND} ./multiseq_factory_mn.x -l 6200 ../testdata/sw175.fna | diff - --strip-trailing-cr ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x -n 20 ../testdata/sw175.fna | diff - --strip-trailing-cr ${TMPFILE}
	@echo "# number of parts	4" > ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x -n 2 ../testdata/varlen_paired_1.fastq ../testdata/varlen_paired_2.fastq | diff --strip-trailing-cr - ${TMPFILE}
	@echo "# number of parts	10" > ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x -n 200 ../testdata/SRR19536726_1_1000.fastq.gz ../testdata/SRR19536726_1_1000.fastq.gz | diff --strip-trailing-cr - ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x --width 70 -n 2 ../testdata/varlen_paired_1.fastq ../testdata/varlen_paired_2.fastq > ${TMPFILE}
	@${VALGRIND} ./fastq_mn.x --width 70 --paired --fasta_output ../testdata/varlen_paired_1.fastq ../testdata/varlen_paired_2.fastq | diff -I '^#' --strip-trailing-cr - ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x --width 70 -n 100 ../testdata/SRR19536726_1_1000.fastq.gz ../testdata/SRR19536726_1_1000.fastq.gz > ${TMPFILE}
	@${VALGRIND} ./fastq_mn.x --width 70 --paired --fasta_output ../testdata/SRR19536726_1_1000.fastq.gz ../testdata/SRR19536726_1_1000.fastq.gz | diff -I '^#' --strip-trailing-cr - ${TMPFILE}
	@echo "# number of parts	25" > ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x -l 8080 ../testdata/SRR19536726_1_1000.fastq.gz ../testdata/SRR19536726_1_1000.fastq.gz | diff --strip-trailing-cr - ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x -n 80 ../testdata/SRR19536726_1_1000.fastq.gz ../testdata/SRR19536726_1_1000.fastq.gz | diff --strip-trailing-cr - ${TMPFILE}
	@${VALGRIND} ./multiseq_factory_mn.x -p 25 ../testdata/SRR19536726_1_1000.fastq.gz ../testdata/SRR19536726_1_1000.fastq.gz | diff --strip-trailing-cr - ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "Congratulations. $@ passed"

.PHONY:test_untar
test_untar:untar_zipped_mn.x
	@for suffix in gz bz2; do \
	  for norapidgzip in "" "--no_rapidgzip"; do \
	    for store_in_memory  in "" "--store_in_memory"; do \
	      ${VALGRIND} ./untar_zipped_mn.x $${norapidgzip} $${store_in_memory} ../testdata/pubmed_small.tar.$${suffix} | diff --strip-trailing-cr - ../testdata/pubmed_small.txt || exit 1;\
	    done \
	  done \
	done
	@echo "$@ passed"

.PHONY:test_read_vector
test_read_vector:read_vector_mn.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@for filename in ${AT1MB} ${SW175} ${Y3} ${VAC}; do \
	  gzip -9 -c $${filename} > ${TMPFILE}.gz || exit 1; \
	  ${VALGRIND} ./read_vector_mn.x $${filename} ${TMPFILE}.gz || exit 1; \
	done
	@${RM} ${TMPFILE} ${TMPFILE}.gz
	@echo "$@ passed"

.PHONY:test_string_of_digits
test_string_of_digits:string_of_digits.x
	@${VALGRIND} ./string_of_digits.x
	@echo "Congratulations. $@ passed"

.PHONY:test_binary_iterator
test_binary_iterator:binary_iterator.x
	@${VALGRIND} ./binary_iterator.x ${AT1MB} | diff --strip-trailing-cr - ${AT1MB}
	@${VALGRIND} ./binary_iterator.x ${SW175} | diff --strip-trailing-cr - ${SW175}
	@while IFS= read -r file; do \
       file=$$(echo "$$file" | tr -d '\r'); \
       ${VALGRIND} ./binary_iterator.x "../testdata/$$file" | diff --strip-trailing-cr - "../testdata/$$file"; \
     done < ../testdata/corrupt_files.txt
	@echo "Congratulations. $@ passed"

.PHONY:test_multiple_options
test_multiple_options:multiple_options_mn.x
	@./multiple_options_mn.x 'p_value+k_mer' | sed -E '/^[A-Z]:\\\\/ s/^([^:]*:){2}/.\/multiple_options_mn.x:/' | diff --strip-trailing-cr - ../testdata/p_value_and_k_mer.txt
	@echo "$@ passed"

.PHONY:test_multibitvector
test_multibitvector:multibitvector.x
	@./multibitvector.x
	@echo "Congratulations. $@ passed"

fastq_mn.x:fastq_mn.o seq_reader_options.o
	$(CXX) ${LDFLAGS} fastq_mn.o seq_reader_options.o -o $@ ${LDLIBS}

fasta_mn.x:fasta_mn.o seq_reader_options.o
	$(CXX) ${LDFLAGS} fasta_mn.o seq_reader_options.o -o $@ ${LDLIBS}

minimizer_mn.x:minimizer_opt.o minimizer_mn.o
	$(CXX) ${LDFLAGS} minimizer_opt.o minimizer_mn.o -o $@ ${LDLIBS}

untar_zipped_mn.x:untar_zipped_mn.o untar_zipped_op.o
	$(CXX) ${LDFLAGS} untar_zipped_mn.o untar_zipped_op.o -o $@ ${LDLIBS}
