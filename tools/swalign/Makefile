GTTL=../..
include ${GTTL}/config.mk

SWAGAOBJ=sw_all_against_all.o

SSW_SIMD_DEBUG?=0
ifneq ($(SSW_SIMD_DEBUG),0)
  CPPFLAGS += -DSSW_SIMD_DEBUG=${SSW_SIMD_DEBUG}
endif

ifeq (${MACHINE},x86_64)
ifneq (${avx},no)
  avxflag=${shell echo "int main(){ return 0; }" | $(CXX) -x c - -DAVX512 -o /dev/null >/dev/null 2>&1 && echo yes || echo no}
  ifeq (${avxflag},yes)
    # we assume avx512 implies avx2 availability
    CFLAGS += -mavx2
    # the following does not work yet
    # CFLAGS += -DAVX512
  else
    avxflag=${shell echo "int main(){ return 0; }" | $(CXX) -x c - -mavx2 -o /dev/null >/dev/null 2>&1 && echo yes || echo no}
    ifeq (${avxflag},yes)
      CFLAGS += -mavx2
    endif
  endif
endif
endif

ifeq ($(SAVETEMPS),yes)
  CFLAGS+=-save-temps
endif

all:sw_all_against_all.x

sw_all_against_all.x:${SWAGAOBJ}
	$(LD) ${LDFLAGS} ${SWAGAOBJ} -o $@ ${LDLIBS}

../../src/alignment/sw_simd_uint8.hpp:ssw_gen.py
	./ssw_gen.py 8 | clang-format --assume-filename=$@ | \
		fold -s -w 80 | \
		sed 's/[ \t]*$$//' > $@

../../src/alignment/sw_simd_uint16.hpp:ssw_gen.py
	./ssw_gen.py 16 | clang-format --assume-filename=$@ | \
		fold -s -w 80 | \
		sed 's/[ \t]*$$//' > $@

../../src/alignment/sw_simd_uint32.hpp:ssw_gen.py
	./ssw_gen.py 32 | clang-format --assume-filename=$@ | \
		fold -s -w 80 | \
		sed 's/[ \t]*$$//' > $@

../../src/alignment/blosum62.hpp:blosum62.txt scorematrix.py
	@./scorematrix.py --constexpr blosum62 $< | clang-format --assume-filename=$@ > $@

../../src/alignment/unit_score_aa.hpp:unit_score_aa.txt scorematrix.py
	@./scorematrix.py --constexpr unit_score_aa $< | clang-format --assume-filename=$@ > $@

../../src/alignment/unit_score_nuc_lower.hpp:unit_score_nuc_lower.txt scorematrix.py
	@./scorematrix.py --constexpr unit_score_nuc_lower $< | clang-format --assume-filename=$@ > $@

../../src/alignment/unit_score_nuc_upper.hpp:unit_score_nuc_upper.txt scorematrix.py
	@./scorematrix.py --constexpr unit_score_nuc_upper $< | clang-format --assume-filename=$@ > $@

../../src/alignment/score_matrix_name.hpp:../../scripts/enum_with_string.py
	@../../scripts/enum_with_string.py --prefix Score_matrix_ \
        ScoreMatrixName \
        undefined \
        blosum62 \
        unit_score_aa \
        unit_score_nuc \
        unit_score_nuc_2_2 \
        unit_score_nuc_lower \
        unit_score_nuc_upper | \
        clang-format --assume-filename=$@ > $@

../../src/alignment/scoring_info_and_seq_trans.hpp:score_class_choices.py
	@./score_class_choices.py --scoring | clang-format -style="{ColumnLimit: 80}" --assume-filename=$@ > $@

alignment_output_function.hpp:score_class_choices.py
	@./score_class_choices.py --alignment_output | clang-format --style="{ColumnLimit: 80}" -assume-filename=$@ > $@

sequence_decode_function.hpp:score_class_choices.py
	@./score_class_choices.py --sequence_decode | clang-format --assume-filename=$@ > $@

.PHONY:score_matrix_list_update
score_matrix_list_update:../../src/alignment/scoring_info_and_seq_trans.hpp \
	../../src/alignment/score_matrix_name.hpp \
        ../../src/alignment/blosum62.hpp \
	../../src/alignment/unit_score_aa.hpp \
	../../src/alignment/unit_score_nuc_lower.hpp \
	../../src/alignment/unit_score_nuc_upper.hpp \
	../../src/utilities/runtime_class.hpp
	make alignment_output_function.hpp
	make sequence_decode_function.hpp

.PHONY:test-agc
test-agc: test-agc-proteins test-agc-DNA test-best-threads test-against-0 

.PHONY:test-agc-proteins
test-agc-proteins:sw_all_against_all.x
	@for prog in ./sw_all_against_all.x; do \
	  ${VALGRIND} $${prog} -d testdata/ref-protein.fsa -q testdata/query-protein.fasta -s blosum62 -g 11 1 -v 1 | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/ref-query-protein-v1.tsv || exit 1; \
	  ${VALGRIND} $${prog} -d testdata/ref-protein.fsa -q testdata/query-protein.fasta -s blosum62 -g 11 1 -v 2 -t 2 -a 1 | ./largerscore.py - testdata/ref-query-protein.tsv || exit 1; \
	  for moption in "" "-m"; do \
	    ${VALGRIND} $${prog} -d testdata/Q6GZS4.fsa -q testdata/query-protein.fasta -s blosum62 -g 11 1 -v 2 -a 5 $${moption} | diff --strip-trailing-cr -I '^#' - testdata/Q6GZS4-query-cigar.tsv || exit 1; \
	    ${VALGRIND} $${prog} -d testdata/Q6GZS4.fsa -q testdata/query-protein.fasta -s blosum62 -g 11 1 -v 2 -a 60 $${moption} | diff --strip-trailing-cr -I '^#' - testdata/Q6GZS4-query-align.txt || exit 1; \
	  done; \
	done
	@echo "tests for protein alignments with affine gap costs (vectorized) passed"

.PHONY:test-agc-DNA
test-agc-DNA:sw_all_against_all.x
	@$(eval TMPFILE1 := $(shell mktemp --tmpdir=.))
	@$(eval TMPFILE2 := $(shell mktemp --tmpdir=.))
	@for prog in ./sw_all_against_all.x; do \
	  ${VALGRIND} $${prog} -v 1 -d testdata/simple_fwd.fna | cut -f 1-5 | diff --strip-trailing-cr -I '^#' - testdata/simple_v1_output.tsv || exit 1; \
	  ${VALGRIND} $${prog} -d testdata/smrt100-2000.fna -s unit_score_nuc -g 4 1 -v 1 | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/smrt100-2000.v1.tsv || exit 1; \
	  ${VALGRIND} $${prog} -d testdata/smrt100-2000.fna -s unit_score_nuc_lower -g 4 1 -v 2 | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/smrt100-2000.v2.tsv || exit 1; \
	  ${VALGRIND} $${prog} -d testdata/smrt100-2000.fna -s unit_score_nuc_lower -g 4 1 -v 2 -a 1 | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/smrt100-2000.v2.tsv || exit 1; \
	  ${VALGRIND} $${prog} -d testdata/smrt100-2000.fna -s unit_score_nuc_lower -g 4 1 -v 2 -t 2 | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/smrt100-2000.v2.tsv || exit 1; \
	  ${VALGRIND} $${prog} -d testdata/smrt100-2000.fna -s unit_score_nuc_lower -g 4 1 -v 2 -a 5 | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/smrt100-2000.a2.tsv || exit 1; \
	  ${VALGRIND} $${prog} -d testdata/smrt100-2000.fna -s unit_score_nuc_lower -g 4 1 -v 2 -a 60 > ${TMPFILE1} || exit 1; \
	  ./alignment_parse.py ${TMPFILE1} || exit 1; \
	  bzip2 -d -c testdata/smrt100-2000.a60_new.tsv.bz2 > ${TMPFILE2} || exit 1; \
	  ./cmp_with_alignment.py -s ${TMPFILE1} ${TMPFILE2} || exit 1; \
	  ${VALGRIND} $${prog} -d testdata/smrt100-2000-first.fna -q testdata/smrt100-2000-last.fna -s unit_score_nuc_lower -g 4 1 -v 2 | grep -v '^#' | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/smrt100-2000-d_vs_q.tsv || exit 1; \
	  $${prog} -d testdata/humhdabcd.fna -q testdata/humdystrop.fna -g 4 1 -v 2 -a 5 | diff --strip-trailing-cr -I '^#' - testdata/humhdabcd_humdystrop.tsv || exit 1; \
	done
	@${RM} ${TMPFILE1} ${TMPFILE2}
	@echo "tests for DNA alignments with affine gap costs (vectorized) passed"

.PHONY:test-best-threads
test-best-threads:sw_all_against_all.x
	@$(eval TMPFILE1 := $(shell mktemp --tmpdir=.))
	@./select-best.py -a 11 testdata/ref-query-protein.tsv | sort -k 1 -n -k 7 -n > ${TMPFILE1}
	@$(eval TMPFILE2 := $(shell mktemp --tmpdir=.))
	@./select-best.py -a 2 testdata/smrt100-2000.v2.best2.tsv | sort -k 1 -n -k 7 -n > ${TMPFILE2}
	@for num_threads in 1 2 4; do \
	  ${VALGRIND} ./sw_all_against_all.x -d testdata/ref-protein.fsa -q testdata/query-protein.fasta -s blosum62 -g 11 1 -v 2 -b 11 -a 1 -t $${num_threads} | sort -k 1 -n -k 7 -n  | diff --strip-trailing-cr -I '^#' - ${TMPFILE1} || exit 1; \
	  ${VALGRIND} ./sw_all_against_all.x -d testdata/smrt100-2000.fna -s unit_score_nuc_lower -g 4 1 -v 2 -b 2 -t $${num_threads} | sort -k 1 -n -k 7 -n | diff --strip-trailing-cr -I '^#' - ${TMPFILE2} || exit 1; \
	done
	@${RM} ${TMPFILE1} ${TMPFILE2}
	@echo "$@ passed"

.PHONY:test-against-0
test-against-0:sw_all_against_all.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@./alignment_parse.py testdata/smrt100-2000-a60-v0.txt
	@for prog in ./sw_all_against_all.x; do \
	  ${VALGRIND} $${prog} -a 60 -v 2 -d testdata/smrt100-2000.fna > ${TMPFILE} || exit 1; \
	  ./alignment_parse.py ${TMPFILE} || exit 1; \
	  ./cmp_with_alignment.py -s testdata/smrt100-2000-a60-v0.txt ${TMPFILE} || exit 1; \
	done
	@${RM} ${TMPFILE}
	@echo "$@ passed"

# tests for linear gaps costs

.PHONY:test-lgc
test-lgc:sw_all_against_all.x testdata/gsa-seqpair.fasta testdata/gsa-align.tsv testdata/seq18.fasta testdata/local-coords18.tsv testdata/seq50.fasta testdata/seq98.fasta
	@$(eval TMPFILE1 := $(shell mktemp --tmpdir=.))
	@paste testdata/local-coords18.tsv testdata/local-cigar18.txt > ${TMPFILE1}
	@$(eval TMPFILE2 := $(shell mktemp --tmpdir=.))
	@gzip -d -c testdata/local-cigar98.tsv.gz > ${TMPFILE2}
	@for prog in ./sw_all_against_all.x; do \
	  ${VALGRIND} $${prog} -g 0 1 -v 2 -s unit_score_aa -d testdata/gsa-seqpair.fasta | diff --strip-trailing-cr -I '^#' - testdata/gsa-align.tsv || exit 1; \
	  ${VALGRIND} $${prog} -g 0 4 -a 1 -v 2 -d testdata/seq18.fasta | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/local-coords18.tsv || exit 1; \
	  ${VALGRIND} $${prog} -g 0 4 -a 1 -v 2 -d testdata/seq18.fasta -t 2 | sort -k1,1n -k2,2n | diff --strip-trailing-cr -I '^#' - testdata/local-coords18.tsv || exit 1; \
	  ${VALGRIND} $${prog} -g 0 4 -a 5 -v 2 -d testdata/seq18.fasta | sort -k1,1n -k2,2n | diff --strip-trailing-cr -I '^#' - ${TMPFILE1} || exit 1; \
	  ${VALGRIND} $${prog} -g 0 4 -a 1 -v 2 -d testdata/seq98.fasta | cut -f 1,2,7 | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/seq98-a1-v0.tsv || exit 1; \
	  ${VALGRIND} $${prog} -g 0 4 -a 5 -v 2 -d testdata/seq98.fasta | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - ${TMPFILE2} || exit 1; \
	done
	@${RM} ${TMPFILE1} ${TMPFILE2}
	@echo "tests for linear gap costs (vectorized) passed"

.PHONY:test_seqoutput
test_seqoutput:sw_all_against_all.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@paste testdata/local-coords18.tsv testdata/local-seqs18.txt > ${TMPFILE}
	@${VALGRIND} ./sw_all_against_all.x -g 0 4 -a 6+7 -d testdata/seq18.fasta |  grep -v '^#' | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -w -I '^#' - ${TMPFILE}
	@cut -f 1 testdata/local-seqs18.txt | paste testdata/local-coords18.tsv - > ${TMPFILE}
	@${VALGRIND} ./sw_all_against_all.x -g 0 4 -a 6 -d testdata/seq18.fasta |  grep -v '^#' | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -w -I '^#' - ${TMPFILE}
	@cut -f 2 testdata/local-seqs18.txt | paste testdata/local-coords18.tsv - > ${TMPFILE}
	@${VALGRIND} ./sw_all_against_all.x -g 0 4 -a 7 -d testdata/seq18.fasta |  grep -v '^#' | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -w -I '^#' - ${TMPFILE}
	@paste testdata/local-coords18.tsv testdata/local-cigar18.txt testdata/local-seqs18.txt > ${TMPFILE}
	@${VALGRIND} ./sw_all_against_all.x -g 0 4 -a 5+6+7 -d testdata/seq18.fasta |  grep -v '^#' | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -w -I '^#' - ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "$@ passed"

.PHONY:test_identity_output
test_identity_output:sw_all_against_all.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@sort -k 1 -n -k 2 -n testdata/local-coords18.tsv | paste - testdata/local-identity18.tsv | cut -f 1-7,10 > ${TMPFILE}
	@${VALGRIND} ./sw_all_against_all.x -g 0 4 -a 4 -d testdata/seq18.fasta | grep -v '^#' | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -w -I '^#' - ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "$@ passed"

.PHONY:test-swaga
test-swaga:test-swaga-proteins-self test-swaga-proteins-query test-swaga-proteins-header test-swaga-DNA test_overflow test_seqoutput test_identity_output
	@echo "$@ passed"

.PHONY:test-swaga-proteins-self
test-swaga-proteins-self:sw_all_against_all.x
	@for fileid in seq18 seq50 query-protein; do \
	  ./sw_all_against_all.x -d testdata/$${fileid}.fasta | sort -k1 -n -k2 -n  | diff --strip-trailing-cr -I '^#' - testdata/$${fileid}-d.tsv || exit 1; \
	  ./sw_all_against_all.x -t 2 -d testdata/$${fileid}.fasta | sort -k1 -n -k2 -n  | diff --strip-trailing-cr -I '^#' - testdata/$${fileid}-d.tsv || exit 1; \
	  ./sw_all_against_all.x -t 4 -d testdata/$${fileid}.fasta | sort -k1 -n -k2 -n  | diff --strip-trailing-cr -I '^#' - testdata/$${fileid}-d.tsv || exit; \
	done
	@echo "$@ passed"

.PHONY:test-swaga-proteins-query
test-swaga-proteins-query:sw_all_against_all.x
	@$(eval TMPFILE1 := $(shell mktemp --tmpdir=.))
	@for num_threads in 1 2 4; do \
	  ./sw_all_against_all.x -t $${num_threads} -d testdata/ref-protein.fsa -q testdata/query-protein.fasta | grep -v '^#' | sort -k1 -n -k2 -n | diff --strip-trailing-cr -I '^#' - testdata/ref-query-protein.tsv || exit 1; \
	done
	@./select-best.py -a 29 testdata/ref-query-protein.tsv | sort -k1 -n -k2 -n > ${TMPFILE1}
	@for num_threads in 1 2 4; do \
	  ./sw_all_against_all.x -c 50 -t $${num_threads} -d testdata/ref-protein.fsa -q testdata/query-protein.fasta | sort -k1 -n -k2 -n | diff --strip-trailing-cr -I '^#' - ${TMPFILE1} || exit 1; \
	done
	@${RM} ${TMPFILE1}
	@echo "$@ passed"

.PHONY:test-swaga-proteins-header
test-swaga-proteins-header:sw_all_against_all.x
	@./sw_all_against_all.x -h -d testdata/ref-protein.fsa -q testdata/query-protein.fasta | sort -k1,2 | diff --strip-trailing-cr -I '^#' - testdata/ref-query-protein-d-h.tsv
	@./sw_all_against_all.x -h -d testdata/seq18.fasta | sed -e 's/seq_//g' | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/seq18-d.tsv
	@echo "$@ passed"

.PHONY:test-swaga-DNA
test-swaga-DNA:sw_all_against_all.x
	@$(eval TMPFILE1 := $(shell mktemp --tmpdir=.))
	@$(eval TMPFILE2 := $(shell mktemp --tmpdir=.))
	@../../scripts/select-seq.py -f 100 ${AT1MB} > ${TMPFILE1}
	@./sw_all_against_all.x -d ${TMPFILE1} | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/at100-d.tsv
	@${RM} ${TMPFILE1} ${TMPFILE2}
	@echo "$@ passed"

.PHONY:test-swaga-b
test-swaga-b:sw_all_against_all.x
	@$(eval TMPFILE1 := $(shell mktemp --tmpdir=.))
	@$(eval TMPFILE2 := $(shell mktemp --tmpdir=.))
	@$(eval TMPFILE3 := $(shell mktemp --tmpdir=.))
	@./sw_all_against_all.x -d testdata/seq50.fasta | ./select-best.py -a 50 - > ${TMPFILE1}
	@../../scripts/select-seq.py -f 100 ${AT1MB} > ${TMPFILE2}
	@./sw_all_against_all.x -d ${TMPFILE2} | ./select-best.py -a 37 - > ${TMPFILE3}
	@for num_threads in 1 2 3 4; do \
	  ./sw_all_against_all.x -t $${num_threads} -b 50 -d testdata/seq50.fasta | diff --strip-trailing-cr -I '^#' - ${TMPFILE1} || exit 1; \
	  ./sw_all_against_all.x -t $${num_threads} -b 37 -d ${TMPFILE2} | diff --strip-trailing-cr -I '^#' - ${TMPFILE3} || exit 1; \
	done
	@${RM} ${TMPFILE1} ${TMPFILE2} ${TMPFILE3}
	@echo "$@ passed"

.PHONY:test-conv
test-conv:sw_all_against_all.x testdata/gsa-seqpair.fasta testdata/gsa-align.tsv testdata/seq18.fasta testdata/local-coords18.tsv testdata/seq50.fasta
	@${VALGRIND} ./sw_all_against_all.x -g 0 1 -v 2 -s unit_score_aa -d testdata/gsa-seqpair.fasta | diff --strip-trailing-cr -I '^#' - testdata/gsa-align.tsv
	@$(eval TMPFILE1 := $(shell mktemp --tmpdir=.))
	@${VALGRIND} ./sw_all_against_all.x -g 0 4 -a 1 -v 2 -d testdata/seq18.fasta | sort -k 1 -n -k 2 -n | cut -f 1,2,7 > ${TMPFILE1}
	@cut -f 1,2,7 testdata/local-coords18.tsv | diff --strip-trailing-cr -I '^#' - ${TMPFILE1}
	@${VALGRIND} ./sw_all_against_all.x -g 11 1 -v 2 -d testdata/seq50.fasta | sort -k 1 -n -k 2 -n | diff --strip-trailing-cr -I '^#' - testdata/seq50-d.tsv
	@${RM} ${TMPFILE1}
	@echo "tests for linear and affine gap costs (not vectorized) passed"

.PHONY:test
test:test_err \
     test-lgc \
     test-agc \
     test-conv \
     test-swaga \
     test-swaga-b \
     test-restrict \
     test_mmseq_FP \
     test-alignment-parse \
     test-calls
	@echo "$@ passed"

# delete sequence 13, 21, 79, 80, 81 from seq98

.PHONY:test_err
test_err:sw_all_against_all.x
	./check_err.py --path ../.. './sw_all_against_all.x -d'

.PHONY:test-calls
txest-calls:
	./gen_calls.py | sh -s

.PHONY:test-alignment-parse
test-alignment-parse:sw_all_against_all.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@../../scripts/select-seq.py -f 250 ${AT1MB} > ${TMPFILE}
	@./sw_all_against_all.x -d ${TMPFILE} -a 60 | ./alignment_parse.py -
	@${RM} ${TMPFILE}
	@echo "$@ passed"

.PHONY:test_overflow
test_overflow:sw_all_against_all.x simd32_seqs
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@touch ${TMPFILE}
	@for filename in `ls simd32_seqs`; do \
	   echo $${filename}; \
           ./sw_all_against_all.x -v 1 -d simd32_seqs/$$filename >> ${TMPFILE} || exit 1; \
           ./sw_all_against_all.x -v 2 -d simd32_seqs/$$filename >> ${TMPFILE} || exit 1; \
        done
	@grep '^# TIME' ${TMPFILE}
	@diff --strip-trailing-cr -I '^# ' ${TMPFILE} testdata/simd32_seqs_results.tsv
	@${RM} ${TMPFILE}
	@echo "$@ passed"

.PHONY:test-restrict
test-restrict:sw_all_against_all.x testdata/seq18.fasta testdata/seq50.fasta testdata/seq50_c1000.tsv testdata/seq50_18_4pairs.tsv testdata/seq50_c1000.tsv  testdata/seq50_18_4pairs.tsv
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@./sw_all_against_all.x -d testdata/seq50.fasta -h -c 1000 > ${TMPFILE}
	@./sw_all_against_all.x -d testdata/seq50.fasta -h -r testdata/seq50_c1000.tsv | diff --strip-trailing-cr -I '^#' - ${TMPFILE}
	@./sw_all_against_all.x -d testdata/seq50.fasta -h -q testdata/seq18.fasta -c 22 > ${TMPFILE}
	@./sw_all_against_all.x -d testdata/seq50.fasta -h -q testdata/seq18.fasta -r testdata/seq50_18_4pairs.tsv | diff --strip-trailing-cr -I '^#' - ${TMPFILE}
	@./sw_all_against_all.x -d testdata/seq50.fasta -c 1000 > ${TMPFILE}
	@$(eval TMPFILE2 := $(shell mktemp --tmpdir=.))
	@sed -e 's/seq_//g' testdata/seq50_c1000.tsv > ${TMPFILE2}
	@./sw_all_against_all.x -d testdata/seq50.fasta -r ${TMPFILE2}  | diff --strip-trailing-cr -I '^#' - ${TMPFILE}
	@${RM} ${TMPFILE} ${TMPFILE2}
	@./sw_all_against_all.x -d testdata/seq50.fasta -q testdata/seq18.fasta -c 22 > ${TMPFILE}
	@sed -e 's/seq_//g' testdata/seq50_18_4pairs.tsv > ${TMPFILE2}
	@./sw_all_against_all.x -d testdata/seq50.fasta -q testdata/seq18.fasta -r ${TMPFILE2} | diff --strip-trailing-cr -I '^#' - ${TMPFILE}
	@${RM} ${TMPFILE} ${TMPFILE2}
	@echo "$@ passed"

.PHONY:test_mmseq_FP
test_mmseq_FP:sw_all_against_all.x
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@./sw_all_against_all.x -h -a 60 -d testdata/A0A059IV01.fna  -q testdata/C0QTH6.fna > ${TMPFILE}
	@./alignment_parse.py ${TMPFILE}
	@diff --strip-trailing-cr -I '^#' ${TMPFILE} testdata/A0A059IV01_C0QTH6.txt
	@${RM} ${TMPFILE}
	@echo "$@ passed"

simd32_seqs:
	@tar -xvjf testdata/simd32_seqs.tar.bz2

.PHONY:include_list
include_list:
	 @gcc -I ../../src -M sw_all_against_all.cpp | grep -v gttl | grep -v 'usr\/include' | grep -v 'usr\/lib'
