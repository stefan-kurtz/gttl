CXX?=g++
LD=${CXX}
CFLAGS=-g -Wall -Werror -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -funroll-loops -O3
GTTL=../..
CPPFLAGS=-I ${GTTL}/src
CXXFLAGS=-std=c++17
LDFLAGS=-g
LDLIBS=-lm -lz -lpthread -lstdc++
ECO?=${ECO29_SINGLE_FILE}

ifneq ($(debug),yes)
  CFLAGS+=-DNDEBUG
endif

OBJ=ntcard_mn.o ntcard_opt.o

all:ntcard_mn.x

ntcard_mn.x:${OBJ}
	$(LD) ${LDFLAGS} ${OBJ} -o $@ ${LDLIBS}

%.o:%.cpp
	$(CXX) ${TIME_OPTION} -o $@ -c $< ${CXXFLAGS} ${CFLAGS} ${CPPFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d)

.PHONY:test
test:test_random test_fastq 
	@echo "$@ passed"

.PHONY:test_large
test_large:test_hiq test_eco test_soil_unique test_hg38chr1
	@echo "$@ passed"

.PHONY:test_hiq
test_hiq:ntcard_mn.x
	@./ntcard_mn.x ${HIQ} | diff --strip-trailing-cr -I "^#" - references/hiq_base.txt
	@for mode in binary fast long ; do \
	   for num_threads in 1 2 3; do \
	     ./ntcard_mn.x --$${mode} --threads $${num_threads} ${HIQ} | diff --strip-trailing-cr -I '^#' - references/hiq_$${mode}.txt; \
	   done \
	done
	@echo "$@ passed"

eco29.fasta:${ECO}
	@gzip -d -c $< > $@

.PHONY:test_eco
test_eco:ntcard_mn.x eco29.fasta
	@./ntcard_mn.x eco29.fasta | diff --strip-trailing-cr -I "^#" - references/eco_base.txt
	@for mode in binary fast long; do \
	   for num_threads in 1 2 3; do \
	     ./ntcard_mn.x --$${mode} --threads $${num_threads} eco29.fasta | diff --strip-trailing-cr -I '^#' - references/eco_$${mode}.txt; \
	   done \
	done
	@echo "$@ passed"

soil_unique.fasta:${SOIL_UNIQUE}
	@gzip -d -c $< > $@

.PHONY:test_soil_unique
test_soil_unique:ntcard_mn.x soil_unique.fasta
	@./ntcard_mn.x soil_unique.fasta | diff --strip-trailing-cr -I "^#" - references/soil_unique_base.txt
	@for mode in binary fast long; do \
	   for num_threads in 1 2 3; do \
	     ./ntcard_mn.x --$${mode} --threads $${num_threads} soil_unique.fasta | diff --strip-trailing-cr -I '^#' - references/soil_unique_$${mode}.txt; \
	   done \
	done
	@echo "$@ passed"

.PHONY:test_hg38chr1
test_hg38chr1:ntcard_mn.x soil_unique.fasta
	@./ntcard_mn.x ${HG38CHR1} | diff --strip-trailing-cr -I "^#" - references/hg38chr1_base.txt
	@for mode in binary fast long; do \
	   ./ntcard_mn.x --$${mode} ${HG38CHR1} | diff --strip-trailing-cr -I '^#' - references/hg38chr1_$${mode}.txt; \
	done
	@echo "$@ passed"

SRR19536726_1_1000.fastq:../../testdata/SRR19536726_1_1000.fastq.gz
	@gzip -d -c $< > $@

70x_161nt_phred64.fastq.gz:../../testdata/70x_161nt_phred64.fastq
	@gzip -9 -c $< > $@

.PHONY:test_fastq
test_fastq:ntcard_mn.x SRR19536726_1_1000.fastq 70x_161nt_phred64.fastq.gz
	@for filename in ../../testdata/70x_161nt_phred64.fastq 70x_161nt_phred64.fastq.gz SRR19536726_1_1000.fastq ../../testdata/SRR19536726_1_1000.fastq.gz; do \
	  bfilename=`basename $${filename}`;\
	  bfilename=`echo $${bfilename} | cut -d '.' -f 1`;\
	  ./ntcard_mn.x $${filename} | diff --strip-trailing-cr -I "^#" - references/$${bfilename}_base.txt; \
	  for mode in binary fast long; do \
	    for num_threads in 1 2 3; do \
	      ./ntcard_mn.x --$${mode} --threads $${num_threads} $${filename} | diff --strip-trailing-cr -I '^#' - references/$${bfilename}_$${mode}.txt; \
	     done \
	   done \
	done
	@echo "$@ passed"

.PHONY:test_random
test_random:ntcard_mn.x
	@rm -rf TMP*
	@./test_cases.py './ntcard_mn.x -b'
	@rm -rf TMP*
	@./test_cases.py -l 100 -r 3 './ntcard_mn.x'
	@rm -rf TMP*
	@echo "$@ passed"

.PHONY:code_check
code_check:
	code_check.py -wt `find . -name '*.[hc]pp'`
	${GTTL}/scripts/check_ifndef.py `find . -name '*.hpp'`

.PHONY:tags
tags:
	@ctags -w --c++-kinds=+p --fields=+iaSKlm --extra=+q `ls *.[hc]pp` \
	`find ${GTTL}/src -name '*.hpp'` > $@

.PHONY:clean
clean:
	@${RM} -r *.[oxd] *.x.dSYM/ gqf/*.[oxd]
	@${RM} -r out TMP*

-include ${wildcard *.d}
