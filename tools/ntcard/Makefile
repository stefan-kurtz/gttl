CXX?=g++
LD=${CXX}
CFLAGS=-g -Wall -Werror -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -funroll-loops -O3
GTTL=../..
CPPFLAGS=-I ${GTTL}/src
CXXFLAGS=-std=c++20
LDFLAGS=-g

ifeq ($(OS),Windows_NT)
	CFLAGS += -D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_WARNINGS
	LDLIBS+=-L C:\\vcpkg\\packages\\zlib_x64-windows\\lib -lzlib
else
	LDLIBS+=-lm -lz -lpthread -lstdc++
endif


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
test:test_random test_fastq test_large
	@echo "$@ passed"

.PHONY:test_large
test_large:test_hiq test_eco test_soil_unique test_hg38chr1
ifeq ($(and $(HIQ),$(ECO29_SINGLE_FILE),$(SOIL_UNIQUE),$(HG38CHR1)),)
	@echo "The following environment vairables for large-file testing are not set:"
endif
ifndef HIQ
	@echo "- HIQ"
endif
ifndef ECO29_SINGLE_FILE
	@echo "- ECO29_SINGLE_FILE"
endif
ifndef SOIL_UNIQUE
	@echo "- SOIL_UNIQUE"
endif
ifndef HG38CHR1
	@echo "- HG38CHR1"
endif
	@echo "$@ passed"

.PHONY:test_hiq
test_hiq:ntcard_mn.x
ifdef HIQ
	@./ntcard_mn.x ${HIQ} | diff --strip-trailing-cr -I "^#" - references/hiq_base.txt
	@for mode in binary fast long ; do \
	   for num_threads in 1 2 3; do \
	     ./ntcard_mn.x --$${mode} --threads $${num_threads} ${HIQ} | diff --strip-trailing-cr -I '^#' - references/hiq_$${mode}.txt; \
	   done \
	done
	@echo "$@ passed"
endif

eco29.fasta:${ECO29_SINGLE_FILE}
ifdef ECO29_SINGLE_FILE
	@gzip -d -c $< > $@
endif

.PHONY:test_eco
test_eco:ntcard_mn.x eco29.fasta
ifdef ECO29_SINGLE_FILE
	@./ntcard_mn.x eco29.fasta | diff --strip-trailing-cr -I "^#" - references/eco_base.txt
	@for mode in binary fast long; do \
	   for num_threads in 1 2 3; do \
	     ./ntcard_mn.x --$${mode} --threads $${num_threads} eco29.fasta | diff --strip-trailing-cr -I '^#' - references/eco_$${mode}.txt; \
	   done \
	done
	@echo "$@ passed"
endif

soil_unique.fasta:${SOIL_UNIQUE}
ifdef SOIL_UNIQUE
	@gzip -d -c $< > $@
endif

.PHONY:test_soil_unique
test_soil_unique:ntcard_mn.x soil_unique.fasta
ifdef SOIL_UNIQUE
	@./ntcard_mn.x soil_unique.fasta | diff --strip-trailing-cr -I "^#" - references/soil_unique_base.txt
	@for mode in binary fast long; do \
	   for num_threads in 1 2 3; do \
	     ./ntcard_mn.x --$${mode} --threads $${num_threads} soil_unique.fasta | diff --strip-trailing-cr -I '^#' - references/soil_unique_$${mode}.txt; \
	   done \
	done
	@echo "$@ passed"
endif

.PHONY:test_hg38chr1
test_hg38chr1:ntcard_mn.x soil_unique.fasta
ifdef HG38CHR1
	@./ntcard_mn.x ${HG38CHR1} | diff --strip-trailing-cr -I "^#" - references/hg38chr1_base.txt
	@for mode in binary fast long; do \
	   ./ntcard_mn.x --$${mode} ${HG38CHR1} | diff --strip-trailing-cr -I '^#' - references/hg38chr1_$${mode}.txt; \
	done
	@echo "$@ passed"
endif

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
	@${RM} -r out TMP* __pycache__

-include ${wildcard *.d}
