GTTL=../..
include ${GTTL}/config.mk

OBJ=ntcard_mn.o ntcard_opt.o

all:ntcard_mn.x

ntcard_mn.x:${OBJ}
	$(LD) ${LDFLAGS} ${OBJ} -o $@ ${LDLIBS}

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
	     ./ntcard_mn.x --$${mode} --threads $${num_threads} ${HIQ} | diff --strip-trailing-cr -I '^#' - references/hiq_$${mode}.txt || exit 1; \
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
	     ./ntcard_mn.x --$${mode} --threads $${num_threads} eco29.fasta | diff --strip-trailing-cr -I '^#' - references/eco_$${mode}.txt || exit 1; \
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
	     ./ntcard_mn.x --$${mode} --threads $${num_threads} soil_unique.fasta | diff --strip-trailing-cr -I '^#' - references/soil_unique_$${mode}.txt || exit 1; \
	   done \
	done
	@echo "$@ passed"
endif

.PHONY:test_hg38chr1
test_hg38chr1:ntcard_mn.x soil_unique.fasta
ifdef HG38CHR1
	@./ntcard_mn.x ${HG38CHR1} | diff --strip-trailing-cr -I "^#" - references/hg38chr1_base.txt
	@for mode in binary fast long; do \
	   ./ntcard_mn.x --$${mode} ${HG38CHR1} | diff --strip-trailing-cr -I '^#' - references/hg38chr1_$${mode}.txt || exit 1; \
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
	  ./ntcard_mn.x $${filename} | diff --strip-trailing-cr -I "^#" - references/$${bfilename}_base.txt || exit 1; \
	  for mode in binary fast long; do \
	    for num_threads in 1 2 3; do \
	      ./ntcard_mn.x --$${mode} --threads $${num_threads} $${filename} | diff --strip-trailing-cr -I '^#' - references/$${bfilename}_$${mode}.txt || exit 1; \
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
