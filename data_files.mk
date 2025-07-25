# Testsuite files
AT1MB=${GTTL}/testdata/at1MB.fna
SW175=${GTTL}/testdata/sw175.fna
Y3=${GTTL}/testdata/ychrIII.fna
VAC=${GTTL}/testdata/vaccg.fna
PROTEIN_FILE=${GTTL}/testdata/protein.fsa

# Chaining files
MATCHES_ND = ${GTTL}/testdata/matches-nd.txt
LECTURE_EXAMPLE = ${GTTL}/testdata/chaining_lecture_example.txt
MATCHDATA = ${MATCHES_ND} \
            ${LECTURE_EXAMPLE}
REPFINDDATA = ${GTTL}/testdata/Duplicate.fna \
              ${GTTL}/testdata/at1MB.fna \
              ${GTTL}/testdata/ychrIII.fna \
              ${GTTL}/testdata/vaccg.fna

# All FastA files (Extensions according to Wikipedia)
CORRUPT_FILES = $(shell cat ${GTTL}/testdata/corrupt_files.txt)
EXCLUDE_PATTERNS = $(foreach f,$(CORRUPT_FILES),! -name $(f))
GTTL_FASTA_FILES = $(shell find ${GTTL}/testdata -type f \
                                                     \( -name "*.fna" \
                                                     -o -name "*.fsa" \
                                                     -o -name "*.fasta" \
                                                     -o -name "*.fas" \
                                                     -o -name "*.fa" \
                                                     -o -name "*.ffn" \
                                                     -o -name "*.faa" \
                                                     -o -name "*.mpfa" \
                                                     -o -name "*.frn" \) \
                                                     $(EXCLUDE_PATTERNS))

# Same, but for FastQ
GTTL_FASTQ_FILES = $(shell find ${GTTL}/testdata -type f -name "*.fastq" \
                                                     -o -name "*.fq")
