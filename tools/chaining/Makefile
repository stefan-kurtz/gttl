GTTL=../..
include ${GTTL}/config.mk

SRC=${wildcard *.cpp}

chaining.x:chaining.o chaining_opt.o
	${LD} ${LDFLAGS} chaining.o chaining_opt.o -o $@ ${LDLIBS}

.PHONY:all
all: chaining.x

.PHONY:test
test: test_small test_large test_repfind

# ecolicmp250.of.gz is from ${GTDIR}/testdata

.PHONY:test_small
test_small:chaining.x
	@echo "$@"
	@./chaining.x ${MATCHES_ND} --silent | grep -q "length 2 score 2"
	@./chaining.x ${LECTURE_EXAMPLE} --silent | grep -q "length 3 score 3"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@for lines in 1 2 4 20 50 283 298 304 665 753 800 900 1500; do \
		gzip -d -c ecolicmp250.of.gz | head -n $$lines | sed -r 's/\t/ /g' > ${TMPFILE}; \
		./chaining.x ${TMPFILE} | diff --strip-trailing-cr - testdata/chain-$$lines-global.txt || exit 1; \
    ./chaining.x --local ${TMPFILE} | diff --strip-trailing-cr - testdata/chain-$$lines-local.txt || exit 1; \
	done
	@${RM} ${TMPFILE}
	@echo "Congratulations, $@ passed."

.PHONY:test_large
test_large:chaining.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@for lines in 8000 10766 58655 62935 64833 733000 733536 733537 800000 900000; do \
		echo '# select first '$$lines' lines of ecoli_vs_ecoli.of'; \
		head -n $$lines testdata/ecoli_vs_ecoli.of > ${TMPFILE} || exit 1; \
    ./chaining.x ${TMPFILE} | diff --strip-trailing-cr - testdata/chain-$$lines-global.txt || exit 1; \
    ./chaining.x --local ${TMPFILE} | diff --strip-trailing-cr - testdata/chain-$$lines-local.txt || exit 1; \
	done
	@${RM} ${TMPFILE}
	@echo "Congratulations, $@ passed."

.PHONY:test_repfind
test_repfind:chaining.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@for file in ${REPFINDDATA}; do \
		echo '# compute matches for '$$file; \
    ./chaining.x testdata/$$(basename $$file)-repdata.txt --silent | diff --strip-trailing-cr - testdata/chain-$$(basename $$file)-global.txt || exit 1; \
    ./chaining.x testdata/$$(basename $$file)-repdata.txt --silent --local | diff --strip-trailing-cr - testdata/chain-$$(basename $$file)-local.txt || exit 1; \
	done
	@${RM} ${TMPFILE}
	@echo "Congratulations, $@ passed."

.PHONY:test_time
test_time:chaining.x
	@echo "# gt"
	@time GT_ENV_OPTIONS= gt chain2dim -m testdata/ecoli_vs_ecoli.of -silent
	@echo "# chaining.x"
	@time ./chaining.x --silent testdata/ecoli_vs_ecoli.of
