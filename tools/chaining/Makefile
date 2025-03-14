CXX?=g++
LD=${CXX}
CFLAGS=-g -m64 -Wall -Werror -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -funroll-loops -O3 -pedantic
GTTL=../..
CPPFLAGS=-I ${GTTL}/src
CXXFLAGS=-std=c++17
LDFLAGS=-g -m64
SYSTEM?=$(shell uname -s)
MACHINE?=$(shell uname -m)
MATCHES_ND = ${GTTL}/testdata/matches-nd.txt
LECTURE_EXAMPLE = ${GTTL}/testdata/chaining_lecture_example.txt
MATCHDATA = ${MATCHES_ND} \
            ${LECTURE_EXAMPLE}
REPFINDDATA = ${GTTL}/testdata/Duplicate.fna \
              ${GTTL}/testdata/at1MB.fna \
              ${GTTL}/testdata/ychrIII.fna \
              ${GTTL}/testdata/vaccg.fna

ifeq ($(OS),Windows_NT)
	CFLAGS += -D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_WARNINGS
	LDLIBS+=-L C:\\vcpkg\\packages\\zlib_x64-windows\\lib -lzlib
else
	LDLIBS+=-lm -lz -lpthread -lstdc++
  SHELL=/bin/bash
endif

ifeq (,$(findstring g++-,$(CXX)))
  CFLAG += -march=native
  TIME_OPTIONS=-time
endif

ifeq ($(debug),yes)
  CFLAGS += -g
  LDFLAGS += -g
  ifneq (${sanitize},no)
    ifneq (${MACHINE},arm64)
      CFLAGS += -fsanitize=address -fno-omit-frame-pointer
      LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
    endif
  endif
else
  CFLAGS += -DNDEBUG
  ifeq ($(SYSTEM),Darwin)
  else
	  ifeq ($(OS),Windows_NT)
	  else
    LDLIBS+=-ldl
	  endif
  endif
endif

ifeq ($(prof),yes)
  CFLAGS += -pg
  LDFLAGS += -pg
endif

SRC=${wildcard *.cpp}

chaining.x:chaining.o chaining_opt.o
	${LD} ${LDFLAGS} chaining.o chaining_opt.o -o $@ ${LDLIBS}

test_chaining.x:test_chaining.o
	${LD} ${LDFLAGS} test_chaining.o -o $@ ${LDLIBS}

%.o:%.cpp
	${CXX} ${TIME_OPTION} -o $@ -c $< ${CXXFLAGS} ${CFLAGS} ${CPPFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d)

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
		./chaining.x ${TMPFILE} | diff --strip-trailing-cr - testdata/chain-$$lines-global.txt; \
    ./chaining.x --local ${TMPFILE} | diff --strip-trailing-cr - testdata/chain-$$lines-local.txt; \
	done
	@${RM} ${TMPFILE}
	@echo "Congratulations, $@ passed."

.PHONY:test_large
test_large:chaining.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@for lines in 8000 10766 58655 62935 64833 733000 733536 733537 800000 900000; do \
		echo '# select first '$$lines' lines of ecoli_vs_ecoli.of'; \
		head -n $$lines testdata/ecoli_vs_ecoli.of > ${TMPFILE}; \
    ./chaining.x ${TMPFILE} | diff --strip-trailing-cr - testdata/chain-$$lines-global.txt; \
    ./chaining.x ${TMPFILE} | diff --strip-trailing-cr - testdata/chain-$$lines-local.txt; \
	done
	@${RM} ${TMPFILE}
	@echo "Congratulations, $@ passed."

.PHONY:test_repfind
test_repfind:chaining.x
	@echo "$@"
	@$(eval TMPFILE := $(shell mktemp --tmpdir=.))
	@for file in ${REPFINDDATA}; do \
		echo '# compute matches for '$$file; \
    ./chaining.x testdata/$$(basename $$file)-repdata.txt --silent | diff --strip-trailing-cr - testdata/chain-$$(basename $$file)-global.txt; \
    ./chaining.x testdata/$$(basename $$file)-repdata.txt --silent --local | diff --strip-trailing-cr - testdata/chain-$$(basename $$file)-local.txt; \
	done
	@${RM} ${TMPFILE}
	@echo "Congratulations, $@ passed."

.PHONY:test_time
test_time:chaining.x
	@echo "# gt"
	@time GT_ENV_OPTIONS= gt chain2dim -m testdata/ecoli_vs_ecoli.of -silent
	@echo "# chaining.x"
	@time ./chaining.x --silent testdata/ecoli_vs_ecoli.of

.PHONY:tags
tags:
	@ctags -w --c++-kinds=+p --fields=+iaSKlm --extra=+q `ls *.[hc]pp` \
        `find ${GTTL}/src -name '*.hpp'` > $@

.PHONY:code_check
code_check:
	@code_check.py -wt `find . -maxdepth 1 -name '*.[hc]pp'`
	@code_check.py -wt `find . -maxdepth 1 -name '*.py'`

.PHONY:clean
clean:
	@${RM} *.[oxd] tmp.*
	@${RM} -r sfx.* *.x.dSYM/ __pycache__

-include ${wildcard *.d}
