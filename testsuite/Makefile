CXX?=g++
CFLAGS=-g -m64 -Wall -Werror -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -funroll-loops -O3 -march=native -pedantic
CPPFLAGS=-I../src/
CXXFLAGS=-std=c++17
LDFLAGS=-g -m64
LDLIBS=-lm -lz -lpthread
SYSTEM?=$(shell uname -s)

ifeq ($(debug),yes)
  CFLAGS += -g
else
  CFLAGS += -DNDEBUG -march=native
endif

ifeq ($(SYSTEM),Darwin)
  ifeq ($(debug),yes)
    CFLAGS += -fsanitize=address -fno-omit-frame-pointer
    LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
  endif
endif

SRC=${wildcard  *.cpp}
EXECS=${subst .cpp,.x,${SRC}}

all:${EXECS}

test:${EXECS}
	./check_err.py --fail_for_protein_sequence enum_nthash.x
	./check_err.py --fail_for_protein_sequence count_non_wc_ranges.x
	./check_err.py multiseq_mn.x

%.o:%.cpp
	$(CXX) $(CXXFLAGS) ${CFLAGS} ${CPPFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d) $< -c -o $@ 

%.x:%.o
	$(CXX) ${LDFLAGS} $< -o $@ ${LDLIBS}

.PHONY:code_check
code_check:
	code_check.py -wt `find . -name '*.cpp'`

.PHONY:clean
clean:
	${RM} -r *.[odx] *.x.dSYM/

-include ${wildcard *.d}
