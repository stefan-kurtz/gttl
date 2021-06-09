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
	@for prog in ${EXECS}; do \
	  echo "$$prog";\
	  ./$$prog;\
	done

%.x:%.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) ${CFLAGS} ${CPPFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d) ${LDLIBS}

.PHONY:clean
clean:
	${RM} -r *.x *.x.dSYM/
