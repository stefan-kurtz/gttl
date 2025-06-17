# What environment are we working in?
SYSTEM ?= $(shell uname -s)
MACHINE ?= $(shell uname -m)

# Set compiler,linker
CXX ?= g++
LD = ${CXX}

# We pass these to the preprocessor
CXX_PURE := $(shell echo ${CXX} | sed -e 's/^ccache //')
CXX_VERSION := $(shell basename ${CXX_PURE})

# Preprocessor flags
CPPFLAGS = -I ${GTTL}/src -DLARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DCXXOPTS_NO_REGEX -DCXX_VERSION=\"${CXX_VERSION}\"

# Compiler flags
CFLAGS = -Wall -Werror -Wextra -pedantic -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -funroll-loops -g -m64 -march=native
CXXFLAGS = -std=c++20

# Linker flags
LDFLAGS = -g -m64

# Show compilation time
ifeq (,$(findstring g++-,$(CXX)))
  ifneq ($(SYSTEM),Darwin)
    ifneq ($(OS,Windows_NT)
      TIME_OPTION=-time
    endif
  endif
endif

# Windows specific fixes, include ZLlib and disable some CRT-warnings from Microsoft
ifeq ($(OS),Windows_NT)
	CPPFLAGS += -D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_WARNINGS
	LDLIBS+=-L C:\\vcpkg\\packages\\zlib_x64-windows\\lib -lzlib
else
	LDLIBS+=-lm -lz -lpthread -lstdc++
  SHELL=/bin/bash
endif

# Call-graph profiling
ifeq ($(prof),yes)
  CFLAGS += -pg
  LDFLAGS += -pg
endif

# Conditionally enable valgrind
ifeq ($(valgrind),yes)
  VALGRIND=${GTTL}/scripts/valgrind.sh
endif

# Undefined-behavior sanitizer
ifeq ($(ubsan),yes)
  CFLAGS += -fsanitize=undefined -fno-omit-frame-pointer
  LDFLAGS += -fsanitize=undefined -fno-omit-frame-pointer
  debug := yes
endif

# Thread sanitizer
ifeq ($(tsan),yes)
  CFLAGS += -fsanitize=thread -fno-omit-frame-pointer
  LDFLAGS += -fsanitize=thread -fno-omit-frame-pointer
  debug := yes
endif

# Memory sanitizer
ifeq ($(msan),yes)
  CFLAGS += -fsanitize=memory -fno-omit-frame-pointer
  LDFLAGS += -fsanitize=memory -fno-omit-frame-pointer
  debug := yes
endif

# Address sanitizer
ifeq ($(asan),yes)
  CFLAGS += -fsanitize=address -fno-omit-frame-pointer
  LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
  debug := yes
endif

# We usually want this enabled. So we explicitly have to set debug=yes to disable it.
ifneq ($(debug),yes)
  CXXFLAGS += -Og
  CPPFLAGS += -DNDEBUG
  ifneq ($(SYSTEM),Darwin)
	 ifneq ($(OS),Windows_NT)
      LDLIBS+=-ldl
	 endif
  endif
else
  CXXFLAGS += -O3
endif

# Do not treat linking against outdated ABI-incompatible objects as an error.
ifeq (,$(findstring g++,$(CXX)))
  CFLAGS += -Wno-error=psabi
endif

# Compile without ZLib?
ifeq ($(without_zlib),yes)
  CPPFLAGS += -DGTTL_WITHOUT_ZLIB
endif

# This includes all auto-generated dependency files. .d files are generated in default_targets.mk
-include ${wildcard *.d}
