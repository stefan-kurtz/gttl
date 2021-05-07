#ifndef GTTL_FILE_OPEN_HPP
#define GTTL_FILE_OPEN_HPP

#ifndef QLI_WITHOUT_ZLIB
#include <zlib.h>
typedef gzFile QliFpType;
#define qli_fp_type_open(FP, MODE) gzopen(FP, MODE)
#define qli_fp_type_close(FP) gzclose(FP)
#else
typedef FILE * QliFpType;
#include <cstdio>
#define qli_fp_type_open(FP, MODE) fopen(FP, MODE)
#define qli_fp_type_close(FP) fclose(FP)
#endif
#endif
