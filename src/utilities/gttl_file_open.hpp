#ifndef GTTL_FILE_OPEN_HPP
#define GTTL_FILE_OPEN_HPP

#ifndef GTTL_WITHOUT_ZLIB
#include <zlib.h>
using GttlFpType = gzFile;
#define gttl_fp_type_open(FP, MODE) gzopen(FP, MODE)
#define gttl_fp_type_close(FP)      gzclose(FP)
#else
using GttlFpType = FILE *;
#include <cstdio>
#define gttl_fp_type_open(FP, MODE) fopen(FP, MODE)
#define gttl_fp_type_close(FP)      fclose(FP)
#endif
#endif
