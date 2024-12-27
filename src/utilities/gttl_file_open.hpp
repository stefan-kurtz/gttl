#ifndef GTTL_FILE_OPEN_HPP
#define GTTL_FILE_OPEN_HPP

#ifndef GTTL_WITHOUT_ZLIB
#include <zlib.h>
using GttlFpType = gzFile;
#define gttl_fp_type_open(FILENAME, MODE)   gzopen(FILENAME, MODE)
#define gttl_fp_type_close(FP)              gzclose(FP)
#define gttl_fp_type_reset(FP)              gzseek(FP, 0, SEEK_SET)
#define gttl_fp_type_gets(FP, BUFFER, SIZE) gzgets(FP, BUFFER, SIZE)
#define gttl_fp_type_is_eof(FP)             gzeof(FP)
#else
using GttlFpType = FILE *;
#include <cstdio>
#define gttl_fp_type_open(FILENAME, MODE)   fopen(FILENAME, MODE)
#define gttl_fp_type_close(FP)              fclose(FP)
#define gttl_fp_type_reset(FP)              fseek(FP, 0, SEEK_SET)
#define gttl_fp_type_gets(FP, BUFFER, SIZE) (fgets(BUFFER, SIZE, FP) != nullptr\
                                             ? BUFFER : nullptr)
#define gttl_fp_type_is_eof(FP)             feof(FP)
#endif
#endif
