#ifndef LCS_LCP_LEN_TYPE_HPP
#define LCS_LCP_LEN_TYPE_HPP
#include <cstddef>

using LcsLcpLenType = size_t (*)(const char *,size_t,
                                 const char *,size_t,
                                 size_t,size_t,
                                 size_t,size_t);

#endif
