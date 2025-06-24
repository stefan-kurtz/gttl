/*
  Copyright (c) 2021 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2021 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
#ifndef QGRAMS_HASH_NTHASH_HPP
#define QGRAMS_HASH_NTHASH_HPP
#include <cstdint>
#include "sequences/alphabet.hpp"
#include "sequences/nthash_fwd.hpp"
#include "sequences/nthash_fwd_aminoacids.hpp"
#include "sequences/qgrams_rec_hash_value_fwd_iter.hpp"
#include "sequences/qgrams_rec_hash_value_iter.hpp"

template<uint8_t undefined_rank>
using QgramNtHashFwdIteratorGeneric
  = QgramRecHashValueFwdIterator<alphabet::nucleotides_upper_lower,
                                 undefined_rank,
                                 NThashTransformer,
                                 char>;

using QgramNtHashFwdIterator4 = QgramNtHashFwdIteratorGeneric<4>;

template<uint8_t undefined_rank>
using QgramNtHashFwdIteratorGenericNoTransform
  = QgramRecHashValueFwdIterator<alphabet::nucleotides_upper_lower,
                                 undefined_rank,
                                 NThashTransformer,
                                 uint8_t>;

using QgramNtHashIterator4
  = QgramRecHashValueIterator<alphabet::nucleotides_upper_lower,
                              4,
                              NThashTransformer>;


template <uint8_t undefined_rank>
using QgramNtHashAAFwdIteratorGeneric =
  QgramRecHashValueFwdIterator<alphabet::amino_acids,
                               undefined_rank,
                               NtHashAminoacidsTransformer,
                               char>;

using QgramNtHashAAFwdIterator20 = QgramNtHashAAFwdIteratorGeneric<20>;

template <uint8_t undefined_rank>
using QgramNtHashAAFwdIteratorGenericNoTransform
= QgramRecHashValueFwdIterator<alphabet::amino_acids,
                               undefined_rank,
                               NtHashAminoacidsTransformer,
                               uint8_t>;

#endif
