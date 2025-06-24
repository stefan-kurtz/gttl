#include "sequences/alphabet.hpp"
#include <cstdlib>

int main(void)
{
  const alphabet::GttlAlphabet_UL_4 dna_alpha{};
  static_assert(dna_alpha.char_to_rank('A') == 0);
  static_assert(dna_alpha.char_to_rank('a') == 0);
  static_assert(dna_alpha.char_to_rank('C') == 1);
  static_assert(dna_alpha.char_to_rank('c') == 1);
  static_assert(dna_alpha.char_to_rank('G') == 2);
  static_assert(dna_alpha.char_to_rank('g') == 2);
  static_assert(dna_alpha.char_to_rank('T') == 3);
  static_assert(dna_alpha.char_to_rank('t') == 3);
  static_assert(dna_alpha.char_to_rank('X') == dna_alpha.undefined_rank());
  static_assert(dna_alpha.size() == 4);
  static_assert(dna_alpha.rank_to_char(0) == 'A');
  static_assert(dna_alpha.rank_to_char(1) == 'C');
  static_assert(dna_alpha.rank_to_char(2) == 'G');
  static_assert(dna_alpha.rank_to_char(3) == 'T');
  static_assert((dna_alpha.characters_get())[0] == 'A');
  static_assert((dna_alpha.characters_get())[1] == 'C');
  static_assert((dna_alpha.characters_get())[2] == 'G');
  static_assert((dna_alpha.characters_get())[3] == 'T');
  dna_alpha.pretty_print();
  const alphabet::GttlAlphabet_20 protein_alpha{};
  static_assert(protein_alpha.char_to_rank('A') == 0);
  static_assert(protein_alpha.char_to_rank('W') == 18);
  static_assert(protein_alpha.char_to_rank('Y') == 19);
  static_assert(protein_alpha.char_to_rank('X') == 20);
  protein_alpha.pretty_print();
  return EXIT_SUCCESS;
}
