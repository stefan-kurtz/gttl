#include <iostream>
#include "sequences/alphabet.hpp"

int main(void)
{
  static constexpr const char nucleotides_upper_lower[] = "Aa|Cc|Gg|Tt";
  const GttlAlphabet<nucleotides_upper_lower,4> dna_alpha{};
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
  for (const char *s = nucleotides_upper_lower; *s != '\0'; s++)
  {
    if (*s != '|')
    {
      std::cout << *s << "\t" << static_cast<int>(dna_alpha.char_to_rank(*s))
                << std::endl;
    }
  }
  for (size_t idx = 0; idx < dna_alpha.size(); idx++)
  {
    std::cout << idx << "\t" << dna_alpha.rank_to_char(idx) << std::endl;
  }
  static constexpr const char amino_acids[]
    = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
  const GttlAlphabet<amino_acids,20> protein_alpha{};
  static_assert(protein_alpha.char_to_rank('A') == 0);
  static_assert(protein_alpha.char_to_rank('W') == 18);
  static_assert(protein_alpha.char_to_rank('Y') == 19);
  static_assert(protein_alpha.char_to_rank('X') == 20);
  for (const char *s = amino_acids ; *s != '\0'; s++)
  {
    if (*s != '|')
    {
      std::cout << *s << "\t" << static_cast<int>(dna_alpha.char_to_rank(*s))
                << std::endl;
    }
  }
  for (size_t idx = 0; idx < protein_alpha.size(); idx++)
  {
    std::cout << idx << "\t" << protein_alpha.rank_to_char(idx) << std::endl;
  }
  return EXIT_SUCCESS;
}
