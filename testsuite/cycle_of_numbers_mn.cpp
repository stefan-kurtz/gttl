#include <cstdio>
#include <cstdlib>
#include <string>
#include "utilities/cycle_of_numbers.hpp"

int main(void)
{
  /* character need to be provided in alphabetical order.
     This does not imply an order. */
  const std::string characters{"ACGNTacgnt"};
  std::vector<uint8_t> forbidden_characters{};
  for (auto &&cc : characters)
  {
    forbidden_characters.push_back(static_cast<uint8_t>(cc));
  }
  CycleOfNumbers cycle(forbidden_characters);
  for (size_t idx = 0; idx < 512; idx++)
  {
    printf("%d\n",cycle.next());
  }
  return EXIT_SUCCESS;
}
