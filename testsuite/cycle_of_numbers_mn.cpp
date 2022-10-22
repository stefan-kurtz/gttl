#include "utilities/cycle_of_numbers.hpp"
#include <cstdio>
#include <cstdlib>

int main(void)
{
  const char *characters = "ACGTNacgtn";
  CycleOfNumbers cycle(characters);
  for (size_t idx = 0; idx < 512; idx++)
  {
    printf("%d\n",cycle.next());
  }
  return EXIT_SUCCESS;
}
