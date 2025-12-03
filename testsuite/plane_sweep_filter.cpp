#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "sequences/stored_match.hpp"
#include "sequences/non_redundant_matches.hpp"

int main(void)
{
  std::vector<GttlStoredMatch> matches{};
  matches.emplace_back(1,3,2,6,0,3.6);
  matches.emplace_back(2,3,4,6,0,3.3);
  matches.emplace_back(3,3,5,7,0,5.0);
  matches.emplace_back(7,3,4,6,0,4.6);
  matches.emplace_back(8,6,5,7,0,2.3);
  matches.emplace_back(9,2,4,6,0,6.0);
  matches.emplace_back(12,1,5,7,0,15.0);

  std::cout << "# query_start,query_end,ref_start,ref_end,identity\n";
  std::cout << "# before filtering\n";
  for (const auto & match : matches)
  {
    std::cout << match.to_string() << '\n';
  }
  const NonRedundantMatches<GttlStoredMatch> non_redundant_matches(matches);
  std::cout << "# after filtering\n";
  for (size_t idx = 0; idx < matches.size(); idx++)
  {
    if (non_redundant_matches[idx])
    {
      std::cout << matches[idx].to_string() << '\n';
    }
  }
  return EXIT_SUCCESS;
}
