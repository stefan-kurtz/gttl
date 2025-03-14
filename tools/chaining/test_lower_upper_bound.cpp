#include <iostream>
#include <vector>
#include <algorithm>

int main(void)
{
  std::vector<int> vec{2,3,4,6,7};

  // precursor
  /*
   std::lower_bound return iterator to the lower bound of val in the range.
   If all the elements in the range compare less than val, the function returns
   last.
   If all the elements in the range are larger than val, the function returns
   a pointer to the first element.
  */
  std::cout << *(std::lower_bound(vec.begin(), vec.end(), 5,
                                  [&](const int& a, const int& b)
                                     { return a < b; }) - 1) << std::endl;
  std::cout << *std::prev((std::lower_bound(vec.begin(), vec.end(), 5,
                                            [&](const int& a, const int& b)
                                               { return a < b; })))
            << std::endl;
  std::cout << *std::prev((std::lower_bound(vec.begin(), vec.end(), 4,
                                            [&](const int& a, const int& b)
                                               { return a < b; })))
            << std::endl;

  // postcursor
  /*
   std::upper_bound return an iterator to the upper bound of val in the range.
   If all the element in the range compare less than val, the function returns
   last.
  */
  std::cout << *std::upper_bound(vec.begin(),vec.end(), 5,
                                 [&](const int& b, const int& a)
                                    { return a > b; }) << std::endl;
  std::cout << *std::upper_bound(vec.begin(),vec.end(), 6,
                                 [&](const int& b, const int& a)
                                    { return a > b; }) << std::endl;
  return EXIT_SUCCESS;
}
