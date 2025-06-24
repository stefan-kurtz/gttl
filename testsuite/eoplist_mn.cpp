#include <cassert>
#include <cmath>
#include <algorithm>
#include <exception>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include "utilities/uniform_random_double.hpp"
#include "sequences/eoplist.hpp"
#include "sequences/alignment_output.hpp"

static bool simple_matching_characters(char a,char b) { return a == b; }

static char to_char_identity(char cc) { return cc; }

class IdentityEncodingHandler
{
  public:
  using source_type = char;
  static bool match_method(source_type a,source_type b) { return a == b; }
  static char to_char(source_type a) { return a; }
};

static void display_alignment(const Eoplist &eoplist)
{
  const size_t ulen = eoplist.count_deletions_get() +
                      eoplist.count_mismatches_get() +
                      eoplist.count_matches_get();
  const size_t vlen = eoplist.count_insertions_get() +
                      eoplist.count_mismatches_get() +
                      eoplist.count_matches_get();
  char *useq = new char [ulen];
  char *vseq = new char [vlen];
  size_t i = 0;
  size_t j = 0;
  const char characters[] = "ACGT";
  size_t char_idx = 0;
  for (auto &&co : eoplist)
  {
    for (size_t iter = 0; iter < co.iteration; iter++)
    {
      switch (co.edit_operation)
      {
        case MatchOp:
          char_idx++;
          assert(i < ulen && j < vlen);
          useq[i++] = characters[char_idx % 4];
          vseq[j++] = characters[char_idx % 4];
          break;
        case MismatchOp:
          char_idx++;
          assert(i < ulen);
          useq[i++] = characters[char_idx % 4];
          char_idx++;
          assert(j < vlen);
          vseq[j++] = characters[char_idx % 4];
          break;
        case DeletionOp:
          char_idx++;
          assert(i < ulen);
          useq[i++] = characters[char_idx % 4];
          break;
        case InsertionOp:
          char_idx++;
          assert(j < vlen);
          vseq[j++] = characters[char_idx % 4];
          break;
        default:
          assert(false);
      }
    }
  }
  if (i != ulen)
  {
    std::cout << "i=" << i << "!=" << ulen << "=ulen\n";
  }
  assert(i == ulen);
  assert(j == vlen);
  const std::string ustring(useq,ulen);
  const std::string vstring(vseq,vlen);
  AlignmentSequenceInfo<std::string, std::string> asi(&ustring,0,&vstring,0,
                                                      ulen,vlen);
  static constexpr const size_t top_seqlength = 0;
  static constexpr const size_t low_reference = 0;
  static constexpr const size_t one_off = 0;
  static constexpr const bool subject_first = true;
  static constexpr const bool distinguish_mismatch_match = true;
  static constexpr const bool forward_strand = true;
  alignment_output<std::string,std::string,char,
                   simple_matching_characters,to_char_identity>
                  (asi,
                   eoplist,
                   top_seqlength,
                   low_reference,
                   one_off,
                   subject_first,
                   distinguish_mismatch_match,
                   60,
                   forward_strand,
                   vlen,
                   vlen,
                   stdout);
  delete[] useq;
  delete[] vseq;
}

int main(int argc,char *argv[])
{
  constexpr const bool debug = false;
  long trials;
  if (argc != 3 || (strcmp(argv[1],"silent") != 0 &&
                    strcmp(argv[1],"display") != 0) ||
      sscanf(argv[2],"%ld",&trials) != 1 || trials < 1)
  {
    std::cerr << "Usage: " << argv[0] << " silent|display <trials>\n";
    return EXIT_FAILURE;
  }
  const bool display = strcmp(argv[1],"display") == 0;
  UniformRandomDouble rgen(0,1.0,12324);
  constexpr const size_t match_length_array_size = 50;
  size_t match_length_array[match_length_array_size];
  for (size_t idx = 0; idx < match_length_array_size; idx++)
  {
    match_length_array[idx]
      = static_cast<size_t>(ceil(10 * exp(static_cast<double>((idx+1)/10))));
  }
  const double boundary_distinguish_mismatch_match[] = {0.25, 0.5, 0.75};
  const double boundary_identify_mismatch_match[] = {0.33, 0.67, 1.0};
  const double *boundary = nullptr;
  for (size_t trial = 0; trial < static_cast<size_t>(trials); trial++)
  {
    for (int mode = 0; mode < 2; mode++)
    {
      bool distinguish_mismatch_match;
      if (mode == 0)
      {
        distinguish_mismatch_match = true;
        boundary = boundary_distinguish_mismatch_match;
      } else
      {
        distinguish_mismatch_match = false;
        boundary = boundary_identify_mismatch_match;
      }
      Eoplist eoplist(distinguish_mismatch_match);
      const size_t num_eops = 1 + static_cast<size_t>(10 * rgen.get());
      bool last_was_match = false;
      size_t count_edit_operations = 0;
      for (size_t idx = 0; idx < num_eops; idx++)
      {
        const double r = rgen.get();
        if (r <= boundary[0])
        {
          if (!last_was_match)
          {
            size_t ml_idx = static_cast<size_t>(match_length_array_size *
                                                    rgen.get());
            ml_idx = std::min(match_length_array_size-1,ml_idx);
            assert(ml_idx < match_length_array_size);
            const size_t match_length = match_length_array[ml_idx];
            if constexpr (debug)
            {
              std::cout << "match_add(" << match_length << ")\n";
            }
            eoplist.match_add(match_length);
            last_was_match = true;
            count_edit_operations++;
          }
        } else
        {
          if (r <= boundary[1])
          {
            if constexpr (debug)
            {
              std::cout << "deletion_add\n";
            }
            eoplist.deletion_add();
            last_was_match = false;
            count_edit_operations++;
          } else
          {
            if (r <= boundary[2])
            {
              if constexpr (debug)
              {
                std::cout << "insertion_add\n";
              }
              eoplist.insertion_add();
              last_was_match = false;
              count_edit_operations++;
            } else
            {
              if (distinguish_mismatch_match)
              {
                if constexpr (debug)
                {
                  std::cout << "mismatch_add\n";
                }
                eoplist.mismatch_add();
                last_was_match = false;
                count_edit_operations++;
              }
            }
          }
        }
      }
      if constexpr (debug)
      {
        std::cout << "# number of editoperations\t" << count_edit_operations
                  << '\n';
      }
      std::string cigar_string
        = eoplist.cigar_string_get(distinguish_mismatch_match);
      Eoplist *eoplist_decode = nullptr;
      try
      {
        eoplist_decode = new Eoplist(distinguish_mismatch_match,cigar_string);
      }
      catch (const std::exception &err) /* check_err.py */
      {
        std::cerr << argv[0] << " " << trials << " failed: " << err.what()
                  << '\n';
        delete eoplist_decode;
        return EXIT_FAILURE;
      }
      if (eoplist != *eoplist_decode)
      {
        std::cerr << "distinguish_mismatch_match\t"
                  << (distinguish_mismatch_match ? "true" : "false") << '\n';
        std::cerr << eoplist.to_string() << '\n';
        std::cerr << cigar_string << '\n';
        std::cerr << eoplist_decode->to_string() << '\n';
        return EXIT_FAILURE;
      }
      if (display)
      {
        std::cout << cigar_string << '\n';
        display_alignment(eoplist);
      }
      delete eoplist_decode;
    }
  }
}
