/*
  Copyright (c) 2021-2022 Stefan Kurtz <stefan.kurtz@uni-hamburg.de>
  Copyright (c) 2021-2022 Center for Bioinformatics, University of Hamburg

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

#ifndef SA_INDUCED_OPTIONS_HPP
#define SA_INDUCED_OPTIONS_HPP

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

class SainOptions
{
  public:
  using LcptabMethod = enum {Lcptab_no,
                             Lcptab_kasai13n,
                             Lcptab_kasai9n,
                             Lcptab_plcp5n};
 private:
  std::vector<std::string> inputfiles;
  std::string show_options_spec;
  std::string indexname;
  std::string lcptab_method_string;
  bool help_opt = false,
       verbose_opt = false,
       plain_input_format_opt = false,
       check_suftab_opt = false,
       abs_suftab_out_opt = false,
       abs_suftab_show_opt = false,
       rel_suftab_out_opt = false,
       rel_suftab_show_opt = false,
       lcptab_show_opt = false,
       tistab_out_opt = false,
       tistab_show_opt = false,
       buffered_option = false,
       reverse_complement_option = false,
       succinct_option = false;

  LcptabMethod lcptab_method = Lcptab_no;
  int intset_sizeof = -1;
 public:
  SainOptions(void) { }
  void parse(int argc, char **argv);

  bool help_opt_is_set(void) const noexcept;
  bool verbose_opt_is_set(void) const noexcept;
  bool plain_input_format_opt_is_set(void) const noexcept;
  bool check_suftab_opt_is_set(void) const noexcept;
  bool abs_suftab_out_opt_is_set(void) const noexcept;
  bool abs_suftab_show_opt_is_set(void) const noexcept;
  bool rel_suftab_out_opt_is_set(void) const noexcept;
  bool rel_suftab_show_opt_is_set(void) const noexcept;
  bool tistab_out_opt_is_set(void) const noexcept;
  bool tistab_show_opt_is_set(void) const noexcept;
  bool reverse_complement_option_is_set(void) const noexcept;
  bool buffered_option_is_set(void) const noexcept;
  bool lcptab_show_opt_is_set(void) const noexcept;
  bool succinct_option_is_set(void) const noexcept;
  LcptabMethod lcptab_method_get(void) const noexcept;
  std::string indexname_get(void) const noexcept;
  const std::vector<std::string> &inputfiles_get(void) const noexcept;
  int intset_sizeof_get(void) const noexcept;
};

#endif
