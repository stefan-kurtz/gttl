#!/usr/bin/env python3

import sys, argparse
from jinja2 import Environment, Template
# https://realpython.com/primer-on-jinja-templating/#install-jinja

environment = Environment()
template \
  = environment.from_string('''/* created by {{ program_call }} DO NOT EDIT */
#ifndef {{ unique_file_key }}
#define {{ unique_file_key }}
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <format>
#include "alignment/blosum62.hpp"
#include "alignment/unit_score_aa.hpp"
#include "alignment/unit_score_nuc.hpp"
#include "alignment/unit_score_nuc_2_2.hpp"
#include "alignment/unit_score_nuc_lower.hpp"
#include "alignment/unit_score_nuc_upper.hpp"
#include "alignment/score_class_base.hpp"
#include "alignment/score_matrix_name.hpp"
#include "sequences/gttl_multiseq.hpp"
{{ include_extra }}

static {{ return_type }} {{ function_name }}
  (const char *score_matrix_id,
   const ScoreMatrixName &score_matrix_name,
   bool dna_alphabet {{ additional_arguments }})
{
  if (!dna_alphabet)
  {
    if (score_matrix_name.is(Score_matrix_undefined) ||
        score_matrix_name.is(Score_matrix_blosum62))
    {
      {{ do_replacement('Blosum62') }}
    } else
    {
      if (score_matrix_name.is(Score_matrix_unit_score_aa))
      {
        {{ do_replacement('Unit_score_aa') }}
      } else
      {
        const ScoreMatrixName smn_instance{};
        auto score_matrix_names = smn_instance.string_values_joined(", ");
        throw std::runtime_error(
                std::format(": score matrix {} is not possible for "
                            "protein sequences; the following "
                            "choices are available: {}",
                            score_matrix_id,
                            score_matrix_names));
      }
    }
  } else
  {
    if (score_matrix_name.is(Score_matrix_undefined) ||
        score_matrix_name.is(Score_matrix_unit_score_nuc))
    {
      {{ do_replacement('Unit_score_nuc') }}
    } else
    {
      if (score_matrix_name.is(Score_matrix_undefined) ||
          score_matrix_name.is(Score_matrix_unit_score_nuc_2_2))
      {
        {{ do_replacement('Unit_score_nuc_2_2') }}
      } else
      {
        if (score_matrix_name.is(Score_matrix_unit_score_nuc_lower))
        {
          {{ do_replacement('Unit_score_nuc_lower') }}
        } else
        {
          if (score_matrix_name.is(Score_matrix_unit_score_nuc_upper))
          {
            {{ do_replacement('Unit_score_nuc_upper') }}
          } else
          {
            const ScoreMatrixName smn_instance{};
            auto score_matrix_names = smn_instance.string_values_joined(", ");
            throw std::runtime_error(
                    std::format(": score matrix {} is not possible "
                                "for DNA sequences; the following "
                                "choices are available: {}",
                                score_matrix_id,
                                score_matrix_names));
          }
        }
      }
    }
  }
  {{ final_return }}
}
#endif /* {{ unique_file_key }} */''')

def parse_arguments(argv):
  p = argparse.ArgumentParser(description=('generate code involving case '
                                           'distinctions on score matrices'))
  choices = p.add_mutually_exclusive_group(required=True)
  choices.add_argument('--scoring',action='store_true',default=False,
                       help=('generate code for '
                             'scoring_info_and_sequence_transformation'))
  choices.add_argument('--alignment_output',action='store_true',default=False,
                       help=('generate code for '
                             'alignment_output_function_get'))
  choices.add_argument('--sequence_decode',action='store_true',default=False,
                       help=('generate code for sequence_decode_function_get'))
  return p.parse_args(argv)

this_program_call=' '.join(sys.argv)
args = parse_arguments(sys.argv[1:])

if args.scoring:
  def do_replacement(x):
    return ('''return {{scorematrix2D_get<{}::num_of_chars>({}::score_matrix),
                      {}::smallest_score,
                      literate_multiseqs<{}>(db_multiseq,query_multiseq)}};''')\
            .format(x,x,x,x)

  template.globals['do_replacement'] = do_replacement

  print(template.render(
          program_call=this_program_call,
          unique_file_key='SCORING_INFO_AND_SEQ_TRANS_HPP',
          return_type='std::tuple<int8_t **,int8_t,size_t>',
          include_extra='',
          function_name='scoring_info_and_seq_trans',
          additional_arguments=(',GttlMultiseq *db_multiseq,'
                                'GttlMultiseq *query_multiseq'),
          final_return=('return std::tuple<int8_t **,int8_t,size_t>'
                        '(nullptr,INT8_MAX,0);')))
elif args.alignment_output:
  def do_replacement(x):
    return ('''return alignment_output<GttlSubstring<uint8_t>,
                                       GttlSubstring<uint8_t>,
                                       uint8_t,
                                       encoded_matching_characters<{}>,
                                       to_char_map<{}>>;'''.format(x,x))
  template.globals['do_replacement'] = do_replacement

  print(template.render(
          program_call=this_program_call,
          unique_file_key='ALIGNMENT_OUTPUT_FUNCTION_HPP',
          include_extra=('#include "sequences/alignment_output.hpp"\n'
                         '#include "sequences/gttl_substring.hpp"'),
          return_type='auto',
          function_name='alignment_output_function_get',
          additional_arguments='',
          final_return=''))
elif args.sequence_decode:
  def do_replacement(x):
    return ('''return sequence_decode<GttlSubstring<uint8_t>,
                                      uint8_t,
                                      to_char_map<{}>>;'''.format(x))
  template.globals['do_replacement'] = do_replacement

  print(template.render(
          program_call=this_program_call,
          unique_file_key='SEQUENCE_DECODE_FUNCTION_HPP',
          include_extra=('''#include "sequences/gttl_substring.hpp"

template<class SeqClass,
         typename CharType,
         char (&to_char)(CharType)>
static std::string sequence_decode(const SeqClass &seq)
{
  std::string s{};
  for (size_t idx = 0; idx < seq.size(); idx++)
  {
    s += to_char(seq[idx]);
  }
  return s;
}'''),
          return_type='auto',
          function_name='sequence_decode_function_get',
          additional_arguments='',
          final_return=''))
