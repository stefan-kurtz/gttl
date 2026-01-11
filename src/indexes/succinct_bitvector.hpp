#ifndef SUCCINCT_BITVECTOR_HPP
#define SUCCINCT_BITVECTOR_HPP

#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <stdexcept>

struct Superblock {
  uint64_t first;
  uint64_t second;
};

class SuccinctBitvector
{
  private:
  std::vector<uint64_t> data_vector;
  size_t length;
  std::vector<Superblock> rank;
  std::vector<size_t> select_1;

  static constexpr const uint64_t INDEXINTERVAL = 8192;
  static constexpr const uint64_t SUPERBLOCKBITS = 4096;
  static constexpr const uint64_t BLOCKBITS = 512;

  public:

  SuccinctBitvector(void)
    : length(0)
  { }

  SuccinctBitvector(const std::string &lls_filename)
    : length(0)
  {
    //deserialize
    FILE *const in_fp = std::fopen(lls_filename.c_str(), "rb");
    if (in_fp == nullptr)
    {
      throw std::runtime_error(std::string("cannot open file ") + lls_filename);
    }
    if (std::fread(&length, sizeof(size_t), 1, in_fp) != 1)
    {
      throw std::runtime_error(std::string("cannot read length of type size_t "
                                           "from file ") + lls_filename);
    }

    size_t number_of_superblocks;
    if (std::fread(&number_of_superblocks, sizeof(size_t), 1, in_fp) != 1)
    {
      throw std::runtime_error(std::string("cannot read number_of_superblocks "
                                           " of type size_t from file ")
                               + lls_filename);
    }
    rank.resize(number_of_superblocks);
    if (std::fread(rank.data(), sizeof(Superblock),
                   number_of_superblocks, in_fp) != number_of_superblocks)
    {
      throw std::runtime_error(std::string("cannot read expected number of ") +
                               std::to_string(number_of_superblocks) +
                               std::string(" superblocks from file ")
                               + lls_filename);
    }
    size_t length_of_select1;
    if (std::fread(&length_of_select1, sizeof(size_t), 1, in_fp) != 1)
    {
      throw std::runtime_error(std::string("cannot read length_of_select1 of "
                                           "type size_t from file ")
                               + lls_filename);
    }
    select_1.resize(length_of_select1);
    if (std::fread(select_1.data(), sizeof(size_t),
                   length_of_select1, in_fp) != length_of_select1)
    {
      throw std::runtime_error(std::string("cannot read expected number of ") +
                               std::to_string(length_of_select1) +
                               std::string(" select_1-value from file ")
                               + lls_filename);
    }
    const size_t size_of_data = (length + 63) / 64;
    data_vector.resize(size_of_data);
    if (std::fread(data_vector.data(), sizeof(uint64_t), size_of_data, in_fp)
         != size_of_data)
    {
      throw std::runtime_error(std::string("cannot read expected number of ") +
                               std::to_string(size_of_data) +
                               std::string(" data value from file ")
                               + lls_filename);
    }
    std::fclose(in_fp);
  }
  void push(bool value) noexcept
  {
    clear();
    if (length % 64 == 0)
    {
      data_vector.push_back(0);
    }
    set(length, value);
    length++;
  }

  void push_false_n(size_t n) noexcept
  {
    clear();
    length += n;
    data_vector.resize((length + 63) / 64);
  }

  void set(size_t index, bool value) noexcept
  {
    clear();
    if (value)
    {
      data_vector[index / 64] |= (uint64_t)1 << (index % 64);
    } else
    {
      data_vector[index / 64] &= ~((uint64_t)1 << (index % 64));
    }
  }

  [[nodiscard]] bool get(size_t index) const noexcept
  {
    return (data_vector[index / 64] >> (index % 64)) & uint64_t(1);
  }

  void buildAccelerationStructures(void) noexcept
  {
    clear();
    size_t global_count = 0;
    size_t index_1 = 0;
    for (size_t outer = 0; outer < (data_vector.size() + 63) / 64; outer++)
    {
      uint64_t first = 0;
      uint64_t second = 0;
      size_t local_count = 0;
      for (size_t inner = 0; inner < 8; inner++)
      {
        for (size_t i = 0; i < 8; i++)
        {
          size_t pop;
          if (i + inner * 8 + outer * 64 >= data_vector.size())
          {
            pop = 0;
          } else
          {
            pop = std::popcount(data_vector[i + inner * 8 + outer * 64]);
            index_1 += pop;

            if (index_1 >= INDEXINTERVAL)
            {
              const size_t rem = index_1 - INDEXINTERVAL;
              index_1 -= INDEXINTERVAL;
              const size_t index = (i + inner * 8 + outer * 64) * 64 +
                                    select_uint64_t(data_vector[i + inner * 8 +
                                                         outer * 64],
                                                    rem, true);
              select_1.push_back(index);
            }
          }
          local_count += pop;
        }
        if (inner < 5)
        {
          first |= local_count << (inner * 12);
        } else
        {
          if (inner == 5)
          {
            first |= local_count << (inner * 12);
            second |= local_count >> 4;
          } else
          {
            if (inner == 6)
            {
              second |= local_count << 8;
            }
          }
        }
      }

      global_count += local_count;
      second |= global_count << 20;
      rank.push_back({first, second});
    }
    select_1.push_back(length);
  }

  [[nodiscard]] size_t get_rank(size_t index, bool value) const
  {
    assert (not rank.empty());
    if (not value)
    {
      return index - get_rank(index, true);
    }

    size_t r = 0;

    const size_t superblock_index = index / SUPERBLOCKBITS;

    if (superblock_index > 0)
    {
      r += get_superblock_count(superblock_index - 1);
    }
    const size_t block_index = (index % SUPERBLOCKBITS) / BLOCKBITS;

    r += get_block_count(superblock_index, block_index);
    for (size_t i = superblock_index * 64 + block_index * 8; i < index / 64;
         i++)
    {
      int pop;
      if (i >= data_vector.size())
      {
        pop = 0;
      } else {
        pop = std::popcount(data_vector[i]);
      }
      r += pop;
    }

    if (index % 64 != 0)
    {
       r += std::popcount((data_vector[index / 64] << (64 - (index % 64))));
    }
    return r;
  }

  [[nodiscard]] size_t get_select(size_t count, bool value) const noexcept
  {
    if (not value)
    {
      return 0;
    }

    assert (not select_1.empty());
    const size_t select_index = count / INDEXINTERVAL;

    const size_t min = select_index > 0 ? select_1[select_index - 1] : 0;
    const size_t max = select_index < select_1.size() ? select_1[select_index]
                                                      : (length - 1);
    size_t superblock_l = min / SUPERBLOCKBITS;
    size_t superblock_r = (max + SUPERBLOCKBITS - 1) / SUPERBLOCKBITS;

    while (superblock_l < superblock_r)
    {
      const size_t m = (superblock_l + superblock_r) / 2;
      const size_t superblock_count = get_superblock_count(m);
      if (superblock_count >= count)
      {
        superblock_r = m;
      } else
      {
        superblock_l = m + 1;
      }
    }

    size_t block_l = 0;
    size_t block_r = 7;


    size_t superblock_count;
    if (superblock_l > 0)
    {
      superblock_count = get_superblock_count(superblock_l - 1);
    } else
    {
      superblock_count = 0;
    }

    while (block_l < block_r)
    {
      const size_t m = (block_l + block_r) / 2;

      if (get_block_count(superblock_l, m + 1) + superblock_count >= count)
      {
        block_r = m;
      } else
      {
        block_l = m + 1;
      }
    }

    const size_t block_count = get_block_count(superblock_l, block_l);
    size_t local_count = 0;
    for (int i = 0; i < 8; i++)
    {
      const size_t index = superblock_l * 64 + block_l * 8 + i;
      const size_t pop
       = index >= data_vector.size()
           ? 0 : std::popcount(data_vector[index]);
      local_count += pop;
      if (superblock_count + block_count + local_count >= count)
      {
        const size_t rem = count - (superblock_count + block_count +
                                    local_count - pop);

        return index * 64 + select_uint64_t(data_vector[index], rem, true);
      }
    }
    assert(false);
    return 0;
  }

  void serialize(const std::string &lls_filename) const
  {
    FILE *const out_fp = std::fopen(lls_filename.c_str(), "wb");

    if (out_fp == nullptr)
    {
      throw std::runtime_error(std::string("cannot write file ") +
                               lls_filename);
    }
    // File format:
    // 8 byte - (length) number of bits
    // 8 byte - number of superblocks
    // (number of superblocks * 16) byte - superblocks
    // 8 byte - length of select1 index
    // (length of selec1 * 8) byte - select1
    // (ceil(length / 64) * 8) byte - data_vector

    std::fwrite(&length, 8, 1, out_fp);

    const size_t number_of_superblocks = rank.size();
    std::fwrite(&number_of_superblocks, sizeof(size_t), 1, out_fp);
    std::fwrite(rank.data(), sizeof(Superblock), number_of_superblocks, out_fp);
    const size_t length_of_select1 = select_1.size();
    std::fwrite(&length_of_select1, sizeof(size_t), 1, out_fp);
    std::fwrite(select_1.data(), sizeof(size_t), length_of_select1, out_fp);
    std::fwrite(data_vector.data(), sizeof(uint64_t), data_vector.size(),
                out_fp);
    std::fclose(out_fp);
  }

  [[nodiscard]] size_t length_get() const noexcept
  {
    return length;
  }

  private:

  void clear(void) noexcept
  {
    rank.clear();
    select_1.clear();
  }

  [[nodiscard]]
  size_t get_superblock_count(size_t superblock_index) const noexcept
  {
    assert(superblock_index < rank.size());
    return rank[superblock_index].second >> 20;
  }

  [[nodiscard]]
  size_t get_block_count(size_t superblock_index,
                         size_t block_index) const noexcept
  {
    assert(superblock_index < rank.size());
    if (block_index == 7)
    {
      return (rank[superblock_index].second >> 8) & 0xFFF;
    }
    if (block_index == 6)
    {
      size_t t = (rank[superblock_index].second & 0xFF) << 4;
      t |= (rank[superblock_index].first >> 60) & 0xF;
      return t;
    }
    if (block_index > 0)
    {
      return (rank[superblock_index].first >> ((block_index - 1) * 12)) & 0xFFF;
    }
    return 0;
  }

  size_t select_uint64_t(uint64_t d, size_t count, bool value) const noexcept
  {
    size_t l = 0;
    size_t r = 63;
    while (l < r)
    {
      const size_t m = (l + r) / 2;
      assert(m <= 63);
      const size_t shift = 63 - m;
      const uint64_t shifted = d << static_cast<int>(shift);
      const size_t pop = value ? std::popcount(shifted)
                               : (std::popcount(shifted) - shift);
      if (pop < count)
      {
        l = m + 1;
      } else
      {
        r = m;
      }
    }
    return l;
  }
};
#endif
