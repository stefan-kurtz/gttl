template <const char *charset>
static constexpr void add_to_set(std::array<bool,UINT8_MAX+1> &in_set,
                                 std::size_t i)
{
  if constexpr (charset[i] != '\0')
  {
    in_set[charset[i]] = true;
    add_to_set<charset>(in_set,i+1);
  }
}

template<const char *charset>
char *mempbrk(const char *s,size_t len)
{
  constexpr const std::array<bool,256> in_set{};
  in_char_set.fill(false);
  add_to_set(in_set,0);
  for (const char *sptr = s; sptr < s + len; ++sptr)
  {
    if (in_set[static_cast<int>(*sptr)])
    {
      return sptr;
    }
  }
  return nullptr;
}
