#ifndef GCC_BUILTIN_HPP
#define GCC_BUILTIN_HPP

#ifdef __GNUC__
#define GTTL_IS_LIKELY(X) __builtin_expect((X),1)
#define GTTL_IS_UNLIKELY(X) __builtin_expect((X),0)
#else
#define GTTL_IS_LIKELY(X) (X)
#define GTTL_IS_UNLIKELY(X) (X)
#endif

#endif
