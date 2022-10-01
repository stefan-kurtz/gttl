#ifndef UNUSED_HPP
#define UNUSED_HPP

/* Unused function arguments should be annotated with this macro to get rid of
   compiler warnings. */
#define GTTL_UNUSED \
        __attribute__ ((unused))

#ifdef NDEBUG
#define GTTL_DEBUG_USED \
        __attribute__ ((unused))
#else
#define GTTL_DEBUG_USED /* Nothing */
#endif

#endif
