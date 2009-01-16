#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.61)
AC_INIT(einspline, 0.8.9, kesler@ciw.edu, einspline)
AM_INIT_AUTOMAKE(einspline, 0.8.9)
AC_CONFIG_SRCDIR([src/bspline.h])
AC_CONFIG_HEADER([src/config.h])


# Checks for programs.
AC_PROG_CC_STDC
AC_PROG_CXX
AC_PROG_LIBTOOL
PKG_PROG_PKG_CONFIG

# Optimal compiler flags
AX_CC_MAXOPT
AX_CXX_MAXOPT
AX_F77_MAXOPT
AC_CXX_RESTRICT



AC_CHECK_LIB([m], [sqrt])
AC_MSG_CHECKING([LIBS = $LIBS])
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])
AC_C_VARARRAYS

AC_MSG_CHECKING([optimal compiler flags])
CFLAGS="${CFLAGS} -g"
CXXFLAGS="${CXXFLAGS} -g"
AC_MSG_RESULT([${CFLAGS}])

# Compile CUDA versions
AC_ARG_ENABLE(cuda, AC_HELP_STRING([--enable-cuda], [compile CUDA routines]))
AM_CONDITIONAL([HAVE_CUDA], [test "$enable_cuda" = "yes"])
if test "x$enable_cuda" = "xyes"; then
   NVCC="nvcc"
   NVCCFLAGS="-DNO_CUDA_MAIN -g -O3"
   CUDA_LIBS="-L/usr/local/cuda/lib -lcuda -lcudart"
   LIBS="$LIBS -L/usr/local/cuda/lib -lcuda -lcudart"
fi
########################
# Check for CUDA
########################
AC_ARG_WITH([cuda],
    [  --with-cuda=PATH        prefix where cuda is installed [default=auto]])
if test "x$with_cuda" != x; then
  CUDA_CFLAGS="-I$with_cuda/include"
  CUDA_LIBS="-L$with_cuda/lib -lcuda -lcudart"
  NVCC="$with_cuda/bin/nvcc"
else
  CUDA_CFLAGS="-I/usr/local/cuda/include"
  CUDA_LIBS="-L/usr/local/cuda/lib -lcuda -lcudart"
  NVCC="nvcc"
fi
NVCCFLAGS="-DNO_CUDA_MAIN -O3 -use_fast_math -g"
AC_SUBST(CUDA_CFLAGS)
AC_SUBST(CUDA_LIBS)
AC_SUBST(NVCC)
AC_SUBST(NVCCFLAGS)



AC_SUBST(NVCC)
AC_SUBST(NVCCFLAGS)
AC_SUBST(CUDA_LIBS)

# Compile with -pthread
AC_ARG_ENABLE(pthread, AC_HELP_STRING([--enable-pthread], [compile with -pthread]))
if test "x$enable_pthread" = "xyes"; then
   CFLAGS="$CFLAGS -pthread"
   CXXFLAGS="$CXXFLAGS -pthread"
   PTHREAD_FLAG="-pthread"
else
   PTHREAD_FLAG=""
fi
AC_SUBST(PTHREAD_FLAG)

# Compile with -openmp
AC_ARG_ENABLE(openmp, AC_HELP_STRING([--enable-openmp], [compile with -openmp]))
if test "x$enable_openmp" = "xyes"; then
   CFLAGS="$CFLAGS -openmp"
   CXXFLAGS="$CXXFLAGS -openmp"
   OPENMP_FLAG="-openmp"
else
   OPENMP_FLAG=""
fi
AC_SUBST(OPENMP_FLAG)



# Check if we want static-only binaries
AC_ARG_ENABLE(prefetch,   AC_HELP_STRING([--enable-prefetch[=N]],[use software prefetch instructions (default=no)]))
if test "$enable_prefetch" = "yes"; then
   AC_DEFINE([USE_PREFETCH],[],[Use SSE prefetch])
   AC_DEFINE([PREFETCH_AHEAD],[12],[Prefetch loop lead distance])
elif test "x$enable_prefetch" != "x"; then
   AC_DEFINE([USE_PREFETCH],[],[Use SSE prefetch])
   AC_DEFINE_UNQUOTED([PREFETCH_AHEAD],$enable_prefetch,[Prefetch loop lead distance])
fi

AC_ARG_ENABLE(all-static, AC_HELP_STRING([--enable-all-static], [build static binaries (default=no)]))
if test "$enable_all_static" = "yes"; then
  ALL_STATIC="-all-static"
else
  ALL_STATIC=""
fi
AC_SUBST(ALL_STATIC)

AC_ARG_ENABLE(altivec, AC_HELP_STRING([--enable-altivec],[enable PowerPC SIMD extensions (default=no)]))
if test "x$enable_altivec" = "xyes"; then
   CFLAGS="$CFLAGS -faltivec -DUSE_ALTIVEC"
   FFLAGS="$FFLAGS -DUSE_ALTIVEC"
fi


AC_ARG_ENABLE(sse, AC_HELP_STRING([--enable-sse],[enable SSE SIMD instructions]))
if test "x$enable_sse" = "xyes"; then
  AC_CHECK_HEADER([xmmintrin.h],
                  [AC_DEFINE([HAVE_XMMINTRIN_H],[1],
                   [xmmintrin.h is available.])],
                  [AC_MSG_ERROR([Cannot utilize SSE2 extensions without xmmintrin.h])])
  AC_DEFINE([HAVE_SSE], [1], [SSE instructions are available])
fi

if test "x$enable_sse" = "xyes"; then
   AX_EXT
fi

AM_CONDITIONAL([HAVE_SSE],    [test "$ax_have_sse_ext"    = "yes"])
AM_CONDITIONAL([HAVE_SSE2],   [test "$ax_have_sse2_ext"   = "yes"])
AM_CONDITIONAL([HAVE_SSE3],   [test "$ax_have_sse3_ext"   = "yes"])
AM_CONDITIONAL([HAVE_SSSE3],  [test "$ax_have_ssse3_ext"  = "yes"])
AM_CONDITIONAL([HAVE_SSE4_1], [test "$ax_have_sse4_1_ext" = "yes"])
AM_CONDITIONAL([HAVE_SSE4_2], [test "$ax_have_sse4_2_ext" = "yes"])

AC_ARG_ENABLE(blips, AC_HELP_STRING([--enable-blips],[enable routines for creating BLIPS (default=no)]))
AM_CONDITIONAL([WANT_BLIPS], [test "$enable_blips" = "yes"])
if test "x$enable_blips" = "xyes"; then
  PKG_CHECK_MODULES(FFTW3,  fftw3,  fftw3_ok=yes,  fftw3_ok=no)
  PKG_CHECK_MODULES(FFTW3F, fftw3f, fftw3f_ok=yes, fftw3f_ok=no)
  if test "x$fftw3_ok" = "xyes"; then
    AC_DEFINE([HAVE_FFTW3], [1], [FFTW3 library is available])
  else
    AC_MSG_ERROR([Blip construction routines require that FFTW3 be installed.])
  fi
  if test "x$fftw3f_ok" = "xyes"; then
    AC_DEFINE([HAVE_FFTW3F], [1], [FFTW3F library is available])
  else
    AC_MSG_ERROR([Blip construction routines require that FFTW3 be installed.])
  fi
fi

AC_ARG_ENABLE(precision, AC_HELP_STRING([--enable-precision],[use double-precision solve for single-precision splines]))
if test "x$enable_precision" = "xyes"; then
   AC_DEFINE([HIGH_PRECISION], [1], [Use double-precision to solve for single-precision splines])
fi

AC_ARG_ENABLE(profile, AC_HELP_STRING([--enable-profile],[instrument code with profiling information]))
if test "x$enable_profile" = "xyes"; then
   CFLAGS="-pg $CFLAGS"
   FFLAGS="-pg $FFLAGS"
   LDFLAGS="-pg $LDFLAGS"
fi

AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug],[enable code for debugging checks]))
if test "x$enable_debug" = "xyes"; then
   CFLAGS="-DDEBUG $CFLAGS"
   FFLAGS="-DDEBUG $FFLAGS"
fi


AC_F77_WRAPPERS
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
AC_CHECK_FUNCS([clock_gettime floor pow sqrt strtol posix_memalign])

AC_SUBST(PKGDATADEF)

AC_CONFIG_FILES([ Makefile                  \
		  src/Makefile              \
                  doc/Makefile              \
		  www/Makefile              \
                  einspline.pc  ])
AC_OUTPUT
