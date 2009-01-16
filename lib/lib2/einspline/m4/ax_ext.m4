dnl Copyright © 2007 Christophe Tournayre <turn3r@users.sourceforge.net>
dnl Copying and distribution of this file, with or without modification, 
dnl are permitted in any medium without royalty provided the copyright 
dnl notice and this notice are preserved.


AC_DEFUN([AX_EXT],
[
  AC_REQUIRE([AX_GCC_X86_CPUID])

  AX_GCC_X86_CPUID(0x00000001)
  ecx=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 3`
  edx=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 4`

 AC_CACHE_CHECK([whether mmx is supported], [ax_have_mmx_ext],
  [
    ax_have_mmx_ext=no
    if test "$((0x$edx>>23&0x01))" = 1; then
      ax_have_mmx_ext=yes
    fi
  ])

 AC_CACHE_CHECK([whether sse is supported], [ax_have_sse_ext],
  [
    ax_have_sse_ext=no
    if test "$((0x$edx>>25&0x01))" = 1; then
      ax_have_sse_ext=yes
    fi
  ])

 AC_CACHE_CHECK([whether sse2 is supported], [ax_have_sse2_ext],
  [
    ax_have_sse2_ext=no
    if test "$((0x$edx>>26&0x01))" = 1; then
      ax_have_sse2_ext=yes
    fi
  ])

 AC_CACHE_CHECK([whether sse3 is supported], [ax_have_sse3_ext],
  [
    ax_have_sse3_ext=no
    if test "$((0x$ecx&0x01))" = 1; then
      ax_have_sse3_ext=yes
    fi
  ])

 AC_CACHE_CHECK([whether ssse3 is supported], [ax_have_ssse3_ext],
  [
    ax_have_ssse3_ext=no
    if test "$((0x$ecx>>9&0x01))" = 1; then
      ax_have_ssse3_ext=yes
    fi
  ])

 AC_CACHE_CHECK([whether sse4.1 is supported], [ax_have_sse4_1_ext],
  [
    ax_have_sse4_1_ext=no
    if test "$((0x$ecx>>19&0x01))" = 1; then
      ax_have_sse4_1_ext=yes
    fi
  ])

 AC_CACHE_CHECK([whether sse4.2 is supported], [ax_have_sse4_2_ext],
  [
    ax_have_sse4_2_ext=no
    if test "$((0x$ecx>>20&0x01))" = 1; then
      ax_have_sse4_2_ext=yes
    fi
  ])

  if test "$ax_have_mmx_ext" = yes; then
    AC_DEFINE(HAVE_MMX,,[Support mmx instructions])
    AX_CHECK_COMPILER_FLAGS(-mmmx, SIMD_FLAGS="$SIMD_FLAGS -mmmx", [])
  fi

  if test "$ax_have_sse_ext" = yes; then
    AC_DEFINE(HAVE_SSE,,[Support SSE (Streaming SIMD Extensions) instructions])
    AX_CHECK_COMPILER_FLAGS(-msse, SIMD_FLAGS="$SIMD_FLAGS -msse", [])
  fi

  if test "$ax_have_sse2_ext" = yes; then
    AC_DEFINE(HAVE_SSE2,,[Support SSE2 (Streaming SIMD Extensions 2) instructions])
    AX_CHECK_COMPILER_FLAGS(-msse2, SIMD_FLAGS="$SIMD_FLAGS -msse2", [])
  fi

  if test "$ax_have_sse3_ext" = yes; then
    AC_DEFINE(HAVE_SSE3,,[Support SSE3 (Streaming SIMD Extensions 3) instructions])
    AX_CHECK_COMPILER_FLAGS(-msse3, SIMD_FLAGS="$SIMD_FLAGS -msse3", [])
  fi

  if test "$ax_have_ssse3_ext" = yes; then
    AC_DEFINE(HAVE_SSSE3,,[Support SSSE3 (Supplemental Streaming SIMD Extensions 3) instructions])
  fi

  if test "$ax_have_sse4_1_ext" = yes; then
    AC_DEFINE(HAVE_SSE4_1,,[Support SSE4.1 (Streaming SIMD Extensions 4.1) instructions])
  fi

  if test "$ax_have_sse4_2_ext" = yes; then
    AC_DEFINE(HAVE_SSE4_2,,[Support SSE4.2 (Streaming SIMD Extensions 4.2) instructions])
  fi

  AC_SUBST(SIMD_FLAGS)
])
