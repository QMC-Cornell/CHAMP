dnl Check routine for "restrict" keyword which was introduced 
dnl in ANSI C99. Some C++ compiler like g++ or kcc does support
dnl the keyword inside C++ as well.
dnl Does nothing if the compiler accepts the keyword.  Otherwise, if 
dnl the compiler supports an equivalent, like gcc's __restrict__ or  
dnl SGI's __restrict, define "restrict" to be that.
dnl Otherwise, define "restrict" to be empty.

AC_DEFUN([AC_CXX_RESTRICT],
[AC_CACHE_CHECK([for restrict], ac_cxx_restrict,
[ac_cxx_restrict=no
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
  for ac_kw in restrict __restrict__ __restrict; do
   AC_TRY_COMPILE(, [void* $ac_kw bar], [ac_cxx_restrict=$ac_kw; break])
  done
 AC_LANG_RESTORE
])
if test "$ac_cxx_restrict" != "restrict"; then
 ac_kw="$ac_cxx_restrict"
 if test "$ac_kw" = unsupported; then ac_kw=""; fi
 AC_DEFINE_UNQUOTED(restrict, $ac_cxx_restrict, [Define to empty if the C99 keyword for C++ does not work.])
fi 
])
