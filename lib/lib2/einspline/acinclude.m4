AC_DEFUN([AX_CC_OPTION], [
AC_REQUIRE([AC_PROG_CC])
AC_MSG_CHECKING([if ${CC-cc} accepts $2 option])
echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} $2 -c conftest.c 2>&1`"; then
      $1=$3
      AC_MSG_RESULT([yes])
else
      $1=$4
      AC_MSG_RESULT([no])
fi
rm -f conftest*
])

AC_DEFUN([AX_F77_OPTION], [
AC_REQUIRE([AC_PROG_F77])
AC_MSG_CHECKING([if ${F77-f77} accepts $2 option])
echo 'void f(){}' > conftest.c
if test -z "`${F77-f77} $2 -c conftest.c 2>&1`"; then
      $1=$3
      AC_MSG_RESULT([yes])
else
      $1=$4
      AC_MSG_RESULT([no])
fi
rm -f conftest*
])

m4_include([m4/acx_pthread.m4])
m4_include([m4/ax_cc_maxopt.m4])
m4_include([m4/ax_cxx_maxopt.m4])
m4_include([m4/ax_f77_maxopt.m4])
m4_include([m4/ax_check_compiler_flags.m4])
m4_include([m4/ax_compiler_vendor.m4])
m4_include([m4/ax_cxx_compiler_vendor.m4])
m4_include([m4/ax_c_compiler_vendor.m4])
m4_include([m4/ax_f77_compiler_vendor.m4])
m4_include([m4/ax_gcc_aligns_stack.m4])
m4_include([m4/ax_gcc_archflag.m4])
m4_include([m4/ax_gxx_archflag.m4])
m4_include([m4/ax_gcc_version.m4])
m4_include([m4/ax_gcc_x86_cpuid.m4])
m4_include([m4/ax_ext.m4])
m4_include([m4/ac_cxx_restrict.m4])
