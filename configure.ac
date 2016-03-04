#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(netabc, 0.1, rmccloskey@cfenet.ubc.ca)
AC_CONFIG_SRCDIR([src/netabc.c])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects foreign 1.9.6])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AM_PROG_AR
AC_PROG_YACC
AC_PROG_LEX

# Doxygen.
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(netabc, Doxyfile, doc)

# Compiler and linker stuff
LT_INIT
AX_COMPILER_FLAGS()
AX_GCC_FUNC_ATTRIBUTE([noreturn])

# Checks for libraries.
AX_PTHREAD([], [AC_MSG_ERROR([pthtreads is required])])
AC_CHECK_LIB([m], [cos], [], [AC_MSG_ERROR([math.h is required])])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [AC_MSG_ERROR([GSL CBLAS is required])])
AC_CHECK_LIB([igraph], [igraph_empty], [], [AC_MSG_ERROR([igraph is required])])
AC_CHECK_LIB([Judy], [JudyLIns], [], [AC_MSG_ERROR([Judy is required])])
AC_CHECK_LIB([yaml], [yaml_parser_initialize], [], [AC_MSG_ERROR([libyaml is required])])
PKG_CHECK_MODULES([GSL], [gsl >= 1.15])
PKG_CHECK_MODULES([CHECK], [check >= 0.9.4])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h stddef.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC

# Output files
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
