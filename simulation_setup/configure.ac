AC_INIT([kh], [0.1], [nicco.francini@libero.it])


AM_INIT_AUTOMAKE(subdir-objects)


AC_CONFIG_HEADERS([config.h])


AC_PROG_CC


AC_CHECK_HEADERS([math.h], [], [AC_MSG_ERROR([Required header math.h not found.])])

AC_CHECK_HEADERS([stdio.h], [], [AC_MSG_ERROR([Required header stdio.h not found.])])

AC_CHECK_HEADERS([stdlib.h], [], [AC_MSG_ERROR([Required header stdlib.h not found.])])

AC_CHECK_HEADERS([time.h], [], [AC_MSG_ERROR([Required header time.h not found.])])

AC_CHECK_HEADERS([stddef.h], [], [AC_MSG_ERROR([Required header stddef.h not found.])])

AC_CHECK_HEADERS([complex.h], [], [AC_MSG_ERROR([Required header complex.h not found.])])

AC_CHECK_HEADERS([errno.h], [], [AC_MSG_ERROR([Required header errno.h not found.])])

AC_CHECK_HEADERS([fenv.h], [], [AC_MSG_ERROR([Required header fenv.h not found.])])

AC_CHECK_HEADERS([sys/types.h], [], [AC_MSG_ERROR([Required header sys/types.h not found.])])

AC_CHECK_HEADERS([sys/stat.h], [], [AC_MSG_ERROR([Required header sys/stat.h not found.])])


AC_SEARCH_LIBS([sqrt], [m], [], [AC_MSG_ERROR([Required c mathematical library not found.])])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
