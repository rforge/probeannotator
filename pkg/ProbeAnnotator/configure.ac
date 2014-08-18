##INIT CONFIG
AC_PREREQ(2.57)
AC_INIT([ProbeAnnotator], 1.0.1)

## CHECK PCRE
#AC_CHECK_HEADER([pcre.h], [HAVE_PCRE=" -DHAVE_PCRE_H" , HAVE_PCRE_BOOL=TRUE], [HAVE_PCRE="" , HAVE_PCRE_BOOL=FALSE])
AC_CHECK_HEADER([pcre.h], 
	[HAVE_PCRE=" -DHAVE_PCRE_H" HAVE_PCRE_BOOL="TRUE"], 
	[HAVE_PCRE="" HAVE_PCRE_BOOL="FALSE"])
## since no warning while compiling C with r-forge, and 'pcre_get_substring' not found. -> check 'pcre_get_substring'
AC_CHECK_LIB(pcre, pcre_get_substring)

