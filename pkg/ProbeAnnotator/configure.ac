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

## EDIT Makevars
AC_SUBST(HAVE_PCRE)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
## EDIT *.R files
AC_SUBST(HAVE_PCRE_BOOL)
AC_CONFIG_FILES([R/probeAnnot.R])
AC_CONFIG_FILES([R/transform_UCSC.R])
AC_OUTPUT
