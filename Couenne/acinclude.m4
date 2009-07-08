# libtool.m4 - Configure libtool for the host system. -*-Autoconf-*-
## Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2003, 2004, 2005
## Free Software Foundation, Inc.
## Originally by Gordon Matzigkeit <gord@gnu.ai.mit.edu>, 1996
##
## This file is free software; the Free Software Foundation gives
## unlimited permission to copy and/or distribute it, with or without
## modifications, as long as this notice is preserved.

# serial 48 AC_PROG_LIBTOOL


# AC_PROVIDE_IFELSE(MACRO-NAME, IF-PROVIDED, IF-NOT-PROVIDED)
# -----------------------------------------------------------
# If this macro is not defined by Autoconf, define it here.
m4_ifdef([AC_PROVIDE_IFELSE],
         [],
         [m4_define([AC_PROVIDE_IFELSE],
	         [m4_ifdef([AC_PROVIDE_$1],
		           [$2], [$3])])])


# AC_PROG_LIBTOOL
# ---------------
AC_DEFUN([AC_PROG_LIBTOOL],
[AC_REQUIRE([_AC_PROG_LIBTOOL])dnl
dnl If AC_PROG_CXX has already been expanded, run AC_LIBTOOL_CXX
dnl immediately, otherwise, hook it in at the end of AC_PROG_CXX.
  AC_PROVIDE_IFELSE([AC_PROG_CXX],
    [AC_LIBTOOL_CXX],
    [define([AC_PROG_CXX], defn([AC_PROG_CXX])[AC_LIBTOOL_CXX
  ])])
dnl And a similar setup for Fortran 77 support
  AC_PROVIDE_IFELSE([AC_PROG_F77],
    [AC_LIBTOOL_F77],
    [define([AC_PROG_F77], defn([AC_PROG_F77])[AC_LIBTOOL_F77
])])

dnl Quote A][M_PROG_GCJ so that aclocal doesn't bring it in needlessly.
dnl If either AC_PROG_GCJ or A][M_PROG_GCJ have already been expanded, run
dnl AC_LIBTOOL_GCJ immediately, otherwise, hook it in at the end of both.
  AC_PROVIDE_IFELSE([AC_PROG_GCJ],
    [AC_LIBTOOL_GCJ],
    [AC_PROVIDE_IFELSE([A][M_PROG_GCJ],
      [AC_LIBTOOL_GCJ],
      [AC_PROVIDE_IFELSE([LT_AC_PROG_GCJ],
	[AC_LIBTOOL_GCJ],
      [ifdef([AC_PROG_GCJ],
	     [define([AC_PROG_GCJ], defn([AC_PROG_GCJ])[AC_LIBTOOL_GCJ])])
       ifdef([A][M_PROG_GCJ],
	     [define([A][M_PROG_GCJ], defn([A][M_PROG_GCJ])[AC_LIBTOOL_GCJ])])
       ifdef([LT_AC_PROG_GCJ],
	     [define([LT_AC_PROG_GCJ],
		defn([LT_AC_PROG_GCJ])[AC_LIBTOOL_GCJ])])])])
])])# AC_PROG_LIBTOOL


# _AC_PROG_LIBTOOL
# ----------------
AC_DEFUN([_AC_PROG_LIBTOOL],
[AC_REQUIRE([AC_LIBTOOL_SETUP])dnl
AC_BEFORE([$0],[AC_LIBTOOL_CXX])dnl
AC_BEFORE([$0],[AC_LIBTOOL_F77])dnl
AC_BEFORE([$0],[AC_LIBTOOL_GCJ])dnl

# This can be used to rebuild libtool when needed
LIBTOOL_DEPS="$ac_aux_dir/ltmain.sh"

# Always use our own libtool.
LIBTOOL='$(SHELL) $(top_builddir)/libtool'
AC_SUBST(LIBTOOL)dnl

# Prevent multiple expansion
define([AC_PROG_LIBTOOL], [])
])# _AC_PROG_LIBTOOL


# AC_LIBTOOL_SETUP
# ----------------
AC_DEFUN([AC_LIBTOOL_SETUP],
[AC_PREREQ(2.50)dnl
AC_REQUIRE([AC_ENABLE_SHARED])dnl
AC_REQUIRE([AC_ENABLE_STATIC])dnl
AC_REQUIRE([AC_ENABLE_FAST_INSTALL])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_PROG_LD])dnl
AC_REQUIRE([AC_PROG_LD_RELOAD_FLAG])dnl
AC_REQUIRE([AC_PROG_NM])dnl

AC_REQUIRE([AC_PROG_LN_S])dnl
AC_REQUIRE([AC_DEPLIBS_CHECK_METHOD])dnl
# Autoconf 2.13's AC_OBJEXT and AC_EXEEXT macros only works for C compilers!
AC_REQUIRE([AC_OBJEXT])dnl
AC_REQUIRE([AC_EXEEXT])dnl
dnl

AC_LIBTOOL_SYS_MAX_CMD_LEN
AC_LIBTOOL_SYS_GLOBAL_SYMBOL_PIPE
AC_LIBTOOL_OBJDIR

AC_REQUIRE([_LT_AC_SYS_COMPILER])dnl
_LT_AC_PROG_ECHO_BACKSLASH

case $host_os in
aix3*)
  # AIX sometimes has problems with the GCC collect2 program.  For some
  # reason, if we set the COLLECT_NAMES environment variable, the problems
  # vanish in a puff of smoke.
  if test "X${COLLECT_NAMES+set}" != Xset; then
    COLLECT_NAMES=
    export COLLECT_NAMES
  fi
  ;;
esac

# Sed substitution that helps us do robust quoting.  It backslashifies
# metacharacters that are still active within double-quoted strings.
Xsed='sed -e 1s/^X//'
[sed_quote_subst='s/\([\\"\\`$\\\\]\)/\\\1/g']

# Same as above, but do not quote variable references.
[double_quote_subst='s/\([\\"\\`\\\\]\)/\\\1/g']

# Sed substitution to delay expansion of an escaped shell variable in a
# double_quote_subst'ed string.
delay_variable_subst='s/\\\\\\\\\\\$/\\\\\\$/g'

# Sed substitution to avoid accidental globbing in evaled expressions
no_glob_subst='s/\*/\\\*/g'

# Constants:
rm="rm -f"

# Global variables:
default_ofile=libtool
can_build_shared=yes

# All known linkers require a `.a' archive for static linking (except MSVC,
# which needs '.lib').
libext=a
ltmain="$ac_aux_dir/ltmain.sh"
ofile="$default_ofile"
with_gnu_ld="$lt_cv_prog_gnu_ld"

AC_CHECK_TOOL(AR, ar, false)
AC_CHECK_TOOL(RANLIB, ranlib, :)
AC_CHECK_TOOL(STRIP, strip, :)

old_CC="$CC"
old_CFLAGS="$CFLAGS"

# Set sane defaults for various variables
test -z "$AR" && AR=ar
test -z "$AR_FLAGS" && AR_FLAGS=cru
test -z "$AS" && AS=as
test -z "$CC" && CC=cc
test -z "$LTCC" && LTCC=$CC
test -z "$LTCFLAGS" && LTCFLAGS=$CFLAGS
test -z "$DLLTOOL" && DLLTOOL=dlltool
test -z "$LD" && LD=ld
test -z "$LN_S" && LN_S="ln -s"
test -z "$MAGIC_CMD" && MAGIC_CMD=file
test -z "$NM" && NM=nm
test -z "$SED" && SED=sed
test -z "$OBJDUMP" && OBJDUMP=objdump
test -z "$RANLIB" && RANLIB=:
test -z "$STRIP" && STRIP=:
test -z "$ac_objext" && ac_objext=o

# Determine commands to create old-style static archives.
old_archive_cmds='$AR $AR_FLAGS $oldlib$oldobjs$old_deplibs'
old_postinstall_cmds='chmod 644 $oldlib'
old_postuninstall_cmds=

if test -n "$RANLIB"; then
  case $host_os in
  openbsd*)
    old_postinstall_cmds="$old_postinstall_cmds~\$RANLIB -t \$oldlib"
    ;;
  *)
    old_postinstall_cmds="$old_postinstall_cmds~\$RANLIB \$oldlib"
    ;;
  esac
  old_archive_cmds="$old_archive_cmds~\$RANLIB \$oldlib"
fi

_LT_CC_BASENAME([$compiler])

# Only perform the check for file, if the check method requires it
case $deplibs_check_method in
file_magic*)
  if test "$file_magic_cmd" = '$MAGIC_CMD'; then
    AC_PATH_MAGIC
  fi
  ;;
esac

AC_PROVIDE_IFELSE([AC_LIBTOOL_DLOPEN], enable_dlopen=yes, enable_dlopen=no)
AC_PROVIDE_IFELSE([AC_LIBTOOL_WIN32_DLL],
enable_win32_dll=yes, enable_win32_dll=no)

AC_ARG_ENABLE([libtool-lock],
    [AC_HELP_STRING([--disable-libtool-lock],
	[avoid locking (might break parallel builds)])])
test "x$enable_libtool_lock" != xno && enable_libtool_lock=yes

AC_ARG_WITH([pic],
    [AC_HELP_STRING([--with-pic],
	[try to use only PIC/non-PIC objects @<:@default=use both@:>@])],
    [pic_mode="$withval"],
    [pic_mode=default])
test -z "$pic_mode" && pic_mode=default

# Use C for the default configuration in the libtool script
tagname=
AC_LIBTOOL_LANG_C_CONFIG
_LT_AC_TAGCONFIG
])# AC_LIBTOOL_SETUP


# _LT_AC_SYS_COMPILER
# -------------------
AC_DEFUN([_LT_AC_SYS_COMPILER],
[AC_REQUIRE([AC_PROG_CC])dnl

# If no C compiler was specified, use CC.
LTCC=${LTCC-"$CC"}

# If no C compiler flags were specified, use CFLAGS.
LTCFLAGS=${LTCFLAGS-"$CFLAGS"}

# Allow CC to be a program name with arguments.
compiler=$CC
])# _LT_AC_SYS_COMPILER


# _LT_CC_BASENAME(CC)
# -------------------
# Calculate cc_basename.  Skip known compiler wrappers and cross-prefix.
AC_DEFUN([_LT_CC_BASENAME],
[for cc_temp in $1""; do
  case $cc_temp in
    compile | *[[\\/]]compile | ccache | *[[\\/]]ccache ) ;;
    distcc | *[[\\/]]distcc | purify | *[[\\/]]purify ) ;;
    \-*) ;;
    *) break;;
  esac
done
cc_basename=`$echo "X$cc_temp" | $Xsed -e 's%.*/%%' -e "s%^$host_alias-%%"`
])


# _LT_COMPILER_BOILERPLATE
# ------------------------
# Check for compiler boilerplate output or warnings with
# the simple compiler test code.
AC_DEFUN([_LT_COMPILER_BOILERPLATE],
[ac_outfile=conftest.$ac_objext
printf "$lt_simple_compile_test_code" >conftest.$ac_ext
eval "$ac_compile" 2>&1 >/dev/null | $SED '/^$/d; /^ *+/d' >conftest.err
_lt_compiler_boilerplate=`cat conftest.err`
$rm conftest*
])# _LT_COMPILER_BOILERPLATE


# _LT_LINKER_BOILERPLATE
# ----------------------
# Check for linker boilerplate output or warnings with
# the simple link test code.
AC_DEFUN([_LT_LINKER_BOILERPLATE],
[ac_outfile=conftest.$ac_objext
printf "$lt_simple_link_test_code" >conftest.$ac_ext
eval "$ac_link" 2>&1 >/dev/null | $SED '/^$/d; /^ *+/d' >conftest.err
_lt_linker_boilerplate=`cat conftest.err`
$rm conftest*
])# _LT_LINKER_BOILERPLATE


# _LT_AC_SYS_LIBPATH_AIX
# ----------------------
# Links a minimal program and checks the executable
# for the system default hardcoded library path. In most cases,
# this is /usr/lib:/lib, but when the MPI compilers are used
# the location of the communication and MPI libs are included too.
# If we don't find anything, use the default library path according
# to the aix ld manual.
AC_DEFUN([_LT_AC_SYS_LIBPATH_AIX],
[AC_LINK_IFELSE(AC_LANG_PROGRAM,[
aix_libpath=`dump -H conftest$ac_exeext 2>/dev/null | $SED -n -e '/Import File Strings/,/^$/ { /^0/ { s/^0  *\(.*\)$/\1/; p; }
}'`
# Check for a 64-bit object if we didn't find anything.
if test -z "$aix_libpath"; then aix_libpath=`dump -HX64 conftest$ac_exeext 2>/dev/null | $SED -n -e '/Import File Strings/,/^$/ { /^0/ { s/^0  *\(.*\)$/\1/; p; }
}'`; fi],[])
if test -z "$aix_libpath"; then aix_libpath="/usr/lib:/lib"; fi
])# _LT_AC_SYS_LIBPATH_AIX


# _LT_AC_SHELL_INIT(ARG)
# ----------------------
AC_DEFUN([_LT_AC_SHELL_INIT],
[ifdef([AC_DIVERSION_NOTICE],
	     [AC_DIVERT_PUSH(AC_DIVERSION_NOTICE)],
	 [AC_DIVERT_PUSH(NOTICE)])
$1
AC_DIVERT_POP
])# _LT_AC_SHELL_INIT


# _LT_AC_PROG_ECHO_BACKSLASH
# --------------------------
# Add some code to the start of the generated configure script which
# will find an echo command which doesn't interpret backslashes.
AC_DEFUN([_LT_AC_PROG_ECHO_BACKSLASH],
[_LT_AC_SHELL_INIT([
# Check that we are running under the correct shell.
SHELL=${CONFIG_SHELL-/bin/sh}

case X$ECHO in
X*--fallback-echo)
  # Remove one level of quotation (which was required for Make).
  ECHO=`echo "$ECHO" | sed 's,\\\\\[$]\\[$]0,'[$]0','`
  ;;
esac

echo=${ECHO-echo}
if test "X[$]1" = X--no-reexec; then
  # Discard the --no-reexec flag, and continue.
  shift
elif test "X[$]1" = X--fallback-echo; then
  # Avoid inline document here, it may be left over
  :
elif test "X`($echo '\t') 2>/dev/null`" = 'X\t' ; then
  # Yippee, $echo works!
  :
else
  # Restart under the correct shell.
  exec $SHELL "[$]0" --no-reexec ${1+"[$]@"}
fi

if test "X[$]1" = X--fallback-echo; then
  # used as fallback echo
  shift
  cat <<EOF
[$]*
EOF
  exit 0
fi

# The HP-UX ksh and POSIX shell print the target directory to stdout
# if CDPATH is set.
(unset CDPATH) >/dev/null 2>&1 && unset CDPATH

if test -z "$ECHO"; then
if test "X${echo_test_string+set}" != Xset; then
# find a string as large as possible, as long as the shell can cope with it
  for cmd in 'sed 50q "[$]0"' 'sed 20q "[$]0"' 'sed 10q "[$]0"' 'sed 2q "[$]0"' 'echo test'; do
    # expected sizes: less than 2Kb, 1Kb, 512 bytes, 16 bytes, ...
    if (echo_test_string=`eval $cmd`) 2>/dev/null &&
       echo_test_string=`eval $cmd` &&
       (test "X$echo_test_string" = "X$echo_test_string") 2>/dev/null
    then
      break
    fi
  done
fi

if test "X`($echo '\t') 2>/dev/null`" = 'X\t' &&
   echo_testing_string=`($echo "$echo_test_string") 2>/dev/null` &&
   test "X$echo_testing_string" = "X$echo_test_string"; then
  :
else
  # The Solaris, AIX, and Digital Unix default echo programs unquote
  # backslashes.  This makes it impossible to quote backslashes using
  #   echo "$something" | sed 's/\\/\\\\/g'
  #
  # So, first we look for a working echo in the user's PATH.

  lt_save_ifs="$IFS"; IFS=$PATH_SEPARATOR
  for dir in $PATH /usr/ucb; do
    IFS="$lt_save_ifs"
    if (test -f $dir/echo || test -f $dir/echo$ac_exeext) &&
       test "X`($dir/echo '\t') 2>/dev/null`" = 'X\t' &&
       echo_testing_string=`($dir/echo "$echo_test_string") 2>/dev/null` &&
       test "X$echo_testing_string" = "X$echo_test_string"; then
      echo="$dir/echo"
      break
    fi
  done
  IFS="$lt_save_ifs"

  if test "X$echo" = Xecho; then
    # We didn't find a better echo, so look for alternatives.
    if test "X`(print -r '\t') 2>/dev/null`" = 'X\t' &&
       echo_testing_string=`(print -r "$echo_test_string") 2>/dev/null` &&
       test "X$echo_testing_string" = "X$echo_test_string"; then
      # This shell has a builtin print -r that does the trick.
      echo='print -r'
    elif (test -f /bin/ksh || test -f /bin/ksh$ac_exeext) &&
	 test "X$CONFIG_SHELL" != X/bin/ksh; then
      # If we have ksh, try running configure again with it.
      ORIGINAL_CONFIG_SHELL=${CONFIG_SHELL-/bin/sh}
      export ORIGINAL_CONFIG_SHELL
      CONFIG_SHELL=/bin/ksh
      export CONFIG_SHELL
      exec $CONFIG_SHELL "[$]0" --no-reexec ${1+"[$]@"}
    else
      # Try using printf.
      echo='printf %s\n'
      if test "X`($echo '\t') 2>/dev/null`" = 'X\t' &&
	 echo_testing_string=`($echo "$echo_test_string") 2>/dev/null` &&
	 test "X$echo_testing_string" = "X$echo_test_string"; then
	# Cool, printf works
	:
      elif echo_testing_string=`($ORIGINAL_CONFIG_SHELL "[$]0" --fallback-echo '\t') 2>/dev/null` &&
	   test "X$echo_testing_string" = 'X\t' &&
	   echo_testing_string=`($ORIGINAL_CONFIG_SHELL "[$]0" --fallback-echo "$echo_test_string") 2>/dev/null` &&
	   test "X$echo_testing_string" = "X$echo_test_string"; then
	CONFIG_SHELL=$ORIGINAL_CONFIG_SHELL
	export CONFIG_SHELL
	SHELL="$CONFIG_SHELL"
	export SHELL
	echo="$CONFIG_SHELL [$]0 --fallback-echo"
      elif echo_testing_string=`($CONFIG_SHELL "[$]0" --fallback-echo '\t') 2>/dev/null` &&
	   test "X$echo_testing_string" = 'X\t' &&
	   echo_testing_string=`($CONFIG_SHELL "[$]0" --fallback-echo "$echo_test_string") 2>/dev/null` &&
	   test "X$echo_testing_string" = "X$echo_test_string"; then
	echo="$CONFIG_SHELL [$]0 --fallback-echo"
      else
	# maybe with a smaller string...
	prev=:

	for cmd in 'echo test' 'sed 2q "[$]0"' 'sed 10q "[$]0"' 'sed 20q "[$]0"' 'sed 50q "[$]0"'; do
	  if (test "X$echo_test_string" = "X`eval $cmd`") 2>/dev/null
	  then
	    break
	  fi
	  prev="$cmd"
	done

	if test "$prev" != 'sed 50q "[$]0"'; then
	  echo_test_string=`eval $prev`
	  export echo_test_string
	  exec ${ORIGINAL_CONFIG_SHELL-${CONFIG_SHELL-/bin/sh}} "[$]0" ${1+"[$]@"}
	else
	  # Oops.  We lost completely, so just stick with echo.
	  echo=echo
	fi
      fi
    fi
  fi
fi
fi

# Copy echo and quote the copy suitably for passing to libtool from
# the Makefile, instead of quoting the original, which is used later.
ECHO=$echo
if test "X$ECHO" = "X$CONFIG_SHELL [$]0 --fallback-echo"; then
   ECHO="$CONFIG_SHELL \\\$\[$]0 --fallback-echo"
fi

AC_SUBST(ECHO)
])])# _LT_AC_PROG_ECHO_BACKSLASH


# _LT_AC_LOCK
# -----------
AC_DEFUN([_LT_AC_LOCK],
[AC_ARG_ENABLE([libtool-lock],
    [AC_HELP_STRING([--disable-libtool-lock],
	[avoid locking (might break parallel builds)])])
test "x$enable_libtool_lock" != xno && enable_libtool_lock=yes

# Some flags need to be propagated to the compiler or linker for good
# libtool support.
case $host in
ia64-*-hpux*)
  # Find out which ABI we are using.
  echo 'int i;' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    case `/usr/bin/file conftest.$ac_objext` in
    *ELF-32*)
      HPUX_IA64_MODE="32"
      ;;
    *ELF-64*)
      HPUX_IA64_MODE="64"
      ;;
    esac
  fi
  rm -rf conftest*
  ;;
*-*-irix6*)
  # Find out which ABI we are using.
  echo '[#]line __oline__ "configure"' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
   if test "$lt_cv_prog_gnu_ld" = yes; then
    case `/usr/bin/file conftest.$ac_objext` in
    *32-bit*)
      LD="${LD-ld} -melf32bsmip"
      ;;
    *N32*)
      LD="${LD-ld} -melf32bmipn32"
      ;;
    *64-bit*)
      LD="${LD-ld} -melf64bmip"
      ;;
    esac
   else
    case `/usr/bin/file conftest.$ac_objext` in
    *32-bit*)
      LD="${LD-ld} -32"
      ;;
    *N32*)
      LD="${LD-ld} -n32"
      ;;
    *64-bit*)
      LD="${LD-ld} -64"
      ;;
    esac
   fi
  fi
  rm -rf conftest*
  ;;

x86_64-*linux*|ppc*-*linux*|powerpc*-*linux*|s390*-*linux*|sparc*-*linux*)
  # Find out which ABI we are using.
  echo 'int i;' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    case `/usr/bin/file conftest.o` in
    *32-bit*)
      case $host in
        x86_64-*linux*)
          LD="${LD-ld} -m elf_i386"
          ;;
        ppc64-*linux*|powerpc64-*linux*)
          LD="${LD-ld} -m elf32ppclinux"
          ;;
        s390x-*linux*)
          LD="${LD-ld} -m elf_s390"
          ;;
        sparc64-*linux*)
          LD="${LD-ld} -m elf32_sparc"
          ;;
      esac
      ;;
    *64-bit*)
      case $host in
        x86_64-*linux*)
          LD="${LD-ld} -m elf_x86_64"
          ;;
        ppc*-*linux*|powerpc*-*linux*)
          LD="${LD-ld} -m elf64ppc"
          ;;
        s390*-*linux*)
          LD="${LD-ld} -m elf64_s390"
          ;;
        sparc*-*linux*)
          LD="${LD-ld} -m elf64_sparc"
          ;;
      esac
      ;;
    esac
  fi
  rm -rf conftest*
  ;;

*-*-sco3.2v5*)
  # On SCO OpenServer 5, we need -belf to get full-featured binaries.
  SAVE_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -belf"
  AC_CACHE_CHECK([whether the C compiler needs -belf], lt_cv_cc_needs_belf,
    [AC_LANG_PUSH(C)
     AC_TRY_LINK([],[],[lt_cv_cc_needs_belf=yes],[lt_cv_cc_needs_belf=no])
     AC_LANG_POP])
  if test x"$lt_cv_cc_needs_belf" != x"yes"; then
    # this is probably gcc 2.8.0, egcs 1.0 or newer; no need for -belf
    CFLAGS="$SAVE_CFLAGS"
  fi
  ;;
sparc*-*solaris*)
  # Find out which ABI we are using.
  echo 'int i;' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    case `/usr/bin/file conftest.o` in
    *64-bit*)
      case $lt_cv_prog_gnu_ld in
      yes*) LD="${LD-ld} -m elf64_sparc" ;;
      *)    LD="${LD-ld} -64" ;;
      esac
      ;;
    esac
  fi
  rm -rf conftest*
  ;;

AC_PROVIDE_IFELSE([AC_LIBTOOL_WIN32_DLL],
[*-*-cygwin* | *-*-mingw* | *-*-pw32*)
  AC_CHECK_TOOL(DLLTOOL, dlltool, false)
  AC_CHECK_TOOL(AS, as, false)
  AC_CHECK_TOOL(OBJDUMP, objdump, false)
  ;;
  ])
esac

need_locks="$enable_libtool_lock"

])# _LT_AC_LOCK


# AC_LIBTOOL_COMPILER_OPTION(MESSAGE, VARIABLE-NAME, FLAGS,
#		[OUTPUT-FILE], [ACTION-SUCCESS], [ACTION-FAILURE])
# ----------------------------------------------------------------
# Check whether the given compiler option works
AC_DEFUN([AC_LIBTOOL_COMPILER_OPTION],
[AC_REQUIRE([LT_AC_PROG_SED])
AC_CACHE_CHECK([$1], [$2],
  [$2=no
  ifelse([$4], , [ac_outfile=conftest.$ac_objext], [ac_outfile=$4])
   printf "$lt_simple_compile_test_code" > conftest.$ac_ext
   lt_compiler_flag="$3"
   # Insert the option either (1) after the last *FLAGS variable, or
   # (2) before a word containing "conftest.", or (3) at the end.
   # Note that $ac_compile itself does not contain backslashes and begins
   # with a dollar sign (not a hyphen), so the echo should work correctly.
   # The option is referenced via a variable to avoid confusing sed.
   lt_compile=`echo "$ac_compile" | $SED \
   -e 's:.*FLAGS}\{0,1\} :&$lt_compiler_flag :; t' \
   -e 's: [[^ ]]*conftest\.: $lt_compiler_flag&:; t' \
   -e 's:$: $lt_compiler_flag:'`
   (eval echo "\"\$as_me:__oline__: $lt_compile\"" >&AS_MESSAGE_LOG_FD)
   (eval "$lt_compile" 2>conftest.err)
   ac_status=$?
   cat conftest.err >&AS_MESSAGE_LOG_FD
   echo "$as_me:__oline__: \$? = $ac_status" >&AS_MESSAGE_LOG_FD
   if (exit $ac_status) && test -s "$ac_outfile"; then
     # The compiler can only warn and ignore the option if not recognized
     # So say no if there are warnings other than the usual output.
     $echo "X$_lt_compiler_boilerplate" | $Xsed -e '/^$/d' >conftest.exp
     $SED '/^$/d; /^ *+/d' conftest.err >conftest.er2
     if test ! -s conftest.er2 || diff conftest.exp conftest.er2 >/dev/null; then
       $2=yes
     fi
   fi
   $rm conftest*
])

if test x"[$]$2" = xyes; then
    ifelse([$5], , :, [$5])
else
    ifelse([$6], , :, [$6])
fi
])# AC_LIBTOOL_COMPILER_OPTION


# AC_LIBTOOL_LINKER_OPTION(MESSAGE, VARIABLE-NAME, FLAGS,
#                          [ACTION-SUCCESS], [ACTION-FAILURE])
# ------------------------------------------------------------
# Check whether the given compiler option works
AC_DEFUN([AC_LIBTOOL_LINKER_OPTION],
[AC_CACHE_CHECK([$1], [$2],
  [$2=no
   save_LDFLAGS="$LDFLAGS"
   LDFLAGS="$LDFLAGS $3"
   printf "$lt_simple_link_test_code" > conftest.$ac_ext
   if (eval $ac_link 2>conftest.err) && test -s conftest$ac_exeext; then
     # The linker can only warn and ignore the option if not recognized
     # So say no if there are warnings
     if test -s conftest.err; then
       # Append any errors to the config.log.
       cat conftest.err 1>&AS_MESSAGE_LOG_FD
       $echo "X$_lt_linker_boilerplate" | $Xsed -e '/^$/d' > conftest.exp
       $SED '/^$/d; /^ *+/d' conftest.err >conftest.er2
       if diff conftest.exp conftest.er2 >/dev/null; then
         $2=yes
       fi
     else
       $2=yes
     fi
   fi
   $rm conftest*
   LDFLAGS="$save_LDFLAGS"
])

if test x"[$]$2" = xyes; then
    ifelse([$4], , :, [$4])
else
    ifelse([$5], , :, [$5])
fi
])# AC_LIBTOOL_LINKER_OPTION


# AC_LIBTOOL_SYS_MAX_CMD_LEN
# --------------------------
AC_DEFUN([AC_LIBTOOL_SYS_MAX_CMD_LEN],
[# find the maximum length of command line arguments
AC_MSG_CHECKING([the maximum length of command line arguments])
AC_CACHE_VAL([lt_cv_sys_max_cmd_len], [dnl
  i=0
  teststring="ABCD"

  case $build_os in
  msdosdjgpp*)
    # On DJGPP, this test can blow up pretty badly due to problems in libc
    # (any single argument exceeding 2000 bytes causes a buffer overrun
    # during glob expansion).  Even if it were fixed, the result of this
    # check would be larger than it should be.
    lt_cv_sys_max_cmd_len=12288;    # 12K is about right
    ;;

  gnu*)
    # Under GNU Hurd, this test is not required because there is
    # no limit to the length of command line arguments.
    # Libtool will interpret -1 as no limit whatsoever
    lt_cv_sys_max_cmd_len=-1;
    ;;

  cygwin* | mingw*)
    # On Win9x/ME, this test blows up -- it succeeds, but takes
    # about 5 minutes as the teststring grows exponentially.
    # Worse, since 9x/ME are not pre-emptively multitasking,
    # you end up with a "frozen" computer, even though with patience
    # the test eventually succeeds (with a max line length of 256k).
    # Instead, let's just punt: use the minimum linelength reported by
    # all of the supported platforms: 8192 (on NT/2K/XP).
    lt_cv_sys_max_cmd_len=8192;
    ;;

  amigaos*)
    # On AmigaOS with pdksh, this test takes hours, literally.
    # So we just punt and use a minimum line length of 8192.
    lt_cv_sys_max_cmd_len=8192;
    ;;

  netbsd* | freebsd* | openbsd* | darwin* | dragonfly*)
    # This has been around since 386BSD, at least.  Likely further.
    if test -x /sbin/sysctl; then
      lt_cv_sys_max_cmd_len=`/sbin/sysctl -n kern.argmax`
    elif test -x /usr/sbin/sysctl; then
      lt_cv_sys_max_cmd_len=`/usr/sbin/sysctl -n kern.argmax`
    else
      lt_cv_sys_max_cmd_len=65536	# usable default for all BSDs
    fi
    # And add a safety zone
    lt_cv_sys_max_cmd_len=`expr $lt_cv_sys_max_cmd_len \/ 4`
    lt_cv_sys_max_cmd_len=`expr $lt_cv_sys_max_cmd_len \* 3`
    ;;

  interix*)
    # We know the value 262144 and hardcode it with a safety zone (like BSD)
    lt_cv_sys_max_cmd_len=196608
    ;;

  osf*)
    # Dr. Hans Ekkehard Plesser reports seeing a kernel panic running configure
    # due to this test when exec_disable_arg_limit is 1 on Tru64. It is not
    # nice to cause kernel panics so lets avoid the loop below.
    # First set a reasonable default.
    lt_cv_sys_max_cmd_len=16384
    #
    if test -x /sbin/sysconfig; then
      case `/sbin/sysconfig -q proc exec_disable_arg_limit` in
        *1*) lt_cv_sys_max_cmd_len=-1 ;;
      esac
    fi
    ;;
  sco3.2v5*)
    lt_cv_sys_max_cmd_len=102400
    ;;
  sysv5* | sco5v6* | sysv4.2uw2*)
    kargmax=`grep ARG_MAX /etc/conf/cf.d/stune 2>/dev/null`
    if test -n "$kargmax"; then
      lt_cv_sys_max_cmd_len=`echo $kargmax | sed 's/.*[[ 	]]//'`
    else
      lt_cv_sys_max_cmd_len=32768
    fi
    ;;
  *)
    # If test is not a shell built-in, we'll probably end up computing a
    # maximum length that is only half of the actual maximum length, but
    # we can't tell.
    SHELL=${SHELL-${CONFIG_SHELL-/bin/sh}}
    while (test "X"`$SHELL [$]0 --fallback-echo "X$teststring" 2>/dev/null` \
	       = "XX$teststring") >/dev/null 2>&1 &&
	    new_result=`expr "X$teststring" : ".*" 2>&1` &&
	    lt_cv_sys_max_cmd_len=$new_result &&
	    test $i != 17 # 1/2 MB should be enough
    do
      i=`expr $i + 1`
      teststring=$teststring$teststring
    done
    teststring=
    # Add a significant safety factor because C++ compilers can tack on massive
    # amounts of additional arguments before passing them to the linker.
    # It appears as though 1/2 is a usable value.
    lt_cv_sys_max_cmd_len=`expr $lt_cv_sys_max_cmd_len \/ 2`
    ;;
  esac
])
if test -n $lt_cv_sys_max_cmd_len ; then
  AC_MSG_RESULT($lt_cv_sys_max_cmd_len)
else
  AC_MSG_RESULT(none)
fi
])# AC_LIBTOOL_SYS_MAX_CMD_LEN


# _LT_AC_CHECK_DLFCN
# ------------------
AC_DEFUN([_LT_AC_CHECK_DLFCN],
[AC_CHECK_HEADERS(dlfcn.h)dnl
])# _LT_AC_CHECK_DLFCN


# _LT_AC_TRY_DLOPEN_SELF (ACTION-IF-TRUE, ACTION-IF-TRUE-W-USCORE,
#                           ACTION-IF-FALSE, ACTION-IF-CROSS-COMPILING)
# ---------------------------------------------------------------------
AC_DEFUN([_LT_AC_TRY_DLOPEN_SELF],
[AC_REQUIRE([_LT_AC_CHECK_DLFCN])dnl
if test "$cross_compiling" = yes; then :
  [$4]
else
  lt_dlunknown=0; lt_dlno_uscore=1; lt_dlneed_uscore=2
  lt_status=$lt_dlunknown
  cat > conftest.$ac_ext <<EOF
[#line __oline__ "configure"
#include "confdefs.h"

#if HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#include <stdio.h>

#ifdef RTLD_GLOBAL
#  define LT_DLGLOBAL		RTLD_GLOBAL
#else
#  ifdef DL_GLOBAL
#    define LT_DLGLOBAL		DL_GLOBAL
#  else
#    define LT_DLGLOBAL		0
#  endif
#endif

/* We may have to define LT_DLLAZY_OR_NOW in the command line if we
   find out it does not work in some platform. */
#ifndef LT_DLLAZY_OR_NOW
#  ifdef RTLD_LAZY
#    define LT_DLLAZY_OR_NOW		RTLD_LAZY
#  else
#    ifdef DL_LAZY
#      define LT_DLLAZY_OR_NOW		DL_LAZY
#    else
#      ifdef RTLD_NOW
#        define LT_DLLAZY_OR_NOW	RTLD_NOW
#      else
#        ifdef DL_NOW
#          define LT_DLLAZY_OR_NOW	DL_NOW
#        else
#          define LT_DLLAZY_OR_NOW	0
#        endif
#      endif
#    endif
#  endif
#endif

#ifdef __cplusplus
extern "C" void exit (int);
#endif

void fnord() { int i=42;}
int main ()
{
  void *self = dlopen (0, LT_DLGLOBAL|LT_DLLAZY_OR_NOW);
  int status = $lt_dlunknown;

  if (self)
    {
      if (dlsym (self,"fnord"))       status = $lt_dlno_uscore;
      else if (dlsym( self,"_fnord")) status = $lt_dlneed_uscore;
      /* dlclose (self); */
    }
  else
    puts (dlerror ());

    exit (status);
}]
EOF
  if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext} 2>/dev/null; then
    (./conftest; exit; ) >&AS_MESSAGE_LOG_FD 2>/dev/null
    lt_status=$?
    case x$lt_status in
      x$lt_dlno_uscore) $1 ;;
      x$lt_dlneed_uscore) $2 ;;
      x$lt_dlunknown|x*) $3 ;;
    esac
  else :
    # compilation failed
    $3
  fi
fi
rm -fr conftest*
])# _LT_AC_TRY_DLOPEN_SELF


# AC_LIBTOOL_DLOPEN_SELF
# ----------------------
AC_DEFUN([AC_LIBTOOL_DLOPEN_SELF],
[AC_REQUIRE([_LT_AC_CHECK_DLFCN])dnl
if test "x$enable_dlopen" != xyes; then
  enable_dlopen=unknown
  enable_dlopen_self=unknown
  enable_dlopen_self_static=unknown
else
  lt_cv_dlopen=no
  lt_cv_dlopen_libs=

  case $host_os in
  beos*)
    lt_cv_dlopen="load_add_on"
    lt_cv_dlopen_libs=
    lt_cv_dlopen_self=yes
    ;;

  mingw* | pw32*)
    lt_cv_dlopen="LoadLibrary"
    lt_cv_dlopen_libs=
   ;;

  cygwin*)
    lt_cv_dlopen="dlopen"
    lt_cv_dlopen_libs=
   ;;

  darwin*)
  # if libdl is installed we need to link against it
    AC_CHECK_LIB([dl], [dlopen],
		[lt_cv_dlopen="dlopen" lt_cv_dlopen_libs="-ldl"],[
    lt_cv_dlopen="dyld"
    lt_cv_dlopen_libs=
    lt_cv_dlopen_self=yes
    ])
   ;;

  *)
    AC_CHECK_FUNC([shl_load],
	  [lt_cv_dlopen="shl_load"],
      [AC_CHECK_LIB([dld], [shl_load],
	    [lt_cv_dlopen="shl_load" lt_cv_dlopen_libs="-dld"],
	[AC_CHECK_FUNC([dlopen],
	      [lt_cv_dlopen="dlopen"],
	  [AC_CHECK_LIB([dl], [dlopen],
		[lt_cv_dlopen="dlopen" lt_cv_dlopen_libs="-ldl"],
	    [AC_CHECK_LIB([svld], [dlopen],
		  [lt_cv_dlopen="dlopen" lt_cv_dlopen_libs="-lsvld"],
	      [AC_CHECK_LIB([dld], [dld_link],
		    [lt_cv_dlopen="dld_link" lt_cv_dlopen_libs="-dld"])
	      ])
	    ])
	  ])
	])
      ])
    ;;
  esac

  if test "x$lt_cv_dlopen" != xno; then
    enable_dlopen=yes
  else
    enable_dlopen=no
  fi

  case $lt_cv_dlopen in
  dlopen)
    save_CPPFLAGS="$CPPFLAGS"
    test "x$ac_cv_header_dlfcn_h" = xyes && CPPFLAGS="$CPPFLAGS -DHAVE_DLFCN_H"

    save_LDFLAGS="$LDFLAGS"
    wl=$lt_prog_compiler_wl eval LDFLAGS=\"\$LDFLAGS $export_dynamic_flag_spec\"

    save_LIBS="$LIBS"
    LIBS="$lt_cv_dlopen_libs $LIBS"

    AC_CACHE_CHECK([whether a program can dlopen itself],
	  lt_cv_dlopen_self, [dnl
	  _LT_AC_TRY_DLOPEN_SELF(
	    lt_cv_dlopen_self=yes, lt_cv_dlopen_self=yes,
	    lt_cv_dlopen_self=no, lt_cv_dlopen_self=cross)
    ])

    if test "x$lt_cv_dlopen_self" = xyes; then
      wl=$lt_prog_compiler_wl eval LDFLAGS=\"\$LDFLAGS $lt_prog_compiler_static\"
      AC_CACHE_CHECK([whether a statically linked program can dlopen itself],
    	  lt_cv_dlopen_self_static, [dnl
	  _LT_AC_TRY_DLOPEN_SELF(
	    lt_cv_dlopen_self_static=yes, lt_cv_dlopen_self_static=yes,
	    lt_cv_dlopen_self_static=no,  lt_cv_dlopen_self_static=cross)
      ])
    fi

    CPPFLAGS="$save_CPPFLAGS"
    LDFLAGS="$save_LDFLAGS"
    LIBS="$save_LIBS"
    ;;
  esac

  case $lt_cv_dlopen_self in
  yes|no) enable_dlopen_self=$lt_cv_dlopen_self ;;
  *) enable_dlopen_self=unknown ;;
  esac

  case $lt_cv_dlopen_self_static in
  yes|no) enable_dlopen_self_static=$lt_cv_dlopen_self_static ;;
  *) enable_dlopen_self_static=unknown ;;
  esac
fi
])# AC_LIBTOOL_DLOPEN_SELF


# AC_LIBTOOL_PROG_CC_C_O([TAGNAME])
# ---------------------------------
# Check to see if options -c and -o are simultaneously supported by compiler
AC_DEFUN([AC_LIBTOOL_PROG_CC_C_O],
[AC_REQUIRE([_LT_AC_SYS_COMPILER])dnl
AC_CACHE_CHECK([if $compiler supports -c -o file.$ac_objext],
  [_LT_AC_TAGVAR(lt_cv_prog_compiler_c_o, $1)],
  [_LT_AC_TAGVAR(lt_cv_prog_compiler_c_o, $1)=no
   $rm -r conftest 2>/dev/null
   mkdir conftest
   cd conftest
   mkdir out
   printf "$lt_simple_compile_test_code" > conftest.$ac_ext

   lt_compiler_flag="-o out/conftest2.$ac_objext"
   # Insert the option either (1) after the last *FLAGS variable, or
   # (2) before a word containing "conftest.", or (3) at the end.
   # Note that $ac_compile itself does not contain backslashes and begins
   # with a dollar sign (not a hyphen), so the echo should work correctly.
   lt_compile=`echo "$ac_compile" | $SED \
   -e 's:.*FLAGS}\{0,1\} :&$lt_compiler_flag :; t' \
   -e 's: [[^ ]]*conftest\.: $lt_compiler_flag&:; t' \
   -e 's:$: $lt_compiler_flag:'`
   (eval echo "\"\$as_me:__oline__: $lt_compile\"" >&AS_MESSAGE_LOG_FD)
   (eval "$lt_compile" 2>out/conftest.err)
   ac_status=$?
   cat out/conftest.err >&AS_MESSAGE_LOG_FD
   echo "$as_me:__oline__: \$? = $ac_status" >&AS_MESSAGE_LOG_FD
   if (exit $ac_status) && test -s out/conftest2.$ac_objext
   then
     # The compiler can only warn and ignore the option if not recognized
     # So say no if there are warnings
     $echo "X$_lt_compiler_boilerplate" | $Xsed -e '/^$/d' > out/conftest.exp
     $SED '/^$/d; /^ *+/d' out/conftest.err >out/conftest.er2
     if test ! -s out/conftest.er2 || diff out/conftest.exp out/conftest.er2 >/dev/null; then
       _LT_AC_TAGVAR(lt_cv_prog_compiler_c_o, $1)=yes
     fi
   fi
   chmod u+w . 2>&AS_MESSAGE_LOG_FD
   $rm conftest*
   # SGI C++ compiler will create directory out/ii_files/ for
   # template instantiation
   test -d out/ii_files && $rm out/ii_files/* && rmdir out/ii_files
   $rm out/* && rmdir out
   cd ..
   rmdir conftest
   $rm conftest*
])
])# AC_LIBTOOL_PROG_CC_C_O


# AC_LIBTOOL_SYS_HARD_LINK_LOCKS([TAGNAME])
# -----------------------------------------
# Check to see if we can do hard links to lock some files if needed
AC_DEFUN([AC_LIBTOOL_SYS_HARD_LINK_LOCKS],
[AC_REQUIRE([_LT_AC_LOCK])dnl

hard_links="nottested"
if test "$_LT_AC_TAGVAR(lt_cv_prog_compiler_c_o, $1)" = no && test "$need_locks" != no; then
  # do not overwrite the value of need_locks provided by the user
  AC_MSG_CHECKING([if we can lock with hard links])
  hard_links=yes
  $rm conftest*
  ln conftest.a conftest.b 2>/dev/null && hard_links=no
  touch conftest.a
  ln conftest.a conftest.b 2>&5 || hard_links=no
  ln conftest.a conftest.b 2>/dev/null && hard_links=no
  AC_MSG_RESULT([$hard_links])
  if test "$hard_links" = no; then
    AC_MSG_WARN([`$CC' does not support `-c -o', so `make -j' may be unsafe])
    need_locks=warn
  fi
else
  need_locks=no
fi
])# AC_LIBTOOL_SYS_HARD_LINK_LOCKS


# AC_LIBTOOL_OBJDIR
# -----------------
AC_DEFUN([AC_LIBTOOL_OBJDIR],
[AC_CACHE_CHECK([for objdir], [lt_cv_objdir],
[rm -f .libs 2>/dev/null
mkdir .libs 2>/dev/null
if test -d .libs; then
  lt_cv_objdir=.libs
else
  # MS-DOS does not allow filenames that begin with a dot.
  lt_cv_objdir=_libs
fi
rmdir .libs 2>/dev/null])
objdir=$lt_cv_objdir
])# AC_LIBTOOL_OBJDIR


# AC_LIBTOOL_PROG_LD_HARDCODE_LIBPATH([TAGNAME])
# ----------------------------------------------
# Check hardcoding attributes.
AC_DEFUN([AC_LIBTOOL_PROG_LD_HARDCODE_LIBPATH],
[AC_MSG_CHECKING([how to hardcode library paths into programs])
_LT_AC_TAGVAR(hardcode_action, $1)=
if test -n "$_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)" || \
   test -n "$_LT_AC_TAGVAR(runpath_var, $1)" || \
   test "X$_LT_AC_TAGVAR(hardcode_automatic, $1)" = "Xyes" ; then

  # We can hardcode non-existant directories.
  if test "$_LT_AC_TAGVAR(hardcode_direct, $1)" != no &&
     # If the only mechanism to avoid hardcoding is shlibpath_var, we
     # have to relink, otherwise we might link with an installed library
     # when we should be linking with a yet-to-be-installed one
     ## test "$_LT_AC_TAGVAR(hardcode_shlibpath_var, $1)" != no &&
     test "$_LT_AC_TAGVAR(hardcode_minus_L, $1)" != no; then
    # Linking always hardcodes the temporary library directory.
    _LT_AC_TAGVAR(hardcode_action, $1)=relink
  else
    # We can link without hardcoding, and we can hardcode nonexisting dirs.
    _LT_AC_TAGVAR(hardcode_action, $1)=immediate
  fi
else
  # We cannot hardcode anything, or else we can only hardcode existing
  # directories.
  _LT_AC_TAGVAR(hardcode_action, $1)=unsupported
fi
AC_MSG_RESULT([$_LT_AC_TAGVAR(hardcode_action, $1)])

if test "$_LT_AC_TAGVAR(hardcode_action, $1)" = relink; then
  # Fast installation is not supported
  enable_fast_install=no
elif test "$shlibpath_overrides_runpath" = yes ||
     test "$enable_shared" = no; then
  # Fast installation is not necessary
  enable_fast_install=needless
fi
])# AC_LIBTOOL_PROG_LD_HARDCODE_LIBPATH


# AC_LIBTOOL_SYS_LIB_STRIP
# ------------------------
AC_DEFUN([AC_LIBTOOL_SYS_LIB_STRIP],
[striplib=
old_striplib=
AC_MSG_CHECKING([whether stripping libraries is possible])
if test -n "$STRIP" && $STRIP -V 2>&1 | grep "GNU strip" >/dev/null; then
  test -z "$old_striplib" && old_striplib="$STRIP --strip-debug"
  test -z "$striplib" && striplib="$STRIP --strip-unneeded"
  AC_MSG_RESULT([yes])
else
# FIXME - insert some real tests, host_os isn't really good enough
  case $host_os in
   darwin*)
       if test -n "$STRIP" ; then
         striplib="$STRIP -x"
         AC_MSG_RESULT([yes])
       else
  AC_MSG_RESULT([no])
fi
       ;;
   *)
  AC_MSG_RESULT([no])
    ;;
  esac
fi
])# AC_LIBTOOL_SYS_LIB_STRIP


# AC_LIBTOOL_SYS_DYNAMIC_LINKER
# -----------------------------
# PORTME Fill in your ld.so characteristics
AC_DEFUN([AC_LIBTOOL_SYS_DYNAMIC_LINKER],
[AC_MSG_CHECKING([dynamic linker characteristics])
library_names_spec=
libname_spec='lib$name'
soname_spec=
shrext_cmds=".so"
postinstall_cmds=
postuninstall_cmds=
finish_cmds=
finish_eval=
shlibpath_var=
shlibpath_overrides_runpath=unknown
version_type=none
dynamic_linker="$host_os ld.so"
sys_lib_dlsearch_path_spec="/lib /usr/lib"
if test "$GCC" = yes; then
  sys_lib_search_path_spec=`$CC -print-search-dirs | grep "^libraries:" | $SED -e "s/^libraries://" -e "s,=/,/,g"`
  if echo "$sys_lib_search_path_spec" | grep ';' >/dev/null ; then
    # if the path contains ";" then we assume it to be the separator
    # otherwise default to the standard path separator (i.e. ":") - it is
    # assumed that no part of a normal pathname contains ";" but that should
    # okay in the real world where ";" in dirpaths is itself problematic.
    sys_lib_search_path_spec=`echo "$sys_lib_search_path_spec" | $SED -e 's/;/ /g'`
  else
    sys_lib_search_path_spec=`echo "$sys_lib_search_path_spec" | $SED  -e "s/$PATH_SEPARATOR/ /g"`
  fi
else
  sys_lib_search_path_spec="/lib /usr/lib /usr/local/lib"
fi
need_lib_prefix=unknown
hardcode_into_libs=no

# when you set need_version to no, make sure it does not cause -set_version
# flags to be left without arguments
need_version=unknown

case $host_os in
aix3*)
  version_type=linux
  library_names_spec='${libname}${release}${shared_ext}$versuffix $libname.a'
  shlibpath_var=LIBPATH

  # AIX 3 has no versioning support, so we append a major version to the name.
  soname_spec='${libname}${release}${shared_ext}$major'
  ;;

aix4* | aix5*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  hardcode_into_libs=yes
  if test "$host_cpu" = ia64; then
    # AIX 5 supports IA64
    library_names_spec='${libname}${release}${shared_ext}$major ${libname}${release}${shared_ext}$versuffix $libname${shared_ext}'
    shlibpath_var=LD_LIBRARY_PATH
  else
    # With GCC up to 2.95.x, collect2 would create an import file
    # for dependence libraries.  The import file would start with
    # the line `#! .'.  This would cause the generated library to
    # depend on `.', always an invalid library.  This was fixed in
    # development snapshots of GCC prior to 3.0.
    case $host_os in
      aix4 | aix4.[[01]] | aix4.[[01]].*)
      if { echo '#if __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ >= 97)'
	   echo ' yes '
	   echo '#endif'; } | ${CC} -E - | grep yes > /dev/null; then
	:
      else
	can_build_shared=no
      fi
      ;;
    esac
    # AIX (on Power*) has no versioning support, so currently we can not hardcode correct
    # soname into executable. Probably we can add versioning support to
    # collect2, so additional links can be useful in future.
    if test "$aix_use_runtimelinking" = yes; then
      # If using run time linking (on AIX 4.2 or later) use lib<name>.so
      # instead of lib<name>.a to let people know that these are not
      # typical AIX shared libraries.
      library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
    else
      # We preserve .a as extension for shared libraries through AIX4.2
      # and later when we are not doing run time linking.
      library_names_spec='${libname}${release}.a $libname.a'
      soname_spec='${libname}${release}${shared_ext}$major'
    fi
    shlibpath_var=LIBPATH
  fi
  ;;

amigaos*)
  library_names_spec='$libname.ixlibrary $libname.a'
  # Create ${libname}_ixlibrary.a entries in /sys/libs.
  finish_eval='for lib in `ls $libdir/*.ixlibrary 2>/dev/null`; do libname=`$echo "X$lib" | $Xsed -e '\''s%^.*/\([[^/]]*\)\.ixlibrary$%\1%'\''`; test $rm /sys/libs/${libname}_ixlibrary.a; $show "cd /sys/libs && $LN_S $lib ${libname}_ixlibrary.a"; cd /sys/libs && $LN_S $lib ${libname}_ixlibrary.a || exit 1; done'
  ;;

beos*)
  library_names_spec='${libname}${shared_ext}'
  dynamic_linker="$host_os ld.so"
  shlibpath_var=LIBRARY_PATH
  ;;

bsdi[[45]]*)
  version_type=linux
  need_version=no
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  finish_cmds='PATH="\$PATH:/sbin" ldconfig $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  sys_lib_search_path_spec="/shlib /usr/lib /usr/X11/lib /usr/contrib/lib /lib /usr/local/lib"
  sys_lib_dlsearch_path_spec="/shlib /usr/lib /usr/local/lib"
  # the default ld.so.conf also contains /usr/contrib/lib and
  # /usr/X11R6/lib (/usr/X11 is a link to /usr/X11R6), but let us allow
  # libtool to hard-code these into programs
  ;;

cygwin* | mingw* | pw32*)
  version_type=windows
  shrext_cmds=".dll"
  need_version=no
  need_lib_prefix=no

  case $GCC,$host_os in
  yes,cygwin* | yes,mingw* | yes,pw32*)
    library_names_spec='$libname.dll.a'
    # DLL is installed to $(libdir)/../bin by postinstall_cmds
    postinstall_cmds='base_file=`basename \${file}`~
      dlpath=`$SHELL 2>&1 -c '\''. $dir/'\''\${base_file}'\''i;echo \$dlname'\''`~
      dldir=$destdir/`dirname \$dlpath`~
      test -d \$dldir || mkdir -p \$dldir~
      $install_prog $dir/$dlname \$dldir/$dlname~
      chmod a+x \$dldir/$dlname'
    postuninstall_cmds='dldll=`$SHELL 2>&1 -c '\''. $file; echo \$dlname'\''`~
      dlpath=$dir/\$dldll~
       $rm \$dlpath'
    shlibpath_overrides_runpath=yes

    case $host_os in
    cygwin*)
      # Cygwin DLLs use 'cyg' prefix rather than 'lib'
      soname_spec='`echo ${libname} | sed -e 's/^lib/cyg/'``echo ${release} | $SED -e 's/[[.]]/-/g'`${versuffix}${shared_ext}'
      sys_lib_search_path_spec="/usr/lib /lib/w32api /lib /usr/local/lib"
      ;;
    mingw*)
      # MinGW DLLs use traditional 'lib' prefix
      soname_spec='${libname}`echo ${release} | $SED -e 's/[[.]]/-/g'`${versuffix}${shared_ext}'
      sys_lib_search_path_spec=`$CC -print-search-dirs | grep "^libraries:" | $SED -e "s/^libraries://" -e "s,=/,/,g"`
      if echo "$sys_lib_search_path_spec" | [grep ';[c-zC-Z]:/' >/dev/null]; then
        # It is most probably a Windows format PATH printed by
        # mingw gcc, but we are running on Cygwin. Gcc prints its search
        # path with ; separators, and with drive letters. We can handle the
        # drive letters (cygwin fileutils understands them), so leave them,
        # especially as we might pass files found there to a mingw objdump,
        # which wouldn't understand a cygwinified path. Ahh.
        sys_lib_search_path_spec=`echo "$sys_lib_search_path_spec" | $SED -e 's/;/ /g'`
      else
        sys_lib_search_path_spec=`echo "$sys_lib_search_path_spec" | $SED  -e "s/$PATH_SEPARATOR/ /g"`
      fi
      ;;
    pw32*)
      # pw32 DLLs use 'pw' prefix rather than 'lib'
      library_names_spec='`echo ${libname} | sed -e 's/^lib/pw/'``echo ${release} | $SED -e 's/[[.]]/-/g'`${versuffix}${shared_ext}'
      ;;
    esac
    ;;

  *)
    library_names_spec='${libname}`echo ${release} | $SED -e 's/[[.]]/-/g'`${versuffix}${shared_ext} $libname.lib'
    ;;
  esac
  dynamic_linker='Win32 ld.exe'
  # FIXME: first we should search . and the directory the executable is in
  shlibpath_var=PATH
  ;;

darwin* | rhapsody*)
  dynamic_linker="$host_os dyld"
  version_type=darwin
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}${versuffix}$shared_ext ${libname}${release}${major}$shared_ext ${libname}$shared_ext'
  soname_spec='${libname}${release}${major}$shared_ext'
  shlibpath_overrides_runpath=yes
  shlibpath_var=DYLD_LIBRARY_PATH
  shrext_cmds='`test .$module = .yes && echo .so || echo .dylib`'
  # Apple's gcc prints 'gcc -print-search-dirs' doesn't operate the same.
  if test "$GCC" = yes; then
    sys_lib_search_path_spec=`$CC -print-search-dirs | tr "\n" "$PATH_SEPARATOR" | sed -e 's/libraries:/@libraries:/' | tr "@" "\n" | grep "^libraries:" | sed -e "s/^libraries://" -e "s,=/,/,g" -e "s,$PATH_SEPARATOR, ,g" -e "s,.*,& /lib /usr/lib /usr/local/lib,g"`
  else
    sys_lib_search_path_spec='/lib /usr/lib /usr/local/lib'
  fi
  sys_lib_dlsearch_path_spec='/usr/local/lib /lib /usr/lib'
  ;;

dgux*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname$shared_ext'
  soname_spec='${libname}${release}${shared_ext}$major'
  shlibpath_var=LD_LIBRARY_PATH
  ;;

freebsd1*)
  dynamic_linker=no
  ;;

kfreebsd*-gnu)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major ${libname}${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=no
  hardcode_into_libs=yes
  dynamic_linker='GNU ld.so'
  ;;

freebsd* | dragonfly*)
  # DragonFly does not have aout.  When/if they implement a new
  # versioning mechanism, adjust this.
  if test -x /usr/bin/objformat; then
    objformat=`/usr/bin/objformat`
  else
    case $host_os in
    freebsd[[123]]*) objformat=aout ;;
    *) objformat=elf ;;
    esac
  fi
  version_type=freebsd-$objformat
  case $version_type in
    freebsd-elf*)
      library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext} $libname${shared_ext}'
      need_version=no
      need_lib_prefix=no
      ;;
    freebsd-*)
      library_names_spec='${libname}${release}${shared_ext}$versuffix $libname${shared_ext}$versuffix'
      need_version=yes
      ;;
  esac
  shlibpath_var=LD_LIBRARY_PATH
  case $host_os in
  freebsd2*)
    shlibpath_overrides_runpath=yes
    ;;
  freebsd3.[[01]]* | freebsdelf3.[[01]]*)
    shlibpath_overrides_runpath=yes
    hardcode_into_libs=yes
    ;;
  freebsd3.[[2-9]]* | freebsdelf3.[[2-9]]* | \
  freebsd4.[[0-5]] | freebsdelf4.[[0-5]] | freebsd4.1.1 | freebsdelf4.1.1)
    shlibpath_overrides_runpath=no
    hardcode_into_libs=yes
    ;;
  freebsd*) # from 4.6 on
    shlibpath_overrides_runpath=yes
    hardcode_into_libs=yes
    ;;
  esac
  ;;

gnu*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}${major} ${libname}${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  shlibpath_var=LD_LIBRARY_PATH
  hardcode_into_libs=yes
  ;;

hpux9* | hpux10* | hpux11*)
  # Give a soname corresponding to the major version so that dld.sl refuses to
  # link against other versions.
  version_type=sunos
  need_lib_prefix=no
  need_version=no
  case $host_cpu in
  ia64*)
    shrext_cmds='.so'
    hardcode_into_libs=yes
    dynamic_linker="$host_os dld.so"
    shlibpath_var=LD_LIBRARY_PATH
    shlibpath_overrides_runpath=yes # Unless +noenvvar is specified.
    library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
    soname_spec='${libname}${release}${shared_ext}$major'
    if test "X$HPUX_IA64_MODE" = X32; then
      sys_lib_search_path_spec="/usr/lib/hpux32 /usr/local/lib/hpux32 /usr/local/lib"
    else
      sys_lib_search_path_spec="/usr/lib/hpux64 /usr/local/lib/hpux64"
    fi
    sys_lib_dlsearch_path_spec=$sys_lib_search_path_spec
    ;;
   hppa*64*)
     shrext_cmds='.sl'
     hardcode_into_libs=yes
     dynamic_linker="$host_os dld.sl"
     shlibpath_var=LD_LIBRARY_PATH # How should we handle SHLIB_PATH
     shlibpath_overrides_runpath=yes # Unless +noenvvar is specified.
     library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
     soname_spec='${libname}${release}${shared_ext}$major'
     sys_lib_search_path_spec="/usr/lib/pa20_64 /usr/ccs/lib/pa20_64"
     sys_lib_dlsearch_path_spec=$sys_lib_search_path_spec
     ;;
   *)
    shrext_cmds='.sl'
    dynamic_linker="$host_os dld.sl"
    shlibpath_var=SHLIB_PATH
    shlibpath_overrides_runpath=no # +s is required to enable SHLIB_PATH
    library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
    soname_spec='${libname}${release}${shared_ext}$major'
    ;;
  esac
  # HP-UX runs *really* slowly unless shared libraries are mode 555.
  postinstall_cmds='chmod 555 $lib'
  ;;

interix3*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major ${libname}${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  dynamic_linker='Interix 3.x ld.so.1 (PE, like ELF)'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=no
  hardcode_into_libs=yes
  ;;

irix5* | irix6* | nonstopux*)
  case $host_os in
    nonstopux*) version_type=nonstopux ;;
    *)
	if test "$lt_cv_prog_gnu_ld" = yes; then
		version_type=linux
	else
		version_type=irix
	fi ;;
  esac
  need_lib_prefix=no
  need_version=no
  soname_spec='${libname}${release}${shared_ext}$major'
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major ${libname}${release}${shared_ext} $libname${shared_ext}'
  case $host_os in
  irix5* | nonstopux*)
    libsuff= shlibsuff=
    ;;
  *)
    case $LD in # libtool.m4 will add one of these switches to LD
    *-32|*"-32 "|*-melf32bsmip|*"-melf32bsmip ")
      libsuff= shlibsuff= libmagic=32-bit;;
    *-n32|*"-n32 "|*-melf32bmipn32|*"-melf32bmipn32 ")
      libsuff=32 shlibsuff=N32 libmagic=N32;;
    *-64|*"-64 "|*-melf64bmip|*"-melf64bmip ")
      libsuff=64 shlibsuff=64 libmagic=64-bit;;
    *) libsuff= shlibsuff= libmagic=never-match;;
    esac
    ;;
  esac
  shlibpath_var=LD_LIBRARY${shlibsuff}_PATH
  shlibpath_overrides_runpath=no
  sys_lib_search_path_spec="/usr/lib${libsuff} /lib${libsuff} /usr/local/lib${libsuff}"
  sys_lib_dlsearch_path_spec="/usr/lib${libsuff} /lib${libsuff}"
  hardcode_into_libs=yes
  ;;

# No shared lib support for Linux oldld, aout, or coff.
linux*oldld* | linux*aout* | linux*coff*)
  dynamic_linker=no
  ;;

# This must be Linux ELF.
linux*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  finish_cmds='PATH="\$PATH:/sbin" ldconfig -n $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=no
  # This implies no fast_install, which is unacceptable.
  # Some rework will be needed to allow for fast_install
  # before this can be enabled.
  hardcode_into_libs=yes

  # Append ld.so.conf contents to the search path
  if test -f /etc/ld.so.conf; then
    lt_ld_extra=`awk '/^include / { system(sprintf("cd /etc; cat %s", \[$]2)); skip = 1; } { if (!skip) print \[$]0; skip = 0; }' < /etc/ld.so.conf | $SED -e 's/#.*//;s/[:,	]/ /g;s/=[^=]*$//;s/=[^= ]* / /g;/^$/d' | tr '\n' ' '`
    sys_lib_dlsearch_path_spec="/lib /usr/lib $lt_ld_extra"
  fi

  # We used to test for /lib/ld.so.1 and disable shared libraries on
  # powerpc, because MkLinux only supported shared libraries with the
  # GNU dynamic linker.  Since this was broken with cross compilers,
  # most powerpc-linux boxes support dynamic linking these days and
  # people can always --disable-shared, the test was removed, and we
  # assume the GNU/Linux dynamic linker is in use.
  dynamic_linker='GNU/Linux ld.so'
  ;;

knetbsd*-gnu)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major ${libname}${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=no
  hardcode_into_libs=yes
  dynamic_linker='GNU ld.so'
  ;;

netbsd*)
  version_type=sunos
  need_lib_prefix=no
  need_version=no
  if echo __ELF__ | $CC -E - | grep __ELF__ >/dev/null; then
    library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${shared_ext}$versuffix'
    finish_cmds='PATH="\$PATH:/sbin" ldconfig -m $libdir'
    dynamic_linker='NetBSD (a.out) ld.so'
  else
    library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major ${libname}${shared_ext}'
    soname_spec='${libname}${release}${shared_ext}$major'
    dynamic_linker='NetBSD ld.elf_so'
  fi
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  hardcode_into_libs=yes
  ;;

newsos6)
  version_type=linux
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  ;;

nto-qnx*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  ;;

openbsd*)
  version_type=sunos
  sys_lib_dlsearch_path_spec="/usr/lib"
  need_lib_prefix=no
  # Some older versions of OpenBSD (3.3 at least) *do* need versioned libs.
  case $host_os in
    openbsd3.3 | openbsd3.3.*) need_version=yes ;;
    *)                         need_version=no  ;;
  esac
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${shared_ext}$versuffix'
  finish_cmds='PATH="\$PATH:/sbin" ldconfig -m $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  if test -z "`echo __ELF__ | $CC -E - | grep __ELF__`" || test "$host_os-$host_cpu" = "openbsd2.8-powerpc"; then
    case $host_os in
      openbsd2.[[89]] | openbsd2.[[89]].*)
	shlibpath_overrides_runpath=no
	;;
      *)
	shlibpath_overrides_runpath=yes
	;;
      esac
  else
    shlibpath_overrides_runpath=yes
  fi
  ;;

os2*)
  libname_spec='$name'
  shrext_cmds=".dll"
  need_lib_prefix=no
  library_names_spec='$libname${shared_ext} $libname.a'
  dynamic_linker='OS/2 ld.exe'
  shlibpath_var=LIBPATH
  ;;

osf3* | osf4* | osf5*)
  version_type=osf
  need_lib_prefix=no
  need_version=no
  soname_spec='${libname}${release}${shared_ext}$major'
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
  shlibpath_var=LD_LIBRARY_PATH
  sys_lib_search_path_spec="/usr/shlib /usr/ccs/lib /usr/lib/cmplrs/cc /usr/lib /usr/local/lib /var/shlib"
  sys_lib_dlsearch_path_spec="$sys_lib_search_path_spec"
  ;;

solaris*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  hardcode_into_libs=yes
  # ldd complains unless libraries are executable
  postinstall_cmds='chmod +x $lib'
  ;;

sunos4*)
  version_type=sunos
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${shared_ext}$versuffix'
  finish_cmds='PATH="\$PATH:/usr/etc" ldconfig $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  if test "$with_gnu_ld" = yes; then
    need_lib_prefix=no
  fi
  need_version=yes
  ;;

sysv4 | sysv4.3*)
  version_type=linux
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  shlibpath_var=LD_LIBRARY_PATH
  case $host_vendor in
    sni)
      shlibpath_overrides_runpath=no
      need_lib_prefix=no
      export_dynamic_flag_spec='${wl}-Blargedynsym'
      runpath_var=LD_RUN_PATH
      ;;
    siemens)
      need_lib_prefix=no
      ;;
    motorola)
      need_lib_prefix=no
      need_version=no
      shlibpath_overrides_runpath=no
      sys_lib_search_path_spec='/lib /usr/lib /usr/ccs/lib'
      ;;
  esac
  ;;

sysv4*MP*)
  if test -d /usr/nec ;then
    version_type=linux
    library_names_spec='$libname${shared_ext}.$versuffix $libname${shared_ext}.$major $libname${shared_ext}'
    soname_spec='$libname${shared_ext}.$major'
    shlibpath_var=LD_LIBRARY_PATH
  fi
  ;;

sysv5* | sco3.2v5* | sco5v6* | unixware* | OpenUNIX* | sysv4*uw2*)
  version_type=freebsd-elf
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext} $libname${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  shlibpath_var=LD_LIBRARY_PATH
  hardcode_into_libs=yes
  if test "$with_gnu_ld" = yes; then
    sys_lib_search_path_spec='/usr/local/lib /usr/gnu/lib /usr/ccs/lib /usr/lib /lib'
    shlibpath_overrides_runpath=no
  else
    sys_lib_search_path_spec='/usr/ccs/lib /usr/lib'
    shlibpath_overrides_runpath=yes
    case $host_os in
      sco3.2v5*)
        sys_lib_search_path_spec="$sys_lib_search_path_spec /lib"
	;;
    esac
  fi
  sys_lib_dlsearch_path_spec='/usr/lib'
  ;;

uts4*)
  version_type=linux
  library_names_spec='${libname}${release}${shared_ext}$versuffix ${libname}${release}${shared_ext}$major $libname${shared_ext}'
  soname_spec='${libname}${release}${shared_ext}$major'
  shlibpath_var=LD_LIBRARY_PATH
  ;;

*)
  dynamic_linker=no
  ;;
esac
AC_MSG_RESULT([$dynamic_linker])
test "$dynamic_linker" = no && can_build_shared=no

variables_saved_for_relink="PATH $shlibpath_var $runpath_var"
if test "$GCC" = yes; then
  variables_saved_for_relink="$variables_saved_for_relink GCC_EXEC_PREFIX COMPILER_PATH LIBRARY_PATH"
fi
])# AC_LIBTOOL_SYS_DYNAMIC_LINKER


# _LT_AC_TAGCONFIG
# ----------------
AC_DEFUN([_LT_AC_TAGCONFIG],
[AC_ARG_WITH([tags],
    [AC_HELP_STRING([--with-tags@<:@=TAGS@:>@],
        [include additional configurations @<:@automatic@:>@])],
    [tagnames="$withval"])

if test -f "$ltmain" && test -n "$tagnames"; then
  if test ! -f "${ofile}"; then
    AC_MSG_WARN([output file `$ofile' does not exist])
  fi

  if test -z "$LTCC"; then
    eval "`$SHELL ${ofile} --config | grep '^LTCC='`"
    if test -z "$LTCC"; then
      AC_MSG_WARN([output file `$ofile' does not look like a libtool script])
    else
      AC_MSG_WARN([using `LTCC=$LTCC', extracted from `$ofile'])
    fi
  fi
  if test -z "$LTCFLAGS"; then
    eval "`$SHELL ${ofile} --config | grep '^LTCFLAGS='`"
  fi

  # Extract list of available tagged configurations in $ofile.
  # Note that this assumes the entire list is on one line.
  available_tags=`grep "^available_tags=" "${ofile}" | $SED -e 's/available_tags=\(.*$\)/\1/' -e 's/\"//g'`

  lt_save_ifs="$IFS"; IFS="${IFS}$PATH_SEPARATOR,"
  for tagname in $tagnames; do
    IFS="$lt_save_ifs"
    # Check whether tagname contains only valid characters
    case `$echo "X$tagname" | $Xsed -e 's:[[-_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890,/]]::g'` in
    "") ;;
    *)  AC_MSG_ERROR([invalid tag name: $tagname])
	;;
    esac

    if grep "^# ### BEGIN LIBTOOL TAG CONFIG: $tagname$" < "${ofile}" > /dev/null
    then
      AC_MSG_ERROR([tag name \"$tagname\" already exists])
    fi

    # Update the list of available tags.
    if test -n "$tagname"; then
      echo appending configuration tag \"$tagname\" to $ofile

      case $tagname in
      CXX)
	if test -n "$CXX" && ( test "X$CXX" != "Xno" &&
	    ( (test "X$CXX" = "Xg++" && `g++ -v >/dev/null 2>&1` ) ||
	    (test "X$CXX" != "Xg++"))) ; then
	  AC_LIBTOOL_LANG_CXX_CONFIG
	else
	  tagname=""
	fi
	;;

      F77)
	if test -n "$F77" && test "X$F77" != "Xno"; then
	  AC_LIBTOOL_LANG_F77_CONFIG
	else
	  tagname=""
	fi
	;;

      GCJ)
	if test -n "$GCJ" && test "X$GCJ" != "Xno"; then
	  AC_LIBTOOL_LANG_GCJ_CONFIG
	else
	  tagname=""
	fi
	;;

      RC)
	AC_LIBTOOL_LANG_RC_CONFIG
	;;

      *)
	AC_MSG_ERROR([Unsupported tag name: $tagname])
	;;
      esac

      # Append the new tag name to the list of available tags.
      if test -n "$tagname" ; then
      available_tags="$available_tags $tagname"
    fi
    fi
  done
  IFS="$lt_save_ifs"

  # Now substitute the updated list of available tags.
  if eval "sed -e 's/^available_tags=.*\$/available_tags=\"$available_tags\"/' \"$ofile\" > \"${ofile}T\""; then
    mv "${ofile}T" "$ofile"
    chmod +x "$ofile"
  else
    rm -f "${ofile}T"
    AC_MSG_ERROR([unable to update list of available tagged configurations.])
  fi
fi
])# _LT_AC_TAGCONFIG


# AC_LIBTOOL_DLOPEN
# -----------------
# enable checks for dlopen support
AC_DEFUN([AC_LIBTOOL_DLOPEN],
 [AC_BEFORE([$0],[AC_LIBTOOL_SETUP])
])# AC_LIBTOOL_DLOPEN


# AC_LIBTOOL_WIN32_DLL
# --------------------
# declare package support for building win32 DLLs
AC_DEFUN([AC_LIBTOOL_WIN32_DLL],
[AC_BEFORE([$0], [AC_LIBTOOL_SETUP])
])# AC_LIBTOOL_WIN32_DLL


# AC_ENABLE_SHARED([DEFAULT])
# ---------------------------
# implement the --enable-shared flag
# DEFAULT is either `yes' or `no'.  If omitted, it defaults to `yes'.
AC_DEFUN([AC_ENABLE_SHARED],
[define([AC_ENABLE_SHARED_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE([shared],
    [AC_HELP_STRING([--enable-shared@<:@=PKGS@:>@],
	[build shared libraries @<:@default=]AC_ENABLE_SHARED_DEFAULT[@:>@])],
    [p=${PACKAGE-default}
    case $enableval in
    yes) enable_shared=yes ;;
    no) enable_shared=no ;;
    *)
      enable_shared=no
      # Look at the argument we got.  We use all the common list separators.
      lt_save_ifs="$IFS"; IFS="${IFS}$PATH_SEPARATOR,"
      for pkg in $enableval; do
	IFS="$lt_save_ifs"
	if test "X$pkg" = "X$p"; then
	  enable_shared=yes
	fi
      done
      IFS="$lt_save_ifs"
      ;;
    esac],
    [enable_shared=]AC_ENABLE_SHARED_DEFAULT)
])# AC_ENABLE_SHARED


# AC_DISABLE_SHARED
# -----------------
# set the default shared flag to --disable-shared
AC_DEFUN([AC_DISABLE_SHARED],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_SHARED(no)
])# AC_DISABLE_SHARED


# AC_ENABLE_STATIC([DEFAULT])
# ---------------------------
# implement the --enable-static flag
# DEFAULT is either `yes' or `no'.  If omitted, it defaults to `yes'.
AC_DEFUN([AC_ENABLE_STATIC],
[define([AC_ENABLE_STATIC_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE([static],
    [AC_HELP_STRING([--enable-static@<:@=PKGS@:>@],
	[build static libraries @<:@default=]AC_ENABLE_STATIC_DEFAULT[@:>@])],
    [p=${PACKAGE-default}
    case $enableval in
    yes) enable_static=yes ;;
    no) enable_static=no ;;
    *)
     enable_static=no
      # Look at the argument we got.  We use all the common list separators.
      lt_save_ifs="$IFS"; IFS="${IFS}$PATH_SEPARATOR,"
      for pkg in $enableval; do
	IFS="$lt_save_ifs"
	if test "X$pkg" = "X$p"; then
	  enable_static=yes
	fi
      done
      IFS="$lt_save_ifs"
      ;;
    esac],
    [enable_static=]AC_ENABLE_STATIC_DEFAULT)
])# AC_ENABLE_STATIC


# AC_DISABLE_STATIC
# -----------------
# set the default static flag to --disable-static
AC_DEFUN([AC_DISABLE_STATIC],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_STATIC(no)
])# AC_DISABLE_STATIC


# AC_ENABLE_FAST_INSTALL([DEFAULT])
# ---------------------------------
# implement the --enable-fast-install flag
# DEFAULT is either `yes' or `no'.  If omitted, it defaults to `yes'.
AC_DEFUN([AC_ENABLE_FAST_INSTALL],
[define([AC_ENABLE_FAST_INSTALL_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE([fast-install],
    [AC_HELP_STRING([--enable-fast-install@<:@=PKGS@:>@],
    [optimize for fast installation @<:@default=]AC_ENABLE_FAST_INSTALL_DEFAULT[@:>@])],
    [p=${PACKAGE-default}
    case $enableval in
    yes) enable_fast_install=yes ;;
    no) enable_fast_install=no ;;
    *)
      enable_fast_install=no
      # Look at the argument we got.  We use all the common list separators.
      lt_save_ifs="$IFS"; IFS="${IFS}$PATH_SEPARATOR,"
      for pkg in $enableval; do
	IFS="$lt_save_ifs"
	if test "X$pkg" = "X$p"; then
	  enable_fast_install=yes
	fi
      done
      IFS="$lt_save_ifs"
      ;;
    esac],
    [enable_fast_install=]AC_ENABLE_FAST_INSTALL_DEFAULT)
])# AC_ENABLE_FAST_INSTALL


# AC_DISABLE_FAST_INSTALL
# -----------------------
# set the default to --disable-fast-install
AC_DEFUN([AC_DISABLE_FAST_INSTALL],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_FAST_INSTALL(no)
])# AC_DISABLE_FAST_INSTALL


# AC_LIBTOOL_PICMODE([MODE])
# --------------------------
# implement the --with-pic flag
# MODE is either `yes' or `no'.  If omitted, it defaults to `both'.
AC_DEFUN([AC_LIBTOOL_PICMODE],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
pic_mode=ifelse($#,1,$1,default)
])# AC_LIBTOOL_PICMODE


# AC_PROG_EGREP
# -------------
# This is predefined starting with Autoconf 2.54, so this conditional
# definition can be removed once we require Autoconf 2.54 or later.
m4_ifndef([AC_PROG_EGREP], [AC_DEFUN([AC_PROG_EGREP],
[AC_CACHE_CHECK([for egrep], [ac_cv_prog_egrep],
   [if echo a | (grep -E '(a|b)') >/dev/null 2>&1
    then ac_cv_prog_egrep='grep -E'
    else ac_cv_prog_egrep='egrep'
    fi])
 EGREP=$ac_cv_prog_egrep
 AC_SUBST([EGREP])
])])


# AC_PATH_TOOL_PREFIX
# -------------------
# find a file program which can recognise shared library
AC_DEFUN([AC_PATH_TOOL_PREFIX],
[AC_REQUIRE([AC_PROG_EGREP])dnl
AC_MSG_CHECKING([for $1])
AC_CACHE_VAL(lt_cv_path_MAGIC_CMD,
[case $MAGIC_CMD in
[[\\/*] |  ?:[\\/]*])
  lt_cv_path_MAGIC_CMD="$MAGIC_CMD" # Let the user override the test with a path.
  ;;
*)
  lt_save_MAGIC_CMD="$MAGIC_CMD"
  lt_save_ifs="$IFS"; IFS=$PATH_SEPARATOR
dnl $ac_dummy forces splitting on constant user-supplied paths.
dnl POSIX.2 word splitting is done only on the output of word expansions,
dnl not every word.  This closes a longstanding sh security hole.
  ac_dummy="ifelse([$2], , $PATH, [$2])"
  for ac_dir in $ac_dummy; do
    IFS="$lt_save_ifs"
    test -z "$ac_dir" && ac_dir=.
    if test -f $ac_dir/$1; then
      lt_cv_path_MAGIC_CMD="$ac_dir/$1"
      if test -n "$file_magic_test_file"; then
	case $deplibs_check_method in
	"file_magic "*)
	  file_magic_regex=`expr "$deplibs_check_method" : "file_magic \(.*\)"`
	  MAGIC_CMD="$lt_cv_path_MAGIC_CMD"
	  if eval $file_magic_cmd \$file_magic_test_file 2> /dev/null |
	    $EGREP "$file_magic_regex" > /dev/null; then
	    :
	  else
	    cat <<EOF 1>&2

*** Warning: the command libtool uses to detect shared libraries,
*** $file_magic_cmd, produces output that libtool cannot recognize.
*** The result is that libtool may fail to recognize shared libraries
*** as such.  This will affect the creation of libtool libraries that
*** depend on shared libraries, but programs linked with such libtool
*** libraries will work regardless of this problem.  Nevertheless, you
*** may want to report the problem to your system manager and/or to
*** bug-libtool@gnu.org

EOF
	  fi ;;
	esac
      fi
      break
    fi
  done
  IFS="$lt_save_ifs"
  MAGIC_CMD="$lt_save_MAGIC_CMD"
  ;;
esac])
MAGIC_CMD="$lt_cv_path_MAGIC_CMD"
if test -n "$MAGIC_CMD"; then
  AC_MSG_RESULT($MAGIC_CMD)
else
  AC_MSG_RESULT(no)
fi
])# AC_PATH_TOOL_PREFIX


# AC_PATH_MAGIC
# -------------
# find a file program which can recognise a shared library
AC_DEFUN([AC_PATH_MAGIC],
[AC_PATH_TOOL_PREFIX(${ac_tool_prefix}file, /usr/bin$PATH_SEPARATOR$PATH)
if test -z "$lt_cv_path_MAGIC_CMD"; then
  if test -n "$ac_tool_prefix"; then
    AC_PATH_TOOL_PREFIX(file, /usr/bin$PATH_SEPARATOR$PATH)
  else
    MAGIC_CMD=:
  fi
fi
])# AC_PATH_MAGIC


# AC_PROG_LD
# ----------
# find the pathname to the GNU or non-GNU linker
AC_DEFUN([AC_PROG_LD],
[AC_ARG_WITH([gnu-ld],
    [AC_HELP_STRING([--with-gnu-ld],
	[assume the C compiler uses GNU ld @<:@default=no@:>@])],
    [test "$withval" = no || with_gnu_ld=yes],
    [with_gnu_ld=no])
AC_REQUIRE([LT_AC_PROG_SED])dnl
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
ac_prog=ld
if test "$GCC" = yes; then
  # Check if gcc -print-prog-name=ld gives a path.
  AC_MSG_CHECKING([for ld used by $CC])
  case $host in
  *-*-mingw*)
    # gcc leaves a trailing carriage return which upsets mingw
    ac_prog=`($CC -print-prog-name=ld) 2>&5 | tr -d '\015'` ;;
  *)
    ac_prog=`($CC -print-prog-name=ld) 2>&5` ;;
  esac
  case $ac_prog in
    # Accept absolute paths.
    [[\\/]]* | ?:[[\\/]]*)
      re_direlt='/[[^/]][[^/]]*/\.\./'
      # Canonicalize the pathname of ld
      ac_prog=`echo $ac_prog| $SED 's%\\\\%/%g'`
      while echo $ac_prog | grep "$re_direlt" > /dev/null 2>&1; do
	ac_prog=`echo $ac_prog| $SED "s%$re_direlt%/%"`
      done
      test -z "$LD" && LD="$ac_prog"
      ;;
  "")
    # If it fails, then pretend we aren't using GCC.
    ac_prog=ld
    ;;
  *)
    # If it is relative, then search for the first ld in PATH.
    with_gnu_ld=unknown
    ;;
  esac
elif test "$with_gnu_ld" = yes; then
  AC_MSG_CHECKING([for GNU ld])
else
  AC_MSG_CHECKING([for non-GNU ld])
fi
AC_CACHE_VAL(lt_cv_path_LD,
[if test -z "$LD"; then
  lt_save_ifs="$IFS"; IFS=$PATH_SEPARATOR
  for ac_dir in $PATH; do
    IFS="$lt_save_ifs"
    test -z "$ac_dir" && ac_dir=.
    if test -f "$ac_dir/$ac_prog" || test -f "$ac_dir/$ac_prog$ac_exeext"; then
      lt_cv_path_LD="$ac_dir/$ac_prog"
      # Check to see if the program is GNU ld.  I'd rather use --version,
      # but apparently some variants of GNU ld only accept -v.
      # Break only if it was the GNU/non-GNU ld that we prefer.
      case `"$lt_cv_path_LD" -v 2>&1 </dev/null` in
      *GNU* | *'with BFD'*)
	test "$with_gnu_ld" != no && break
	;;
      *)
	test "$with_gnu_ld" != yes && break
	;;
      esac
    fi
  done
  IFS="$lt_save_ifs"
else
  lt_cv_path_LD="$LD" # Let the user override the test with a path.
fi])
LD="$lt_cv_path_LD"
if test -n "$LD"; then
  AC_MSG_RESULT($LD)
else
  AC_MSG_RESULT(no)
fi
test -z "$LD" && AC_MSG_ERROR([no acceptable ld found in \$PATH])
AC_PROG_LD_GNU
])# AC_PROG_LD


# AC_PROG_LD_GNU
# --------------
AC_DEFUN([AC_PROG_LD_GNU],
[AC_REQUIRE([AC_PROG_EGREP])dnl
AC_CACHE_CHECK([if the linker ($LD) is GNU ld], lt_cv_prog_gnu_ld,
[# I'd rather use --version here, but apparently some GNU lds only accept -v.
case `$LD -v 2>&1 </dev/null` in
*GNU* | *'with BFD'*)
  lt_cv_prog_gnu_ld=yes
  ;;
*)
  lt_cv_prog_gnu_ld=no
  ;;
esac])
with_gnu_ld=$lt_cv_prog_gnu_ld
])# AC_PROG_LD_GNU


# AC_PROG_LD_RELOAD_FLAG
# ----------------------
# find reload flag for linker
#   -- PORTME Some linkers may need a different reload flag.
AC_DEFUN([AC_PROG_LD_RELOAD_FLAG],
[AC_CACHE_CHECK([for $LD option to reload object files],
  lt_cv_ld_reload_flag,
  [lt_cv_ld_reload_flag='-r'])
reload_flag=$lt_cv_ld_reload_flag
case $reload_flag in
"" | " "*) ;;
*) reload_flag=" $reload_flag" ;;
esac
reload_cmds='$LD$reload_flag -o $output$reload_objs'
case $host_os in
  darwin*)
    if test "$GCC" = yes; then
      reload_cmds='$LTCC $LTCFLAGS -nostdlib ${wl}-r -o $output$reload_objs'
    else
      reload_cmds='$LD$reload_flag -o $output$reload_objs'
    fi
    ;;
esac
])# AC_PROG_LD_RELOAD_FLAG


# AC_DEPLIBS_CHECK_METHOD
# -----------------------
# how to check for library dependencies
#  -- PORTME fill in with the dynamic library characteristics
AC_DEFUN([AC_DEPLIBS_CHECK_METHOD],
[AC_CACHE_CHECK([how to recognise dependent libraries],
lt_cv_deplibs_check_method,
[lt_cv_file_magic_cmd='$MAGIC_CMD'
lt_cv_file_magic_test_file=
lt_cv_deplibs_check_method='unknown'
# Need to set the preceding variable on all platforms that support
# interlibrary dependencies.
# 'none' -- dependencies not supported.
# `unknown' -- same as none, but documents that we really don't know.
# 'pass_all' -- all dependencies passed with no checks.
# 'test_compile' -- check by making test program.
# 'file_magic [[regex]]' -- check by looking for files in library path
# which responds to the $file_magic_cmd with a given extended regex.
# If you have `file' or equivalent on your system and you're not sure
# whether `pass_all' will *always* work, you probably want this one.

case $host_os in
aix4* | aix5*)
  lt_cv_deplibs_check_method=pass_all
  ;;

beos*)
  lt_cv_deplibs_check_method=pass_all
  ;;

bsdi[[45]]*)
  lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[ML]]SB (shared object|dynamic lib)'
  lt_cv_file_magic_cmd='/usr/bin/file -L'
  lt_cv_file_magic_test_file=/shlib/libc.so
  ;;

cygwin*)
  # func_win32_libid is a shell function defined in ltmain.sh
  lt_cv_deplibs_check_method='file_magic ^x86 archive import|^x86 DLL'
  lt_cv_file_magic_cmd='func_win32_libid'
  ;;

mingw* | pw32*)
  # Base MSYS/MinGW do not provide the 'file' command needed by
  # func_win32_libid shell function, so use a weaker test based on 'objdump'.
  lt_cv_deplibs_check_method='file_magic file format pei*-i386(.*architecture: i386)?'
  lt_cv_file_magic_cmd='$OBJDUMP -f'
  ;;

darwin* | rhapsody*)
  lt_cv_deplibs_check_method=pass_all
  ;;

freebsd* | kfreebsd*-gnu | dragonfly*)
  if echo __ELF__ | $CC -E - | grep __ELF__ > /dev/null; then
    case $host_cpu in
    i*86 )
      # Not sure whether the presence of OpenBSD here was a mistake.
      # Let's accept both of them until this is cleared up.
      lt_cv_deplibs_check_method='file_magic (FreeBSD|OpenBSD|DragonFly)/i[[3-9]]86 (compact )?demand paged shared library'
      lt_cv_file_magic_cmd=/usr/bin/file
      lt_cv_file_magic_test_file=`echo /usr/lib/libc.so.*`
      ;;
    esac
  else
    lt_cv_deplibs_check_method=pass_all
  fi
  ;;

gnu*)
  lt_cv_deplibs_check_method=pass_all
  ;;

hpux10.20* | hpux11*)
  lt_cv_file_magic_cmd=/usr/bin/file
  case $host_cpu in
  ia64*)
    lt_cv_deplibs_check_method='file_magic (s[[0-9]][[0-9]][[0-9]]|ELF-[[0-9]][[0-9]]) shared object file - IA64'
    lt_cv_file_magic_test_file=/usr/lib/hpux32/libc.so
    ;;
  hppa*64*)
    [lt_cv_deplibs_check_method='file_magic (s[0-9][0-9][0-9]|ELF-[0-9][0-9]) shared object file - PA-RISC [0-9].[0-9]']
    lt_cv_file_magic_test_file=/usr/lib/pa20_64/libc.sl
    ;;
  *)
    lt_cv_deplibs_check_method='file_magic (s[[0-9]][[0-9]][[0-9]]|PA-RISC[[0-9]].[[0-9]]) shared library'
    lt_cv_file_magic_test_file=/usr/lib/libc.sl
    ;;
  esac
  ;;

interix3*)
  # PIC code is broken on Interix 3.x, that's why |\.a not |_pic\.a here
  lt_cv_deplibs_check_method='match_pattern /lib[[^/]]+(\.so|\.a)$'
  ;;

irix5* | irix6* | nonstopux*)
  case $LD in
  *-32|*"-32 ") libmagic=32-bit;;
  *-n32|*"-n32 ") libmagic=N32;;
  *-64|*"-64 ") libmagic=64-bit;;
  *) libmagic=never-match;;
  esac
  lt_cv_deplibs_check_method=pass_all
  ;;

# This must be Linux ELF.
linux*)
  lt_cv_deplibs_check_method=pass_all
  ;;

netbsd*)
  if echo __ELF__ | $CC -E - | grep __ELF__ > /dev/null; then
    lt_cv_deplibs_check_method='match_pattern /lib[[^/]]+(\.so\.[[0-9]]+\.[[0-9]]+|_pic\.a)$'
  else
    lt_cv_deplibs_check_method='match_pattern /lib[[^/]]+(\.so|_pic\.a)$'
  fi
  ;;

newos6*)
  lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[ML]]SB (executable|dynamic lib)'
  lt_cv_file_magic_cmd=/usr/bin/file
  lt_cv_file_magic_test_file=/usr/lib/libnls.so
  ;;

nto-qnx*)
  lt_cv_deplibs_check_method=unknown
  ;;

openbsd*)
  if test -z "`echo __ELF__ | $CC -E - | grep __ELF__`" || test "$host_os-$host_cpu" = "openbsd2.8-powerpc"; then
    lt_cv_deplibs_check_method='match_pattern /lib[[^/]]+(\.so\.[[0-9]]+\.[[0-9]]+|\.so|_pic\.a)$'
  else
    lt_cv_deplibs_check_method='match_pattern /lib[[^/]]+(\.so\.[[0-9]]+\.[[0-9]]+|_pic\.a)$'
  fi
  ;;

osf3* | osf4* | osf5*)
  lt_cv_deplibs_check_method=pass_all
  ;;

solaris*)
  lt_cv_deplibs_check_method=pass_all
  ;;

sysv4 | sysv4.3*)
  case $host_vendor in
  motorola)
    lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[ML]]SB (shared object|dynamic lib) M[[0-9]][[0-9]]* Version [[0-9]]'
    lt_cv_file_magic_test_file=`echo /usr/lib/libc.so*`
    ;;
  ncr)
    lt_cv_deplibs_check_method=pass_all
    ;;
  sequent)
    lt_cv_file_magic_cmd='/bin/file'
    lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[LM]]SB (shared object|dynamic lib )'
    ;;
  sni)
    lt_cv_file_magic_cmd='/bin/file'
    lt_cv_deplibs_check_method="file_magic ELF [[0-9]][[0-9]]*-bit [[LM]]SB dynamic lib"
    lt_cv_file_magic_test_file=/lib/libc.so
    ;;
  siemens)
    lt_cv_deplibs_check_method=pass_all
    ;;
  pc)
    lt_cv_deplibs_check_method=pass_all
    ;;
  esac
  ;;

sysv5* | sco3.2v5* | sco5v6* | unixware* | OpenUNIX* | sysv4*uw2*)
  lt_cv_deplibs_check_method=pass_all
  ;;
esac
])
file_magic_cmd=$lt_cv_file_magic_cmd
deplibs_check_method=$lt_cv_deplibs_check_method
test -z "$deplibs_check_method" && deplibs_check_method=unknown
])# AC_DEPLIBS_CHECK_METHOD


# AC_PROG_NM
# ----------
# find the pathname to a BSD-compatible name lister
AC_DEFUN([AC_PROG_NM],
[AC_CACHE_CHECK([for BSD-compatible nm], lt_cv_path_NM,
[if test -n "$NM"; then
  # Let the user override the test.
  lt_cv_path_NM="$NM"
else
  lt_nm_to_check="${ac_tool_prefix}nm"
  if test -n "$ac_tool_prefix" && test "$build" = "$host"; then 
    lt_nm_to_check="$lt_nm_to_check nm"
  fi
  for lt_tmp_nm in $lt_nm_to_check; do
    lt_save_ifs="$IFS"; IFS=$PATH_SEPARATOR
    for ac_dir in $PATH /usr/ccs/bin/elf /usr/ccs/bin /usr/ucb /bin; do
      IFS="$lt_save_ifs"
      test -z "$ac_dir" && ac_dir=.
      tmp_nm="$ac_dir/$lt_tmp_nm"
      if test -f "$tmp_nm" || test -f "$tmp_nm$ac_exeext" ; then
	# Check to see if the nm accepts a BSD-compat flag.
	# Adding the `sed 1q' prevents false positives on HP-UX, which says:
	#   nm: unknown option "B" ignored
	# Tru64's nm complains that /dev/null is an invalid object file
	case `"$tmp_nm" -B /dev/null 2>&1 | sed '1q'` in
	*/dev/null* | *'Invalid file or object type'*)
	  lt_cv_path_NM="$tmp_nm -B"
	  break
	  ;;
	*)
	  case `"$tmp_nm" -p /dev/null 2>&1 | sed '1q'` in
	  */dev/null*)
	    lt_cv_path_NM="$tmp_nm -p"
	    break
	    ;;
	  *)
	    lt_cv_path_NM=${lt_cv_path_NM="$tmp_nm"} # keep the first match, but
	    continue # so that we can try to find one that supports BSD flags
	    ;;
	  esac
	  ;;
	esac
      fi
    done
    IFS="$lt_save_ifs"
  done
  test -z "$lt_cv_path_NM" && lt_cv_path_NM=nm
fi])
NM="$lt_cv_path_NM"
])# AC_PROG_NM


# AC_CHECK_LIBM
# -------------
# check for math library
AC_DEFUN([AC_CHECK_LIBM],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
LIBM=
case $host in
*-*-beos* | *-*-cygwin* | *-*-pw32* | *-*-darwin*)
  # These system don't have libm, or don't need it
  ;;
*-ncr-sysv4.3*)
  AC_CHECK_LIB(mw, _mwvalidcheckl, LIBM="-lmw")
  AC_CHECK_LIB(m, cos, LIBM="$LIBM -lm")
  ;;
*)
  AC_CHECK_LIB(m, cos, LIBM="-lm")
  ;;
esac
])# AC_CHECK_LIBM


# AC_LIBLTDL_CONVENIENCE([DIRECTORY])
# -----------------------------------
# sets LIBLTDL to the link flags for the libltdl convenience library and
# LTDLINCL to the include flags for the libltdl header and adds
# --enable-ltdl-convenience to the configure arguments.  Note that
# AC_CONFIG_SUBDIRS is not called here.  If DIRECTORY is not provided,
# it is assumed to be `libltdl'.  LIBLTDL will be prefixed with
# '${top_builddir}/' and LTDLINCL will be prefixed with '${top_srcdir}/'
# (note the single quotes!).  If your package is not flat and you're not
# using automake, define top_builddir and top_srcdir appropriately in
# the Makefiles.
AC_DEFUN([AC_LIBLTDL_CONVENIENCE],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
  case $enable_ltdl_convenience in
  no) AC_MSG_ERROR([this package needs a convenience libltdl]) ;;
  "") enable_ltdl_convenience=yes
      ac_configure_args="$ac_configure_args --enable-ltdl-convenience" ;;
  esac
  LIBLTDL='${top_builddir}/'ifelse($#,1,[$1],['libltdl'])/libltdlc.la
  LTDLINCL='-I${top_srcdir}/'ifelse($#,1,[$1],['libltdl'])
  # For backwards non-gettext consistent compatibility...
  INCLTDL="$LTDLINCL"
])# AC_LIBLTDL_CONVENIENCE


# AC_LIBLTDL_INSTALLABLE([DIRECTORY])
# -----------------------------------
# sets LIBLTDL to the link flags for the libltdl installable library and
# LTDLINCL to the include flags for the libltdl header and adds
# --enable-ltdl-install to the configure arguments.  Note that
# AC_CONFIG_SUBDIRS is not called here.  If DIRECTORY is not provided,
# and an installed libltdl is not found, it is assumed to be `libltdl'.
# LIBLTDL will be prefixed with '${top_builddir}/'# and LTDLINCL with
# '${top_srcdir}/' (note the single quotes!).  If your package is not
# flat and you're not using automake, define top_builddir and top_srcdir
# appropriately in the Makefiles.
# In the future, this macro may have to be called after AC_PROG_LIBTOOL.
AC_DEFUN([AC_LIBLTDL_INSTALLABLE],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
  AC_CHECK_LIB(ltdl, lt_dlinit,
  [test x"$enable_ltdl_install" != xyes && enable_ltdl_install=no],
  [if test x"$enable_ltdl_install" = xno; then
     AC_MSG_WARN([libltdl not installed, but installation disabled])
   else
     enable_ltdl_install=yes
   fi
  ])
  if test x"$enable_ltdl_install" = x"yes"; then
    ac_configure_args="$ac_configure_args --enable-ltdl-install"
    LIBLTDL='${top_builddir}/'ifelse($#,1,[$1],['libltdl'])/libltdl.la
    LTDLINCL='-I${top_srcdir}/'ifelse($#,1,[$1],['libltdl'])
  else
    ac_configure_args="$ac_configure_args --enable-ltdl-install=no"
    LIBLTDL="-lltdl"
    LTDLINCL=
  fi
  # For backwards non-gettext consistent compatibility...
  INCLTDL="$LTDLINCL"
])# AC_LIBLTDL_INSTALLABLE


# AC_LIBTOOL_CXX
# --------------
# enable support for C++ libraries
AC_DEFUN([AC_LIBTOOL_CXX],
[AC_REQUIRE([_LT_AC_LANG_CXX])
])# AC_LIBTOOL_CXX


# _LT_AC_LANG_CXX
# ---------------
AC_DEFUN([_LT_AC_LANG_CXX],
[AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([_LT_AC_PROG_CXXCPP])
_LT_AC_SHELL_INIT([tagnames=${tagnames+${tagnames},}CXX])
])# _LT_AC_LANG_CXX

# _LT_AC_PROG_CXXCPP
# ------------------
AC_DEFUN([_LT_AC_PROG_CXXCPP],
[
AC_REQUIRE([AC_PROG_CXX])
if test -n "$CXX" && ( test "X$CXX" != "Xno" &&
    ( (test "X$CXX" = "Xg++" && `g++ -v >/dev/null 2>&1` ) ||
    (test "X$CXX" != "Xg++"))) ; then
  AC_PROG_CXXCPP
fi
])# _LT_AC_PROG_CXXCPP

# AC_LIBTOOL_F77
# --------------
# enable support for Fortran 77 libraries
AC_DEFUN([AC_LIBTOOL_F77],
[AC_REQUIRE([_LT_AC_LANG_F77])
])# AC_LIBTOOL_F77


# _LT_AC_LANG_F77
# ---------------
AC_DEFUN([_LT_AC_LANG_F77],
[AC_REQUIRE([AC_PROG_F77])
_LT_AC_SHELL_INIT([tagnames=${tagnames+${tagnames},}F77])
])# _LT_AC_LANG_F77


# AC_LIBTOOL_GCJ
# --------------
# enable support for GCJ libraries
AC_DEFUN([AC_LIBTOOL_GCJ],
[AC_REQUIRE([_LT_AC_LANG_GCJ])
])# AC_LIBTOOL_GCJ


# _LT_AC_LANG_GCJ
# ---------------
AC_DEFUN([_LT_AC_LANG_GCJ],
[AC_PROVIDE_IFELSE([AC_PROG_GCJ],[],
  [AC_PROVIDE_IFELSE([A][M_PROG_GCJ],[],
    [AC_PROVIDE_IFELSE([LT_AC_PROG_GCJ],[],
      [ifdef([AC_PROG_GCJ],[AC_REQUIRE([AC_PROG_GCJ])],
	 [ifdef([A][M_PROG_GCJ],[AC_REQUIRE([A][M_PROG_GCJ])],
	   [AC_REQUIRE([A][C_PROG_GCJ_OR_A][M_PROG_GCJ])])])])])])
_LT_AC_SHELL_INIT([tagnames=${tagnames+${tagnames},}GCJ])
])# _LT_AC_LANG_GCJ


# AC_LIBTOOL_RC
# -------------
# enable support for Windows resource files
AC_DEFUN([AC_LIBTOOL_RC],
[AC_REQUIRE([LT_AC_PROG_RC])
_LT_AC_SHELL_INIT([tagnames=${tagnames+${tagnames},}RC])
])# AC_LIBTOOL_RC


# AC_LIBTOOL_LANG_C_CONFIG
# ------------------------
# Ensure that the configuration vars for the C compiler are
# suitably defined.  Those variables are subsequently used by
# AC_LIBTOOL_CONFIG to write the compiler configuration to `libtool'.
AC_DEFUN([AC_LIBTOOL_LANG_C_CONFIG], [_LT_AC_LANG_C_CONFIG])
AC_DEFUN([_LT_AC_LANG_C_CONFIG],
[lt_save_CC="$CC"
AC_LANG_PUSH(C)

# Source file extension for C test sources.
ac_ext=c

# Object file extension for compiled C test sources.
objext=o
_LT_AC_TAGVAR(objext, $1)=$objext

# Code to be used in simple compile tests
lt_simple_compile_test_code="int some_variable = 0;\n"

# Code to be used in simple link tests
lt_simple_link_test_code='int main(){return(0);}\n'

_LT_AC_SYS_COMPILER

# save warnings/boilerplate of simple test code
_LT_COMPILER_BOILERPLATE
_LT_LINKER_BOILERPLATE

## CAVEAT EMPTOR:
## There is no encapsulation within the following macros, do not change
## the running order or otherwise move them around unless you know exactly
## what you are doing...
AC_LIBTOOL_PROG_COMPILER_NO_RTTI($1)
AC_LIBTOOL_PROG_COMPILER_PIC($1)
AC_LIBTOOL_PROG_CC_C_O($1)
AC_LIBTOOL_SYS_HARD_LINK_LOCKS($1)
AC_LIBTOOL_PROG_LD_SHLIBS($1)
AC_LIBTOOL_SYS_DYNAMIC_LINKER($1)
AC_LIBTOOL_PROG_LD_HARDCODE_LIBPATH($1)
AC_LIBTOOL_SYS_LIB_STRIP
AC_LIBTOOL_DLOPEN_SELF

# Report which library types will actually be built
AC_MSG_CHECKING([if libtool supports shared libraries])
AC_MSG_RESULT([$can_build_shared])

AC_MSG_CHECKING([whether to build shared libraries])
test "$can_build_shared" = "no" && enable_shared=no

# On AIX, shared libraries and static libraries use the same namespace, and
# are all built from PIC.
case $host_os in
aix3*)
  test "$enable_shared" = yes && enable_static=no
  if test -n "$RANLIB"; then
    archive_cmds="$archive_cmds~\$RANLIB \$lib"
    postinstall_cmds='$RANLIB $lib'
  fi
  ;;

aix4* | aix5*)
  if test "$host_cpu" != ia64 && test "$aix_use_runtimelinking" = no ; then
    test "$enable_shared" = yes && enable_static=no
  fi
    ;;
esac
AC_MSG_RESULT([$enable_shared])

AC_MSG_CHECKING([whether to build static libraries])
# Make sure either enable_shared or enable_static is yes.
test "$enable_shared" = yes || enable_static=yes
AC_MSG_RESULT([$enable_static])

AC_LIBTOOL_CONFIG($1)

AC_LANG_POP
CC="$lt_save_CC"
])# AC_LIBTOOL_LANG_C_CONFIG


# AC_LIBTOOL_LANG_CXX_CONFIG
# --------------------------
# Ensure that the configuration vars for the C compiler are
# suitably defined.  Those variables are subsequently used by
# AC_LIBTOOL_CONFIG to write the compiler configuration to `libtool'.
AC_DEFUN([AC_LIBTOOL_LANG_CXX_CONFIG], [_LT_AC_LANG_CXX_CONFIG(CXX)])
AC_DEFUN([_LT_AC_LANG_CXX_CONFIG],
[AC_LANG_PUSH(C++)
AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([_LT_AC_PROG_CXXCPP])

_LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
_LT_AC_TAGVAR(allow_undefined_flag, $1)=
_LT_AC_TAGVAR(always_export_symbols, $1)=no
_LT_AC_TAGVAR(archive_expsym_cmds, $1)=
_LT_AC_TAGVAR(export_dynamic_flag_spec, $1)=
_LT_AC_TAGVAR(hardcode_direct, $1)=no
_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)=
_LT_AC_TAGVAR(hardcode_libdir_flag_spec_ld, $1)=
_LT_AC_TAGVAR(hardcode_libdir_separator, $1)=
_LT_AC_TAGVAR(hardcode_minus_L, $1)=no
_LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=unsupported
_LT_AC_TAGVAR(hardcode_automatic, $1)=no
_LT_AC_TAGVAR(module_cmds, $1)=
_LT_AC_TAGVAR(module_expsym_cmds, $1)=
_LT_AC_TAGVAR(link_all_deplibs, $1)=unknown
_LT_AC_TAGVAR(old_archive_cmds, $1)=$old_archive_cmds
_LT_AC_TAGVAR(no_undefined_flag, $1)=
_LT_AC_TAGVAR(whole_archive_flag_spec, $1)=
_LT_AC_TAGVAR(enable_shared_with_static_runtimes, $1)=no

# Dependencies to place before and after the object being linked:
_LT_AC_TAGVAR(predep_objects, $1)=
_LT_AC_TAGVAR(postdep_objects, $1)=
_LT_AC_TAGVAR(predeps, $1)=
_LT_AC_TAGVAR(postdeps, $1)=
_LT_AC_TAGVAR(compiler_lib_search_path, $1)=

# Source file extension for C++ test sources.
ac_ext=cpp

# Object file extension for compiled C++ test sources.
objext=o
_LT_AC_TAGVAR(objext, $1)=$objext

# Code to be used in simple compile tests
lt_simple_compile_test_code="int some_variable = 0;\n"

# Code to be used in simple link tests
lt_simple_link_test_code='int main(int, char *[[]]) { return(0); }\n'

# ltmain only uses $CC for tagged configurations so make sure $CC is set.
_LT_AC_SYS_COMPILER

# save warnings/boilerplate of simple test code
_LT_COMPILER_BOILERPLATE
_LT_LINKER_BOILERPLATE

# Allow CC to be a program name with arguments.
lt_save_CC=$CC
lt_save_LD=$LD
lt_save_GCC=$GCC
GCC=$GXX
lt_save_with_gnu_ld=$with_gnu_ld
lt_save_path_LD=$lt_cv_path_LD
if test -n "${lt_cv_prog_gnu_ldcxx+set}"; then
  lt_cv_prog_gnu_ld=$lt_cv_prog_gnu_ldcxx
else
  $as_unset lt_cv_prog_gnu_ld
fi
if test -n "${lt_cv_path_LDCXX+set}"; then
  lt_cv_path_LD=$lt_cv_path_LDCXX
else
  $as_unset lt_cv_path_LD
fi
test -z "${LDCXX+set}" || LD=$LDCXX
CC=${CXX-"c++"}
compiler=$CC
_LT_AC_TAGVAR(compiler, $1)=$CC
_LT_CC_BASENAME([$compiler])

# We don't want -fno-exception wen compiling C++ code, so set the
# no_builtin_flag separately
if test "$GXX" = yes; then
  _LT_AC_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)=' -fno-builtin'
else
  _LT_AC_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)=
fi

if test "$GXX" = yes; then
  # Set up default GNU C++ configuration

  AC_PROG_LD

  # Check if GNU C++ uses GNU ld as the underlying linker, since the
  # archiving commands below assume that GNU ld is being used.
  if test "$with_gnu_ld" = yes; then
    _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname $wl$soname -o $lib'
    _LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname $wl$soname ${wl}-retain-symbols-file $wl$export_symbols -o $lib'

    _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}--rpath ${wl}$libdir'
    _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}--export-dynamic'

    # If archive_cmds runs LD, not CC, wlarc should be empty
    # XXX I think wlarc can be eliminated in ltcf-cxx, but I need to
    #     investigate it a little bit more. (MM)
    wlarc='${wl}'

    # ancient GNU ld didn't support --whole-archive et. al.
    if eval "`$CC -print-prog-name=ld` --help 2>&1" | \
	grep 'no-whole-archive' > /dev/null; then
      _LT_AC_TAGVAR(whole_archive_flag_spec, $1)="$wlarc"'--whole-archive$convenience '"$wlarc"'--no-whole-archive'
    else
      _LT_AC_TAGVAR(whole_archive_flag_spec, $1)=
    fi
  else
    with_gnu_ld=no
    wlarc=

    # A generic and very simple default shared library creation
    # command for GNU C++ for the case where it uses the native
    # linker, instead of GNU ld.  If possible, this setting should
    # overridden to take advantage of the native linker features on
    # the platform it is being used on.
    _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -o $lib'
  fi

  # Commands to make compiler produce verbose output that lists
  # what "hidden" libraries, object files and flags are used when
  # linking a shared library.
  output_verbose_link_cmd='$CC -shared $CFLAGS -v conftest.$objext 2>&1 | grep "\-L"'

else
  GXX=no
  with_gnu_ld=no
  wlarc=
fi

# PORTME: fill in a description of your system's C++ link characteristics
AC_MSG_CHECKING([whether the $compiler linker ($LD) supports shared libraries])
_LT_AC_TAGVAR(ld_shlibs, $1)=yes
case $host_os in
  aix3*)
    # FIXME: insert proper C++ library support
    _LT_AC_TAGVAR(ld_shlibs, $1)=no
    ;;
  aix4* | aix5*)
    if test "$host_cpu" = ia64; then
      # On IA64, the linker does run time linking by default, so we don't
      # have to do anything special.
      aix_use_runtimelinking=no
      exp_sym_flag='-Bexport'
      no_entry_flag=""
    else
      aix_use_runtimelinking=no

      # Test if we are trying to use run time linking or normal
      # AIX style linking. If -brtl is somewhere in LDFLAGS, we
      # need to do runtime linking.
      case $host_os in aix4.[[23]]|aix4.[[23]].*|aix5*)
	for ld_flag in $LDFLAGS; do
	  case $ld_flag in
	  *-brtl*)
	    aix_use_runtimelinking=yes
	    break
	    ;;
	  esac
	done
	;;
      esac

      exp_sym_flag='-bexport'
      no_entry_flag='-bnoentry'
    fi

    # When large executables or shared objects are built, AIX ld can
    # have problems creating the table of contents.  If linking a library
    # or program results in "error TOC overflow" add -mminimal-toc to
    # CXXFLAGS/CFLAGS for g++/gcc.  In the cases where that is not
    # enough to fix the problem, add -Wl,-bbigtoc to LDFLAGS.

    _LT_AC_TAGVAR(archive_cmds, $1)=''
    _LT_AC_TAGVAR(hardcode_direct, $1)=yes
    _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=':'
    _LT_AC_TAGVAR(link_all_deplibs, $1)=yes

    if test "$GXX" = yes; then
      case $host_os in aix4.[[012]]|aix4.[[012]].*)
      # We only want to do this on AIX 4.2 and lower, the check
      # below for broken collect2 doesn't work under 4.3+
	collect2name=`${CC} -print-prog-name=collect2`
	if test -f "$collect2name" && \
	   strings "$collect2name" | grep resolve_lib_name >/dev/null
	then
	  # We have reworked collect2
	  _LT_AC_TAGVAR(hardcode_direct, $1)=yes
	else
	  # We have old collect2
	  _LT_AC_TAGVAR(hardcode_direct, $1)=unsupported
	  # It fails to find uninstalled libraries when the uninstalled
	  # path is not listed in the libpath.  Setting hardcode_minus_L
	  # to unsupported forces relinking
	  _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes
	  _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
	  _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=
	fi
	;;
      esac
      shared_flag='-shared'
      if test "$aix_use_runtimelinking" = yes; then
	shared_flag="$shared_flag "'${wl}-G'
      fi
    else
      # not using gcc
      if test "$host_cpu" = ia64; then
	# VisualAge C++, Version 5.5 for AIX 5L for IA-64, Beta 3 Release
	# chokes on -Wl,-G. The following line is correct:
	shared_flag='-G'
      else
	if test "$aix_use_runtimelinking" = yes; then
	  shared_flag='${wl}-G'
	else
	  shared_flag='${wl}-bM:SRE'
	fi
      fi
    fi

    # It seems that -bexpall does not export symbols beginning with
    # underscore (_), so it is better to generate a list of symbols to export.
    _LT_AC_TAGVAR(always_export_symbols, $1)=yes
    if test "$aix_use_runtimelinking" = yes; then
      # Warning - without using the other runtime loading flags (-brtl),
      # -berok will link without error, but may produce a broken library.
      _LT_AC_TAGVAR(allow_undefined_flag, $1)='-berok'
      # Determine the default libpath from the value encoded in an empty executable.
      _LT_AC_SYS_LIBPATH_AIX
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-blibpath:$libdir:'"$aix_libpath"

      _LT_AC_TAGVAR(archive_expsym_cmds, $1)="\$CC"' -o $output_objdir/$soname $libobjs $deplibs '"\${wl}$no_entry_flag"' $compiler_flags `if test "x${allow_undefined_flag}" != "x"; then echo "${wl}${allow_undefined_flag}"; else :; fi` '"\${wl}$exp_sym_flag:\$export_symbols $shared_flag"
     else
      if test "$host_cpu" = ia64; then
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-R $libdir:/usr/lib:/lib'
	_LT_AC_TAGVAR(allow_undefined_flag, $1)="-z nodefs"
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)="\$CC $shared_flag"' -o $output_objdir/$soname $libobjs $deplibs '"\${wl}$no_entry_flag"' $compiler_flags ${wl}${allow_undefined_flag} '"\${wl}$exp_sym_flag:\$export_symbols"
      else
	# Determine the default libpath from the value encoded in an empty executable.
	_LT_AC_SYS_LIBPATH_AIX
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-blibpath:$libdir:'"$aix_libpath"
	# Warning - without using the other run time loading flags,
	# -berok will link without error, but may produce a broken library.
	_LT_AC_TAGVAR(no_undefined_flag, $1)=' ${wl}-bernotok'
	_LT_AC_TAGVAR(allow_undefined_flag, $1)=' ${wl}-berok'
	# Exported symbols can be pulled into shared objects from archives
	_LT_AC_TAGVAR(whole_archive_flag_spec, $1)='$convenience'
	_LT_AC_TAGVAR(archive_cmds_need_lc, $1)=yes
	# This is similar to how AIX traditionally builds its shared libraries.
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)="\$CC $shared_flag"' -o $output_objdir/$soname $libobjs $deplibs ${wl}-bnoentry $compiler_flags ${wl}-bE:$export_symbols${allow_undefined_flag}~$AR $AR_FLAGS $output_objdir/$libname$release.a $output_objdir/$soname'
      fi
    fi
    ;;

  beos*)
    if $LD --help 2>&1 | grep ': supported targets:.* elf' > /dev/null; then
      _LT_AC_TAGVAR(allow_undefined_flag, $1)=unsupported
      # Joseph Beckenbach <jrb3@best.com> says some releases of gcc
      # support --undefined.  This deserves some investigation.  FIXME
      _LT_AC_TAGVAR(archive_cmds, $1)='$CC -nostart $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'
    else
      _LT_AC_TAGVAR(ld_shlibs, $1)=no
    fi
    ;;

  chorus*)
    case $cc_basename in
      *)
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
    esac
    ;;

  cygwin* | mingw* | pw32*)
    # _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1) is actually meaningless,
    # as there is no search path for DLLs.
    _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
    _LT_AC_TAGVAR(allow_undefined_flag, $1)=unsupported
    _LT_AC_TAGVAR(always_export_symbols, $1)=no
    _LT_AC_TAGVAR(enable_shared_with_static_runtimes, $1)=yes

    if $LD --help 2>&1 | grep 'auto-import' > /dev/null; then
      _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -o $output_objdir/$soname ${wl}--enable-auto-image-base -Xlinker --out-implib -Xlinker $lib'
      # If the export-symbols file already is a .def file (1st line
      # is EXPORTS), use it as is; otherwise, prepend...
      _LT_AC_TAGVAR(archive_expsym_cmds, $1)='if test "x`$SED 1q $export_symbols`" = xEXPORTS; then
	cp $export_symbols $output_objdir/$soname.def;
      else
	echo EXPORTS > $output_objdir/$soname.def;
	cat $export_symbols >> $output_objdir/$soname.def;
      fi~
      $CC -shared -nostdlib $output_objdir/$soname.def $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -o $output_objdir/$soname ${wl}--enable-auto-image-base -Xlinker --out-implib -Xlinker $lib'
    else
      _LT_AC_TAGVAR(ld_shlibs, $1)=no
    fi
  ;;
      darwin* | rhapsody*)
        case $host_os in
        rhapsody* | darwin1.[[012]])
         _LT_AC_TAGVAR(allow_undefined_flag, $1)='${wl}-undefined ${wl}suppress'
         ;;
       *) # Darwin 1.3 on
         if test -z ${MACOSX_DEPLOYMENT_TARGET} ; then
           _LT_AC_TAGVAR(allow_undefined_flag, $1)='${wl}-flat_namespace ${wl}-undefined ${wl}suppress'
         else
           case ${MACOSX_DEPLOYMENT_TARGET} in
             10.[[012]])
               _LT_AC_TAGVAR(allow_undefined_flag, $1)='${wl}-flat_namespace ${wl}-undefined ${wl}suppress'
               ;;
             10.*)
               _LT_AC_TAGVAR(allow_undefined_flag, $1)='${wl}-undefined ${wl}dynamic_lookup'
               ;;
           esac
         fi
         ;;
        esac
      _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
      _LT_AC_TAGVAR(hardcode_direct, $1)=no
      _LT_AC_TAGVAR(hardcode_automatic, $1)=yes
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=unsupported
      _LT_AC_TAGVAR(whole_archive_flag_spec, $1)=''
      _LT_AC_TAGVAR(link_all_deplibs, $1)=yes

    if test "$GXX" = yes ; then
      lt_int_apple_cc_single_mod=no
      output_verbose_link_cmd='echo'
      if $CC -dumpspecs 2>&1 | $EGREP 'single_module' >/dev/null ; then
       lt_int_apple_cc_single_mod=yes
      fi
      if test "X$lt_int_apple_cc_single_mod" = Xyes ; then
       _LT_AC_TAGVAR(archive_cmds, $1)='$CC -dynamiclib -single_module $allow_undefined_flag -o $lib $libobjs $deplibs $compiler_flags -install_name $rpath/$soname $verstring'
      else
          _LT_AC_TAGVAR(archive_cmds, $1)='$CC -r -keep_private_externs -nostdlib -o ${lib}-master.o $libobjs~$CC -dynamiclib $allow_undefined_flag -o $lib ${lib}-master.o $deplibs $compiler_flags -install_name $rpath/$soname $verstring'
        fi
        _LT_AC_TAGVAR(module_cmds, $1)='$CC $allow_undefined_flag -o $lib -bundle $libobjs $deplibs$compiler_flags'
        # Don't fix this by using the ld -exported_symbols_list flag, it doesn't exist in older darwin lds
          if test "X$lt_int_apple_cc_single_mod" = Xyes ; then
            _LT_AC_TAGVAR(archive_expsym_cmds, $1)='sed -e "s,#.*,," -e "s,^[    ]*,," -e "s,^\(..*\),_&," < $export_symbols > $output_objdir/${libname}-symbols.expsym~$CC -dynamiclib -single_module $allow_undefined_flag -o $lib $libobjs $deplibs $compiler_flags -install_name $rpath/$soname $verstring~nmedit -s $output_objdir/${libname}-symbols.expsym ${lib}'
          else
            _LT_AC_TAGVAR(archive_expsym_cmds, $1)='sed -e "s,#.*,," -e "s,^[    ]*,," -e "s,^\(..*\),_&," < $export_symbols > $output_objdir/${libname}-symbols.expsym~$CC -r -keep_private_externs -nostdlib -o ${lib}-master.o $libobjs~$CC -dynamiclib $allow_undefined_flag -o $lib ${lib}-master.o $deplibs $compiler_flags -install_name $rpath/$soname $verstring~nmedit -s $output_objdir/${libname}-symbols.expsym ${lib}'
          fi
            _LT_AC_TAGVAR(module_expsym_cmds, $1)='sed -e "s,#.*,," -e "s,^[    ]*,," -e "s,^\(..*\),_&," < $export_symbols > $output_objdir/${libname}-symbols.expsym~$CC $allow_undefined_flag  -o $lib -bundle $libobjs $deplibs$compiler_flags~nmedit -s $output_objdir/${libname}-symbols.expsym ${lib}'
      else
      case $cc_basename in
        xlc*)
         output_verbose_link_cmd='echo'
          _LT_AC_TAGVAR(archive_cmds, $1)='$CC -qmkshrobj ${wl}-single_module $allow_undefined_flag -o $lib $libobjs $deplibs $compiler_flags ${wl}-install_name ${wl}`echo $rpath/$soname` $verstring'
          _LT_AC_TAGVAR(module_cmds, $1)='$CC $allow_undefined_flag -o $lib -bundle $libobjs $deplibs$compiler_flags'
          # Don't fix this by using the ld -exported_symbols_list flag, it doesn't exist in older darwin lds
          _LT_AC_TAGVAR(archive_expsym_cmds, $1)='sed -e "s,#.*,," -e "s,^[    ]*,," -e "s,^\(..*\),_&," < $export_symbols > $output_objdir/${libname}-symbols.expsym~$CC -qmkshrobj ${wl}-single_module $allow_undefined_flag -o $lib $libobjs $deplibs $compiler_flags ${wl}-install_name ${wl}$rpath/$soname $verstring~nmedit -s $output_objdir/${libname}-symbols.expsym ${lib}'
          _LT_AC_TAGVAR(module_expsym_cmds, $1)='sed -e "s,#.*,," -e "s,^[    ]*,," -e "s,^\(..*\),_&," < $export_symbols > $output_objdir/${libname}-symbols.expsym~$CC $allow_undefined_flag  -o $lib -bundle $libobjs $deplibs$compiler_flags~nmedit -s $output_objdir/${libname}-symbols.expsym ${lib}'
          ;;
       *)
         _LT_AC_TAGVAR(ld_shlibs, $1)=no
          ;;
      esac
      fi
        ;;

  dgux*)
    case $cc_basename in
      ec++*)
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
      ghcx*)
	# Green Hills C++ Compiler
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
      *)
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
    esac
    ;;
  freebsd[[12]]*)
    # C++ shared libraries reported to be fairly broken before switch to ELF
    _LT_AC_TAGVAR(ld_shlibs, $1)=no
    ;;
  freebsd-elf*)
    _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
    ;;
  freebsd* | kfreebsd*-gnu | dragonfly*)
    # FreeBSD 3 and later use GNU C++ and GNU ld with standard ELF
    # conventions
    _LT_AC_TAGVAR(ld_shlibs, $1)=yes
    ;;
  gnu*)
    ;;
  hpux9*)
    _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}+b ${wl}$libdir'
    _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:
    _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-E'
    _LT_AC_TAGVAR(hardcode_direct, $1)=yes
    _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes # Not in the search PATH,
				# but as the default
				# location of the library.

    case $cc_basename in
    CC*)
      # FIXME: insert proper C++ library support
      _LT_AC_TAGVAR(ld_shlibs, $1)=no
      ;;
    aCC*)
      _LT_AC_TAGVAR(archive_cmds, $1)='$rm $output_objdir/$soname~$CC -b ${wl}+b ${wl}$install_libdir -o $output_objdir/$soname $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags~test $output_objdir/$soname = $lib || mv $output_objdir/$soname $lib'
      # Commands to make compiler produce verbose output that lists
      # what "hidden" libraries, object files and flags are used when
      # linking a shared library.
      #
      # There doesn't appear to be a way to prevent this compiler from
      # explicitly linking system object files so we need to strip them
      # from the output so that they don't get included in the library
      # dependencies.
      output_verbose_link_cmd='templist=`($CC -b $CFLAGS -v conftest.$objext 2>&1) | grep "[[-]]L"`; list=""; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; echo $list'
      ;;
    *)
      if test "$GXX" = yes; then
        _LT_AC_TAGVAR(archive_cmds, $1)='$rm $output_objdir/$soname~$CC -shared -nostdlib -fPIC ${wl}+b ${wl}$install_libdir -o $output_objdir/$soname $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags~test $output_objdir/$soname = $lib || mv $output_objdir/$soname $lib'
      else
        # FIXME: insert proper C++ library support
        _LT_AC_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;
    esac
    ;;
  hpux10*|hpux11*)
    if test $with_gnu_ld = no; then
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}+b ${wl}$libdir'
      _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:

      case $host_cpu in
      hppa*64*|ia64*)
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec_ld, $1)='+b $libdir'
        ;;
      *)
	_LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-E'
        ;;
      esac
    fi
    case $host_cpu in
    hppa*64*|ia64*)
      _LT_AC_TAGVAR(hardcode_direct, $1)=no
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;
    *)
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes # Not in the search PATH,
					      # but as the default
					      # location of the library.
      ;;
    esac

    case $cc_basename in
      CC*)
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
      aCC*)
	case $host_cpu in
	hppa*64*)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -b ${wl}+h ${wl}$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	  ;;
	ia64*)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -b ${wl}+h ${wl}$soname ${wl}+nodefaultrpath -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	  ;;
	*)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -b ${wl}+h ${wl}$soname ${wl}+b ${wl}$install_libdir -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	  ;;
	esac
	# Commands to make compiler produce verbose output that lists
	# what "hidden" libraries, object files and flags are used when
	# linking a shared library.
	#
	# There doesn't appear to be a way to prevent this compiler from
	# explicitly linking system object files so we need to strip them
	# from the output so that they don't get included in the library
	# dependencies.
	output_verbose_link_cmd='templist=`($CC -b $CFLAGS -v conftest.$objext 2>&1) | grep "\-L"`; list=""; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; echo $list'
	;;
      *)
	if test "$GXX" = yes; then
	  if test $with_gnu_ld = no; then
	    case $host_cpu in
	    hppa*64*)
	      _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib -fPIC ${wl}+h ${wl}$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	      ;;
	    ia64*)
	      _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib -fPIC ${wl}+h ${wl}$soname ${wl}+nodefaultrpath -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	      ;;
	    *)
	      _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib -fPIC ${wl}+h ${wl}$soname ${wl}+b ${wl}$install_libdir -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	      ;;
	    esac
	  fi
	else
	  # FIXME: insert proper C++ library support
	  _LT_AC_TAGVAR(ld_shlibs, $1)=no
	fi
	;;
    esac
    ;;
  interix3*)
    _LT_AC_TAGVAR(hardcode_direct, $1)=no
    _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
    _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath,$libdir'
    _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-E'
    # Hack: On Interix 3.x, we cannot compile PIC because of a broken gcc.
    # Instead, shared libraries are loaded at an image base (0x10000000 by
    # default) and relocated if they conflict, which is a slow very memory
    # consuming and fragmenting process.  To avoid this, we pick a random,
    # 256 KiB-aligned image base between 0x50000000 and 0x6FFC0000 at link
    # time.  Moving up from 0x10000000 also allows more sbrk(2) space.
    _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags ${wl}-h,$soname ${wl}--image-base,`expr ${RANDOM-$$} % 4096 / 2 \* 262144 + 1342177280` -o $lib'
    _LT_AC_TAGVAR(archive_expsym_cmds, $1)='sed "s,^,_," $export_symbols >$output_objdir/$soname.expsym~$CC -shared $pic_flag $libobjs $deplibs $compiler_flags ${wl}-h,$soname ${wl}--retain-symbols-file,$output_objdir/$soname.expsym ${wl}--image-base,`expr ${RANDOM-$$} % 4096 / 2 \* 262144 + 1342177280` -o $lib'
    ;;
  irix5* | irix6*)
    case $cc_basename in
      CC*)
	# SGI C++
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -all -multigot $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -soname $soname `test -n "$verstring" && echo -set_version $verstring` -update_registry ${output_objdir}/so_locations -o $lib'

	# Archives containing C++ object files must be created using
	# "CC -ar", where "CC" is the IRIX C++ compiler.  This is
	# necessary to make sure instantiated templates are included
	# in the archive.
	_LT_AC_TAGVAR(old_archive_cmds, $1)='$CC -ar -WR,-u -o $oldlib $oldobjs'
	;;
      *)
	if test "$GXX" = yes; then
	  if test "$with_gnu_ld" = no; then
	    _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname ${wl}$soname `test -n "$verstring" && echo ${wl}-set_version ${wl}$verstring` ${wl}-update_registry ${wl}${output_objdir}/so_locations -o $lib'
	  else
	    _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname ${wl}$soname `test -n "$verstring" && echo ${wl}-set_version ${wl}$verstring` -o $lib'
	  fi
	fi
	_LT_AC_TAGVAR(link_all_deplibs, $1)=yes
	;;
    esac
    _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath ${wl}$libdir'
    _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:
    ;;
  linux*)
    case $cc_basename in
      KCC*)
	# Kuck and Associates, Inc. (KAI) C++ Compiler

	# KCC will only create a shared library if the output file
	# ends with ".so" (or ".sl" for HP-UX), so rename the library
	# to its proper name (with version) after linking.
	_LT_AC_TAGVAR(archive_cmds, $1)='tempext=`echo $shared_ext | $SED -e '\''s/\([[^()0-9A-Za-z{}]]\)/\\\\\1/g'\''`; templib=`echo $lib | $SED -e "s/\${tempext}\..*/.so/"`; $CC $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags --soname $soname -o \$templib; mv \$templib $lib'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='tempext=`echo $shared_ext | $SED -e '\''s/\([[^()0-9A-Za-z{}]]\)/\\\\\1/g'\''`; templib=`echo $lib | $SED -e "s/\${tempext}\..*/.so/"`; $CC $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags --soname $soname -o \$templib ${wl}-retain-symbols-file,$export_symbols; mv \$templib $lib'
	# Commands to make compiler produce verbose output that lists
	# what "hidden" libraries, object files and flags are used when
	# linking a shared library.
	#
	# There doesn't appear to be a way to prevent this compiler from
	# explicitly linking system object files so we need to strip them
	# from the output so that they don't get included in the library
	# dependencies.
	output_verbose_link_cmd='templist=`$CC $CFLAGS -v conftest.$objext -o libconftest$shared_ext 2>&1 | grep "ld"`; rm -f libconftest$shared_ext; list=""; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; echo $list'

	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}--rpath,$libdir'
	_LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}--export-dynamic'

	# Archives containing C++ object files must be created using
	# "CC -Bstatic", where "CC" is the KAI C++ compiler.
	_LT_AC_TAGVAR(old_archive_cmds, $1)='$CC -Bstatic -o $oldlib $oldobjs'
	;;
      icpc*)
	# Intel C++
	with_gnu_ld=yes
	# version 8.0 and above of icpc choke on multiply defined symbols
	# if we add $predep_objects and $postdep_objects, however 7.1 and
	# earlier do not add the objects themselves.
	case `$CC -V 2>&1` in
	*"Version 7."*)
  	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname $wl$soname -o $lib'
  	  _LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname $wl$soname ${wl}-retain-symbols-file $wl$export_symbols -o $lib'
	  ;;
	*)  # Version 8.0 or newer
	  tmp_idyn=
	  case $host_cpu in
	    ia64*) tmp_idyn=' -i_dynamic';;
	  esac
  	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared'"$tmp_idyn"' $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'
	  _LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared'"$tmp_idyn"' $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname ${wl}-retain-symbols-file $wl$export_symbols -o $lib'
	  ;;
	esac
	_LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath,$libdir'
	_LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}--export-dynamic'
	_LT_AC_TAGVAR(whole_archive_flag_spec, $1)='${wl}--whole-archive$convenience ${wl}--no-whole-archive'
	;;
      pgCC*)
        # Portland Group C++ compiler
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname ${wl}$soname -o $lib'
  	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $pic_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname ${wl}$soname ${wl}-retain-symbols-file ${wl}$export_symbols -o $lib'

	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}--rpath ${wl}$libdir'
	_LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}--export-dynamic'
	_LT_AC_TAGVAR(whole_archive_flag_spec, $1)='${wl}--whole-archive`for conv in $convenience\"\"; do test  -n \"$conv\" && new_convenience=\"$new_convenience,$conv\"; done; $echo \"$new_convenience\"` ${wl}--no-whole-archive'
        ;;
      cxx*)
	# Compaq C++
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname $wl$soname -o $lib'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname $wl$soname  -o $lib ${wl}-retain-symbols-file $wl$export_symbols'

	runpath_var=LD_RUN_PATH
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-rpath $libdir'
	_LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:

	# Commands to make compiler produce verbose output that lists
	# what "hidden" libraries, object files and flags are used when
	# linking a shared library.
	#
	# There doesn't appear to be a way to prevent this compiler from
	# explicitly linking system object files so we need to strip them
	# from the output so that they don't get included in the library
	# dependencies.
	output_verbose_link_cmd='templist=`$CC -shared $CFLAGS -v conftest.$objext 2>&1 | grep "ld"`; templist=`echo $templist | $SED "s/\(^.*ld.*\)\( .*ld .*$\)/\1/"`; list=""; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; echo $list'
	;;
    esac
    ;;
  lynxos*)
    # FIXME: insert proper C++ library support
    _LT_AC_TAGVAR(ld_shlibs, $1)=no
    ;;
  m88k*)
    # FIXME: insert proper C++ library support
    _LT_AC_TAGVAR(ld_shlibs, $1)=no
    ;;
  mvs*)
    case $cc_basename in
      cxx*)
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
      *)
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
    esac
    ;;
  netbsd*)
    if echo __ELF__ | $CC -E - | grep __ELF__ >/dev/null; then
      _LT_AC_TAGVAR(archive_cmds, $1)='$LD -Bshareable  -o $lib $predep_objects $libobjs $deplibs $postdep_objects $linker_flags'
      wlarc=
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
    fi
    # Workaround some broken pre-1.5 toolchains
    output_verbose_link_cmd='$CC -shared $CFLAGS -v conftest.$objext 2>&1 | grep conftest.$objext | $SED -e "s:-lgcc -lc -lgcc::"'
    ;;
  openbsd2*)
    # C++ shared libraries are fairly broken
    _LT_AC_TAGVAR(ld_shlibs, $1)=no
    ;;
  openbsd*)
    _LT_AC_TAGVAR(hardcode_direct, $1)=yes
    _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
    _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -o $lib'
    _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath,$libdir'
    if test -z "`echo __ELF__ | $CC -E - | grep __ELF__`" || test "$host_os-$host_cpu" = "openbsd2.8-powerpc"; then
      _LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $pic_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-retain-symbols-file,$export_symbols -o $lib'
      _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-E'
      _LT_AC_TAGVAR(whole_archive_flag_spec, $1)="$wlarc"'--whole-archive$convenience '"$wlarc"'--no-whole-archive'
    fi
    output_verbose_link_cmd='echo'
    ;;
  osf3*)
    case $cc_basename in
      KCC*)
	# Kuck and Associates, Inc. (KAI) C++ Compiler

	# KCC will only create a shared library if the output file
	# ends with ".so" (or ".sl" for HP-UX), so rename the library
	# to its proper name (with version) after linking.
	_LT_AC_TAGVAR(archive_cmds, $1)='tempext=`echo $shared_ext | $SED -e '\''s/\([[^()0-9A-Za-z{}]]\)/\\\\\1/g'\''`; templib=`echo $lib | $SED -e "s/\${tempext}\..*/.so/"`; $CC $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags --soname $soname -o \$templib; mv \$templib $lib'

	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath,$libdir'
	_LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:

	# Archives containing C++ object files must be created using
	# "CC -Bstatic", where "CC" is the KAI C++ compiler.
	_LT_AC_TAGVAR(old_archive_cmds, $1)='$CC -Bstatic -o $oldlib $oldobjs'

	;;
      RCC*)
	# Rational C++ 2.4.1
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
      cxx*)
	_LT_AC_TAGVAR(allow_undefined_flag, $1)=' ${wl}-expect_unresolved ${wl}\*'
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared${allow_undefined_flag} $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname $soname `test -n "$verstring" && echo ${wl}-set_version $verstring` -update_registry ${output_objdir}/so_locations -o $lib'

	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath ${wl}$libdir'
	_LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:

	# Commands to make compiler produce verbose output that lists
	# what "hidden" libraries, object files and flags are used when
	# linking a shared library.
	#
	# There doesn't appear to be a way to prevent this compiler from
	# explicitly linking system object files so we need to strip them
	# from the output so that they don't get included in the library
	# dependencies.
	output_verbose_link_cmd='templist=`$CC -shared $CFLAGS -v conftest.$objext 2>&1 | grep "ld" | grep -v "ld:"`; templist=`echo $templist | $SED "s/\(^.*ld.*\)\( .*ld.*$\)/\1/"`; list=""; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; echo $list'
	;;
      *)
	if test "$GXX" = yes && test "$with_gnu_ld" = no; then
	  _LT_AC_TAGVAR(allow_undefined_flag, $1)=' ${wl}-expect_unresolved ${wl}\*'
	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib ${allow_undefined_flag} $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-soname ${wl}$soname `test -n "$verstring" && echo ${wl}-set_version ${wl}$verstring` ${wl}-update_registry ${wl}${output_objdir}/so_locations -o $lib'

	  _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath ${wl}$libdir'
	  _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:

	  # Commands to make compiler produce verbose output that lists
	  # what "hidden" libraries, object files and flags are used when
	  # linking a shared library.
	  output_verbose_link_cmd='$CC -shared $CFLAGS -v conftest.$objext 2>&1 | grep "\-L"'

	else
	  # FIXME: insert proper C++ library support
	  _LT_AC_TAGVAR(ld_shlibs, $1)=no
	fi
	;;
    esac
    ;;
  osf4* | osf5*)
    case $cc_basename in
      KCC*)
	# Kuck and Associates, Inc. (KAI) C++ Compiler

	# KCC will only create a shared library if the output file
	# ends with ".so" (or ".sl" for HP-UX), so rename the library
	# to its proper name (with version) after linking.
	_LT_AC_TAGVAR(archive_cmds, $1)='tempext=`echo $shared_ext | $SED -e '\''s/\([[^()0-9A-Za-z{}]]\)/\\\\\1/g'\''`; templib=`echo $lib | $SED -e "s/\${tempext}\..*/.so/"`; $CC $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags --soname $soname -o \$templib; mv \$templib $lib'

	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath,$libdir'
	_LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:

	# Archives containing C++ object files must be created using
	# the KAI C++ compiler.
	_LT_AC_TAGVAR(old_archive_cmds, $1)='$CC -o $oldlib $oldobjs'
	;;
      RCC*)
	# Rational C++ 2.4.1
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
      cxx*)
	_LT_AC_TAGVAR(allow_undefined_flag, $1)=' -expect_unresolved \*'
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared${allow_undefined_flag} $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -msym -soname $soname `test -n "$verstring" && echo -set_version $verstring` -update_registry ${output_objdir}/so_locations -o $lib'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='for i in `cat $export_symbols`; do printf "%s %s\\n" -exported_symbol "\$i" >> $lib.exp; done~
	  echo "-hidden">> $lib.exp~
	  $CC -shared$allow_undefined_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -msym -soname $soname -Wl,-input -Wl,$lib.exp  `test -n "$verstring" && echo -set_version	$verstring` -update_registry ${output_objdir}/so_locations -o $lib~
	  $rm $lib.exp'

	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-rpath $libdir'
	_LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:

	# Commands to make compiler produce verbose output that lists
	# what "hidden" libraries, object files and flags are used when
	# linking a shared library.
	#
	# There doesn't appear to be a way to prevent this compiler from
	# explicitly linking system object files so we need to strip them
	# from the output so that they don't get included in the library
	# dependencies.
	output_verbose_link_cmd='templist=`$CC -shared $CFLAGS -v conftest.$objext 2>&1 | grep "ld" | grep -v "ld:"`; templist=`echo $templist | $SED "s/\(^.*ld.*\)\( .*ld.*$\)/\1/"`; list=""; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; echo $list'
	;;
      *)
	if test "$GXX" = yes && test "$with_gnu_ld" = no; then
	  _LT_AC_TAGVAR(allow_undefined_flag, $1)=' ${wl}-expect_unresolved ${wl}\*'
	 _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib ${allow_undefined_flag} $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-msym ${wl}-soname ${wl}$soname `test -n "$verstring" && echo ${wl}-set_version ${wl}$verstring` ${wl}-update_registry ${wl}${output_objdir}/so_locations -o $lib'

	  _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath ${wl}$libdir'
	  _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:

	  # Commands to make compiler produce verbose output that lists
	  # what "hidden" libraries, object files and flags are used when
	  # linking a shared library.
	  output_verbose_link_cmd='$CC -shared $CFLAGS -v conftest.$objext 2>&1 | grep "\-L"'

	else
	  # FIXME: insert proper C++ library support
	  _LT_AC_TAGVAR(ld_shlibs, $1)=no
	fi
	;;
    esac
    ;;
  psos*)
    # FIXME: insert proper C++ library support
    _LT_AC_TAGVAR(ld_shlibs, $1)=no
    ;;
  sunos4*)
    case $cc_basename in
      CC*)
	# Sun C++ 4.x
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
      lcc*)
	# Lucid
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
      *)
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
    esac
    ;;
  solaris*)
    case $cc_basename in
      CC*)
	# Sun C++ 4.2, 5.x and Centerline C++
        _LT_AC_TAGVAR(archive_cmds_need_lc,$1)=yes
	_LT_AC_TAGVAR(no_undefined_flag, $1)=' -zdefs'
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -G${allow_undefined_flag}  -h$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~$echo "local: *; };" >> $lib.exp~
	$CC -G${allow_undefined_flag}  ${wl}-M ${wl}$lib.exp -h$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags~$rm $lib.exp'

	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
	_LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
	case $host_os in
	  solaris2.[[0-5]] | solaris2.[[0-5]].*) ;;
	  *)
	    # The C++ compiler is used as linker so we must use $wl
	    # flag to pass the commands to the underlying system
	    # linker. We must also pass each convience library through
	    # to the system linker between allextract/defaultextract.
	    # The C++ compiler will combine linker options so we
	    # cannot just pass the convience library names through
	    # without $wl.
	    # Supported since Solaris 2.6 (maybe 2.5.1?)
	    _LT_AC_TAGVAR(whole_archive_flag_spec, $1)='${wl}-z ${wl}allextract`for conv in $convenience\"\"; do test -n \"$conv\" && new_convenience=\"$new_convenience,$conv\"; done; $echo \"$new_convenience\"` ${wl}-z ${wl}defaultextract'
	    ;;
	esac
	_LT_AC_TAGVAR(link_all_deplibs, $1)=yes

	output_verbose_link_cmd='echo'

	# Archives containing C++ object files must be created using
	# "CC -xar", where "CC" is the Sun C++ compiler.  This is
	# necessary to make sure instantiated templates are included
	# in the archive.
	_LT_AC_TAGVAR(old_archive_cmds, $1)='$CC -xar -o $oldlib $oldobjs'
	;;
      gcx*)
	# Green Hills C++ Compiler
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-h $wl$soname -o $lib'

	# The C++ compiler must be used to create the archive.
	_LT_AC_TAGVAR(old_archive_cmds, $1)='$CC $LDFLAGS -archive -o $oldlib $oldobjs'
	;;
      *)
	# GNU C++ compiler with Solaris linker
	if test "$GXX" = yes && test "$with_gnu_ld" = no; then
	  _LT_AC_TAGVAR(no_undefined_flag, $1)=' ${wl}-z ${wl}defs'
	  if $CC --version | grep -v '^2\.7' > /dev/null; then
	    _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $LDFLAGS $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-h $wl$soname -o $lib'
	    _LT_AC_TAGVAR(archive_expsym_cmds, $1)='$echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~$echo "local: *; };" >> $lib.exp~
		$CC -shared -nostdlib ${wl}-M $wl$lib.exp -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags~$rm $lib.exp'

	    # Commands to make compiler produce verbose output that lists
	    # what "hidden" libraries, object files and flags are used when
	    # linking a shared library.
	    output_verbose_link_cmd="$CC -shared $CFLAGS -v conftest.$objext 2>&1 | grep \"\-L\""
	  else
	    # g++ 2.7 appears to require `-G' NOT `-shared' on this
	    # platform.
	    _LT_AC_TAGVAR(archive_cmds, $1)='$CC -G -nostdlib $LDFLAGS $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags ${wl}-h $wl$soname -o $lib'
	    _LT_AC_TAGVAR(archive_expsym_cmds, $1)='$echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~$echo "local: *; };" >> $lib.exp~
		$CC -G -nostdlib ${wl}-M $wl$lib.exp -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags~$rm $lib.exp'

	    # Commands to make compiler produce verbose output that lists
	    # what "hidden" libraries, object files and flags are used when
	    # linking a shared library.
	    output_verbose_link_cmd="$CC -G $CFLAGS -v conftest.$objext 2>&1 | grep \"\-L\""
	  fi

	  _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-R $wl$libdir'
	fi
	;;
    esac
    ;;
  sysv4*uw2* | sysv5OpenUNIX* | sysv5UnixWare7.[[01]].[[10]]* | unixware7* | sco3.2v5.0.[[024]]*)
    _LT_AC_TAGVAR(no_undefined_flag, $1)='${wl}-z,text'
    _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
    _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
    runpath_var='LD_RUN_PATH'

    case $cc_basename in
      CC*)
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -G ${wl}-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -G ${wl}-Bexport:$export_symbols ${wl}-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	;;
      *)
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared ${wl}-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared ${wl}-Bexport:$export_symbols ${wl}-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	;;
    esac
    ;;
  sysv5* | sco3.2v5* | sco5v6*)
    # Note: We can NOT use -z defs as we might desire, because we do not
    # link with -lc, and that would cause any symbols used from libc to
    # always be unresolved, which means just about no library would
    # ever link correctly.  If we're not using GNU ld we use -z text
    # though, which does catch some bad symbols but isn't as heavy-handed
    # as -z defs.
    # For security reasons, it is highly recommended that you always
    # use absolute paths for naming shared libraries, and exclude the
    # DT_RUNPATH tag from executables and libraries.  But doing so
    # requires that you compile everything twice, which is a pain.
    # So that behaviour is only enabled if SCOABSPATH is set to a
    # non-empty value in the environment.  Most likely only useful for
    # creating official distributions of packages.
    # This is a hack until libtool officially supports absolute path
    # names for shared libraries.
    _LT_AC_TAGVAR(no_undefined_flag, $1)='${wl}-z,text'
    _LT_AC_TAGVAR(allow_undefined_flag, $1)='${wl}-z,nodefs'
    _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
    _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
    _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='`test -z "$SCOABSPATH" && echo ${wl}-R,$libdir`'
    _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=':'
    _LT_AC_TAGVAR(link_all_deplibs, $1)=yes
    _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-Bexport'
    runpath_var='LD_RUN_PATH'

    case $cc_basename in
      CC*)
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -G ${wl}-h,\${SCOABSPATH:+${install_libdir}/}$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -G ${wl}-Bexport:$export_symbols ${wl}-h,\${SCOABSPATH:+${install_libdir}/}$soname -o $lib $libobjs $deplibs $compiler_flags'
	;;
      *)
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared ${wl}-h,\${SCOABSPATH:+${install_libdir}/}$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared ${wl}-Bexport:$export_symbols ${wl}-h,\${SCOABSPATH:+${install_libdir}/}$soname -o $lib $libobjs $deplibs $compiler_flags'
	;;
    esac
    ;;
  tandem*)
    case $cc_basename in
      NCC*)
	# NonStop-UX NCC 3.20
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
      *)
	# FIXME: insert proper C++ library support
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	;;
    esac
    ;;
  vxworks*)
    # FIXME: insert proper C++ library support
    _LT_AC_TAGVAR(ld_shlibs, $1)=no
    ;;
  *)
    # FIXME: insert proper C++ library support
    _LT_AC_TAGVAR(ld_shlibs, $1)=no
    ;;
esac
AC_MSG_RESULT([$_LT_AC_TAGVAR(ld_shlibs, $1)])
test "$_LT_AC_TAGVAR(ld_shlibs, $1)" = no && can_build_shared=no

_LT_AC_TAGVAR(GCC, $1)="$GXX"
_LT_AC_TAGVAR(LD, $1)="$LD"

## CAVEAT EMPTOR:
## There is no encapsulation within the following macros, do not change
## the running order or otherwise move them around unless you know exactly
## what you are doing...
AC_LIBTOOL_POSTDEP_PREDEP($1)
AC_LIBTOOL_PROG_COMPILER_PIC($1)
AC_LIBTOOL_PROG_CC_C_O($1)
AC_LIBTOOL_SYS_HARD_LINK_LOCKS($1)
AC_LIBTOOL_PROG_LD_SHLIBS($1)
AC_LIBTOOL_SYS_DYNAMIC_LINKER($1)
AC_LIBTOOL_PROG_LD_HARDCODE_LIBPATH($1)

AC_LIBTOOL_CONFIG($1)

AC_LANG_POP
CC=$lt_save_CC
LDCXX=$LD
LD=$lt_save_LD
GCC=$lt_save_GCC
with_gnu_ldcxx=$with_gnu_ld
with_gnu_ld=$lt_save_with_gnu_ld
lt_cv_path_LDCXX=$lt_cv_path_LD
lt_cv_path_LD=$lt_save_path_LD
lt_cv_prog_gnu_ldcxx=$lt_cv_prog_gnu_ld
lt_cv_prog_gnu_ld=$lt_save_with_gnu_ld
])# AC_LIBTOOL_LANG_CXX_CONFIG

# AC_LIBTOOL_POSTDEP_PREDEP([TAGNAME])
# ------------------------------------
# Figure out "hidden" library dependencies from verbose
# compiler output when linking a shared library.
# Parse the compiler output and extract the necessary
# objects, libraries and library flags.
AC_DEFUN([AC_LIBTOOL_POSTDEP_PREDEP],[
dnl we can't use the lt_simple_compile_test_code here,
dnl because it contains code intended for an executable,
dnl not a library.  It's possible we should let each
dnl tag define a new lt_????_link_test_code variable,
dnl but it's only used here...
ifelse([$1],[],[cat > conftest.$ac_ext <<EOF
int a;
void foo (void) { a = 0; }
EOF
],[$1],[CXX],[cat > conftest.$ac_ext <<EOF
class Foo
{
public:
  Foo (void) { a = 0; }
private:
  int a;
};
EOF
],[$1],[F77],[cat > conftest.$ac_ext <<EOF
      subroutine foo
      implicit none
      integer*4 a
      a=0
      return
      end
EOF
],[$1],[GCJ],[cat > conftest.$ac_ext <<EOF
public class foo {
  private int a;
  public void bar (void) {
    a = 0;
  }
};
EOF
])
dnl Parse the compiler output and extract the necessary
dnl objects, libraries and library flags.
if AC_TRY_EVAL(ac_compile); then
  # Parse the compiler output and extract the necessary
  # objects, libraries and library flags.

  # Sentinel used to keep track of whether or not we are before
  # the conftest object file.
  pre_test_object_deps_done=no

  # The `*' in the case matches for architectures that use `case' in
  # $output_verbose_cmd can trigger glob expansion during the loop
  # eval without this substitution.
  output_verbose_link_cmd=`$echo "X$output_verbose_link_cmd" | $Xsed -e "$no_glob_subst"`

  for p in `eval $output_verbose_link_cmd`; do
    case $p in

    -L* | -R* | -l*)
       # Some compilers place space between "-{L,R}" and the path.
       # Remove the space.
       if test $p = "-L" \
	  || test $p = "-R"; then
	 prev=$p
	 continue
       else
	 prev=
       fi

       if test "$pre_test_object_deps_done" = no; then
	 case $p in
	 -L* | -R*)
	   # Internal compiler library paths should come after those
	   # provided the user.  The postdeps already come after the
	   # user supplied libs so there is no need to process them.
	   if test -z "$_LT_AC_TAGVAR(compiler_lib_search_path, $1)"; then
	     _LT_AC_TAGVAR(compiler_lib_search_path, $1)="${prev}${p}"
	   else
	     _LT_AC_TAGVAR(compiler_lib_search_path, $1)="${_LT_AC_TAGVAR(compiler_lib_search_path, $1)} ${prev}${p}"
	   fi
	   ;;
	 # The "-l" case would never come before the object being
	 # linked, so don't bother handling this case.
	 esac
       else
	 if test -z "$_LT_AC_TAGVAR(postdeps, $1)"; then
	   _LT_AC_TAGVAR(postdeps, $1)="${prev}${p}"
	 else
	   _LT_AC_TAGVAR(postdeps, $1)="${_LT_AC_TAGVAR(postdeps, $1)} ${prev}${p}"
	 fi
       fi
       ;;

    *.$objext)
       # This assumes that the test object file only shows up
       # once in the compiler output.
       if test "$p" = "conftest.$objext"; then
	 pre_test_object_deps_done=yes
	 continue
       fi

       if test "$pre_test_object_deps_done" = no; then
	 if test -z "$_LT_AC_TAGVAR(predep_objects, $1)"; then
	   _LT_AC_TAGVAR(predep_objects, $1)="$p"
	 else
	   _LT_AC_TAGVAR(predep_objects, $1)="$_LT_AC_TAGVAR(predep_objects, $1) $p"
	 fi
       else
	 if test -z "$_LT_AC_TAGVAR(postdep_objects, $1)"; then
	   _LT_AC_TAGVAR(postdep_objects, $1)="$p"
	 else
	   _LT_AC_TAGVAR(postdep_objects, $1)="$_LT_AC_TAGVAR(postdep_objects, $1) $p"
	 fi
       fi
       ;;

    *) ;; # Ignore the rest.

    esac
  done

  # Clean up.
  rm -f a.out a.exe
else
  echo "libtool.m4: error: problem compiling $1 test program"
fi

$rm -f confest.$objext

# PORTME: override above test on systems where it is broken
ifelse([$1],[CXX],
[case $host_os in
interix3*)
  # Interix 3.5 installs completely hosed .la files for C++, so rather than
  # hack all around it, let's just trust "g++" to DTRT.
  _LT_AC_TAGVAR(predep_objects,$1)=
  _LT_AC_TAGVAR(postdep_objects,$1)=
  _LT_AC_TAGVAR(postdeps,$1)=
  ;;

solaris*)
  case $cc_basename in
  CC*)
    # Adding this requires a known-good setup of shared libraries for
    # Sun compiler versions before 5.6, else PIC objects from an old
    # archive will be linked into the output, leading to subtle bugs.
    _LT_AC_TAGVAR(postdeps,$1)='-lCstd -lCrun'
    ;;
  esac
  ;;
esac
])

case " $_LT_AC_TAGVAR(postdeps, $1) " in
*" -lc "*) _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no ;;
esac
])# AC_LIBTOOL_POSTDEP_PREDEP

# AC_LIBTOOL_LANG_F77_CONFIG
# --------------------------
# Ensure that the configuration vars for the C compiler are
# suitably defined.  Those variables are subsequently used by
# AC_LIBTOOL_CONFIG to write the compiler configuration to `libtool'.
AC_DEFUN([AC_LIBTOOL_LANG_F77_CONFIG], [_LT_AC_LANG_F77_CONFIG(F77)])
AC_DEFUN([_LT_AC_LANG_F77_CONFIG],
[AC_REQUIRE([AC_PROG_F77])
AC_LANG_PUSH(Fortran 77)

_LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
_LT_AC_TAGVAR(allow_undefined_flag, $1)=
_LT_AC_TAGVAR(always_export_symbols, $1)=no
_LT_AC_TAGVAR(archive_expsym_cmds, $1)=
_LT_AC_TAGVAR(export_dynamic_flag_spec, $1)=
_LT_AC_TAGVAR(hardcode_direct, $1)=no
_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)=
_LT_AC_TAGVAR(hardcode_libdir_flag_spec_ld, $1)=
_LT_AC_TAGVAR(hardcode_libdir_separator, $1)=
_LT_AC_TAGVAR(hardcode_minus_L, $1)=no
_LT_AC_TAGVAR(hardcode_automatic, $1)=no
_LT_AC_TAGVAR(module_cmds, $1)=
_LT_AC_TAGVAR(module_expsym_cmds, $1)=
_LT_AC_TAGVAR(link_all_deplibs, $1)=unknown
_LT_AC_TAGVAR(old_archive_cmds, $1)=$old_archive_cmds
_LT_AC_TAGVAR(no_undefined_flag, $1)=
_LT_AC_TAGVAR(whole_archive_flag_spec, $1)=
_LT_AC_TAGVAR(enable_shared_with_static_runtimes, $1)=no

# Source file extension for f77 test sources.
ac_ext=f

# Object file extension for compiled f77 test sources.
objext=o
_LT_AC_TAGVAR(objext, $1)=$objext

# Code to be used in simple compile tests
lt_simple_compile_test_code="      subroutine t\n      return\n      end\n"

# Code to be used in simple link tests
lt_simple_link_test_code="      program t\n      end\n"

# ltmain only uses $CC for tagged configurations so make sure $CC is set.
_LT_AC_SYS_COMPILER

# save warnings/boilerplate of simple test code
_LT_COMPILER_BOILERPLATE
_LT_LINKER_BOILERPLATE

# Allow CC to be a program name with arguments.
lt_save_CC="$CC"
CC=${F77-"f77"}
compiler=$CC
_LT_AC_TAGVAR(compiler, $1)=$CC
_LT_CC_BASENAME([$compiler])

AC_MSG_CHECKING([if libtool supports shared libraries])
AC_MSG_RESULT([$can_build_shared])

AC_MSG_CHECKING([whether to build shared libraries])
test "$can_build_shared" = "no" && enable_shared=no

# On AIX, shared libraries and static libraries use the same namespace, and
# are all built from PIC.
case $host_os in
aix3*)
  test "$enable_shared" = yes && enable_static=no
  if test -n "$RANLIB"; then
    archive_cmds="$archive_cmds~\$RANLIB \$lib"
    postinstall_cmds='$RANLIB $lib'
  fi
  ;;
aix4* | aix5*)
  if test "$host_cpu" != ia64 && test "$aix_use_runtimelinking" = no ; then
    test "$enable_shared" = yes && enable_static=no
  fi
  ;;
esac
AC_MSG_RESULT([$enable_shared])

AC_MSG_CHECKING([whether to build static libraries])
# Make sure either enable_shared or enable_static is yes.
test "$enable_shared" = yes || enable_static=yes
AC_MSG_RESULT([$enable_static])

_LT_AC_TAGVAR(GCC, $1)="$G77"
_LT_AC_TAGVAR(LD, $1)="$LD"

AC_LIBTOOL_PROG_COMPILER_PIC($1)
AC_LIBTOOL_PROG_CC_C_O($1)
AC_LIBTOOL_SYS_HARD_LINK_LOCKS($1)
AC_LIBTOOL_PROG_LD_SHLIBS($1)
AC_LIBTOOL_SYS_DYNAMIC_LINKER($1)
AC_LIBTOOL_PROG_LD_HARDCODE_LIBPATH($1)

AC_LIBTOOL_CONFIG($1)

AC_LANG_POP
CC="$lt_save_CC"
])# AC_LIBTOOL_LANG_F77_CONFIG


# AC_LIBTOOL_LANG_GCJ_CONFIG
# --------------------------
# Ensure that the configuration vars for the C compiler are
# suitably defined.  Those variables are subsequently used by
# AC_LIBTOOL_CONFIG to write the compiler configuration to `libtool'.
AC_DEFUN([AC_LIBTOOL_LANG_GCJ_CONFIG], [_LT_AC_LANG_GCJ_CONFIG(GCJ)])
AC_DEFUN([_LT_AC_LANG_GCJ_CONFIG],
[AC_LANG_SAVE

# Source file extension for Java test sources.
ac_ext=java

# Object file extension for compiled Java test sources.
objext=o
_LT_AC_TAGVAR(objext, $1)=$objext

# Code to be used in simple compile tests
lt_simple_compile_test_code="class foo {}\n"

# Code to be used in simple link tests
lt_simple_link_test_code='public class conftest { public static void main(String[[]] argv) {}; }\n'

# ltmain only uses $CC for tagged configurations so make sure $CC is set.
_LT_AC_SYS_COMPILER

# save warnings/boilerplate of simple test code
_LT_COMPILER_BOILERPLATE
_LT_LINKER_BOILERPLATE

# Allow CC to be a program name with arguments.
lt_save_CC="$CC"
CC=${GCJ-"gcj"}
compiler=$CC
_LT_AC_TAGVAR(compiler, $1)=$CC
_LT_CC_BASENAME([$compiler])

# GCJ did not exist at the time GCC didn't implicitly link libc in.
_LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no

_LT_AC_TAGVAR(old_archive_cmds, $1)=$old_archive_cmds

## CAVEAT EMPTOR:
## There is no encapsulation within the following macros, do not change
## the running order or otherwise move them around unless you know exactly
## what you are doing...
AC_LIBTOOL_PROG_COMPILER_NO_RTTI($1)
AC_LIBTOOL_PROG_COMPILER_PIC($1)
AC_LIBTOOL_PROG_CC_C_O($1)
AC_LIBTOOL_SYS_HARD_LINK_LOCKS($1)
AC_LIBTOOL_PROG_LD_SHLIBS($1)
AC_LIBTOOL_SYS_DYNAMIC_LINKER($1)
AC_LIBTOOL_PROG_LD_HARDCODE_LIBPATH($1)

AC_LIBTOOL_CONFIG($1)

AC_LANG_RESTORE
CC="$lt_save_CC"
])# AC_LIBTOOL_LANG_GCJ_CONFIG


# AC_LIBTOOL_LANG_RC_CONFIG
# -------------------------
# Ensure that the configuration vars for the Windows resource compiler are
# suitably defined.  Those variables are subsequently used by
# AC_LIBTOOL_CONFIG to write the compiler configuration to `libtool'.
AC_DEFUN([AC_LIBTOOL_LANG_RC_CONFIG], [_LT_AC_LANG_RC_CONFIG(RC)])
AC_DEFUN([_LT_AC_LANG_RC_CONFIG],
[AC_LANG_SAVE

# Source file extension for RC test sources.
ac_ext=rc

# Object file extension for compiled RC test sources.
objext=o
_LT_AC_TAGVAR(objext, $1)=$objext

# Code to be used in simple compile tests
lt_simple_compile_test_code='sample MENU { MENUITEM "&Soup", 100, CHECKED }\n'

# Code to be used in simple link tests
lt_simple_link_test_code="$lt_simple_compile_test_code"

# ltmain only uses $CC for tagged configurations so make sure $CC is set.
_LT_AC_SYS_COMPILER

# save warnings/boilerplate of simple test code
_LT_COMPILER_BOILERPLATE
_LT_LINKER_BOILERPLATE

# Allow CC to be a program name with arguments.
lt_save_CC="$CC"
CC=${RC-"windres"}
compiler=$CC
_LT_AC_TAGVAR(compiler, $1)=$CC
_LT_CC_BASENAME([$compiler])
_LT_AC_TAGVAR(lt_cv_prog_compiler_c_o, $1)=yes

AC_LIBTOOL_CONFIG($1)

AC_LANG_RESTORE
CC="$lt_save_CC"
])# AC_LIBTOOL_LANG_RC_CONFIG


# AC_LIBTOOL_CONFIG([TAGNAME])
# ----------------------------
# If TAGNAME is not passed, then create an initial libtool script
# with a default configuration from the untagged config vars.  Otherwise
# add code to config.status for appending the configuration named by
# TAGNAME from the matching tagged config vars.
AC_DEFUN([AC_LIBTOOL_CONFIG],
[# The else clause should only fire when bootstrapping the
# libtool distribution, otherwise you forgot to ship ltmain.sh
# with your package, and you will get complaints that there are
# no rules to generate ltmain.sh.
if test -f "$ltmain"; then
  # See if we are running on zsh, and set the options which allow our commands through
  # without removal of \ escapes.
  if test -n "${ZSH_VERSION+set}" ; then
    setopt NO_GLOB_SUBST
  fi
  # Now quote all the things that may contain metacharacters while being
  # careful not to overquote the AC_SUBSTed values.  We take copies of the
  # variables and quote the copies for generation of the libtool script.
  for var in echo old_CC old_CFLAGS AR AR_FLAGS EGREP RANLIB LN_S LTCC LTCFLAGS NM \
    SED SHELL STRIP \
    libname_spec library_names_spec soname_spec extract_expsyms_cmds \
    old_striplib striplib file_magic_cmd finish_cmds finish_eval \
    deplibs_check_method reload_flag reload_cmds need_locks \
    lt_cv_sys_global_symbol_pipe lt_cv_sys_global_symbol_to_cdecl \
    lt_cv_sys_global_symbol_to_c_name_address \
    sys_lib_search_path_spec sys_lib_dlsearch_path_spec \
    old_postinstall_cmds old_postuninstall_cmds \
    _LT_AC_TAGVAR(compiler, $1) \
    _LT_AC_TAGVAR(CC, $1) \
    _LT_AC_TAGVAR(LD, $1) \
    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1) \
    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1) \
    _LT_AC_TAGVAR(lt_prog_compiler_static, $1) \
    _LT_AC_TAGVAR(lt_prog_compiler_no_builtin_flag, $1) \
    _LT_AC_TAGVAR(export_dynamic_flag_spec, $1) \
    _LT_AC_TAGVAR(thread_safe_flag_spec, $1) \
    _LT_AC_TAGVAR(whole_archive_flag_spec, $1) \
    _LT_AC_TAGVAR(enable_shared_with_static_runtimes, $1) \
    _LT_AC_TAGVAR(old_archive_cmds, $1) \
    _LT_AC_TAGVAR(old_archive_from_new_cmds, $1) \
    _LT_AC_TAGVAR(predep_objects, $1) \
    _LT_AC_TAGVAR(postdep_objects, $1) \
    _LT_AC_TAGVAR(predeps, $1) \
    _LT_AC_TAGVAR(postdeps, $1) \
    _LT_AC_TAGVAR(compiler_lib_search_path, $1) \
    _LT_AC_TAGVAR(archive_cmds, $1) \
    _LT_AC_TAGVAR(archive_expsym_cmds, $1) \
    _LT_AC_TAGVAR(postinstall_cmds, $1) \
    _LT_AC_TAGVAR(postuninstall_cmds, $1) \
    _LT_AC_TAGVAR(old_archive_from_expsyms_cmds, $1) \
    _LT_AC_TAGVAR(allow_undefined_flag, $1) \
    _LT_AC_TAGVAR(no_undefined_flag, $1) \
    _LT_AC_TAGVAR(export_symbols_cmds, $1) \
    _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1) \
    _LT_AC_TAGVAR(hardcode_libdir_flag_spec_ld, $1) \
    _LT_AC_TAGVAR(hardcode_libdir_separator, $1) \
    _LT_AC_TAGVAR(hardcode_automatic, $1) \
    _LT_AC_TAGVAR(module_cmds, $1) \
    _LT_AC_TAGVAR(module_expsym_cmds, $1) \
    _LT_AC_TAGVAR(lt_cv_prog_compiler_c_o, $1) \
    _LT_AC_TAGVAR(exclude_expsyms, $1) \
    _LT_AC_TAGVAR(include_expsyms, $1); do

    case $var in
    _LT_AC_TAGVAR(old_archive_cmds, $1) | \
    _LT_AC_TAGVAR(old_archive_from_new_cmds, $1) | \
    _LT_AC_TAGVAR(archive_cmds, $1) | \
    _LT_AC_TAGVAR(archive_expsym_cmds, $1) | \
    _LT_AC_TAGVAR(module_cmds, $1) | \
    _LT_AC_TAGVAR(module_expsym_cmds, $1) | \
    _LT_AC_TAGVAR(old_archive_from_expsyms_cmds, $1) | \
    _LT_AC_TAGVAR(export_symbols_cmds, $1) | \
    extract_expsyms_cmds | reload_cmds | finish_cmds | \
    postinstall_cmds | postuninstall_cmds | \
    old_postinstall_cmds | old_postuninstall_cmds | \
    sys_lib_search_path_spec | sys_lib_dlsearch_path_spec)
      # Double-quote double-evaled strings.
      eval "lt_$var=\\\"\`\$echo \"X\$$var\" | \$Xsed -e \"\$double_quote_subst\" -e \"\$sed_quote_subst\" -e \"\$delay_variable_subst\"\`\\\""
      ;;
    *)
      eval "lt_$var=\\\"\`\$echo \"X\$$var\" | \$Xsed -e \"\$sed_quote_subst\"\`\\\""
      ;;
    esac
  done

  case $lt_echo in
  *'\[$]0 --fallback-echo"')
    lt_echo=`$echo "X$lt_echo" | $Xsed -e 's/\\\\\\\[$]0 --fallback-echo"[$]/[$]0 --fallback-echo"/'`
    ;;
  esac

ifelse([$1], [],
  [cfgfile="${ofile}T"
  trap "$rm \"$cfgfile\"; exit 1" 1 2 15
  $rm -f "$cfgfile"
  AC_MSG_NOTICE([creating $ofile])],
  [cfgfile="$ofile"])

  cat <<__EOF__ >> "$cfgfile"
ifelse([$1], [],
[#! $SHELL

# `$echo "$cfgfile" | sed 's%^.*/%%'` - Provide generalized library-building support services.
# Generated automatically by $PROGRAM (GNU $PACKAGE $VERSION$TIMESTAMP)
# NOTE: Changes made to this file will be lost: look at ltmain.sh.
#
# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001
# Free Software Foundation, Inc.
#
# This file is part of GNU Libtool:
# Originally by Gordon Matzigkeit <gord@gnu.ai.mit.edu>, 1996
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# As a special exception to the GNU General Public License, if you
# distribute this file as part of a program that contains a
# configuration script generated by Autoconf, you may include it under
# the same distribution terms that you use for the rest of that program.

# A sed program that does not truncate output.
SED=$lt_SED

# Sed that helps us avoid accidentally triggering echo(1) options like -n.
Xsed="$SED -e 1s/^X//"

# The HP-UX ksh and POSIX shell print the target directory to stdout
# if CDPATH is set.
(unset CDPATH) >/dev/null 2>&1 && unset CDPATH

# The names of the tagged configurations supported by this script.
available_tags=

# ### BEGIN LIBTOOL CONFIG],
[# ### BEGIN LIBTOOL TAG CONFIG: $tagname])

# Libtool was configured on host `(hostname || uname -n) 2>/dev/null | sed 1q`:

# Shell to use when invoking shell scripts.
SHELL=$lt_SHELL

# Whether or not to build shared libraries.
build_libtool_libs=$enable_shared

# Whether or not to build static libraries.
build_old_libs=$enable_static

# Whether or not to add -lc for building shared libraries.
build_libtool_need_lc=$_LT_AC_TAGVAR(archive_cmds_need_lc, $1)

# Whether or not to disallow shared libs when runtime libs are static
allow_libtool_libs_with_static_runtimes=$_LT_AC_TAGVAR(enable_shared_with_static_runtimes, $1)

# Whether or not to optimize for fast installation.
fast_install=$enable_fast_install

# The host system.
host_alias=$host_alias
host=$host
host_os=$host_os

# The build system.
build_alias=$build_alias
build=$build
build_os=$build_os

# An echo program that does not interpret backslashes.
echo=$lt_echo

# The archiver.
AR=$lt_AR
AR_FLAGS=$lt_AR_FLAGS

# A C compiler.
LTCC=$lt_LTCC

# LTCC compiler flags.
LTCFLAGS=$lt_LTCFLAGS

# A language-specific compiler.
CC=$lt_[]_LT_AC_TAGVAR(compiler, $1)

# Is the compiler the GNU C compiler?
with_gcc=$_LT_AC_TAGVAR(GCC, $1)

# An ERE matcher.
EGREP=$lt_EGREP

# The linker used to build libraries.
LD=$lt_[]_LT_AC_TAGVAR(LD, $1)

# Whether we need hard or soft links.
LN_S=$lt_LN_S

# A BSD-compatible nm program.
NM=$lt_NM

# A symbol stripping program
STRIP=$lt_STRIP

# Used to examine libraries when file_magic_cmd begins "file"
MAGIC_CMD=$MAGIC_CMD

# Used on cygwin: DLL creation program.
DLLTOOL="$DLLTOOL"

# Used on cygwin: object dumper.
OBJDUMP="$OBJDUMP"

# Used on cygwin: assembler.
AS="$AS"

# The name of the directory that contains temporary libtool files.
objdir=$objdir

# How to create reloadable object files.
reload_flag=$lt_reload_flag
reload_cmds=$lt_reload_cmds

# How to pass a linker flag through the compiler.
wl=$lt_[]_LT_AC_TAGVAR(lt_prog_compiler_wl, $1)

# Object file suffix (normally "o").
objext="$ac_objext"

# Old archive suffix (normally "a").
libext="$libext"

# Shared library suffix (normally ".so").
shrext_cmds='$shrext_cmds'

# Executable file suffix (normally "").
exeext="$exeext"

# Additional compiler flags for building library objects.
pic_flag=$lt_[]_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)
pic_mode=$pic_mode

# What is the maximum length of a command?
max_cmd_len=$lt_cv_sys_max_cmd_len

# Does compiler simultaneously support -c and -o options?
compiler_c_o=$lt_[]_LT_AC_TAGVAR(lt_cv_prog_compiler_c_o, $1)

# Must we lock files when doing compilation?
need_locks=$lt_need_locks

# Do we need the lib prefix for modules?
need_lib_prefix=$need_lib_prefix

# Do we need a version for libraries?
need_version=$need_version

# Whether dlopen is supported.
dlopen_support=$enable_dlopen

# Whether dlopen of programs is supported.
dlopen_self=$enable_dlopen_self

# Whether dlopen of statically linked programs is supported.
dlopen_self_static=$enable_dlopen_self_static

# Compiler flag to prevent dynamic linking.
link_static_flag=$lt_[]_LT_AC_TAGVAR(lt_prog_compiler_static, $1)

# Compiler flag to turn off builtin functions.
no_builtin_flag=$lt_[]_LT_AC_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)

# Compiler flag to allow reflexive dlopens.
export_dynamic_flag_spec=$lt_[]_LT_AC_TAGVAR(export_dynamic_flag_spec, $1)

# Compiler flag to generate shared objects directly from archives.
whole_archive_flag_spec=$lt_[]_LT_AC_TAGVAR(whole_archive_flag_spec, $1)

# Compiler flag to generate thread-safe objects.
thread_safe_flag_spec=$lt_[]_LT_AC_TAGVAR(thread_safe_flag_spec, $1)

# Library versioning type.
version_type=$version_type

# Format of library name prefix.
libname_spec=$lt_libname_spec

# List of archive names.  First name is the real one, the rest are links.
# The last name is the one that the linker finds with -lNAME.
library_names_spec=$lt_library_names_spec

# The coded name of the library, if different from the real name.
soname_spec=$lt_soname_spec

# Commands used to build and install an old-style archive.
RANLIB=$lt_RANLIB
old_archive_cmds=$lt_[]_LT_AC_TAGVAR(old_archive_cmds, $1)
old_postinstall_cmds=$lt_old_postinstall_cmds
old_postuninstall_cmds=$lt_old_postuninstall_cmds

# Create an old-style archive from a shared archive.
old_archive_from_new_cmds=$lt_[]_LT_AC_TAGVAR(old_archive_from_new_cmds, $1)

# Create a temporary old-style archive to link instead of a shared archive.
old_archive_from_expsyms_cmds=$lt_[]_LT_AC_TAGVAR(old_archive_from_expsyms_cmds, $1)

# Commands used to build and install a shared archive.
archive_cmds=$lt_[]_LT_AC_TAGVAR(archive_cmds, $1)
archive_expsym_cmds=$lt_[]_LT_AC_TAGVAR(archive_expsym_cmds, $1)
postinstall_cmds=$lt_postinstall_cmds
postuninstall_cmds=$lt_postuninstall_cmds

# Commands used to build a loadable module (assumed same as above if empty)
module_cmds=$lt_[]_LT_AC_TAGVAR(module_cmds, $1)
module_expsym_cmds=$lt_[]_LT_AC_TAGVAR(module_expsym_cmds, $1)

# Commands to strip libraries.
old_striplib=$lt_old_striplib
striplib=$lt_striplib

# Dependencies to place before the objects being linked to create a
# shared library.
predep_objects=$lt_[]_LT_AC_TAGVAR(predep_objects, $1)

# Dependencies to place after the objects being linked to create a
# shared library.
postdep_objects=$lt_[]_LT_AC_TAGVAR(postdep_objects, $1)

# Dependencies to place before the objects being linked to create a
# shared library.
predeps=$lt_[]_LT_AC_TAGVAR(predeps, $1)

# Dependencies to place after the objects being linked to create a
# shared library.
postdeps=$lt_[]_LT_AC_TAGVAR(postdeps, $1)

# The library search path used internally by the compiler when linking
# a shared library.
compiler_lib_search_path=$lt_[]_LT_AC_TAGVAR(compiler_lib_search_path, $1)

# Method to check whether dependent libraries are shared objects.
deplibs_check_method=$lt_deplibs_check_method

# Command to use when deplibs_check_method == file_magic.
file_magic_cmd=$lt_file_magic_cmd

# Flag that allows shared libraries with undefined symbols to be built.
allow_undefined_flag=$lt_[]_LT_AC_TAGVAR(allow_undefined_flag, $1)

# Flag that forces no undefined symbols.
no_undefined_flag=$lt_[]_LT_AC_TAGVAR(no_undefined_flag, $1)

# Commands used to finish a libtool library installation in a directory.
finish_cmds=$lt_finish_cmds

# Same as above, but a single script fragment to be evaled but not shown.
finish_eval=$lt_finish_eval

# Take the output of nm and produce a listing of raw symbols and C names.
global_symbol_pipe=$lt_lt_cv_sys_global_symbol_pipe

# Transform the output of nm in a proper C declaration
global_symbol_to_cdecl=$lt_lt_cv_sys_global_symbol_to_cdecl

# Transform the output of nm in a C name address pair
global_symbol_to_c_name_address=$lt_lt_cv_sys_global_symbol_to_c_name_address

# This is the shared library runtime path variable.
runpath_var=$runpath_var

# This is the shared library path variable.
shlibpath_var=$shlibpath_var

# Is shlibpath searched before the hard-coded library search path?
shlibpath_overrides_runpath=$shlibpath_overrides_runpath

# How to hardcode a shared library path into an executable.
hardcode_action=$_LT_AC_TAGVAR(hardcode_action, $1)

# Whether we should hardcode library paths into libraries.
hardcode_into_libs=$hardcode_into_libs

# Flag to hardcode \$libdir into a binary during linking.
# This must work even if \$libdir does not exist.
hardcode_libdir_flag_spec=$lt_[]_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)

# If ld is used when linking, flag to hardcode \$libdir into
# a binary during linking. This must work even if \$libdir does
# not exist.
hardcode_libdir_flag_spec_ld=$lt_[]_LT_AC_TAGVAR(hardcode_libdir_flag_spec_ld, $1)

# Whether we need a single -rpath flag with a separated argument.
hardcode_libdir_separator=$lt_[]_LT_AC_TAGVAR(hardcode_libdir_separator, $1)

# Set to yes if using DIR/libNAME${shared_ext} during linking hardcodes DIR into the
# resulting binary.
hardcode_direct=$_LT_AC_TAGVAR(hardcode_direct, $1)

# Set to yes if using the -LDIR flag during linking hardcodes DIR into the
# resulting binary.
hardcode_minus_L=$_LT_AC_TAGVAR(hardcode_minus_L, $1)

# Set to yes if using SHLIBPATH_VAR=DIR during linking hardcodes DIR into
# the resulting binary.
hardcode_shlibpath_var=$_LT_AC_TAGVAR(hardcode_shlibpath_var, $1)

# Set to yes if building a shared library automatically hardcodes DIR into the library
# and all subsequent libraries and executables linked against it.
hardcode_automatic=$_LT_AC_TAGVAR(hardcode_automatic, $1)

# Variables whose values should be saved in libtool wrapper scripts and
# restored at relink time.
variables_saved_for_relink="$variables_saved_for_relink"

# Whether libtool must link a program against all its dependency libraries.
link_all_deplibs=$_LT_AC_TAGVAR(link_all_deplibs, $1)

# Compile-time system search path for libraries
sys_lib_search_path_spec=$lt_sys_lib_search_path_spec

# Run-time system search path for libraries
sys_lib_dlsearch_path_spec=$lt_sys_lib_dlsearch_path_spec

# Fix the shell variable \$srcfile for the compiler.
fix_srcfile_path="$_LT_AC_TAGVAR(fix_srcfile_path, $1)"

# Set to yes if exported symbols are required.
always_export_symbols=$_LT_AC_TAGVAR(always_export_symbols, $1)

# The commands to list exported symbols.
export_symbols_cmds=$lt_[]_LT_AC_TAGVAR(export_symbols_cmds, $1)

# The commands to extract the exported symbol list from a shared archive.
extract_expsyms_cmds=$lt_extract_expsyms_cmds

# Symbols that should not be listed in the preloaded symbols.
exclude_expsyms=$lt_[]_LT_AC_TAGVAR(exclude_expsyms, $1)

# Symbols that must always be exported.
include_expsyms=$lt_[]_LT_AC_TAGVAR(include_expsyms, $1)

ifelse([$1],[],
[# ### END LIBTOOL CONFIG],
[# ### END LIBTOOL TAG CONFIG: $tagname])

__EOF__

ifelse([$1],[], [
  case $host_os in
  aix3*)
    cat <<\EOF >> "$cfgfile"

# AIX sometimes has problems with the GCC collect2 program.  For some
# reason, if we set the COLLECT_NAMES environment variable, the problems
# vanish in a puff of smoke.
if test "X${COLLECT_NAMES+set}" != Xset; then
  COLLECT_NAMES=
  export COLLECT_NAMES
fi
EOF
    ;;
  esac

  # We use sed instead of cat because bash on DJGPP gets confused if
  # if finds mixed CR/LF and LF-only lines.  Since sed operates in
  # text mode, it properly converts lines to CR/LF.  This bash problem
  # is reportedly fixed, but why not run on old versions too?
  sed '$q' "$ltmain" >> "$cfgfile" || (rm -f "$cfgfile"; exit 1)

  mv -f "$cfgfile" "$ofile" || \
    (rm -f "$ofile" && cp "$cfgfile" "$ofile" && rm -f "$cfgfile")
  chmod +x "$ofile"
])
else
  # If there is no Makefile yet, we rely on a make rule to execute
  # `config.status --recheck' to rerun these tests and create the
  # libtool script then.
  ltmain_in=`echo $ltmain | sed -e 's/\.sh$/.in/'`
  if test -f "$ltmain_in"; then
    test -f Makefile && make "$ltmain"
  fi
fi
])# AC_LIBTOOL_CONFIG


# AC_LIBTOOL_PROG_COMPILER_NO_RTTI([TAGNAME])
# -------------------------------------------
AC_DEFUN([AC_LIBTOOL_PROG_COMPILER_NO_RTTI],
[AC_REQUIRE([_LT_AC_SYS_COMPILER])dnl

_LT_AC_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)=

if test "$GCC" = yes; then
  _LT_AC_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)=' -fno-builtin'

  AC_LIBTOOL_COMPILER_OPTION([if $compiler supports -fno-rtti -fno-exceptions],
    lt_cv_prog_compiler_rtti_exceptions,
    [-fno-rtti -fno-exceptions], [],
    [_LT_AC_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)="$_LT_AC_TAGVAR(lt_prog_compiler_no_builtin_flag, $1) -fno-rtti -fno-exceptions"])
fi
])# AC_LIBTOOL_PROG_COMPILER_NO_RTTI


# AC_LIBTOOL_SYS_GLOBAL_SYMBOL_PIPE
# ---------------------------------
AC_DEFUN([AC_LIBTOOL_SYS_GLOBAL_SYMBOL_PIPE],
[AC_REQUIRE([AC_CANONICAL_HOST])
AC_REQUIRE([AC_PROG_NM])
AC_REQUIRE([AC_OBJEXT])
# Check for command to grab the raw symbol name followed by C symbol from nm.
AC_MSG_CHECKING([command to parse $NM output from $compiler object])
AC_CACHE_VAL([lt_cv_sys_global_symbol_pipe],
[
# These are sane defaults that work on at least a few old systems.
# [They come from Ultrix.  What could be older than Ultrix?!! ;)]

# Character class describing NM global symbol codes.
symcode='[[BCDEGRST]]'

# Regexp to match symbols that can be accessed directly from C.
sympat='\([[_A-Za-z]][[_A-Za-z0-9]]*\)'

# Transform an extracted symbol line into a proper C declaration
lt_cv_sys_global_symbol_to_cdecl="sed -n -e 's/^. .* \(.*\)$/extern int \1;/p'"

# Transform an extracted symbol line into symbol name and symbol address
lt_cv_sys_global_symbol_to_c_name_address="sed -n -e 's/^: \([[^ ]]*\) $/  {\\\"\1\\\", (lt_ptr) 0},/p' -e 's/^$symcode \([[^ ]]*\) \([[^ ]]*\)$/  {\"\2\", (lt_ptr) \&\2},/p'"

# Define system-specific variables.
case $host_os in
aix*)
  symcode='[[BCDT]]'
  ;;
cygwin* | mingw* | pw32*)
  symcode='[[ABCDGISTW]]'
  ;;
hpux*) # Its linker distinguishes data from code symbols
  if test "$host_cpu" = ia64; then
    symcode='[[ABCDEGRST]]'
  fi
  lt_cv_sys_global_symbol_to_cdecl="sed -n -e 's/^T .* \(.*\)$/extern int \1();/p' -e 's/^$symcode* .* \(.*\)$/extern char \1;/p'"
  lt_cv_sys_global_symbol_to_c_name_address="sed -n -e 's/^: \([[^ ]]*\) $/  {\\\"\1\\\", (lt_ptr) 0},/p' -e 's/^$symcode* \([[^ ]]*\) \([[^ ]]*\)$/  {\"\2\", (lt_ptr) \&\2},/p'"
  ;;
linux*)
  if test "$host_cpu" = ia64; then
    symcode='[[ABCDGIRSTW]]'
    lt_cv_sys_global_symbol_to_cdecl="sed -n -e 's/^T .* \(.*\)$/extern int \1();/p' -e 's/^$symcode* .* \(.*\)$/extern char \1;/p'"
    lt_cv_sys_global_symbol_to_c_name_address="sed -n -e 's/^: \([[^ ]]*\) $/  {\\\"\1\\\", (lt_ptr) 0},/p' -e 's/^$symcode* \([[^ ]]*\) \([[^ ]]*\)$/  {\"\2\", (lt_ptr) \&\2},/p'"
  fi
  ;;
irix* | nonstopux*)
  symcode='[[BCDEGRST]]'
  ;;
osf*)
  symcode='[[BCDEGQRST]]'
  ;;
solaris*)
  symcode='[[BDRT]]'
  ;;
sco3.2v5*)
  symcode='[[DT]]'
  ;;
sysv4.2uw2*)
  symcode='[[DT]]'
  ;;
sysv5* | sco5v6* | unixware* | OpenUNIX*)
  symcode='[[ABDT]]'
  ;;
sysv4)
  symcode='[[DFNSTU]]'
  ;;
esac

# Handle CRLF in mingw tool chain
opt_cr=
case $build_os in
mingw*)
  opt_cr=`echo 'x\{0,1\}' | tr x '\015'` # option cr in regexp
  ;;
esac

# If we're using GNU nm, then use its standard symbol codes.
case `$NM -V 2>&1` in
*GNU* | *'with BFD'*)
  symcode='[[ABCDGIRSTW]]' ;;
esac

# Try without a prefix undercore, then with it.
for ac_symprfx in "" "_"; do

  # Transform symcode, sympat, and symprfx into a raw symbol and a C symbol.
  symxfrm="\\1 $ac_symprfx\\2 \\2"

  # Write the raw and C identifiers.
  lt_cv_sys_global_symbol_pipe="sed -n -e 's/^.*[[ 	]]\($symcode$symcode*\)[[ 	]][[ 	]]*$ac_symprfx$sympat$opt_cr$/$symxfrm/p'"

  # Check to see that the pipe works correctly.
  pipe_works=no

  rm -f conftest*
  cat > conftest.$ac_ext <<EOF
#ifdef __cplusplus
extern "C" {
#endif
char nm_test_var;
void nm_test_func(){}
#ifdef __cplusplus
}
#endif
int main(){nm_test_var='a';nm_test_func();return(0);}
EOF

  if AC_TRY_EVAL(ac_compile); then
    # Now try to grab the symbols.
    nlist=conftest.nm
    if AC_TRY_EVAL(NM conftest.$ac_objext \| $lt_cv_sys_global_symbol_pipe \> $nlist) && test -s "$nlist"; then
      # Try sorting and uniquifying the output.
      if sort "$nlist" | uniq > "$nlist"T; then
	mv -f "$nlist"T "$nlist"
      else
	rm -f "$nlist"T
      fi

      # Make sure that we snagged all the symbols we need.
      if grep ' nm_test_var$' "$nlist" >/dev/null; then
	if grep ' nm_test_func$' "$nlist" >/dev/null; then
	  cat <<EOF > conftest.$ac_ext
#ifdef __cplusplus
extern "C" {
#endif

EOF
	  # Now generate the symbol file.
	  eval "$lt_cv_sys_global_symbol_to_cdecl"' < "$nlist" | grep -v main >> conftest.$ac_ext'

	  cat <<EOF >> conftest.$ac_ext
#if defined (__STDC__) && __STDC__
# define lt_ptr_t void *
#else
# define lt_ptr_t char *
# define const
#endif

/* The mapping between symbol names and symbols. */
const struct {
  const char *name;
  lt_ptr_t address;
}
lt_preloaded_symbols[[]] =
{
EOF
	  $SED "s/^$symcode$symcode* \(.*\) \(.*\)$/  {\"\2\", (lt_ptr_t) \&\2},/" < "$nlist" | grep -v main >> conftest.$ac_ext
	  cat <<\EOF >> conftest.$ac_ext
  {0, (lt_ptr_t) 0}
};

#ifdef __cplusplus
}
#endif
EOF
	  # Now try linking the two files.
	  mv conftest.$ac_objext conftstm.$ac_objext
	  lt_save_LIBS="$LIBS"
	  lt_save_CFLAGS="$CFLAGS"
	  LIBS="conftstm.$ac_objext"
	  CFLAGS="$CFLAGS$_LT_AC_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)"
	  if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
	    pipe_works=yes
	  fi
	  LIBS="$lt_save_LIBS"
	  CFLAGS="$lt_save_CFLAGS"
	else
	  echo "cannot find nm_test_func in $nlist" >&AS_MESSAGE_LOG_FD
	fi
      else
	echo "cannot find nm_test_var in $nlist" >&AS_MESSAGE_LOG_FD
      fi
    else
      echo "cannot run $lt_cv_sys_global_symbol_pipe" >&AS_MESSAGE_LOG_FD
    fi
  else
    echo "$progname: failed program was:" >&AS_MESSAGE_LOG_FD
    cat conftest.$ac_ext >&5
  fi
  rm -f conftest* conftst*

  # Do not use the global_symbol_pipe unless it works.
  if test "$pipe_works" = yes; then
    break
  else
    lt_cv_sys_global_symbol_pipe=
  fi
done
])
if test -z "$lt_cv_sys_global_symbol_pipe"; then
  lt_cv_sys_global_symbol_to_cdecl=
fi
if test -z "$lt_cv_sys_global_symbol_pipe$lt_cv_sys_global_symbol_to_cdecl"; then
  AC_MSG_RESULT(failed)
else
  AC_MSG_RESULT(ok)
fi
]) # AC_LIBTOOL_SYS_GLOBAL_SYMBOL_PIPE


# AC_LIBTOOL_PROG_COMPILER_PIC([TAGNAME])
# ---------------------------------------
AC_DEFUN([AC_LIBTOOL_PROG_COMPILER_PIC],
[_LT_AC_TAGVAR(lt_prog_compiler_wl, $1)=
_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)=
_LT_AC_TAGVAR(lt_prog_compiler_static, $1)=

AC_MSG_CHECKING([for $compiler option to produce PIC])
 ifelse([$1],[CXX],[
  # C++ specific cases for pic, static, wl, etc.
  if test "$GXX" = yes; then
    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-static'

    case $host_os in
    aix*)
      # All AIX code is PIC.
      if test "$host_cpu" = ia64; then
	# AIX 5 now supports IA64 processor
	_LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      fi
      ;;
    amigaos*)
      # FIXME: we need at least 68020 code to build shared libraries, but
      # adding the `-m68020' flag to GCC prevents building anything better,
      # like `-m68040'.
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-m68020 -resident32 -malways-restore-a4'
      ;;
    beos* | cygwin* | irix5* | irix6* | nonstopux* | osf3* | osf4* | osf5*)
      # PIC is the default for these OSes.
      ;;
    mingw* | os2* | pw32*)
      # This hack is so that the source file can tell whether it is being
      # built for inclusion in a dll (and should export symbols for example).
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-DDLL_EXPORT'
      ;;
    darwin* | rhapsody*)
      # PIC is the default on this platform
      # Common symbols not allowed in MH_DYLIB files
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-fno-common'
      ;;
    *djgpp*)
      # DJGPP does not support shared libraries at all
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)=
      ;;
    interix3*)
      # Interix 3.x gcc -fpic/-fPIC options generate broken code.
      # Instead, we relocate shared libraries at runtime.
      ;;
    sysv4*MP*)
      if test -d /usr/nec; then
	_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)=-Kconform_pic
      fi
      ;;
    hpux*)
      # PIC is the default for IA64 HP-UX and 64-bit HP-UX, but
      # not for PA HP-UX.
      case $host_cpu in
      hppa*64*|ia64*)
	;;
      *)
	_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
	;;
      esac
      ;;
    *)
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
      ;;
    esac
  else
    case $host_os in
      aix4* | aix5*)
	# All AIX code is PIC.
	if test "$host_cpu" = ia64; then
	  # AIX 5 now supports IA64 processor
	  _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	else
	  _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-bnso -bI:/lib/syscalls.exp'
	fi
	;;
      chorus*)
	case $cc_basename in
	cxch68*)
	  # Green Hills C++ Compiler
	  # _LT_AC_TAGVAR(lt_prog_compiler_static, $1)="--no_auto_instantiation -u __main -u __premain -u _abort -r $COOL_DIR/lib/libOrb.a $MVME_DIR/lib/CC/libC.a $MVME_DIR/lib/classix/libcx.s.a"
	  ;;
	esac
	;;
       darwin*)
         # PIC is the default on this platform
         # Common symbols not allowed in MH_DYLIB files
         case $cc_basename in
           xlc*)
           _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-qnocommon'
           _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
           ;;
         esac
       ;;
      dgux*)
	case $cc_basename in
	  ec++*)
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	    ;;
	  ghcx*)
	    # Green Hills C++ Compiler
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-pic'
	    ;;
	  *)
	    ;;
	esac
	;;
      freebsd* | kfreebsd*-gnu | dragonfly*)
	# FreeBSD uses GNU C++
	;;
      hpux9* | hpux10* | hpux11*)
	case $cc_basename in
	  CC*)
	    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='${wl}-a ${wl}archive'
	    if test "$host_cpu" != ia64; then
	      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='+Z'
	    fi
	    ;;
	  aCC*)
	    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='${wl}-a ${wl}archive'
	    case $host_cpu in
	    hppa*64*|ia64*)
	      # +Z the default
	      ;;
	    *)
	      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='+Z'
	      ;;
	    esac
	    ;;
	  *)
	    ;;
	esac
	;;
      interix*)
	# This is c89, which is MS Visual C++ (no shared libs)
	# Anyone wants to do a port?
	;;
      irix5* | irix6* | nonstopux*)
	case $cc_basename in
	  CC*)
	    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
	    # CC pic flag -KPIC is the default.
	    ;;
	  *)
	    ;;
	esac
	;;
      linux*)
	case $cc_basename in
	  KCC*)
	    # KAI C++ Compiler
	    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='--backend -Wl,'
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
	    ;;
	  icpc* | ecpc*)
	    # Intel C++
	    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-static'
	    ;;
	  pgCC*)
	    # Portland Group C++ compiler.
	    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-fpic'
	    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	    ;;
	  cxx*)
	    # Compaq C++
	    # Make sure the PIC flag is empty.  It appears that all Alpha
	    # Linux and Compaq Tru64 Unix objects are PIC.
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)=
	    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
	    ;;
	  *)
	    ;;
	esac
	;;
      lynxos*)
	;;
      m88k*)
	;;
      mvs*)
	case $cc_basename in
	  cxx*)
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-W c,exportall'
	    ;;
	  *)
	    ;;
	esac
	;;
      netbsd*)
	;;
      osf3* | osf4* | osf5*)
	case $cc_basename in
	  KCC*)
	    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='--backend -Wl,'
	    ;;
	  RCC*)
	    # Rational C++ 2.4.1
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-pic'
	    ;;
	  cxx*)
	    # Digital/Compaq C++
	    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    # Make sure the PIC flag is empty.  It appears that all Alpha
	    # Linux and Compaq Tru64 Unix objects are PIC.
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)=
	    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
	    ;;
	  *)
	    ;;
	esac
	;;
      psos*)
	;;
      solaris*)
	case $cc_basename in
	  CC*)
	    # Sun C++ 4.2, 5.x and Centerline C++
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Qoption ld '
	    ;;
	  gcx*)
	    # Green Hills C++ Compiler
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-PIC'
	    ;;
	  *)
	    ;;
	esac
	;;
      sunos4*)
	case $cc_basename in
	  CC*)
	    # Sun C++ 4.x
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-pic'
	    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	    ;;
	  lcc*)
	    # Lucid
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-pic'
	    ;;
	  *)
	    ;;
	esac
	;;
      tandem*)
	case $cc_basename in
	  NCC*)
	    # NonStop-UX NCC 3.20
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	    ;;
	  *)
	    ;;
	esac
	;;
      sysv5* | unixware* | sco3.2v5* | sco5v6* | OpenUNIX*)
	case $cc_basename in
	  CC*)
	    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	    ;;
	esac
	;;
      vxworks*)
	;;
      *)
	_LT_AC_TAGVAR(lt_prog_compiler_can_build_shared, $1)=no
	;;
    esac
  fi
],
[
  if test "$GCC" = yes; then
    _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
    _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-static'

    case $host_os in
      aix*)
      # All AIX code is PIC.
      if test "$host_cpu" = ia64; then
	# AIX 5 now supports IA64 processor
	_LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      fi
      ;;

    amigaos*)
      # FIXME: we need at least 68020 code to build shared libraries, but
      # adding the `-m68020' flag to GCC prevents building anything better,
      # like `-m68040'.
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-m68020 -resident32 -malways-restore-a4'
      ;;

    beos* | cygwin* | irix5* | irix6* | nonstopux* | osf3* | osf4* | osf5*)
      # PIC is the default for these OSes.
      ;;

    mingw* | pw32* | os2*)
      # This hack is so that the source file can tell whether it is being
      # built for inclusion in a dll (and should export symbols for example).
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-DDLL_EXPORT'
      ;;

    darwin* | rhapsody*)
      # PIC is the default on this platform
      # Common symbols not allowed in MH_DYLIB files
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-fno-common'
      ;;

    interix3*)
      # Interix 3.x gcc -fpic/-fPIC options generate broken code.
      # Instead, we relocate shared libraries at runtime.
      ;;

    msdosdjgpp*)
      # Just because we use GCC doesn't mean we suddenly get shared libraries
      # on systems that don't support them.
      _LT_AC_TAGVAR(lt_prog_compiler_can_build_shared, $1)=no
      enable_shared=no
      ;;

    sysv4*MP*)
      if test -d /usr/nec; then
	_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)=-Kconform_pic
      fi
      ;;

    hpux*)
      # PIC is the default for IA64 HP-UX and 64-bit HP-UX, but
      # not for PA HP-UX.
      case $host_cpu in
      hppa*64*|ia64*)
	# +Z the default
	;;
      *)
	_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
	;;
      esac
      ;;

    *)
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
      ;;
    esac
  else
    # PORTME Check for flag to pass linker flags through the system compiler.
    case $host_os in
    aix*)
      _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      if test "$host_cpu" = ia64; then
	# AIX 5 now supports IA64 processor
	_LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      else
	_LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-bnso -bI:/lib/syscalls.exp'
      fi
      ;;
      darwin*)
        # PIC is the default on this platform
        # Common symbols not allowed in MH_DYLIB files
       case $cc_basename in
         xlc*)
         _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-qnocommon'
         _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
         ;;
       esac
       ;;

    mingw* | pw32* | os2*)
      # This hack is so that the source file can tell whether it is being
      # built for inclusion in a dll (and should export symbols for example).
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-DDLL_EXPORT'
      ;;

    hpux9* | hpux10* | hpux11*)
      _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      # PIC is the default for IA64 HP-UX and 64-bit HP-UX, but
      # not for PA HP-UX.
      case $host_cpu in
      hppa*64*|ia64*)
	# +Z the default
	;;
      *)
	_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='+Z'
	;;
      esac
      # Is there a better lt_prog_compiler_static that works with the bundled CC?
      _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='${wl}-a ${wl}archive'
      ;;

    irix5* | irix6* | nonstopux*)
      _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      # PIC (with -KPIC) is the default.
      _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
      ;;

    newsos6)
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
      _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      ;;

    linux*)
      case $cc_basename in
      icc* | ecc*)
	_LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	_LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-static'
        ;;
      pgcc* | pgf77* | pgf90* | pgf95*)
        # Portland Group compilers (*not* the Pentium gcc compiler,
	# which looks to be a dead project)
	_LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-fpic'
	_LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
        ;;
      ccc*)
        _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
        # All Alpha code is PIC.
        _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
        ;;
      esac
      ;;

    osf3* | osf4* | osf5*)
      _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      # All OSF/1 code is PIC.
      _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
      ;;

    solaris*)
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
      _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      case $cc_basename in
      f77* | f90* | f95*)
	_LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Qoption ld ';;
      *)
	_LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,';;
      esac
      ;;

    sunos4*)
      _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Qoption ld '
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-PIC'
      _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      ;;

    sysv4 | sysv4.2uw2* | sysv4.3*)
      _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
      _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      ;;

    sysv4*MP*)
      if test -d /usr/nec ;then
	_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-Kconform_pic'
	_LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      fi
      ;;

    sysv5* | unixware* | sco3.2v5* | sco5v6* | OpenUNIX*)
      _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
      _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      ;;

    unicos*)
      _LT_AC_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      _LT_AC_TAGVAR(lt_prog_compiler_can_build_shared, $1)=no
      ;;

    uts4*)
      _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)='-pic'
      _LT_AC_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      ;;

    *)
      _LT_AC_TAGVAR(lt_prog_compiler_can_build_shared, $1)=no
      ;;
    esac
  fi
])
AC_MSG_RESULT([$_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)])

#
# Check to make sure the PIC flag actually works.
#
if test -n "$_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)"; then
  AC_LIBTOOL_COMPILER_OPTION([if $compiler PIC flag $_LT_AC_TAGVAR(lt_prog_compiler_pic, $1) works],
    _LT_AC_TAGVAR(lt_prog_compiler_pic_works, $1),
    [$_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)ifelse([$1],[],[ -DPIC],[ifelse([$1],[CXX],[ -DPIC],[])])], [],
    [case $_LT_AC_TAGVAR(lt_prog_compiler_pic, $1) in
     "" | " "*) ;;
     *) _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)=" $_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)" ;;
     esac],
    [_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)=
     _LT_AC_TAGVAR(lt_prog_compiler_can_build_shared, $1)=no])
fi
case $host_os in
  # For platforms which do not support PIC, -DPIC is meaningless:
  *djgpp*)
    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)=
    ;;
  *)
    _LT_AC_TAGVAR(lt_prog_compiler_pic, $1)="$_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)ifelse([$1],[],[ -DPIC],[ifelse([$1],[CXX],[ -DPIC],[])])"
    ;;
esac

#
# Check to make sure the static flag actually works.
#
wl=$_LT_AC_TAGVAR(lt_prog_compiler_wl, $1) eval lt_tmp_static_flag=\"$_LT_AC_TAGVAR(lt_prog_compiler_static, $1)\"
AC_LIBTOOL_LINKER_OPTION([if $compiler static flag $lt_tmp_static_flag works],
  _LT_AC_TAGVAR(lt_prog_compiler_static_works, $1),
  $lt_tmp_static_flag,
  [],
  [_LT_AC_TAGVAR(lt_prog_compiler_static, $1)=])
])


# AC_LIBTOOL_PROG_LD_SHLIBS([TAGNAME])
# ------------------------------------
# See if the linker supports building shared libraries.
AC_DEFUN([AC_LIBTOOL_PROG_LD_SHLIBS],
[AC_MSG_CHECKING([whether the $compiler linker ($LD) supports shared libraries])
ifelse([$1],[CXX],[
  _LT_AC_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED '\''s/.* //'\'' | sort | uniq > $export_symbols'
  case $host_os in
  aix4* | aix5*)
    # If we're using GNU nm, then we don't want the "-C" option.
    # -C means demangle to AIX nm, but means don't demangle with GNU nm
    if $NM -V 2>&1 | grep 'GNU' > /dev/null; then
      _LT_AC_TAGVAR(export_symbols_cmds, $1)='$NM -Bpg $libobjs $convenience | awk '\''{ if (((\[$]2 == "T") || (\[$]2 == "D") || (\[$]2 == "B")) && ([substr](\[$]3,1,1) != ".")) { print \[$]3 } }'\'' | sort -u > $export_symbols'
    else
      _LT_AC_TAGVAR(export_symbols_cmds, $1)='$NM -BCpg $libobjs $convenience | awk '\''{ if (((\[$]2 == "T") || (\[$]2 == "D") || (\[$]2 == "B")) && ([substr](\[$]3,1,1) != ".")) { print \[$]3 } }'\'' | sort -u > $export_symbols'
    fi
    ;;
  pw32*)
    _LT_AC_TAGVAR(export_symbols_cmds, $1)="$ltdll_cmds"
  ;;
  cygwin* | mingw*)
    _LT_AC_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED -e '\''/^[[BCDGRS]] /s/.* \([[^ ]]*\)/\1 DATA/;/^.* __nm__/s/^.* __nm__\([[^ ]]*\) [[^ ]]*/\1 DATA/;/^I /d;/^[[AITW]] /s/.* //'\'' | sort | uniq > $export_symbols'
  ;;
  *)
    _LT_AC_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED '\''s/.* //'\'' | sort | uniq > $export_symbols'
  ;;
  esac
],[
  runpath_var=
  _LT_AC_TAGVAR(allow_undefined_flag, $1)=
  _LT_AC_TAGVAR(enable_shared_with_static_runtimes, $1)=no
  _LT_AC_TAGVAR(archive_cmds, $1)=
  _LT_AC_TAGVAR(archive_expsym_cmds, $1)=
  _LT_AC_TAGVAR(old_archive_From_new_cmds, $1)=
  _LT_AC_TAGVAR(old_archive_from_expsyms_cmds, $1)=
  _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)=
  _LT_AC_TAGVAR(whole_archive_flag_spec, $1)=
  _LT_AC_TAGVAR(thread_safe_flag_spec, $1)=
  _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)=
  _LT_AC_TAGVAR(hardcode_libdir_flag_spec_ld, $1)=
  _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=
  _LT_AC_TAGVAR(hardcode_direct, $1)=no
  _LT_AC_TAGVAR(hardcode_minus_L, $1)=no
  _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=unsupported
  _LT_AC_TAGVAR(link_all_deplibs, $1)=unknown
  _LT_AC_TAGVAR(hardcode_automatic, $1)=no
  _LT_AC_TAGVAR(module_cmds, $1)=
  _LT_AC_TAGVAR(module_expsym_cmds, $1)=
  _LT_AC_TAGVAR(always_export_symbols, $1)=no
  _LT_AC_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED '\''s/.* //'\'' | sort | uniq > $export_symbols'
  # include_expsyms should be a list of space-separated symbols to be *always*
  # included in the symbol list
  _LT_AC_TAGVAR(include_expsyms, $1)=
  # exclude_expsyms can be an extended regexp of symbols to exclude
  # it will be wrapped by ` (' and `)$', so one must not match beginning or
  # end of line.  Example: `a|bc|.*d.*' will exclude the symbols `a' and `bc',
  # as well as any symbol that contains `d'.
  _LT_AC_TAGVAR(exclude_expsyms, $1)="_GLOBAL_OFFSET_TABLE_"
  # Although _GLOBAL_OFFSET_TABLE_ is a valid symbol C name, most a.out
  # platforms (ab)use it in PIC code, but their linkers get confused if
  # the symbol is explicitly referenced.  Since portable code cannot
  # rely on this symbol name, it's probably fine to never include it in
  # preloaded symbol tables.
  extract_expsyms_cmds=
  # Just being paranoid about ensuring that cc_basename is set.
  _LT_CC_BASENAME([$compiler])
  case $host_os in
  cygwin* | mingw* | pw32*)
    # FIXME: the MSVC++ port hasn't been tested in a loooong time
    # When not using gcc, we currently assume that we are using
    # Microsoft Visual C++.
    if test "$GCC" != yes; then
      with_gnu_ld=no
    fi
    ;;
  interix*)
    # we just hope/assume this is gcc and not c89 (= MSVC++)
    with_gnu_ld=yes
    ;;
  openbsd*)
    with_gnu_ld=no
    ;;
  esac

  _LT_AC_TAGVAR(ld_shlibs, $1)=yes
  if test "$with_gnu_ld" = yes; then
    # If archive_cmds runs LD, not CC, wlarc should be empty
    wlarc='${wl}'

    # Set some defaults for GNU ld with shared library support. These
    # are reset later if shared libraries are not supported. Putting them
    # here allows them to be overridden if necessary.
    runpath_var=LD_RUN_PATH
    _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}--rpath ${wl}$libdir'
    _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}--export-dynamic'
    # ancient GNU ld didn't support --whole-archive et. al.
    if $LD --help 2>&1 | grep 'no-whole-archive' > /dev/null; then
	_LT_AC_TAGVAR(whole_archive_flag_spec, $1)="$wlarc"'--whole-archive$convenience '"$wlarc"'--no-whole-archive'
      else
  	_LT_AC_TAGVAR(whole_archive_flag_spec, $1)=
    fi
    supports_anon_versioning=no
    case `$LD -v 2>/dev/null` in
      *\ [[01]].* | *\ 2.[[0-9]].* | *\ 2.10.*) ;; # catch versions < 2.11
      *\ 2.11.93.0.2\ *) supports_anon_versioning=yes ;; # RH7.3 ...
      *\ 2.11.92.0.12\ *) supports_anon_versioning=yes ;; # Mandrake 8.2 ...
      *\ 2.11.*) ;; # other 2.11 versions
      *) supports_anon_versioning=yes ;;
    esac

    # See if GNU ld supports shared libraries.
    case $host_os in
    aix3* | aix4* | aix5*)
      # On AIX/PPC, the GNU linker is very broken
      if test "$host_cpu" != ia64; then
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	cat <<EOF 1>&2

*** Warning: the GNU linker, at least up to release 2.9.1, is reported
*** to be unable to reliably create shared libraries on AIX.
*** Therefore, libtool is disabling shared libraries support.  If you
*** really care for shared libraries, you may want to modify your PATH
*** so that a non-GNU linker is found, and then restart.

EOF
      fi
      ;;

    amigaos*)
      _LT_AC_TAGVAR(archive_cmds, $1)='$rm $output_objdir/a2ixlibrary.data~$echo "#define NAME $libname" > $output_objdir/a2ixlibrary.data~$echo "#define LIBRARY_ID 1" >> $output_objdir/a2ixlibrary.data~$echo "#define VERSION $major" >> $output_objdir/a2ixlibrary.data~$echo "#define REVISION $revision" >> $output_objdir/a2ixlibrary.data~$AR $AR_FLAGS $lib $libobjs~$RANLIB $lib~(cd $output_objdir && a2ixlibrary -32)'
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes

      # Samuel A. Falvo II <kc5tja@dolphin.openprojects.net> reports
      # that the semantics of dynamic libraries on AmigaOS, at least up
      # to version 4, is to share data among multiple programs linked
      # with the same dynamic library.  Since this doesn't match the
      # behavior of shared libraries on other platforms, we can't use
      # them.
      _LT_AC_TAGVAR(ld_shlibs, $1)=no
      ;;

    beos*)
      if $LD --help 2>&1 | grep ': supported targets:.* elf' > /dev/null; then
	_LT_AC_TAGVAR(allow_undefined_flag, $1)=unsupported
	# Joseph Beckenbach <jrb3@best.com> says some releases of gcc
	# support --undefined.  This deserves some investigation.  FIXME
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -nostart $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'
      else
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;

    cygwin* | mingw* | pw32*)
      # _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1) is actually meaningless,
      # as there is no search path for DLLs.
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_AC_TAGVAR(allow_undefined_flag, $1)=unsupported
      _LT_AC_TAGVAR(always_export_symbols, $1)=no
      _LT_AC_TAGVAR(enable_shared_with_static_runtimes, $1)=yes
      _LT_AC_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED -e '\''/^[[BCDGRS]] /s/.* \([[^ ]]*\)/\1 DATA/'\'' | $SED -e '\''/^[[AITW]] /s/.* //'\'' | sort | uniq > $export_symbols'

      if $LD --help 2>&1 | grep 'auto-import' > /dev/null; then
        _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags -o $output_objdir/$soname ${wl}--enable-auto-image-base -Xlinker --out-implib -Xlinker $lib'
	# If the export-symbols file already is a .def file (1st line
	# is EXPORTS), use it as is; otherwise, prepend...
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='if test "x`$SED 1q $export_symbols`" = xEXPORTS; then
	  cp $export_symbols $output_objdir/$soname.def;
	else
	  echo EXPORTS > $output_objdir/$soname.def;
	  cat $export_symbols >> $output_objdir/$soname.def;
	fi~
	$CC -shared $output_objdir/$soname.def $libobjs $deplibs $compiler_flags -o $output_objdir/$soname ${wl}--enable-auto-image-base -Xlinker --out-implib -Xlinker $lib'
      else
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;

    interix3*)
      _LT_AC_TAGVAR(hardcode_direct, $1)=no
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath,$libdir'
      _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-E'
      # Hack: On Interix 3.x, we cannot compile PIC because of a broken gcc.
      # Instead, shared libraries are loaded at an image base (0x10000000 by
      # default) and relocated if they conflict, which is a slow very memory
      # consuming and fragmenting process.  To avoid this, we pick a random,
      # 256 KiB-aligned image base between 0x50000000 and 0x6FFC0000 at link
      # time.  Moving up from 0x10000000 also allows more sbrk(2) space.
      _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags ${wl}-h,$soname ${wl}--image-base,`expr ${RANDOM-$$} % 4096 / 2 \* 262144 + 1342177280` -o $lib'
      _LT_AC_TAGVAR(archive_expsym_cmds, $1)='sed "s,^,_," $export_symbols >$output_objdir/$soname.expsym~$CC -shared $pic_flag $libobjs $deplibs $compiler_flags ${wl}-h,$soname ${wl}--retain-symbols-file,$output_objdir/$soname.expsym ${wl}--image-base,`expr ${RANDOM-$$} % 4096 / 2 \* 262144 + 1342177280` -o $lib'
      ;;

    linux*)
      if $LD --help 2>&1 | grep ': supported targets:.* elf' > /dev/null; then
	tmp_addflag=
	case $cc_basename,$host_cpu in
	pgcc*)				# Portland Group C compiler
	  _LT_AC_TAGVAR(whole_archive_flag_spec, $1)='${wl}--whole-archive`for conv in $convenience\"\"; do test  -n \"$conv\" && new_convenience=\"$new_convenience,$conv\"; done; $echo \"$new_convenience\"` ${wl}--no-whole-archive'
	  tmp_addflag=' $pic_flag'
	  ;;
	pgf77* | pgf90* | pgf95*)	# Portland Group f77 and f90 compilers
	  _LT_AC_TAGVAR(whole_archive_flag_spec, $1)='${wl}--whole-archive`for conv in $convenience\"\"; do test  -n \"$conv\" && new_convenience=\"$new_convenience,$conv\"; done; $echo \"$new_convenience\"` ${wl}--no-whole-archive'
	  tmp_addflag=' $pic_flag -Mnomain' ;;
	ecc*,ia64* | icc*,ia64*)		# Intel C compiler on ia64
	  tmp_addflag=' -i_dynamic' ;;
	efc*,ia64* | ifort*,ia64*)	# Intel Fortran compiler on ia64
	  tmp_addflag=' -i_dynamic -nofor_main' ;;
	ifc* | ifort*)			# Intel Fortran compiler
	  tmp_addflag=' -nofor_main' ;;
	esac
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared'"$tmp_addflag"' $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'

	if test $supports_anon_versioning = yes; then
	  _LT_AC_TAGVAR(archive_expsym_cmds, $1)='$echo "{ global:" > $output_objdir/$libname.ver~
  cat $export_symbols | sed -e "s/\(.*\)/\1;/" >> $output_objdir/$libname.ver~
  $echo "local: *; };" >> $output_objdir/$libname.ver~
	  $CC -shared'"$tmp_addflag"' $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname ${wl}-version-script ${wl}$output_objdir/$libname.ver -o $lib'
	fi
      else
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;

    netbsd*)
      if echo __ELF__ | $CC -E - | grep __ELF__ >/dev/null; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$LD -Bshareable $libobjs $deplibs $linker_flags -o $lib'
	wlarc=
      else
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname ${wl}-retain-symbols-file $wl$export_symbols -o $lib'
      fi
      ;;

    solaris*)
      if $LD -v 2>&1 | grep 'BFD 2\.8' > /dev/null; then
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	cat <<EOF 1>&2

*** Warning: The releases 2.8.* of the GNU linker cannot reliably
*** create shared libraries on Solaris systems.  Therefore, libtool
*** is disabling shared libraries support.  We urge you to upgrade GNU
*** binutils to release 2.9.1 or newer.  Another option is to modify
*** your PATH or compiler configuration so that the native linker is
*** used, and then restart.

EOF
      elif $LD --help 2>&1 | grep ': supported targets:.* elf' > /dev/null; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname ${wl}-retain-symbols-file $wl$export_symbols -o $lib'
      else
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;

    sysv5* | sco3.2v5* | sco5v6* | unixware* | OpenUNIX*)
      case `$LD -v 2>&1` in
        *\ [[01]].* | *\ 2.[[0-9]].* | *\ 2.1[[0-5]].*) 
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
	cat <<_LT_EOF 1>&2

*** Warning: Releases of the GNU linker prior to 2.16.91.0.3 can not
*** reliably create shared libraries on SCO systems.  Therefore, libtool
*** is disabling shared libraries support.  We urge you to upgrade GNU
*** binutils to release 2.16.91.0.3 or newer.  Another option is to modify
*** your PATH or compiler configuration so that the native linker is
*** used, and then restart.

_LT_EOF
	;;
	*)
	  if $LD --help 2>&1 | grep ': supported targets:.* elf' > /dev/null; then
	    _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='`test -z "$SCOABSPATH" && echo ${wl}-rpath,$libdir`'
	    _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname,\${SCOABSPATH:+${install_libdir}/}$soname -o $lib'
	    _LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname,\${SCOABSPATH:+${install_libdir}/}$soname,-retain-symbols-file,$export_symbols -o $lib'
	  else
	    _LT_AC_TAGVAR(ld_shlibs, $1)=no
	  fi
	;;
      esac
      ;;

    sunos4*)
      _LT_AC_TAGVAR(archive_cmds, $1)='$LD -assert pure-text -Bshareable -o $lib $libobjs $deplibs $linker_flags'
      wlarc=
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    *)
      if $LD --help 2>&1 | grep ': supported targets:.* elf' > /dev/null; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname ${wl}-retain-symbols-file $wl$export_symbols -o $lib'
      else
	_LT_AC_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;
    esac

    if test "$_LT_AC_TAGVAR(ld_shlibs, $1)" = no; then
      runpath_var=
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)=
      _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)=
      _LT_AC_TAGVAR(whole_archive_flag_spec, $1)=
    fi
  else
    # PORTME fill in a description of your system's linker (not GNU ld)
    case $host_os in
    aix3*)
      _LT_AC_TAGVAR(allow_undefined_flag, $1)=unsupported
      _LT_AC_TAGVAR(always_export_symbols, $1)=yes
      _LT_AC_TAGVAR(archive_expsym_cmds, $1)='$LD -o $output_objdir/$soname $libobjs $deplibs $linker_flags -bE:$export_symbols -T512 -H512 -bM:SRE~$AR $AR_FLAGS $lib $output_objdir/$soname'
      # Note: this linker hardcodes the directories in LIBPATH if there
      # are no directories specified by -L.
      _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes
      if test "$GCC" = yes && test -z "$lt_prog_compiler_static"; then
	# Neither direct hardcoding nor static linking is supported with a
	# broken collect2.
	_LT_AC_TAGVAR(hardcode_direct, $1)=unsupported
      fi
      ;;

    aix4* | aix5*)
      if test "$host_cpu" = ia64; then
	# On IA64, the linker does run time linking by default, so we don't
	# have to do anything special.
	aix_use_runtimelinking=no
	exp_sym_flag='-Bexport'
	no_entry_flag=""
      else
	# If we're using GNU nm, then we don't want the "-C" option.
	# -C means demangle to AIX nm, but means don't demangle with GNU nm
	if $NM -V 2>&1 | grep 'GNU' > /dev/null; then
	  _LT_AC_TAGVAR(export_symbols_cmds, $1)='$NM -Bpg $libobjs $convenience | awk '\''{ if (((\[$]2 == "T") || (\[$]2 == "D") || (\[$]2 == "B")) && ([substr](\[$]3,1,1) != ".")) { print \[$]3 } }'\'' | sort -u > $export_symbols'
	else
	  _LT_AC_TAGVAR(export_symbols_cmds, $1)='$NM -BCpg $libobjs $convenience | awk '\''{ if (((\[$]2 == "T") || (\[$]2 == "D") || (\[$]2 == "B")) && ([substr](\[$]3,1,1) != ".")) { print \[$]3 } }'\'' | sort -u > $export_symbols'
	fi
	aix_use_runtimelinking=no

	# Test if we are trying to use run time linking or normal
	# AIX style linking. If -brtl is somewhere in LDFLAGS, we
	# need to do runtime linking.
	case $host_os in aix4.[[23]]|aix4.[[23]].*|aix5*)
	  for ld_flag in $LDFLAGS; do
  	  if (test $ld_flag = "-brtl" || test $ld_flag = "-Wl,-brtl"); then
  	    aix_use_runtimelinking=yes
  	    break
  	  fi
	  done
	  ;;
	esac

	exp_sym_flag='-bexport'
	no_entry_flag='-bnoentry'
      fi

      # When large executables or shared objects are built, AIX ld can
      # have problems creating the table of contents.  If linking a library
      # or program results in "error TOC overflow" add -mminimal-toc to
      # CXXFLAGS/CFLAGS for g++/gcc.  In the cases where that is not
      # enough to fix the problem, add -Wl,-bbigtoc to LDFLAGS.

      _LT_AC_TAGVAR(archive_cmds, $1)=''
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=':'
      _LT_AC_TAGVAR(link_all_deplibs, $1)=yes

      if test "$GCC" = yes; then
	case $host_os in aix4.[[012]]|aix4.[[012]].*)
	# We only want to do this on AIX 4.2 and lower, the check
	# below for broken collect2 doesn't work under 4.3+
	  collect2name=`${CC} -print-prog-name=collect2`
	  if test -f "$collect2name" && \
  	   strings "$collect2name" | grep resolve_lib_name >/dev/null
	  then
  	  # We have reworked collect2
  	  _LT_AC_TAGVAR(hardcode_direct, $1)=yes
	  else
  	  # We have old collect2
  	  _LT_AC_TAGVAR(hardcode_direct, $1)=unsupported
  	  # It fails to find uninstalled libraries when the uninstalled
  	  # path is not listed in the libpath.  Setting hardcode_minus_L
  	  # to unsupported forces relinking
  	  _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes
  	  _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
  	  _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=
	  fi
	  ;;
	esac
	shared_flag='-shared'
	if test "$aix_use_runtimelinking" = yes; then
	  shared_flag="$shared_flag "'${wl}-G'
	fi
      else
	# not using gcc
	if test "$host_cpu" = ia64; then
  	# VisualAge C++, Version 5.5 for AIX 5L for IA-64, Beta 3 Release
  	# chokes on -Wl,-G. The following line is correct:
	  shared_flag='-G'
	else
	  if test "$aix_use_runtimelinking" = yes; then
	    shared_flag='${wl}-G'
	  else
	    shared_flag='${wl}-bM:SRE'
	  fi
	fi
      fi

      # It seems that -bexpall does not export symbols beginning with
      # underscore (_), so it is better to generate a list of symbols to export.
      _LT_AC_TAGVAR(always_export_symbols, $1)=yes
      if test "$aix_use_runtimelinking" = yes; then
	# Warning - without using the other runtime loading flags (-brtl),
	# -berok will link without error, but may produce a broken library.
	_LT_AC_TAGVAR(allow_undefined_flag, $1)='-berok'
       # Determine the default libpath from the value encoded in an empty executable.
       _LT_AC_SYS_LIBPATH_AIX
       _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-blibpath:$libdir:'"$aix_libpath"
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)="\$CC"' -o $output_objdir/$soname $libobjs $deplibs '"\${wl}$no_entry_flag"' $compiler_flags `if test "x${allow_undefined_flag}" != "x"; then echo "${wl}${allow_undefined_flag}"; else :; fi` '"\${wl}$exp_sym_flag:\$export_symbols $shared_flag"
       else
	if test "$host_cpu" = ia64; then
	  _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-R $libdir:/usr/lib:/lib'
	  _LT_AC_TAGVAR(allow_undefined_flag, $1)="-z nodefs"
	  _LT_AC_TAGVAR(archive_expsym_cmds, $1)="\$CC $shared_flag"' -o $output_objdir/$soname $libobjs $deplibs '"\${wl}$no_entry_flag"' $compiler_flags ${wl}${allow_undefined_flag} '"\${wl}$exp_sym_flag:\$export_symbols"
	else
	 # Determine the default libpath from the value encoded in an empty executable.
	 _LT_AC_SYS_LIBPATH_AIX
	 _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-blibpath:$libdir:'"$aix_libpath"
	  # Warning - without using the other run time loading flags,
	  # -berok will link without error, but may produce a broken library.
	  _LT_AC_TAGVAR(no_undefined_flag, $1)=' ${wl}-bernotok'
	  _LT_AC_TAGVAR(allow_undefined_flag, $1)=' ${wl}-berok'
	  # Exported symbols can be pulled into shared objects from archives
	  _LT_AC_TAGVAR(whole_archive_flag_spec, $1)='$convenience'
	  _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=yes
	  # This is similar to how AIX traditionally builds its shared libraries.
	  _LT_AC_TAGVAR(archive_expsym_cmds, $1)="\$CC $shared_flag"' -o $output_objdir/$soname $libobjs $deplibs ${wl}-bnoentry $compiler_flags ${wl}-bE:$export_symbols${allow_undefined_flag}~$AR $AR_FLAGS $output_objdir/$libname$release.a $output_objdir/$soname'
	fi
      fi
      ;;

    amigaos*)
      _LT_AC_TAGVAR(archive_cmds, $1)='$rm $output_objdir/a2ixlibrary.data~$echo "#define NAME $libname" > $output_objdir/a2ixlibrary.data~$echo "#define LIBRARY_ID 1" >> $output_objdir/a2ixlibrary.data~$echo "#define VERSION $major" >> $output_objdir/a2ixlibrary.data~$echo "#define REVISION $revision" >> $output_objdir/a2ixlibrary.data~$AR $AR_FLAGS $lib $libobjs~$RANLIB $lib~(cd $output_objdir && a2ixlibrary -32)'
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes
      # see comment about different semantics on the GNU ld section
      _LT_AC_TAGVAR(ld_shlibs, $1)=no
      ;;

    bsdi[[45]]*)
      _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)=-rdynamic
      ;;

    cygwin* | mingw* | pw32*)
      # When not using gcc, we currently assume that we are using
      # Microsoft Visual C++.
      # hardcode_libdir_flag_spec is actually meaningless, as there is
      # no search path for DLLs.
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)=' '
      _LT_AC_TAGVAR(allow_undefined_flag, $1)=unsupported
      # Tell ltmain to make .lib files, not .a files.
      libext=lib
      # Tell ltmain to make .dll files, not .so files.
      shrext_cmds=".dll"
      # FIXME: Setting linknames here is a bad hack.
      _LT_AC_TAGVAR(archive_cmds, $1)='$CC -o $lib $libobjs $compiler_flags `echo "$deplibs" | $SED -e '\''s/ -lc$//'\''` -link -dll~linknames='
      # The linker will automatically build a .lib file if we build a DLL.
      _LT_AC_TAGVAR(old_archive_From_new_cmds, $1)='true'
      # FIXME: Should let the user specify the lib program.
      _LT_AC_TAGVAR(old_archive_cmds, $1)='lib /OUT:$oldlib$oldobjs$old_deplibs'
      _LT_AC_TAGVAR(fix_srcfile_path, $1)='`cygpath -w "$srcfile"`'
      _LT_AC_TAGVAR(enable_shared_with_static_runtimes, $1)=yes
      ;;

    darwin* | rhapsody*)
      case $host_os in
        rhapsody* | darwin1.[[012]])
         _LT_AC_TAGVAR(allow_undefined_flag, $1)='${wl}-undefined ${wl}suppress'
         ;;
       *) # Darwin 1.3 on
         if test -z ${MACOSX_DEPLOYMENT_TARGET} ; then
           _LT_AC_TAGVAR(allow_undefined_flag, $1)='${wl}-flat_namespace ${wl}-undefined ${wl}suppress'
         else
           case ${MACOSX_DEPLOYMENT_TARGET} in
             10.[[012]])
               _LT_AC_TAGVAR(allow_undefined_flag, $1)='${wl}-flat_namespace ${wl}-undefined ${wl}suppress'
               ;;
             10.*)
               _LT_AC_TAGVAR(allow_undefined_flag, $1)='${wl}-undefined ${wl}dynamic_lookup'
               ;;
           esac
         fi
         ;;
      esac
      _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
      _LT_AC_TAGVAR(hardcode_direct, $1)=no
      _LT_AC_TAGVAR(hardcode_automatic, $1)=yes
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=unsupported
      _LT_AC_TAGVAR(whole_archive_flag_spec, $1)=''
      _LT_AC_TAGVAR(link_all_deplibs, $1)=yes
    if test "$GCC" = yes ; then
    	output_verbose_link_cmd='echo'
        _LT_AC_TAGVAR(archive_cmds, $1)='$CC -dynamiclib $allow_undefined_flag -o $lib $libobjs $deplibs $compiler_flags -install_name $rpath/$soname $verstring'
      _LT_AC_TAGVAR(module_cmds, $1)='$CC $allow_undefined_flag -o $lib -bundle $libobjs $deplibs$compiler_flags'
      # Don't fix this by using the ld -exported_symbols_list flag, it doesn't exist in older darwin lds
      _LT_AC_TAGVAR(archive_expsym_cmds, $1)='sed -e "s,#.*,," -e "s,^[    ]*,," -e "s,^\(..*\),_&," < $export_symbols > $output_objdir/${libname}-symbols.expsym~$CC -dynamiclib $allow_undefined_flag -o $lib $libobjs $deplibs $compiler_flags -install_name $rpath/$soname $verstring~nmedit -s $output_objdir/${libname}-symbols.expsym ${lib}'
      _LT_AC_TAGVAR(module_expsym_cmds, $1)='sed -e "s,#.*,," -e "s,^[    ]*,," -e "s,^\(..*\),_&," < $export_symbols > $output_objdir/${libname}-symbols.expsym~$CC $allow_undefined_flag  -o $lib -bundle $libobjs $deplibs$compiler_flags~nmedit -s $output_objdir/${libname}-symbols.expsym ${lib}'
    else
      case $cc_basename in
        xlc*)
         output_verbose_link_cmd='echo'
         _LT_AC_TAGVAR(archive_cmds, $1)='$CC -qmkshrobj $allow_undefined_flag -o $lib $libobjs $deplibs $compiler_flags ${wl}-install_name ${wl}`echo $rpath/$soname` $verstring'
         _LT_AC_TAGVAR(module_cmds, $1)='$CC $allow_undefined_flag -o $lib -bundle $libobjs $deplibs$compiler_flags'
          # Don't fix this by using the ld -exported_symbols_list flag, it doesn't exist in older darwin lds
         _LT_AC_TAGVAR(archive_expsym_cmds, $1)='sed -e "s,#.*,," -e "s,^[    ]*,," -e "s,^\(..*\),_&," < $export_symbols > $output_objdir/${libname}-symbols.expsym~$CC -qmkshrobj $allow_undefined_flag -o $lib $libobjs $deplibs $compiler_flags ${wl}-install_name ${wl}$rpath/$soname $verstring~nmedit -s $output_objdir/${libname}-symbols.expsym ${lib}'
          _LT_AC_TAGVAR(module_expsym_cmds, $1)='sed -e "s,#.*,," -e "s,^[    ]*,," -e "s,^\(..*\),_&," < $export_symbols > $output_objdir/${libname}-symbols.expsym~$CC $allow_undefined_flag  -o $lib -bundle $libobjs $deplibs$compiler_flags~nmedit -s $output_objdir/${libname}-symbols.expsym ${lib}'
          ;;
       *)
         _LT_AC_TAGVAR(ld_shlibs, $1)=no
          ;;
      esac
    fi
      ;;

    dgux*)
      _LT_AC_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    freebsd1*)
      _LT_AC_TAGVAR(ld_shlibs, $1)=no
      ;;

    # FreeBSD 2.2.[012] allows us to include c++rt0.o to get C++ constructor
    # support.  Future versions do this automatically, but an explicit c++rt0.o
    # does not break anything, and helps significantly (at the cost of a little
    # extra space).
    freebsd2.2*)
      _LT_AC_TAGVAR(archive_cmds, $1)='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags /usr/lib/c++rt0.o'
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    # Unfortunately, older versions of FreeBSD 2 do not have this feature.
    freebsd2*)
      _LT_AC_TAGVAR(archive_cmds, $1)='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags'
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    # FreeBSD 3 and greater uses gcc -shared to do shared libraries.
    freebsd* | kfreebsd*-gnu | dragonfly*)
      _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -o $lib $libobjs $deplibs $compiler_flags'
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    hpux9*)
      if test "$GCC" = yes; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$rm $output_objdir/$soname~$CC -shared -fPIC ${wl}+b ${wl}$install_libdir -o $output_objdir/$soname $libobjs $deplibs $compiler_flags~test $output_objdir/$soname = $lib || mv $output_objdir/$soname $lib'
      else
	_LT_AC_TAGVAR(archive_cmds, $1)='$rm $output_objdir/$soname~$LD -b +b $install_libdir -o $output_objdir/$soname $libobjs $deplibs $linker_flags~test $output_objdir/$soname = $lib || mv $output_objdir/$soname $lib'
      fi
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}+b ${wl}$libdir'
      _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes

      # hardcode_minus_L: Not really in the search PATH,
      # but as the default location of the library.
      _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes
      _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-E'
      ;;

    hpux10*)
      if test "$GCC" = yes -a "$with_gnu_ld" = no; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -fPIC ${wl}+h ${wl}$soname ${wl}+b ${wl}$install_libdir -o $lib $libobjs $deplibs $compiler_flags'
      else
	_LT_AC_TAGVAR(archive_cmds, $1)='$LD -b +h $soname +b $install_libdir -o $lib $libobjs $deplibs $linker_flags'
      fi
      if test "$with_gnu_ld" = no; then
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}+b ${wl}$libdir'
	_LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:

	_LT_AC_TAGVAR(hardcode_direct, $1)=yes
	_LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-E'

	# hardcode_minus_L: Not really in the search PATH,
	# but as the default location of the library.
	_LT_AC_TAGVAR(hardcode_minus_L, $1)=yes
      fi
      ;;

    hpux11*)
      if test "$GCC" = yes -a "$with_gnu_ld" = no; then
	case $host_cpu in
	hppa*64*)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared ${wl}+h ${wl}$soname -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	ia64*)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared ${wl}+h ${wl}$soname ${wl}+nodefaultrpath -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	*)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared -fPIC ${wl}+h ${wl}$soname ${wl}+b ${wl}$install_libdir -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	esac
      else
	case $host_cpu in
	hppa*64*)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -b ${wl}+h ${wl}$soname -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	ia64*)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -b ${wl}+h ${wl}$soname ${wl}+nodefaultrpath -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	*)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$CC -b ${wl}+h ${wl}$soname ${wl}+b ${wl}$install_libdir -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	esac
      fi
      if test "$with_gnu_ld" = no; then
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}+b ${wl}$libdir'
	_LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:

	case $host_cpu in
	hppa*64*|ia64*)
	  _LT_AC_TAGVAR(hardcode_libdir_flag_spec_ld, $1)='+b $libdir'
	  _LT_AC_TAGVAR(hardcode_direct, $1)=no
	  _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
	  ;;
	*)
	  _LT_AC_TAGVAR(hardcode_direct, $1)=yes
	  _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-E'

	  # hardcode_minus_L: Not really in the search PATH,
	  # but as the default location of the library.
	  _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes
	  ;;
	esac
      fi
      ;;

    irix5* | irix6* | nonstopux*)
      if test "$GCC" = yes; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname ${wl}$soname `test -n "$verstring" && echo ${wl}-set_version ${wl}$verstring` ${wl}-update_registry ${wl}${output_objdir}/so_locations -o $lib'
      else
	_LT_AC_TAGVAR(archive_cmds, $1)='$LD -shared $libobjs $deplibs $linker_flags -soname $soname `test -n "$verstring" && echo -set_version $verstring` -update_registry ${output_objdir}/so_locations -o $lib'
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec_ld, $1)='-rpath $libdir'
      fi
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath ${wl}$libdir'
      _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:
      _LT_AC_TAGVAR(link_all_deplibs, $1)=yes
      ;;

    netbsd*)
      if echo __ELF__ | $CC -E - | grep __ELF__ >/dev/null; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags'  # a.out
      else
	_LT_AC_TAGVAR(archive_cmds, $1)='$LD -shared -o $lib $libobjs $deplibs $linker_flags'      # ELF
      fi
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    newsos6)
      _LT_AC_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath ${wl}$libdir'
      _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    openbsd*)
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      if test -z "`echo __ELF__ | $CC -E - | grep __ELF__`" || test "$host_os-$host_cpu" = "openbsd2.8-powerpc"; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag -o $lib $libobjs $deplibs $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $pic_flag -o $lib $libobjs $deplibs $compiler_flags ${wl}-retain-symbols-file,$export_symbols'
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath,$libdir'
	_LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-E'
      else
       case $host_os in
	 openbsd[[01]].* | openbsd2.[[0-7]] | openbsd2.[[0-7]].*)
	   _LT_AC_TAGVAR(archive_cmds, $1)='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags'
	   _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
	   ;;
	 *)
	   _LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag -o $lib $libobjs $deplibs $compiler_flags'
	   _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath,$libdir'
	   ;;
       esac
      fi
      ;;

    os2*)
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes
      _LT_AC_TAGVAR(allow_undefined_flag, $1)=unsupported
      _LT_AC_TAGVAR(archive_cmds, $1)='$echo "LIBRARY $libname INITINSTANCE" > $output_objdir/$libname.def~$echo "DESCRIPTION \"$libname\"" >> $output_objdir/$libname.def~$echo DATA >> $output_objdir/$libname.def~$echo " SINGLE NONSHARED" >> $output_objdir/$libname.def~$echo EXPORTS >> $output_objdir/$libname.def~emxexp $libobjs >> $output_objdir/$libname.def~$CC -Zdll -Zcrtdll -o $lib $libobjs $deplibs $compiler_flags $output_objdir/$libname.def'
      _LT_AC_TAGVAR(old_archive_From_new_cmds, $1)='emximp -o $output_objdir/$libname.a $output_objdir/$libname.def'
      ;;

    osf3*)
      if test "$GCC" = yes; then
	_LT_AC_TAGVAR(allow_undefined_flag, $1)=' ${wl}-expect_unresolved ${wl}\*'
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared${allow_undefined_flag} $libobjs $deplibs $compiler_flags ${wl}-soname ${wl}$soname `test -n "$verstring" && echo ${wl}-set_version ${wl}$verstring` ${wl}-update_registry ${wl}${output_objdir}/so_locations -o $lib'
      else
	_LT_AC_TAGVAR(allow_undefined_flag, $1)=' -expect_unresolved \*'
	_LT_AC_TAGVAR(archive_cmds, $1)='$LD -shared${allow_undefined_flag} $libobjs $deplibs $linker_flags -soname $soname `test -n "$verstring" && echo -set_version $verstring` -update_registry ${output_objdir}/so_locations -o $lib'
      fi
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath ${wl}$libdir'
      _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:
      ;;

    osf4* | osf5*)	# as osf3* with the addition of -msym flag
      if test "$GCC" = yes; then
	_LT_AC_TAGVAR(allow_undefined_flag, $1)=' ${wl}-expect_unresolved ${wl}\*'
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared${allow_undefined_flag} $libobjs $deplibs $compiler_flags ${wl}-msym ${wl}-soname ${wl}$soname `test -n "$verstring" && echo ${wl}-set_version ${wl}$verstring` ${wl}-update_registry ${wl}${output_objdir}/so_locations -o $lib'
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='${wl}-rpath ${wl}$libdir'
      else
	_LT_AC_TAGVAR(allow_undefined_flag, $1)=' -expect_unresolved \*'
	_LT_AC_TAGVAR(archive_cmds, $1)='$LD -shared${allow_undefined_flag} $libobjs $deplibs $linker_flags -msym -soname $soname `test -n "$verstring" && echo -set_version $verstring` -update_registry ${output_objdir}/so_locations -o $lib'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='for i in `cat $export_symbols`; do printf "%s %s\\n" -exported_symbol "\$i" >> $lib.exp; done; echo "-hidden">> $lib.exp~
	$LD -shared${allow_undefined_flag} -input $lib.exp $linker_flags $libobjs $deplibs -soname $soname `test -n "$verstring" && echo -set_version $verstring` -update_registry ${output_objdir}/so_locations -o $lib~$rm $lib.exp'

	# Both c and cxx compiler support -rpath directly
	_LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-rpath $libdir'
      fi
      _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=:
      ;;

    solaris*)
      _LT_AC_TAGVAR(no_undefined_flag, $1)=' -z text'
      if test "$GCC" = yes; then
	wlarc='${wl}'
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared ${wl}-h ${wl}$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~$echo "local: *; };" >> $lib.exp~
	  $CC -shared ${wl}-M ${wl}$lib.exp ${wl}-h ${wl}$soname -o $lib $libobjs $deplibs $compiler_flags~$rm $lib.exp'
      else
	wlarc=''
	_LT_AC_TAGVAR(archive_cmds, $1)='$LD -G${allow_undefined_flag} -h $soname -o $lib $libobjs $deplibs $linker_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~$echo "local: *; };" >> $lib.exp~
  	$LD -G${allow_undefined_flag} -M $lib.exp -h $soname -o $lib $libobjs $deplibs $linker_flags~$rm $lib.exp'
      fi
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      case $host_os in
      solaris2.[[0-5]] | solaris2.[[0-5]].*) ;;
      *)
 	# The compiler driver will combine linker options so we
 	# cannot just pass the convience library names through
 	# without $wl, iff we do not link with $LD.
 	# Luckily, gcc supports the same syntax we need for Sun Studio.
 	# Supported since Solaris 2.6 (maybe 2.5.1?)
 	case $wlarc in
 	'')
 	  _LT_AC_TAGVAR(whole_archive_flag_spec, $1)='-z allextract$convenience -z defaultextract' ;;
 	*)
 	  _LT_AC_TAGVAR(whole_archive_flag_spec, $1)='${wl}-z ${wl}allextract`for conv in $convenience\"\"; do test -n \"$conv\" && new_convenience=\"$new_convenience,$conv\"; done; $echo \"$new_convenience\"` ${wl}-z ${wl}defaultextract' ;;
 	esac ;;
      esac
      _LT_AC_TAGVAR(link_all_deplibs, $1)=yes
      ;;

    sunos4*)
      if test "x$host_vendor" = xsequent; then
	# Use $CC to link under sequent, because it throws in some extra .o
	# files that make .init and .fini sections work.
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -G ${wl}-h $soname -o $lib $libobjs $deplibs $compiler_flags'
      else
	_LT_AC_TAGVAR(archive_cmds, $1)='$LD -assert pure-text -Bstatic -o $lib $libobjs $deplibs $linker_flags'
      fi
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_AC_TAGVAR(hardcode_direct, $1)=yes
      _LT_AC_TAGVAR(hardcode_minus_L, $1)=yes
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    sysv4)
      case $host_vendor in
	sni)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
	  _LT_AC_TAGVAR(hardcode_direct, $1)=yes # is this really true???
	;;
	siemens)
	  ## LD is ld it makes a PLAMLIB
	  ## CC just makes a GrossModule.
	  _LT_AC_TAGVAR(archive_cmds, $1)='$LD -G -o $lib $libobjs $deplibs $linker_flags'
	  _LT_AC_TAGVAR(reload_cmds, $1)='$CC -r -o $output$reload_objs'
	  _LT_AC_TAGVAR(hardcode_direct, $1)=no
        ;;
	motorola)
	  _LT_AC_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
	  _LT_AC_TAGVAR(hardcode_direct, $1)=no #Motorola manual says yes, but my tests say they lie
	;;
      esac
      runpath_var='LD_RUN_PATH'
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    sysv4.3*)
      _LT_AC_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='-Bexport'
      ;;

    sysv4*MP*)
      if test -d /usr/nec; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
	_LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
	runpath_var=LD_RUN_PATH
	hardcode_runpath_var=yes
	_LT_AC_TAGVAR(ld_shlibs, $1)=yes
      fi
      ;;

    sysv4*uw2* | sysv5OpenUNIX* | sysv5UnixWare7.[[01]].[[10]]* | unixware7*)
      _LT_AC_TAGVAR(no_undefined_flag, $1)='${wl}-z,text'
      _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      runpath_var='LD_RUN_PATH'

      if test "$GCC" = yes; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared ${wl}-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared ${wl}-Bexport:$export_symbols ${wl}-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
      else
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -G ${wl}-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -G ${wl}-Bexport:$export_symbols ${wl}-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
      fi
      ;;

    sysv5* | sco3.2v5* | sco5v6*)
      # Note: We can NOT use -z defs as we might desire, because we do not
      # link with -lc, and that would cause any symbols used from libc to
      # always be unresolved, which means just about no library would
      # ever link correctly.  If we're not using GNU ld we use -z text
      # though, which does catch some bad symbols but isn't as heavy-handed
      # as -z defs.
      _LT_AC_TAGVAR(no_undefined_flag, $1)='${wl}-z,text'
      _LT_AC_TAGVAR(allow_undefined_flag, $1)='${wl}-z,nodefs'
      _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='`test -z "$SCOABSPATH" && echo ${wl}-R,$libdir`'
      _LT_AC_TAGVAR(hardcode_libdir_separator, $1)=':'
      _LT_AC_TAGVAR(link_all_deplibs, $1)=yes
      _LT_AC_TAGVAR(export_dynamic_flag_spec, $1)='${wl}-Bexport'
      runpath_var='LD_RUN_PATH'

      if test "$GCC" = yes; then
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -shared ${wl}-h,\${SCOABSPATH:+${install_libdir}/}$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -shared ${wl}-Bexport:$export_symbols ${wl}-h,\${SCOABSPATH:+${install_libdir}/}$soname -o $lib $libobjs $deplibs $compiler_flags'
      else
	_LT_AC_TAGVAR(archive_cmds, $1)='$CC -G ${wl}-h,\${SCOABSPATH:+${install_libdir}/}$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_AC_TAGVAR(archive_expsym_cmds, $1)='$CC -G ${wl}-Bexport:$export_symbols ${wl}-h,\${SCOABSPATH:+${install_libdir}/}$soname -o $lib $libobjs $deplibs $compiler_flags'
      fi
      ;;

    uts4*)
      _LT_AC_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
      _LT_AC_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_AC_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    *)
      _LT_AC_TAGVAR(ld_shlibs, $1)=no
      ;;
    esac
  fi
])
AC_MSG_RESULT([$_LT_AC_TAGVAR(ld_shlibs, $1)])
test "$_LT_AC_TAGVAR(ld_shlibs, $1)" = no && can_build_shared=no

#
# Do we need to explicitly link libc?
#
case "x$_LT_AC_TAGVAR(archive_cmds_need_lc, $1)" in
x|xyes)
  # Assume -lc should be added
  _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=yes

  if test "$enable_shared" = yes && test "$GCC" = yes; then
    case $_LT_AC_TAGVAR(archive_cmds, $1) in
    *'~'*)
      # FIXME: we may have to deal with multi-command sequences.
      ;;
    '$CC '*)
      # Test whether the compiler implicitly links with -lc since on some
      # systems, -lgcc has to come before -lc. If gcc already passes -lc
      # to ld, don't add -lc before -lgcc.
      AC_MSG_CHECKING([whether -lc should be explicitly linked in])
      $rm conftest*
      printf "$lt_simple_compile_test_code" > conftest.$ac_ext

      if AC_TRY_EVAL(ac_compile) 2>conftest.err; then
        soname=conftest
        lib=conftest
        libobjs=conftest.$ac_objext
        deplibs=
        wl=$_LT_AC_TAGVAR(lt_prog_compiler_wl, $1)
	pic_flag=$_LT_AC_TAGVAR(lt_prog_compiler_pic, $1)
        compiler_flags=-v
        linker_flags=-v
        verstring=
        output_objdir=.
        libname=conftest
        lt_save_allow_undefined_flag=$_LT_AC_TAGVAR(allow_undefined_flag, $1)
        _LT_AC_TAGVAR(allow_undefined_flag, $1)=
        if AC_TRY_EVAL(_LT_AC_TAGVAR(archive_cmds, $1) 2\>\&1 \| grep \" -lc \" \>/dev/null 2\>\&1)
        then
	  _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=no
        else
	  _LT_AC_TAGVAR(archive_cmds_need_lc, $1)=yes
        fi
        _LT_AC_TAGVAR(allow_undefined_flag, $1)=$lt_save_allow_undefined_flag
      else
        cat conftest.err 1>&5
      fi
      $rm conftest*
      AC_MSG_RESULT([$_LT_AC_TAGVAR(archive_cmds_need_lc, $1)])
      ;;
    esac
  fi
  ;;
esac
])# AC_LIBTOOL_PROG_LD_SHLIBS


# _LT_AC_FILE_LTDLL_C
# -------------------
# Be careful that the start marker always follows a newline.
AC_DEFUN([_LT_AC_FILE_LTDLL_C], [
# /* ltdll.c starts here */
# #define WIN32_LEAN_AND_MEAN
# #include <windows.h>
# #undef WIN32_LEAN_AND_MEAN
# #include <stdio.h>
#
# #ifndef __CYGWIN__
# #  ifdef __CYGWIN32__
# #    define __CYGWIN__ __CYGWIN32__
# #  endif
# #endif
#
# #ifdef __cplusplus
# extern "C" {
# #endif
# BOOL APIENTRY DllMain (HINSTANCE hInst, DWORD reason, LPVOID reserved);
# #ifdef __cplusplus
# }
# #endif
#
# #ifdef __CYGWIN__
# #include <cygwin/cygwin_dll.h>
# DECLARE_CYGWIN_DLL( DllMain );
# #endif
# HINSTANCE __hDllInstance_base;
#
# BOOL APIENTRY
# DllMain (HINSTANCE hInst, DWORD reason, LPVOID reserved)
# {
#   __hDllInstance_base = hInst;
#   return TRUE;
# }
# /* ltdll.c ends here */
])# _LT_AC_FILE_LTDLL_C


# _LT_AC_TAGVAR(VARNAME, [TAGNAME])
# ---------------------------------
AC_DEFUN([_LT_AC_TAGVAR], [ifelse([$2], [], [$1], [$1_$2])])


# old names
AC_DEFUN([AM_PROG_LIBTOOL],   [AC_PROG_LIBTOOL])
AC_DEFUN([AM_ENABLE_SHARED],  [AC_ENABLE_SHARED($@)])
AC_DEFUN([AM_ENABLE_STATIC],  [AC_ENABLE_STATIC($@)])
AC_DEFUN([AM_DISABLE_SHARED], [AC_DISABLE_SHARED($@)])
AC_DEFUN([AM_DISABLE_STATIC], [AC_DISABLE_STATIC($@)])
AC_DEFUN([AM_PROG_LD],        [AC_PROG_LD])
AC_DEFUN([AM_PROG_NM],        [AC_PROG_NM])

# This is just to silence aclocal about the macro not being used
ifelse([AC_DISABLE_FAST_INSTALL])

AC_DEFUN([LT_AC_PROG_GCJ],
[AC_CHECK_TOOL(GCJ, gcj, no)
  test "x${GCJFLAGS+set}" = xset || GCJFLAGS="-g -O2"
  AC_SUBST(GCJFLAGS)
])

AC_DEFUN([LT_AC_PROG_RC],
[AC_CHECK_TOOL(RC, windres, no)
])

############################################################
# NOTE: This macro has been submitted for inclusion into   #
#  GNU Autoconf as AC_PROG_SED.  When it is available in   #
#  a released version of Autoconf we should remove this    #
#  macro and use it instead.                               #
############################################################
# LT_AC_PROG_SED
# --------------
# Check for a fully-functional sed program, that truncates
# as few characters as possible.  Prefer GNU sed if found.
AC_DEFUN([LT_AC_PROG_SED],
[AC_MSG_CHECKING([for a sed that does not truncate output])
AC_CACHE_VAL(lt_cv_path_SED,
[# Loop through the user's path and test for sed and gsed.
# Then use that list of sed's as ones to test for truncation.
as_save_IFS=$IFS; IFS=$PATH_SEPARATOR
for as_dir in $PATH
do
  IFS=$as_save_IFS
  test -z "$as_dir" && as_dir=.
  for lt_ac_prog in sed gsed; do
    for ac_exec_ext in '' $ac_executable_extensions; do
      if $as_executable_p "$as_dir/$lt_ac_prog$ac_exec_ext"; then
        lt_ac_sed_list="$lt_ac_sed_list $as_dir/$lt_ac_prog$ac_exec_ext"
      fi
    done
  done
done
lt_ac_max=0
lt_ac_count=0
# Add /usr/xpg4/bin/sed as it is typically found on Solaris
# along with /bin/sed that truncates output.
for lt_ac_sed in $lt_ac_sed_list /usr/xpg4/bin/sed; do
  test ! -f $lt_ac_sed && continue
  cat /dev/null > conftest.in
  lt_ac_count=0
  echo $ECHO_N "0123456789$ECHO_C" >conftest.in
  # Check for GNU sed and select it if it is found.
  if "$lt_ac_sed" --version 2>&1 < /dev/null | grep 'GNU' > /dev/null; then
    lt_cv_path_SED=$lt_ac_sed
    break
  fi
  while true; do
    cat conftest.in conftest.in >conftest.tmp
    mv conftest.tmp conftest.in
    cp conftest.in conftest.nl
    echo >>conftest.nl
    $lt_ac_sed -e 's/a$//' < conftest.nl >conftest.out || break
    cmp -s conftest.out conftest.nl || break
    # 10000 chars as input seems more than enough
    test $lt_ac_count -gt 10 && break
    lt_ac_count=`expr $lt_ac_count + 1`
    if test $lt_ac_count -gt $lt_ac_max; then
      lt_ac_max=$lt_ac_count
      lt_cv_path_SED=$lt_ac_sed
    fi
  done
done
])
SED=$lt_cv_path_SED
AC_MSG_RESULT([$SED])
])
# Copyright (C) 2006, 2009 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id: coin.m4 1313 2009-07-07 22:07:09Z andreasw $
#
# Author: Andreas Wachter    IBM      2006-04-14

# This file defines the common autoconf macros for COIN
#

# Check requirements
AC_PREREQ(2.59)

###########################################################################
#                           COIN_MAIN_SUBDIRS                             #
###########################################################################

# This is the macro for AC_COIN_MAIN_SUBDIRS taking care of ONE argument
AC_DEFUN([AC_COIN_MAIN_SUBDIR],
[AC_ARG_VAR([COIN_SKIP_PROJECTS],[Set to the subdirectories of projects that should be skipped in the configuration])

# This is a hack to find out if there is a "/" in the name, which would
# break the AC_COIN_HAS_PROJECT macro
m4_case($1,m4_bpatsubsts($1,[/],[_]),
[# We have no "/" in the $1 argument
  AC_COIN_HAS_PROJECT($1)
  AC_MSG_CHECKING(whether directory $1 should be recursed into)
  if test "$m4_tolower(coin_has_$1)" != skipped &&
     test "$m4_tolower(coin_has_$1)" != installed; then
    if test -r $srcdir/$1/configure; then
      coin_subdirs="$coin_subdirs $1"
      AC_MSG_RESULT(yes)
      AC_CONFIG_SUBDIRS($1)
    else
      AC_MSG_RESULT(no)
    fi
  else
    AC_MSG_RESULT(no)
  fi],
[# This must be Data/Simple or something else
  AC_MSG_CHECKING(whether directory $1 should be recursed into)
  coin_skip=no
  if test x"$COIN_SKIP_PROJECTS" != x; then
    for dir in $COIN_SKIP_PROJECTS; do
      if test $dir = $1; then
	coin_skip=yes
      fi
    done
  fi
  if test $coin_skip = yes; then
    AC_MSG_RESULT(skipping)
  elif test -r $srcdir/$1/configure; then
    coin_subdirs="$coin_subdirs $1"
    AC_MSG_RESULT(yes)
    AC_CONFIG_SUBDIRS($1)
  else
    AC_MSG_RESULT(no)
  fi])
])

# This macro sets up the recursion into configure scripts into
# subdirectories.  Each possible subdirectory should be provided as a
# new argument to this macro.  The current limit is 10 subdirectories.
# This automatically also checks for the Data subdirectory.

AC_DEFUN([AC_COIN_MAIN_SUBDIRS],
[m4_ifvaln([$1],[AC_COIN_MAIN_SUBDIR($1)])
m4_ifvaln([$2],[AC_COIN_MAIN_SUBDIR($2)])
m4_ifvaln([$3],[AC_COIN_MAIN_SUBDIR($3)])
m4_ifvaln([$4],[AC_COIN_MAIN_SUBDIR($4)])
m4_ifvaln([$5],[AC_COIN_MAIN_SUBDIR($5)])
m4_ifvaln([$6],[AC_COIN_MAIN_SUBDIR($6)])
m4_ifvaln([$7],[AC_COIN_MAIN_SUBDIR($7)])
m4_ifvaln([$8],[AC_COIN_MAIN_SUBDIR($8)])
m4_ifvaln([$9],[AC_COIN_MAIN_SUBDIR($9)])
m4_ifvaln([$10],[AC_COIN_MAIN_SUBDIR($10)])

]) # AC_COIN_MAIN_SUBDIRS

###########################################################################
#                            COIN_CHECK_FILE                              #
###########################################################################

# A simple replacement for AC_CHECK_FILE that works for cross compilation

AC_DEFUN([AC_COIN_CHECK_FILE],
[if test -r $1; then
  $2
  :
else
  $3
  :
fi
])

###########################################################################
#                        COIN_THIRDPARTY_SUBDIRS                          #
###########################################################################

# This macro sets up the recursion into the configure script in a
# subdirectory for compilation of third party code.  The first
# argument is just the string that appears in the configure output.
# The second argument is the directory with the configure script, and
# the third one is a file that should exists in that directory.  If
# this file does not exist, we assume that the user has not downloaded
# the code, and we are not going to compile it

AC_DEFUN([AC_COIN_THIRDPARTY_SUBDIRS],
[AC_MSG_CHECKING(whether code for third party package $1 is available)
coin_skip=no
if test x"$COIN_SKIP_PROJECTS" != x; then
  for dir in $COIN_SKIP_PROJECTS; do
    if test $dir = $2; then
      coin_skip=yes
    fi
  done
fi
if test $coin_skip = yes; then
  AC_MSG_RESULT(skipping)
else
  coin_tmp=`echo $srcdir/$2/$3`
  # There is probably a more elegant way to get the first thing out here...
  for i in $coin_tmp; do
    coin_tmp2=$i
    #break 1
  done
  if test -r $coin_tmp2; then
    coin_subdirs="$coin_subdirs $2"
    AC_MSG_RESULT(yes)
    AC_CONFIG_SUBDIRS($2)
  else
    AC_MSG_RESULT(no)
  fi
fi
]) # AC_COIN_THIRDPARTY_SUBDIRS

###########################################################################
#                           COIN_CHECK_VPATH                              #
###########################################################################

# This macro sets the variable coin_vpath_config to true if this is a
# VPATH configuration, otherwise it sets it to false.
AC_DEFUN([AC_COIN_CHECK_VPATH],
[AC_MSG_CHECKING(whether this is a VPATH configuration)
if test `cd $srcdir; pwd` != `pwd`; then
  coin_vpath_config=yes;
else
  coin_vpath_config=no;
fi
AC_MSG_RESULT($coin_vpath_config)
]) # AC_COIN_CHECK_VPATH

###########################################################################
#                         COIN_PROJECTDIR_INIT                            #
###########################################################################

# This macro does everything that is required in the early part in the
# configure script, such as defining a few variables.  This should only
# be used in the main directory of a project directory (the one under
# which src is)

AC_DEFUN([AC_COIN_PROJECTDIR_INIT],
[# Initialize the ADDLIBS variable
ADDLIBS="-lm $LIBS"
AC_SUBST(ADDLIBS)

# Initialize the FADDLIBS variable (which is to be used with a fortran
# compiler and will not include FLIBS)
FADDLIBS="$LIBS"
AC_SUBST(FADDLIBS)

# A useful makefile conditional that is always false
AM_CONDITIONAL(ALWAYS_FALSE, false)

# We set the following variable so that we know later in AC_COIN_FINALIZE
# that we are in a project main directory
coin_projectdir=yes
]) # AC_COIN_PROJECTDIR_INIT

###########################################################################
#                          COIN_DEBUG_COMPILE                             #
###########################################################################

# enable the configure flags --enable-debug and --enable-debug-prjct
# (where prcjt is the name of the project in lower case) and set the
# variable coin_debug_compile to true or false This is used by
# COIN_PROG_CXX, COIN_PROG_CC and COIN_PROG_F77 to determine the
# compilation flags.  This macro also makes the switches
# --with-prjct-verbosity and --with-prjct-checklevel available, which
# define the preprocessor macros COIN_PRJCT_VERBOSITY and
# COIN_PRJCT_CHECKLEVEL to the specified value (default is 0).
#
# The project specific flags are only made available, if one gives the
# name of the project as first argument to this macro.

AC_DEFUN([AC_COIN_DEBUG_COMPILE],
[AC_BEFORE([$0],[AC_COIN_PROG_CXX])dnl
AC_BEFORE([$0],[AC_COIN_PROG_CC])dnl
AC_BEFORE([$0],[AC_COIN_PROG_F77])dnl

AC_MSG_CHECKING([whether we want to compile in debug mode])

AC_ARG_ENABLE([debug],
[AC_HELP_STRING([--enable-debug],
                [compile all projects with debug options tests])],
[case "${enableval}" in
   yes) coin_debug_compile=true
        if test "${enable_shared+set}" = set; then :; else
          enable_shared=no
        fi 
        ;;
   no)  coin_debug_compile=false
        ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug)
        ;;
esac],
[coin_debug_compile=false])

m4_ifvaln([$1],
[AC_ARG_ENABLE(debug-m4_tolower($1),
 [AC_HELP_STRING([--enable-debug-m4_tolower($1)],
                 [compile this project ($1) with debug options])],
 [case "${enableval}" in
    yes) coin_debug_compile=true
         enable_shared=no
         ;;
    no)  coin_debug_compile=false
         ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug-m4_tolower($1))
         ;;
 esac],[:])
]) # m4_ifvaln([$1],

if test $coin_debug_compile = true; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

m4_ifvaln([$1],
[AC_ARG_WITH(m4_tolower($1)-verbosity,
             AC_HELP_STRING([--with-m4_tolower($1)-verbosity],
                            [specify the debug verbosity level for project $1]),
             [if test "$withval" = yes; then
                withval=1
              fi
              m4_tolower(coin_$1_verbosity)=$withval],
             [m4_tolower(coin_$1_verbosity)=0])
 AC_DEFINE_UNQUOTED(m4_toupper(COIN_$1_VERBOSITY),
                    m4_tolower($coin_$1_verbosity),
                    [Define to the debug verbosity level (0 is no output)])

 AC_ARG_WITH(m4_tolower($1)-checklevel,
             AC_HELP_STRING([--with-m4_tolower($1)-checklevel],
                            [specify the sanity check level for project $1]),
             [if test "$withval" = yes; then
                withval=1
              fi
              m4_tolower(coin_$1_checklevel)=$withval],
             [m4_tolower(coin_$1_checklevel)=0])
 AC_DEFINE_UNQUOTED(m4_toupper(COIN_$1_CHECKLEVEL),
                    m4_tolower($coin_$1_checklevel),
                    [Define to the debug sanity check level (0 is no test)])

# We use the following variable to have a string with the upper case
# version of the project name
COIN_PRJCT=m4_toupper($1)

]) # m4_ifvaln([$1],
 
]) # AC_COIN_DEBUG_COMPILE

###########################################################################
#                          COIN_MINGW_LD_FIX                              #
###########################################################################

# This macro is included by any PROG_compiler macro, to set the LD
# environment variable on MinGW to the correct value (link). But note that
# if we're building in cygwin with -mno-cygwin, we still want ld! If we're
# building from cygwin with MSVC tools (cl/link), then we do want link and
# you'd better have your PATH variable straight, else you'll be doing file
# links instead of code links! Arguably, LDFLAGS should include -mno-cygwin
# but in practice all linking seems to be handled through the compilers, so
# CFLAGS and CXXFLAGS suffice.

AC_DEFUN([AC_COIN_MINGW_LD_FIX],
[case $build in
  *-mingw*)
    if test "${LD+set}" = set; then :; else
      LD=link
    fi
    ;;
 esac
 case $enable_doscompile in
   msvc)
     if test "x${LD+set}" = xset; then :; else
       LD=link
     fi
     ;;
 esac
])

###########################################################################
#                        COIN_ENABLE_DOSCOMPILE                           #
###########################################################################

# This macro is invoked by any PROG_compiler macro to establish the
# --enable-doscompile option, used when one wants to compile an executable
# under Cygwin which also runs directly under DOS (without requiring
# Cygwin1.dll). The job of this macro is to make sure the option is correct and
# to set enable_doscompile. Legal values are mingw, msvc, and no (disabled).
# mingw: Fake mingw under cygwin, using GCC tools and -mno-cygwin. The most
#	 important thing here is that preprocess, compile, and link steps
#	 *all* see -mno-cygwin.
# msvc:  Assume the presence of cl/link. It's the user's responsibility to
#        make sure their PATH is correct. In particular, that MSVC link is
#	 found and not cygwin link (we want to do code linking, not file
#	 linking).
# It's the responsibility of individual PROG_compiler macros to ensure that
# they correctly set the compiler search list and preprocess, compile, and
# link flags. This is tied to compiler setup because in practice invocations
# of the preprocessor and linker are made through the compiler.

AC_DEFUN([AC_COIN_ENABLE_DOSCOMPILE],
[AC_REQUIRE([AC_CANONICAL_BUILD])
 AC_ARG_ENABLE([doscompile],
  [AC_HELP_STRING([--enable-doscompile],
                  [Under Cygwin, compile so that executables run under DOS.
		   Set to mingw to use gcc/g++/ld with -mno-cygwin.
		   Set to msvc to use cl/link (or icl/link).
		   Default when mentioned: mingw.
		   Default when not mentioned: disabled.])],
  [if test "$enable_doscompile" != no; then
     case $build in
       *-cygwin* | *-mingw*) ;;
       *) AC_MSG_ERROR([--enable-doscompile option makes sense only under Cygwin or MinGW]) ;;
     esac
   fi],
  [enable_doscompile=no])
 case "$enable_doscompile" in
   mingw)
     case $build in
       *-mingw*) enable_doscompile=no ;;
     esac
     ;;
   msvc|no) ;;
   yes) enable_doscompile=mingw ;;
   *) AC_MSG_ERROR([Invalid value $enable_doscompile for --enable-doscompile.
		    Try configure --help=recursive.])
      ;;
  esac
  if test "$enable_doscompile" != no ; then
    AC_MSG_NOTICE([DOS compile style is: $enable_doscompile])
  fi
])

###########################################################################
#                             COIN_PROG_CXX                               #
###########################################################################

# Find the compile command by running AC_PROG_CXX (with compiler names for
# different operating systems) and put it into CXX (unless it was given by the
# user). Then find an appropriate value for CXXFLAGS. If either of CXXFLAGS or
# PRJCT_CXXFLAGS is defined, that value is used (replace PRJCT with the upper
# case name of this project).  It is possible to provide additional -D flags
# in the variable CXXDEFS, and additional compilation flags with ADD_CXXFLAGS.

AC_DEFUN([AC_COIN_PROG_CXX],
[AC_REQUIRE([AC_COIN_PROG_CC]) #Let's try if that overcomes configuration problem with VC++ 6.0
AC_REQUIRE([AC_COIN_ENABLE_DOSCOMPILE])
AC_LANG_PUSH(C++)

AC_ARG_VAR(CXXDEFS,[Additional -D flags to be used when compiling C++ code.])
AC_ARG_VAR(ADD_CXXFLAGS,[Additional C++ compiler options])
AC_ARG_VAR(DBG_CXXFLAGS,[Debug C++ compiler options])
AC_ARG_VAR(OPT_CXXFLAGS,[Optimize C++ compiler options])

coin_has_cxx=yes

save_cxxflags="$CXXFLAGS"
# For *-*-solaris*, promote Studio/Workshop compiler to front of list.
case $build in
  *-cygwin* | *-mingw*)
  	     if test "$enable_doscompile" = msvc ; then
	       comps="icl cl"
	     else
	       comps="g++ cl"
	     fi ;;
  *-*-solaris*)
  	     comps="CC xlC_r aCC g++ c++ pgCC icpc gpp cxx cc++ cl FCC KCC RCC" ;;
  *-darwin*) comps="g++ c++ CC" ;;
          *) comps="xlC_r aCC CC g++ c++ pgCC icpc gpp cxx cc++ cl FCC KCC RCC" ;;
esac

# We delete the cached value, since the test might not have been
# performed with our choice of compilers earlier
$as_unset ac_cv_prog_CXX || test "${ac_cv_prog_CXX+set}" != set || { ac_cv_prog_CXX=; export ac_cv_prog_CXX; }
# AC_MSG_NOTICE([C++ compiler candidates: $comps])
AC_PROG_CXX([$comps])

#AC_PROG_CXX sets CXX to g++ if it cannot find a working C++ compiler
#thus, we test here whether $CXX is actually working 
AC_LANG_PUSH(C++)
AC_MSG_CHECKING([whether C++ compiler $CXX works]);
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM(, [int i=0;])],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_RESULT(no)
   AC_MSG_ERROR(failed to find a C++ compiler or C++ compiler $CXX does not work)]
)
AC_LANG_POP(C++)

# It seems that we need to cleanup something here for the Windows 
case "$CXX" in
  cl* | */cl* | CL* | */CL* | icl* | */icl* | ICL* | */ICL*)
    sed -e 's/^void exit (int);//' confdefs.h >> confdefs.hh
    mv confdefs.hh confdefs.h
    ;;
esac

# Autoconf incorrectly concludes that cl recognises -g. It doesn't.
case "$CXX" in
  cl* | */cl* | CL* | */CL* )
    if test "$ac_cv_prog_cxx_g" = yes ; then
      ac_cv_prog_cxx_g=no
      AC_MSG_NOTICE([Overruling autoconf; cl does not recognise -g.])
    fi ;;
esac
CXXFLAGS="$save_cxxflags"

# Check if a project specific CXXFLAGS variable has been set
if test x$COIN_PRJCT != x; then
  eval coin_tmp=\${${COIN_PRJCT}_CXXFLAGS+set}
  if test x$coin_tmp = xset; then
    eval CXXFLAGS=\${${COIN_PRJCT}_CXXFLAGS}
  fi
fi

if test x"$CXXFLAGS" = x; then

# ToDo decide whether we want -DNDEBUG for optimization
  coin_add_cxxflags=
  coin_opt_cxxflags=
  coin_dbg_cxxflags=
  coin_warn_cxxflags=

  if test "$GXX" = "yes"; then
    case "$CXX" in
      icpc* | */icpc*)
        ;;
      *)
# ToDo decide about unroll-loops
        coin_opt_cxxflags="-O3 -fomit-frame-pointer"
        coin_add_cxxflags="-pipe"
        coin_dbg_cxxflags="-g"
        coin_warn_cxxflags="-Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas"
        case $build in
          *-darwin*)
            ;;
          *)
            coin_warn_cxxflags="-pedantic-errors $coin_warn_cxxflags"
            ;;
        esac

        case $enable_doscompile in
          mingw)
            CXXFLAGS="-mno-cygwin"
            AC_TRY_LINK(,[int i=0; i++;],[coin_add_cxxflags="-mno-cygwin $coin_add_cxxflags"])
            CXXFLAGS=
          ;;
        esac
    esac
  fi

# Note that we do not need to cover GCC in the following tests.

  if test -z "$coin_opt_cxxflags"; then
    case $build in
      *-cygwin* | *-mingw*)
        case "$CXX" in
          cl* | */cl* | CL* | */CL*)
	    # The MT and MTd options are mutually exclusive
            coin_opt_cxxflags='-MT -O2'
            coin_add_cxxflags='-nologo -EHsc -GR -wd4996 -D_CRT_SECURE_NO_DEPRECATE'
            coin_dbg_cxxflags='-MTd'
            ;;
          icl* | */icl* | ICL* | */ICL*)
	    # The MT and MTd options are mutually exclusive
            coin_opt_cxxflags='-MT -Ox'
            coin_add_cxxflags='-nologo -EHsc -GR -D_CRT_SECURE_NO_DEPRECATE'
            coin_dbg_cxxflags='-MTd -debug'
            ;;
        esac
        ;;
      *-linux-*)
        case "$CXX" in
          icpc* | */icpc*)
            coin_opt_cxxflags="-O3 -ip -mp1"
            coin_add_cxxflags=""
            coin_dbg_cxxflags="-g"
            # Check if -i_dynamic is necessary (for new glibc library)
            CXXFLAGS=
            AC_TRY_LINK(,[int i=0; i++;],[],
                        [coin_add_cxxflags="-i_dynamic $coin_add_cxxflags"])
            ;;
          pgCC* | */pgCC*)
            coin_opt_cxxflags="-fast"
            coin_add_cxxflags="-Kieee -pc 64"
            coin_dbg_cxxflags="-g"
            ;;
        esac
        ;;
      *-ibm-*)
        case "$CXX" in
          xlC* | */xlC* | mpxlC* | */mpxlC*)
            coin_opt_cxxflags="-O -qarch=auto -qcache=auto -qtune=auto -qmaxmem=-1"
            coin_add_cxxflags="-bmaxdata:0x80000000 -qrtti=dyna -qsuppress=1500-036 -qsuppress=1500-029 -qsourcetype=c++"
            coin_dbg_cxxflags="-g"
            ;;
        esac
        ;;
      *-hp-*)
        case "$CXX" in
          aCC* | */aCC* )
            coin_opt_cxxflags="-O"
            coin_add_cxxflags="-AA"
            coin_dbg_cxxflags="-g"
            ;;
        esac
        ;;
      *-*-solaris*)
          coin_opt_cxxflags="-O4"
          coin_dbg_cxxflags="-g"
        ;;
    esac
  fi

# Generic flag settings. If these don't work, add a case above.

  if test "$ac_cv_prog_cxx_g" = yes && test -z "$coin_dbg_cxxflags" ; then
    coin_dbg_cxxflags="-g"
  fi

  if test -z "$coin_opt_cxxflags"; then
    # Try if -O option works if nothing else is set
    CXXFLAGS=-O
    AC_TRY_LINK([],[int i=0; i++;],[coin_opt_cxxflags="-O"])
  fi

  # if PM doesn't want the warning messages, take them out
  if test x"$coin_skip_warn_cxxflags" = xyes; then
    coin_warn_cxxflags=
  fi

# Do final setup of flags based on values determined above.

  if test x${DBG_CXXFLAGS+set} != xset; then
    DBG_CXXFLAGS="$coin_dbg_cxxflags $coin_add_cxxflags $coin_warn_cxxflags"
  fi
  if test x${OPT_CXXFLAGS+set} != xset; then
    OPT_CXXFLAGS="$coin_opt_cxxflags $coin_add_cxxflags -DNDEBUG $coin_warn_cxxflags"
  fi

  DBG_CXXFLAGS="$DBG_CXXFLAGS $ADD_CXXFLAGS $CXXDEFS"
  OPT_CXXFLAGS="$OPT_CXXFLAGS $ADD_CXXFLAGS $CXXDEFS"

  if test "$coin_debug_compile" = "true"; then
    CXXFLAGS="$DBG_CXXFLAGS"
  else
    CXXFLAGS="$OPT_CXXFLAGS"
  fi

# Handle the case where CXXFLAGS was set externally.
else
  CXXFLAGS="$CXXFLAGS $ADD_CXXFLAGS $CXXDEFS"
  if test x${DBG_CXXFLAGS+set} != xset; then
    DBG_CXXFLAGS="$CXXFLAGS"
  fi
  if test x${OPT_CXXFLAGS+set} != xset; then
    OPT_CXXFLAGS="$CXXFLAGS"
  fi
fi

# If CXXFLAGS contains -mno-cygwin, CPPFLAGS must also have it.

case "$CXXFLAGS" in
  *-mno-cygwin*)
    if test x${CPPFLAGS+set} != xset ; then
      CPPFLAGS="-mno-cygwin"
    else
      case "$CPPFLAGS" in
        *-mno-cygwin*)
	  ;;
	*)
	  CPPFLAGS="$CPPFLAGS -mno-cygwin"
	  ;;
      esac
    fi ;;
esac

# Try if CXXFLAGS works
save_CXXFLAGS="$CXXFLAGS"
AC_TRY_LINK([],[int i=0; i++;],[],[CXXFLAGS=])
if test -z "$CXXFLAGS"; then
  AC_MSG_WARN([The flags CXXFLAGS="$save_CXXFLAGS" do not work.  I will now just try '-O', but you might want to set CXXFLAGS manually.])
  CXXFLAGS='-O'
  AC_TRY_LINK([],[int i=0; i++;],[],[CXXFLAGS=])
  if test -z "$CXXFLAGS"; then
    AC_MSG_WARN([This value for CXXFLAGS does not work.  I will continue with empty CXXFLAGS, but you might want to set CXXFLAGS manually.])
  fi
fi

AC_MSG_NOTICE([C++ compiler options are: $CXXFLAGS])

AC_ARG_VAR(MPICXX,[C++ MPI Compiler])
if test x"$MPICXX" = x; then :; else
  AC_MSG_NOTICE([Will use MPI C++ compiler $MPICXX])
  CXX="$MPICXX"
fi

case "$CXX" in
  cl* | */cl* | CL* | */CL* )
    AC_COIN_MINGW_LD_FIX
    ;;
esac

AC_LANG_POP(C++)
]) # AC_COIN_PROG_CXX


###########################################################################
#                             COIN_CXXLIBS                                #
###########################################################################

# Determine the C++ runtime libraries required for linking a C++ library
# with a Fortran or C compiler.  The result is available in CXXLIBS.

AC_DEFUN([AC_COIN_CXXLIBS],
[AC_REQUIRE([AC_PROG_CXX])dnl
AC_LANG_PUSH(C++)
AC_ARG_VAR(CXXLIBS,[Libraries necessary for linking C++ code with Fortran compiler])
if test -z "$CXXLIBS"; then
  if test "$GXX" = "yes"; then
    case "$CXX" in
      icpc* | */icpc*)
        CXXLIBS="-lstdc++"
        ;;
      *)
        CXXLIBS="-lstdc++ -lm" # -lgcc"
        ;;
    esac
  else
    case $build in
     *-mingw32 | *-cygwin* )
      case "$CXX" in
      cl* | */cl* | CL* | */CL*)
        CXXLIBS=nothing;;
      esac;;
     *-linux-*)
      case "$CXX" in
      icpc* | */icpc*)
        CXXLIBS="-lstdc++"
             ;;
      pgCC* | */pgCC*)
        CXXLIBS="-lstd -lC -lc"
             ;;
      esac;;
    *-ibm-*)
      CXXLIBS="-lC -lc"
      ;;
    *-hp-*)
      CXXLIBS="-L/opt/aCC/lib -l++ -lstd_v2 -lCsup_v2 -lm -lcl -lc"
      ;;
    *-*-solaris*)
      CXXLIBS="-lCstd -lCrun"
    esac
  fi
fi
if test -z "$CXXLIBS"; then
  AC_MSG_WARN([Could not automatically determine CXXLIBS (C++ link libraries; necessary if main program is in Fortran or C).])
else
  AC_MSG_NOTICE([Assuming that CXXLIBS is \"$CXXLIBS\".])
fi
if test x"$CXXLIBS" = xnothing; then
  CXXLIBS=
fi
AC_LANG_POP(C++)
]) # AC_COIN_CXXLIBS

###########################################################################
#                           COIN_CHECK_HEADER                             #
###########################################################################

# This macro checks for a header file, but it does so without the
# standard header.  This avoids warning messages like:
#
# configure: WARNING: dlfcn.h: present but cannot be compiled
# configure: WARNING: dlfcn.h:     check for missing prerequisite headers?
# configure: WARNING: dlfcn.h: see the Autoconf documentation
# configure: WARNING: dlfcn.h:     section "Present But Cannot Be Compiled"
# configure: WARNING: dlfcn.h: proceeding with the preprocessor's result
# configure: WARNING: dlfcn.h: in the future, the compiler will take precedence

# My guess that the problem lay with CPPFLAGS seems to be correct. With this
# macro reduced to a direct call to AC_CHECK_HEADERS, there are no warnings
# now that CPPFLAGS contains -mno-cygwin when it needs it. -- lh, 061214 --

# AW 070609: I restored the previous version, since otherwise the warnings
# still pop up for the cl compiler

AC_DEFUN([AC_COIN_CHECK_HEADER],
[#if test x"$4" = x; then
#  hdr="#include <$1>"
#else
#  hdr="$4"
#fi
#AC_CHECK_HEADERS([$1],[$2],[$3],[$hdr])
AC_CHECK_HEADERS([$1],[$2],[$3],[$4])
]) # AC_COIN_CHECK_HEADER

###########################################################################
#                       COIN_CHECK_CXX_CHEADER                             #
###########################################################################

# This macro checks for C header files that are used from C++.  For a give
# stub (e.g., math), it first checks if the C++ library (cmath) is available.
# If it is, it defines HAVE_CMATH (or whatever the stub is).  If it is not
# available, it checks for the old C head (math.h) and defines HAVE_MATH_H
# if that one exists.

AC_DEFUN([AC_COIN_CHECK_CXX_CHEADER],
[AC_LANG_PUSH(C++)
AC_COIN_CHECK_HEADER([c$1],[$2],[$3],[$4])
if test "$ac_cv_header_c$1" != "yes"; then
  AC_COIN_CHECK_HEADER([$1.h],[$2],[$3],[$4])
fi
AC_LANG_POP(C++)
]) # AC_COIN_CHECK_CXX_CHEADER

###########################################################################
#                             COIN_PROG_CC                                #
###########################################################################

# Find the compile command by running AC_PROG_CC (with compiler names
# for different operating systems) and put it into CC (unless it was
# given my the user), and find an appropriate value for CFLAGS.  It is
# possible to provide additional -D flags in the variable CDEFS.

AC_DEFUN([AC_COIN_PROG_CC],
[AC_REQUIRE([AC_COIN_MINGW_LD_FIX])
AC_REQUIRE([AC_COIN_ENABLE_DOSCOMPILE])
AC_LANG_PUSH(C)

# For consistency, we set the C compiler to the same value of the C++
# compiler, if the C++ is set, but the C compiler isn't (only for CXX=cl)
if test x"$CXX" != x; then
  case "$CXX" in
    cl* | */cl* | CL* | */CL* | icl* | */icl* | ICL* | */ICL*)
      if test x"$CC" = x; then
        CC="$CXX"
        AC_MSG_WARN([C++ compiler name provided as $CXX, but CC is unset. Setting CC to $CXX])
      fi
      ;;
  esac
fi

AC_ARG_VAR(CDEFS,[Additional -D flags to be used when compiling C code.])
AC_ARG_VAR(ADD_CFLAGS,[Additional C compiler options])
AC_ARG_VAR(DBG_CFLAGS,[Debug C compiler options])
AC_ARG_VAR(OPT_CFLAGS,[Optimize C compiler options])

coin_has_cc=yes

save_cflags="$CFLAGS"

# For *-*-solaris*, promote Studio/Workshop cc compiler to front of list.
# Depending on the user's PATH, when Studio/Workshop cc is not present we may
# find /usr/ucb/cc, which is almost certainly not a good choice for the C
# compiler. In this case, put cc after gcc.

case $build in
  *-cygwin* | *-mingw*)
  	     if test "$enable_doscompile" = msvc ; then
	       comps="icl cl"
	     else
	       comps="gcc cl"
	     fi ;;
  *-*-solaris*)
	     AC_CHECK_PROG(sol_cc_compiler,cc,cc,[],[],/usr/ucb/cc)
	     if test "$sol_cc_compiler" = "cc" ; then
	       comps="cc xlc gcc pgcc icc"
	     else
	       comps="xlc gcc pgcc icc cc"
	     fi
	     ;;
  *-linux-*) comps="xlc gcc cc pgcc icc" ;;
  *)         comps="xlc_r xlc cc gcc pgcc icc" ;;
esac

# We delete the cached value, since the test might not have been
# performed with our choice of compilers earlier
$as_unset ac_cv_prog_CC || test "${ac_cv_prog_CC+set}" != set || { ac_cv_prog_CC=; export ac_cv_prog_CC; }
# AC_MSG_NOTICE([C compiler candidates: $comps])
AC_PROG_CC([$comps])
if test -z "$CC" ; then
  AC_MSG_ERROR([Failed to find a C compiler!])
fi
# Autoconf incorrectly concludes that cl recognises -g. It doesn't.
case "$CC" in
  cl* | */cl* | CL* | */CL* )
    if test "$ac_cv_prog_cc_g" = yes ; then
      ac_cv_prog_cc_g=no
      AC_MSG_NOTICE([Overruling autoconf; cl does not recognise -g.])
    fi ;;
esac
CFLAGS="$save_cflags"

# Check if a project specific CFLAGS variable has been set
if test x$COIN_PRJCT != x; then
  eval coin_tmp=\${${COIN_PRJCT}_CFLAGS+set}
  if test x$coin_tmp = xset; then
    eval CFLAGS=\${${COIN_PRJCT}_CFLAGS}
  fi
fi

if test x"$CFLAGS" = x; then

  coin_add_cflags=
  coin_opt_cflags=
  coin_dbg_cflags=
  coin_warn_cflags=

  if test "$GCC" = "yes"; then
    case "$CC" in
      icc* | */icc*)
        ;;
      *)
        coin_opt_cflags="-O3 -fomit-frame-pointer"
        coin_add_cflags="-pipe"
        coin_dbg_cflags="-g"
        coin_warn_cflags="-Wimplicit -Wparentheses -Wsequence-point -Wreturn-type -Wcast-qual -Wall -Wno-unknown-pragmas"
        case $build in
          *-darwin*)
            ;;
          *)
            coin_warn_cflags="-pedantic-errors $coin_warn_cflags"
            ;;
        esac
	case $enable_doscompile in
	  mingw)
	    CFLAGS="-mno-cygwin"
	    AC_TRY_LINK([],[int i=0; i++;],
			[coin_add_cflags="-mno-cygwin $coin_add_cflags"])
	    CFLAGS=
	  ;;
	esac
    esac
  fi
  if test -z "$coin_opt_cflags"; then
    case $build in
      *-cygwin* | *-mingw*)
        case "$CC" in
          cl* | */cl* | CL* | */CL*)
            coin_opt_cflags='-MT -O2'
            coin_add_cflags='-nologo -wd4996 -D_CRT_SECURE_NO_DEPRECATE'
            coin_dbg_cflags='-MTd'
            ;;
          icl* | */icl* | ICL* | */ICL*)
            coin_opt_cflags='-MT -Ox'
            coin_add_cflags='-nologo -D_CRT_SECURE_NO_DEPRECATE'
            coin_dbg_cflags='-MTd -debug'
            ;;
        esac
        ;;
      *-linux-*)
        case "$CC" in
          icc* | */icc*)
            coin_opt_cflags="-O3 -ip -mp1"
            coin_add_cflags=""
            coin_dbg_cflags="-g"
            # Check if -i_dynamic is necessary (for new glibc library)
            CFLAGS=
            AC_TRY_LINK([],[int i=0; i++;],[],
                        [coin_add_cflags="-i_dynamic $coin_add_cflags"])
            ;;
          pgcc* | */pgcc*)
            coin_opt_cflags="-fast"
            coin_add_cflags="-Kieee -pc 64"
            coin_dbg_cflags="-g"
            ;;
        esac
        ;;
      *-ibm-*)
        case "$CC" in
          xlc* | */xlc* | mpxlc* | */mpxlc*)
            coin_opt_cflags="-O -qarch=auto -qcache=auto -qtune=auto -qmaxmem=-1"
            coin_add_cflags="-bmaxdata:0x80000000 -qsuppress=1500-036 -qsuppress=1500-029"
            coin_dbg_cflags="-g"
          ;;
        esac
        ;;
      *-hp-*)
        coin_opt_cflags="-O"
        coin_add_cflags="-Ae"
        coin_dbg_cflags="-g"
        ;;
      *-*-solaris*)
        coin_opt_cflags="-xO4"
        coin_dbg_cflags="-g"
        ;;
      *-sgi-*)
        coin_opt_cflags="-O -OPT:Olimit=0"
        coin_dbg_cflags="-g"
        ;;
    esac
  fi

  if test "$ac_cv_prog_cc_g" = yes && test -z "$coin_dbg_cflags" ; then
    coin_dbg_cflags="-g"
  fi

  if test -z "$coin_opt_cflags"; then
    # Try if -O option works if nothing else is set
    CFLAGS="-O"
    AC_TRY_LINK([],[int i=0; i++;],[coin_opt_cflags="-O"])
  fi

  # if PM doesn't want the warning messages, take them out
  if test x"$coin_skip_warn_cflags" = xyes; then
    coin_warn_cflags=
  fi

  if test x${DBG_CFLAGS+set} != xset; then
    DBG_CFLAGS="$coin_dbg_cflags $coin_add_cflags $coin_warn_cflags"
  fi
  if test x${OPT_CFLAGS+set} != xset; then
    OPT_CFLAGS="$coin_opt_cflags $coin_add_cflags -DNDEBUG $coin_warn_cflags"
  fi

  DBG_CFLAGS="$DBG_CFLAGS $ADD_CFLAGS $CDEFS"
  OPT_CFLAGS="$OPT_CFLAGS $ADD_CFLAGS $CDEFS"

  if test "$coin_debug_compile" = "true"; then
    CFLAGS="$DBG_CFLAGS"
  else
    CFLAGS="$OPT_CFLAGS"
  fi
else
  CFLAGS="$CFLAGS $ADD_CFLAGS $CDEFS"
  if test x${DBG_CFLAGS+set} != xset; then
    DBG_CFLAGS="$CFLAGS"
  fi
  if test x${OPT_CFLAGS+set} != xset; then
    OPT_CFLAGS="$CFLAGS"
  fi
fi

# If CFLAGS contains -mno-cygwin, CPPFLAGS must also have it.

case "$CFLAGS" in
  *-mno-cygwin*)
    if test x${CPPFLAGS+set} != xset ; then
      CPPFLAGS="-mno-cygwin"
    else
      case "$CPPFLAGS" in
        *-mno-cygwin*)
	  ;;
	*)
	  CPPFLAGS="$CPPFLAGS -mno-cygwin"
	  ;;
      esac
    fi ;;
esac

# Try if CFLAGS works
save_CFLAGS="$CFLAGS"
AC_TRY_LINK([],[int i=0; i++;],[],[CFLAGS=])
if test -z "$CFLAGS"; then
  AC_MSG_WARN([The value CFLAGS="$save_CFLAGS" do not work.  I will now just try '-O', but you might want to set CFLAGS manually.])
  CFLAGS='-O'
  AC_TRY_LINK([],[int i=0; i++;],[],[CFLAGS=])
  if test -z "$CFLAGS"; then
    AC_MSG_WARN([This value for CFLAGS does not work.  I will continue with empty CFLAGS, but you might want to set CFLAGS manually.])
  fi
fi

AC_MSG_NOTICE([C compiler options are: $CFLAGS])

AC_ARG_VAR(MPICC,[C MPI Compiler])
if test x"$MPICC" = x; then :; else
  AC_MSG_NOTICE([Will use MPI C compiler $MPICC])
  CC="$MPICC"
fi

# Correct ADDLIBS initialization if we are using the MS compiler
case "$CC" in
  cl* | */cl* | CL* | */CL* | icl* | */icl* | ICL* | */ICL*)
    ADDLIBS=
    AC_COIN_MINGW_LD_FIX
    ;;
esac

AC_LANG_POP(C)
]) # AC_COIN_PROG_CC

###########################################################################
#                             COIN_PROG_F77                               #
###########################################################################

# Find the compile command by running AC_PROG_F77 (with compiler names
# for different operating systems) and put it into F77 (unless it was
# given by the user), and find an appropriate value for FFLAGS

AC_DEFUN([AC_COIN_PROG_F77],
[AC_REQUIRE([AC_COIN_MINGW_LD_FIX])
AC_REQUIRE([AC_COIN_ENABLE_DOSCOMPILE])
AC_REQUIRE([AC_COIN_PROG_CC])
AC_REQUIRE([AC_COIN_F77_COMPS])
AC_LANG_PUSH([Fortran 77])

AC_ARG_VAR(ADD_FFLAGS,[Additional Fortran compiler options])
AC_ARG_VAR(DBG_FFLAGS,[Debug Fortran compiler options])
AC_ARG_VAR(OPT_FFLAGS,[Optimize Fortran compiler options])

coin_has_f77=yes

save_fflags="$FFLAGS"

# We delete the cached value, since the test might not have been
# performed with our choice of compilers earlier
$as_unset ac_cv_prog_F77 || test "${ac_cv_prog_F77+set}" != set || { ac_cv_prog_F77=; export ac_cv_prog_F77; }

# This is a real belt-and-suspenders approach. AC_COIN_FIND_F77 will use
# coin_f77_comps to see if there's a program that matches one of the names.
# If there's no such program, F77 = unavailable. If we match the name,
# feed AC_PROG_F77 the same search list, just to be sure it's a functioning
# compiler.
# AC_MSG_NOTICE([Fortran compiler candidates: $coin_f77_comps])
AC_COIN_FIND_F77
if test "$F77" != "unavailable" ; then
  AC_PROG_F77($coin_f77_comps)
else
  AC_MSG_WARN([Failed to find a Fortran compiler!])
fi

FFLAGS="$save_fflags"

# Check if a project specific FFLAGS variable has been set
if test x$COIN_PRJCT != x; then
  eval coin_tmp=\${${COIN_PRJCT}_FFLAGS+set}
  if test x$coin_tmp = xset; then
    eval FFLAGS=\${${COIN_PRJCT}_FFLAGS}
  fi
fi

if test "$F77" != "unavailable" && test x"$FFLAGS" = x ; then

  coin_add_fflags=
  coin_opt_fflags=
  coin_dbg_fflags=
  coin_warn_fflags=

  if test "$G77" = "yes"; then
    coin_opt_fflags="-O3 -fomit-frame-pointer"
    coin_add_fflags="-pipe"
    coin_dbg_fflags="-g"
    case $enable_doscompile in
      mingw)
	FFLAGS="-mno-cygwin"
	AC_TRY_LINK(,[      write(*,*) 'Hello world'],
		    [coin_add_fflags="-mno-cygwin $coin_add_fflags"])
	FFLAGS=
      ;;
    esac
  else
    case $build in
      *-cygwin* | *-mingw*)
        case $F77 in
          ifort* | */ifort* | IFORT* | */IFORT* )
            coin_opt_fflags='-MT -O3'
            coin_add_fflags='-fpp -nologo'
            coin_dbg_fflags='-MTd -debug'
          ;;
          compile_f2c*)
            coin_opt_fflags='-MT -O2'
            coin_add_fflags='-nologo -wd4996'
            coin_dbg_fflags='-MTd'
          ;;
        esac
        ;;
      *-linux-*)
        case $F77 in
          ifc* | */ifc* | ifort* | */ifort*)
            coin_opt_fflags="-O3 -ip"
            coin_add_fflags="-cm -w90 -w95"
            coin_dbg_fflags="-g -CA -CB -CS"
            # Check if -i_dynamic is necessary (for new glibc library)
            FFLAGS=
            AC_TRY_LINK(,[      write(*,*) 'Hello world'],[],
                        [coin_add_fflags="-i_dynamic $coin_add_fflags"])
            ;;
          pgf77* | */pgf77* | pgf90* | */pgf90*)
            coin_opt_fflags="-fast"
            coin_add_fflags="-Kieee -pc 64"
            coin_dbg_fflags="-g"
          ;;
        esac
        ;;
      *-ibm-*)
        case "$F77" in
          xlf* | */xlf* | mpxlf* | */mpxlf* )
            coin_opt_fflags="-O -qarch=auto -qcache=auto -qtune=auto -qmaxmem=-1"
            coin_add_fflags="-bmaxdata:0x80000000 -qsuppress=1500-036 -qsuppress=1500-029"
            coin_dbg_fflags="-g -C"
            ;;
        esac
        ;;
      *-hp-*)
        coin_opt_fflags="+O3"
        coin_add_fflags="+U77"
        coin_dbg_fflags="-C -g"
        ;;
      *-*-solaris*)
        coin_opt_fflags="-O4"
        coin_dbg_fflags="-g"
        ;;
      *-sgi-*)
        coin_opt_fflags="-O5 -OPT:Olimit=0"
        coin_dbg_fflags="-g"
        ;;
    esac
  fi

  if test "$ac_cv_prog_f77_g" = yes && test -z "$coin_dbg_fflags" ; then
    coin_dbg_fflags="-g"
  fi

  if test -z "$coin_opt_fflags"; then
    # Try if -O option works if nothing else is set
    FFLAGS=-O
    AC_TRY_LINK(,[      integer i], [coin_opt_fflags="-O"])
  fi

  # if PM doesn't want the warning messages, take them out
  if test x"$coin_skip_warn_fflags" = xyes; then
    coin_warn_fflags=
  fi

  if test x${DBG_FFLAGS+set} != xset; then
    DBG_FFLAGS="$coin_dbg_fflags $coin_add_fflags $coin_warn_fflags"
  fi
  if test x${OPT_FFLAGS+set} != xset; then
    OPT_FFLAGS="$coin_opt_fflags $coin_add_fflags $coin_warn_fflags"
  fi

  DBG_FFLAGS="$DBG_FFLAGS $ADD_FFLAGS"
  OPT_FFLAGS="$OPT_FFLAGS $ADD_FFLAGS"

  if test "$coin_debug_compile" = "true"; then
    FFLAGS="$DBG_FFLAGS"
  else
    FFLAGS="$OPT_FFLAGS"
  fi
else
  FFLAGS="$FFLAGS $ADD_FFLAGS"
  if test x${DBG_FFLAGS+set} != xset; then
    DBG_FFLAGS="$FFLAGS"
  fi
  if test x${OPT_FFLAGS+set} != xset; then
    OPT_FFLAGS="$FFLAGS"
  fi
fi

# If FFLAGS contains -mno-cygwin, CPPFLAGS must have it.
case "$FFLAGS" in
  *-mno-cygwin*)
    if test x${CPPFLAGS+set} != xset ; then
      CPPFLAGS="-mno-cygwin"
    else
      case "$CPPFLAGS" in
        *-mno-cygwin*)
	  ;;
	*)
	  CPPFLAGS="$CPPFLAGS -mno-cygwin"
	  ;;
      esac
    fi ;;
esac

# Try if FFLAGS works
if test "$F77" != "unavailable" ; then
  AC_TRY_LINK(,[      integer i],[],[FFLAGS=])
  if test -z "$FFLAGS"; then
    AC_MSG_WARN([The flags FFLAGS="$FFLAGS" do not work.  I will now just try '-O', but you might want to set FFLAGS manually.])
    FFLAGS='-O'
    AC_TRY_LINK(,[      integer i],[],[FFLAGS=])
    if test -z "$FFLAGS"; then
      AC_MSG_WARN([This value for FFLAGS does not work.  I will continue with empty FFLAGS, but you might want to set FFLAGS manually.])
    fi
  fi
fi

AC_MSG_NOTICE([Fortran compiler options are: $FFLAGS])

AC_ARG_VAR(MPIF77,[Fortran MPI Compiler])
if test x"$MPIF77" = x; then :; else
  AC_MSG_NOTICE([Will use MPI Fortran compiler $MPIF77])
  F77="$MPIF77"
fi

case "$F77" in
  ifort* | */ifort* | IFORT* | */IFORT*)
    AC_COIN_MINGW_LD_FIX
    ;;
esac

AC_LANG_POP([Fortran 77])
]) # AC_COIN_PROG_F77

###########################################################################
#                           COIN_F77_WRAPPERS                             #
###########################################################################

# Calls autoconfs AC_F77_WRAPPERS and does additional corrections to FLIBS

AC_DEFUN([AC_COIN_F77_WRAPPERS],
[AC_BEFORE([AC_COIN_PROG_F77],[$0])dnl
AC_BEFORE([AC_PROG_F77],[$0])dnl

AC_LANG_PUSH([Fortran 77])

save_FLIBS="$FLIBS"

AC_F77_WRAPPERS

if test x"$save_FLIBS" != x; then
  FLIBS="$save_FLIBS"
else
  # This is to correct a missing exclusion in autoconf 2.59
  if test x"$FLIBS" != x; then
    my_flibs=
    for flag in $FLIBS; do
      case $flag in
        -lcrt*.o) ;;
               *) my_flibs="$my_flibs $flag" ;;
      esac
    done
    FLIBS="$my_flibs"
  fi

  case $build in
  # The following is a fix to define FLIBS for ifort on Windows
  # In its original version, it linked in libifcorert.lib or libifcorertd.lib on Windows/ifort explicitly.
  # However, this seem to create a dependency on libifcorert.dll (or libifcorertd.dll) in the executables.
  # This is seem to be unnecessary, libifcorert(d).lib has been removed from the link line.
     *-cygwin* | *-mingw*)
       case "$F77" in
         ifort* | */ifort* | IFORT* | */IFORT*)
           FLIBS="-link $LIBS /NODEFAULTLIB:libc.lib"
#           if "$coin_debug_compile" = true ; then
#             FLIBS="-link $LIBS /NODEFAULTLIB:libc.lib /NODEFAULTLIB:libcmt.lib"
#           else
#             FLIBS="-link $LIBS /NODEFAULTLIB:libc.lib /NODEFAULTLIB:libcmtd.lib"
#           fi
           ;;
         compile_f2c*)
           FLIBS=`$F77 -FLIBS` ;;
       esac;;
     *-hp-*)
         FLIBS="$FLIBS -lm";;
     *-ibm-*)
         FLIBS=`echo $FLIBS | sed 's/-lc)/-lc/g'` ;;
     *-linux-*)
       case "$F77" in
         pgf77* | */pgf77* | pgf90* | */pgf90*)
# ask linker to go through the archives multiple times
# (the Fortran compiler seems to do that automatically...
           FLIBS="-Wl,--start-group $FLIBS -Wl,--end-group" ;;
       esac
  esac
  ac_cv_f77_libs="$FLIBS"
fi

AC_LANG_POP([Fortran 77])

]) # AC_COIN_F77_WRAPPERS

###########################################################################
#                             COIN_FIND_F77                               #
###########################################################################

# Attempt to preempt autoconf by locating an appropriate F77 program. This
# macro will not trigger a fatal error if a suitable compiler cannot be
# found. It should be called before COIN_PROG_F77 or COIN_TRY_FLINK.

AC_DEFUN([AC_COIN_FIND_F77],
[AC_REQUIRE([AC_COIN_ENABLE_DOSCOMPILE])
AC_REQUIRE([AC_COIN_F77_COMPS])
AC_MSG_NOTICE([Trying to determine Fortran compiler name])
AC_CHECK_PROGS([F77],[$coin_f77_comps],[unavailable])
])

# Auxilliary macro to make sure both COIN_PROG_F77 and COIN_FIND_F77 use
# the same search lists for compiler names.
# For *-*-solaris*, promote Studio/Workshop compilers to front of list.
AC_DEFUN([AC_COIN_F77_COMPS],
[case $build in
  *-cygwin* | *-mingw*)
     if test "$enable_doscompile" = msvc ; then
       coin_f77_comps="ifort fl32 compile_f2c"
     else
       coin_f77_comps="gfortran g77 ifort fl32 compile_f2c"
     fi ;;
  *-*-solaris*)
     coin_f77_comps="f95 f90 f77 xlf_r fort77 gfortran g77 pgf90 pgf77 ifort ifc frt af77" ;;
  *) coin_f77_comps="xlf_r fort77 gfortran f77 g77 pgf90 pgf77 ifort ifc frt af77" ;;
 esac
])

###########################################################################
#                          COIN_INIT_AUTOMAKE                             #
###########################################################################

# This macro calls the regular INIT_AUTOMAKE and MAINTAINER_MODE
# macros, and defines additional variables and makefile conditionals,
# that are used in the maintainer parts of the makfile.  It also
# checks if the correct versions of the autotools are used.
#
# This also defines the AC_SUBST variables:
# abs_source_dir     absolute path to source code for this package
# abs_bin_dir        absolute path to the directory where binaries are
#                    going to be installed (prefix/bin)
# abs_lib_dir        absolute path to the directory where libraries are
#                    going to be installed (prefix/lib)
# abs_include_dir    absolute path to the directory where the header files
#                    are installed (prefix/include)

AC_DEFUN([AC_COIN_INIT_AUTOMAKE],
[AC_REQUIRE([AC_PROG_EGREP])

# AC_MSG_NOTICE([Beginning automake initialisation.])
# Stuff for automake
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

coin_have_externals=no
if test "$enable_maintainer_mode" = yes; then

  # If maintainer mode is chosen, we make sure that the correct versions
  # of the tools are used, and that we know where libtool.m4 is (to
  # recreate acinclude.m4)

  AC_SUBST(LIBTOOLM4)
  LIBTOOLM4=
  # Normally, $HOME
  AUTOTOOLS_DFLT=$HOME

  AC_CACHE_CHECK([whether we are using the correct autotools],
                 [ac_cv_use_correct_autotools],
                 [ac_cv_use_correct_autotools=check])

  if test $ac_cv_use_correct_autotools = check; then
    ac_cv_use_correct_autotools=yes
    # Check if we have autoconf
    AC_CHECK_PROG([have_autoconf],[autoconf],[yes],[no])
    if test $have_autoconf = no; then
      AC_MSG_ERROR([You specified you want to use maintainer mode, but I cannot find autoconf in your path.])
    fi

    # Check whether autoconf is the correct version
    correct_version='2.59'
    grep_version=`echo  $correct_version | sed -e 's/\\./\\\\\\./g'`
    AC_MSG_CHECKING([whether we are using the correct version ($correct_version) of autoconf])
    autoconf --version > confauto.out 2>&1
    if $EGREP $grep_version confauto.out >/dev/null 2>&1; then
      AC_MSG_RESULT([yes])
    else
      rm -f confauto.out
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([You don't have the correct version of autoconf as the first one in your path.])
    fi
    rm -f confauto.out

    # Check if the executable autoconf is picked up from the correct location
    AC_MSG_CHECKING([whether autoconf is coming from the correct location])
    autoconf_dir=`which autoconf | sed -e 's=/autoconf=='`
    autoconf_dir=`cd $autoconf_dir; pwd`
    if test x$AUTOTOOLS_DIR = x; then
      want_dir=$AUTOTOOLS_DFLT/bin
    else
      want_dir=$AUTOTOOLS_DIR/bin
    fi
    if test $autoconf_dir = `cd $want_dir; pwd`; then
      AC_MSG_RESULT([yes])
    else
      rm -f confauto.out
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([The autoconf executable should be picked up from \$AUTOTOOLS_DFLT/bin or \$AUTOTOOLS_DIR/bin.])
    fi

    # Check if we have automake
    AC_CHECK_PROG([have_automake],[automake],[yes],[no])
    if test $have_automake = no; then
      AC_MSG_ERROR([You specified you want to use maintainer mode, but I cannot find automake in your path.])
    fi
  
    # Check whether automake is the correct version
    correct_version='1.9.6'
    grep_version=`echo  $correct_version | sed -e 's/\\./\\\\\\./g'`
    AC_MSG_CHECKING([whether we are using the correct version ($correct_version) of automake])
    automake --version > confauto.out 2>&1
    if $EGREP $grep_version confauto.out >/dev/null 2>&1; then
      AC_MSG_RESULT([yes])
    else
      rm -f confauto.out
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([You don't have the correct version of automake as the first one in your path.])
    fi
    rm -f confauto.out

    # Check if the executable automake is picked up from the correct location
    AC_MSG_CHECKING([whether automake is coming from the correct location])
    automake_dir=`which automake | sed -e 's=/automake=='`
    automake_dir=`cd $automake_dir; pwd`
    if test x$AUTOTOOLS_DIR = x; then
      want_dir=$AUTOTOOLS_DFLT/bin
    else
      want_dir=$AUTOTOOLS_DIR/bin
    fi
    if test $automake_dir = `cd $want_dir; pwd`; then
      AC_MSG_RESULT([yes])
    else
      rm -f confauto.out
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([The automake executable should be picked up from \$AUTOTOOLS_DFLT/bin or \$AUTOTOOLS_DIR/bin.])
    fi

    # Check if this is the correct version of libtool (with escaped dots)
    if test x$AUTOTOOLS_DIR = x; then
      want_dir=$AUTOTOOLS_DFLT/share
    else
      want_dir=$AUTOTOOLS_DIR/share
    fi
    correct_version='1.5.22'
    grep_version=`echo  $correct_version | sed -e 's/\\./\\\\\\./g'`
    AC_COIN_CHECK_FILE([$want_dir/libtool/ltmain.sh],
	               [have_ltmain=yes],
                       [have_ltmain=no])
    AC_MSG_CHECKING([whether we are using the correct version ($correct_version) of libtool.])
    if test $have_ltmain = yes; then
    if $EGREP $grep_version $want_dir/libtool/ltmain.sh >/dev/null 2>&1; then
        AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([You don't have the correct version of libtool.])
      fi
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([I cannot find the ltmain.sh file.])
    fi
  fi

  # Check if we can find the libtool file
  if test x$AUTOTOOLS_DIR = x; then
    want_dir=$AUTOTOOLS_DFLT/share
  else
    want_dir=$AUTOTOOLS_DIR/share
  fi
  AC_COIN_CHECK_FILE([$want_dir/aclocal/libtool.m4],
                     [LIBTOOLM4="$want_dir/aclocal/libtool.m4"],
                     [AC_MSG_ERROR([I cannot find the libtool.m4 file.])])

  # Check if we have an Externals file
  if test -r $srcdir/Externals; then
    coin_have_externals=yes
  fi
  # Check if subversion is installed and understands https
  AC_CHECK_PROG([have_svn],[svn],[yes],[no])
  if test x$have_svn = xyes; then
    AC_CACHE_CHECK([whether svn understands https],
                   [ac_cv_svn_understands_https],
                   [svn --version > confauto.out 2>&1
                    if $EGREP https confauto.out >/dev/null 2>&1; then
                      ac_cv_svn_understands_https=yes
                    else
                      ac_cv_svn_understands_https=no
                      have_svn=no
                      ac_cv_prog_have_svn=no
                    fi
                    rm -f confauto.out])
  fi

  # Find the location of the BuildTools directory
  BUILDTOOLSDIR=
  if test -r $srcdir/BuildTools/coin.m4; then
    BUILDTOOLSDIR=$srcdir/BuildTools
  else
    if test -r $srcdir/../BuildTools/coin.m4; then
      BUILDTOOLSDIR=$srcdir/../BuildTools
    else
      if test -r $srcdir/../../BuildTools/coin.m4; then
        BUILDTOOLSDIR=$srcdir/../../BuildTools
      else
        AC_MSG_ERROR(Cannot find the BuildTools directory)
      fi
    fi
  fi
  AC_SUBST(BUILDTOOLSDIR)

  # The following variable is set to the name of the directory where
  # the autotool scripts are located
  AC_SUBST(AUX_DIR)
  AUX_DIR=$ac_aux_dir
fi

# helpful variable for the base directory of this package
abs_source_dir=`cd $srcdir; pwd`
AC_SUBST(abs_source_dir)

# Stuff for example Makefiles
if test x$prefix = xNONE; then
  full_prefix=$ac_default_prefix
else
  full_prefix=$prefix
fi
full_prefix=`cd $full_prefix ; pwd`
AC_SUBST(abs_lib_dir)
abs_lib_dir=$full_prefix/lib
AC_SUBST(abs_include_dir)
abs_include_dir=$full_prefix/include
AC_SUBST(abs_bin_dir)
abs_bin_dir=$full_prefix/bin

AM_CONDITIONAL(HAVE_EXTERNALS,
               test $coin_have_externals = yes && test x$have_svn = xyes)

# AC_MSG_NOTICE([End automake initialisation.])

]) # AC_COIN_INIT_AUTOMAKE

###########################################################################
#                          COIN_CREATE_LIBTOOL                            #
###########################################################################

# This does all the tests necessary to create the libtool script in the
# package base directory.  If this is used, then the COIN_INIT_AUTO_TOOLS
# test in the subdirectories will be able to skip the libtool tests and
# just use the one in the package base directory.

m4_define([AC_COIN_CREATE_LIBTOOL],
[AC_CANONICAL_BUILD

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE

# Get the name of the C compiler and appropriate compiler options
AC_COIN_PROG_CC

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Get the name of the Fortran compiler and appropriate compiler options
AC_COIN_PROG_F77

# Initialize automake and libtool
# AC_MSG_NOTICE([Calling INIT_AUTO_TOOLS from CREATE_LIBTOOL.])
AC_COIN_INIT_AUTO_TOOLS
# AC_MSG_NOTICE([Finished INIT_AUTO_TOOLS from CREATE_LIBTOOL.])
])

###########################################################################
#                         COIN_INIT_AUTO_TOOLS                            #
###########################################################################

# Initialize the auto tools automake and libtool, with all
# modifications we want for COIN packages.
#
# RPATH_FLAGS        link flags for hardcoding path to shared objects

# This is a trick to have this code before AC_COIN_PROG_LIBTOOL
AC_DEFUN([AC_COIN_DISABLE_STATIC],
[
# Test if force_shared has been set
if test "x$1" = xforce_shared; then
  if test x$enable_shared = xno; then
    AC_MSG_ERROR([Shared libraries are disabled by user, but this is not feasible with the given options])
  fi
  enable_shared=yes;
else
  # On Cygwin and AIX, building DLLs doesn't work
  case $build in
    *-cygwin*)
      coin_disable_shared=yes
      platform=Cygwin
    ;;
    *-aix*)
      coin_disable_shared=yes
      platform=AIX
    ;;
    *-mingw*)
      coin_disable_shared=yes
      platform="Msys"
#      case "$CXX" in
#        cl*)
#          coin_disable_shared=yes
#          platform="Msys with cl"
#      ;;
#      esac
    ;;
  esac
fi
if test x"$coin_disable_shared" = xyes; then
  if test x"$enable_shared" = xyes; then
    AC_MSG_WARN([On $platform, shared objects are not supported. I'm disabling your choice.])
  fi
  enable_shared=no
fi
# By default, we only want the shared objects to be compiled
AC_DISABLE_STATIC
])

m4_define([AC_COIN_INIT_AUTO_TOOLS],
[{AC_BEFORE([AC_COIN_PROG_CXX],[$0])
AC_BEFORE([AC_COIN_PROG_CC],[$0])
AC_BEFORE([AC_COIN_PROG_F77],[$0])

# START
AC_COIN_DISABLE_STATIC([$1])

# Initialize automake
AC_COIN_INIT_AUTOMAKE

LIBTOOL=
if test -r ../libtool; then
  coin_config_dir=..
  LIBTOOL='$(SHELL) $(top_builddir)/../libtool'
fi
if test "x$LIBTOOL" = x; then
  if test -r ../../libtool; then
    coin_config_dir=../..
    LIBTOOL='$(SHELL) $(top_builddir)/../../libtool'
  fi
fi

if test "x$LIBTOOL" = x; then
# AC_MSG_NOTICE([Creating libtool script (calling COIN_PROG_LIBTOOL).])
  # Stuff for libtool
  AC_COIN_PROG_LIBTOOL
# AC_MSG_NOTICE([Finished COIN_PROG_LIBTOOL.])
  # set RPATH_FLAGS to the compiler link flags required to hardcode location
  # of the shared objects
  AC_COIN_RPATH_FLAGS($abs_lib_dir)

else

  AC_MSG_NOTICE([Using libtool script in directory $coin_config_dir])
  # get all missing information from the config.log file 

  # output variables and defines
  as_save_IFS=$IFS
  IFS='
'
  for oneline in `cat $coin_config_dir/config.status`; do
    case "$oneline" in
         # First some automake conditionals
      s,@am__fastdep* | s,@AR@* | s,@CPP@*  | s,@CPPFLAGS@* | s,@CXXCPP@*  | \
      s,@RANLIB@* | s,@STRIP@* | s,@ac_ct_AR@* | s,@ac_ct_RANLIB@* | \
      s,@ac_ct_STRIP@* | s,@host* | s,@LN_S@* | s,@RPATH_FLAGS@* | \
      s,@ac_c_preproc_warn_flag@* |  s,@ac_cxx_preproc_warn_flag@* ) 
        command=`echo $oneline | sed -e 's/^s,@//' -e 's/@,/="/' -e 's/,;t t/"/'`
#        echo "$command"
        eval "$command"
        ;;
      s,@DEFS@* )
        command=`echo $oneline | sed -e 's/^s,@DEFS@,/defsline="/' -e 's/,;t t/"/'`
#        echo "$command"
        eval "$command"
        ;;
    esac
  done
  IFS=$as_save_IFS

  # And some defines (assuming here that the packages base dir
  # doesn't have a config.h file
  for word in $defsline; do
#    echo word $word 
    case $word in
      -DHAVE_@<:@A-Z_@:>@*_H=1 | -DSTDC_HEADERS=1 )
        i=`echo $word | sed -e 's/-D/#define /' -e 's/=/ /'`
#        echo dd $i
        echo $i >>confdefs.h
        ;;
    esac
  done
fi

# AC_MSG_NOTICE([End of INIT_AUTO_TOOLS.])

# ToDo
# For now, don't use the -no-undefined flag, since the Makefiles are
# not yet set up that way.  But we need to fix this, when we want
# to comile DLLs under Windows.
LT_LDFLAGS=
AC_SUBST(LT_LDFLAGS)

#END
}])


###########################################################################
#                      COIN_PATCH_LIBTOOL_CYGWIN                          #
###########################################################################

# Patches to libtool for cygwin. Lots for cl, a few for GCC.
# For cl:
# - cygpath is not correctly quoted in fix_srcfile_path
# - paths generated for .lib files is not run through cygpath -w


AC_DEFUN([AC_COIN_PATCH_LIBTOOL_CYGWIN],
[ case "$CXX" in
    cl* | */cl* | CL* | */CL* | icl* | */icl* | ICL* | */ICL*) 
      AC_MSG_NOTICE(Applying patches to libtool for cl compiler)
      sed -e 's|fix_srcfile_path=\"`cygpath -w \"\$srcfile\"`\"|fix_srcfile_path=\"\\\`'"$CYGPATH_W"' \\\"\\$srcfile\\\"\\\`\"|' \
	  -e 's|fix_srcfile_path=\"\"|fix_srcfile_path=\"\\\`'"$CYGPATH_W"' \\\"\\$srcfile\\\"\\\`\"|' \
	  -e 's%compile_deplibs=\"\$dir/\$old_library \$compile_deplibs\"%compile_deplibs="'\`"$CYGPATH_W"' \$dir/\$old_library | sed -e '"'"'sY\\\\\\\\Y/Yg'"'"\`' \$compile_deplibs\"'% \
	  -e 's%compile_deplibs=\"\$dir/\$linklib \$compile_deplibs\"%compile_deplibs="'\`"$CYGPATH_W"' \$dir/\$linklib | sed -e '"'"'sY\\\\\\\\Y/Yg'"'"\`' \$compile_deplibs\"'% \
	  -e 's%lib /OUT:%lib -OUT:%' \
	  -e "s%cygpath -w%$CYGPATH_W%" \
	  -e 's%$AR x \\$f_ex_an_ar_oldlib%bla=\\`lib -nologo -list \\$f_ex_an_ar_oldlib | xargs echo '"$mydos2unix"'\\`; echo \\$bla; for i in \\$bla; do lib -nologo -extract:\\$i \\$f_ex_an_ar_oldlib; done%' \
	  -e 's/$AR t/lib -nologo -list/' \
	  -e 's%f_ex_an_ar_oldlib="\($?*1*\)"%f_ex_an_ar_oldlib='\`"$CYGPATH_W"' \1`%' \ 
	  -e 's%^archive_cmds=.*%archive_cmds="\\$CC -o \\$lib \\$libobjs \\$compiler_flags \\\\\\`echo \\\\\\"\\$deplibs\\\\\\" | \\$SED -e '"\'"'s/ -lc\\$//'"\'"'\\\\\\` -link -dll~linknames="%' \
	  -e 's%old_archive_cmds="lib -OUT:\\$oldlib\\$oldobjs\\$old_deplibs"%old_archive_cmds="if test -r \\$oldlib; then bla=\\"\\$oldlib\\"; else bla=; fi; lib -OUT:\\$oldlib \\\\\\$bla\\$oldobjs\\$old_deplibs"%' \
      libtool > conftest.bla

      mv conftest.bla libtool
      chmod 755 libtool
      ;;
    *)
      AC_MSG_NOTICE(Applying patches to libtool for GNU compiler)
      sed -e 's|fix_srcfile_path=\"`cygpath -w \"\$srcfile\"`\"|fix_srcfile_path=\"\\\`'"$CYGPATH_W"' \\\"\\$srcfile\\\"\\\`\"|' \
	  -e 's|"lib /OUT:\\$oldlib\\$oldobjs\\$old_deplibs"|"\\$AR \\$AR_FLAGS \\$oldlib\\$oldobjs\\$old_deplibs~\\$RANLIB \\$oldlib"|' \
	  -e 's|libext="lib"|libext="a"|' \
      libtool > conftest.bla

      mv conftest.bla libtool
      chmod 755 libtool
      ;;
  esac ]) # COIN_PATCH_LIBTOOL_CYGWIN

###########################################################################
#                    COIN_PATCH_LIBTOOL_SOLARIS                           #
###########################################################################
# If we want to do a 64-bit build with GCC on Solaris, the system search
# libraries need to point to 64-bit subdirectories. If they do not already do
# that, fix them. This patch is evolving, as are GCC compilers.  GCC 4.2.1
# reports the correct search list, given the correct call. GCC 4.1.1 does not.
# `Correct call' means -m64 is specified. `Correct search list' seems to amount
# to prepending the list of 64-bit subdirectories to the 32-bit directories.
# Both SPARC and x86 have this issue, but a different hardware id string is
# required depending on the underlying CPU. The macro executes isainfo to get
# the string. This will fail, of course, if we're cross-compiling. The
# alternative is to fail on a regular basis each time a new CPU identifier is
# needed. This macro will also fail if the search list reported with
# -print-search-dirs differs between the C, C++, and Fortran compilers; each
# have their own setting in libtool.  If GCC reports the correct search list
# given the -m64 flag, the best solution is to define CC='gcc -m64', and
# similarly for CXX, F77, so that libtool will make the correct call.
###########################################################################
AC_DEFUN([AC_COIN_PATCH_LIBTOOL_SOLARIS],
[ if test "$GCC" = yes && \
     (echo $CXXFLAGS $CFLAGS $FFLAGS | $EGREP 'm64' >/dev/null 2>&1) ; then
    hdwisa=`isainfo | sed -e 's/\(@<:@^ @:>@*\) .*$/\1/'`
    if `$EGREP 'sys_lib_search_path_spec=' libtool | $EGREP -v $hdwisa >/dev/null 2>&1` ; then
      AC_MSG_NOTICE([Applying patches to libtool for 64-bit GCC compilation])
      fixlibtmp=`$CC -m64 -print-search-dirs | $EGREP '^libraries:'`
      fixlibtmp=`echo $fixlibtmp | sed -e 's/libraries: =//' -e 's/:/ /g'`
      if `echo "$fixlibtmp" | $EGREP -v $hdwisa  >/dev/null 2>&1` ; then
	# AC_MSG_NOTICE(Compensating for broken gcc)
	for lib in $fixlibtmp ; do
	  if test -d "${lib}${hdwisa}" ; then
	    syslibpath64="$syslibpath64 ${lib}${hdwisa}/"
	  fi
	done
	syslibpath64="${syslibpath64} ${fixlibtmp}"
      else
	syslibpath64="$fixlibtmp"
      fi
      sed -e 's|sys_lib_search_path_spec=".*"|sys_lib_search_path_spec="'"$syslibpath64"'"|' libtool > conftest.bla
      mv conftest.bla libtool
      chmod 755 libtool  
    fi
    # AC_MSG_NOTICE(Result is )
    # $EGREP 'sys_lib_search_path_spec=' libtool
  fi ])	# COIN_PATCH_LIBTOOL_SOLARIS

###########################################################################
#                           COIN_PROG_LIBTOOL                             #
###########################################################################

# Setup the libtool stuff together with any modifications to make it
# work on additional platforms

AC_DEFUN([AC_COIN_LIBTOOL_WRAPPER],
[AC_BEFORE([AC_COIN_BLA],[$0])
AC_PROG_LIBTOOL])

AC_DEFUN([AC_COIN_BLA],[echo bla])

AC_DEFUN([AC_COIN_PROG_LIBTOOL],
[# No longer needed now that CPPFLAGS is correctly set -- lh, 061214 --
 # AC_REQUIRE([AC_COIN_DLFCN_H])

# NEW: If libtool exists in the directory higher up, we use that one
#      instead of creating a new one

# It turns out that the code for AC_PROG_LIBTOOL is somehow AC_REQUIRED
# out in front of this macro body. You'll notice that LIBTOOL is already
# defined here.  We'll have to count on this macro not being called if libtool
# already exists, or at least move the libtool fixes outside the conditional.
# AC_MSG_NOTICE([Entering coin_prog_libtool, LIBTOOL = "$LIBTOOL".])
# This test is therefore removed.  -- lh, 061214 --
# if test "x$LIBTOOL" = x; then

# AC_MSG_NOTICE([Calling PROG_LIBTOOL.])
  AC_PROG_LIBTOOL
# AC_MSG_NOTICE([Finished PROG_LIBTOOL.])
  AC_SUBST(ac_c_preproc_warn_flag)
  AC_SUBST(ac_cxx_preproc_warn_flag)

  AC_MSG_NOTICE([Build is "$build".])
  mydos2unix='| dos2unix'
  case $build in
    *-mingw*)
      CYGPATH_W=echo
      mydos2unix=
      ;;
  esac

  case $build in
    # Here we need to check if -m32 is specified.  If so, we need to correct
    # sys_lib_search_path_spec
    *x86_64-*)
      if test "$GCC" = yes && (echo $CXXFLAGS $CFLAGS $FFLAGS | $EGREP 'm32' >& /dev/null); then 
        AC_MSG_NOTICE(Applying patches to libtool for 32bit compilation)
        sed -e 's|sys_lib_search_path_spec=".*"|sys_lib_search_path_spec="/lib /usr/lib"|' libtool > conftest.bla
        mv conftest.bla libtool
        chmod 755 libtool  
      fi
      ;;

    *-solaris*)
      AC_COIN_PATCH_LIBTOOL_SOLARIS
      ;;
    # Cygwin. Ah, cygwin. Too big and ugly to inline; see the macro.
    *-cygwin* | *-mingw*)
      AC_COIN_PATCH_LIBTOOL_CYGWIN
      ;;
    *-darwin*)
      AC_MSG_NOTICE(Applying patches to libtool for Darwin)
      sed -e 's/verstring="${wl}-compatibility_version ${wl}$minor_current ${wl}-current_version ${wl}$minor_current.$revision"/verstring="-compatibility_version $minor_current -current_version $minor_current.$revision"/' \
        -e 's/ -dynamiclib / -dynamiclib -single_module /g' \
      libtool > conftest.bla

      mv conftest.bla libtool
      chmod 755 libtool
      ;;
  esac
# This fi matches the commented `if test "x$LIBTOOL" = x;' up at the head of
# the macro. -- lh, 061214 --
# fi

# AC_MSG_NOTICE([End libtool initialisation.])
]) # AC_COIN_PROG_LIBTOOL

# This is a trick to force the check for the dlfcn header to be done before
# the checks for libtool
# No longer needed now that CPPFLAGS is correctly set.  -- lh, 061214 --
# ACDEFUN([AC_COIN_DLFCN_H],
# [AC_LANG_PUSH(C)
# AC_COIN_CHECK_HEADER([dlfcn.h])
# AC_LANG_POP(C)
# ]) # AC_COIN_DLFCN_H

###########################################################################
#                            COIN_RPATH_FLAGS                             #
###########################################################################

# This macro, in case shared objects are used, defines a variable
# RPATH_FLAGS that can be used by the linker to hardwire the library
# search path for the given directories.  This is useful for example
# Makefiles

AC_DEFUN([AC_COIN_RPATH_FLAGS],
[RPATH_FLAGS=

if test $enable_shared = yes; then
  case $build in
    *-linux-*)
      if test "$GXX" = "yes"; then
        RPATH_FLAGS=
        for dir in $1; do
          RPATH_FLAGS="$RPATH_FLAGS -Wl,--rpath -Wl,$dir"
        done
      fi ;;
    *-darwin*)
        RPATH_FLAGS=nothing ;;
    *-ibm-*)
      case "$CXX" in
      xlC* | */xlC* | mpxlC* | */mpxlC*)
        RPATH_FLAGS=nothing ;;
      esac ;;
    *-hp-*)
        RPATH_FLAGS=nothing ;;
    *-mingw32)
        RPATH_FLAGS=nothing ;;
    *-*-solaris*)
        RPATH_FLAGS=
        for dir in $1; do
          RPATH_FLAGS="$RPATH_FLAGS -R$dir"
        done
  esac

  if test "$RPATH_FLAGS" = ""; then
    AC_MSG_WARN([Could not automatically determine how to tell the linker about automatic inclusion of the path for shared libraries.  The test examples might not work if you link against shared objects.  You will need to set the LD_LIBRARY_PATH, DYLP_LIBRARY_PATH, or LIBDIR variable manually.])
  fi
  if test "$RPATH_FLAGS" = "nothing"; then
    RPATH_FLAGS=
  fi
fi

AC_SUBST(RPATH_FLAGS)
]) # AC_COIN_RPATH_FLAGS

###########################################################################
#                        COIN_LINK_INPUT_CMD                              #
###########################################################################

# This macro determines which command should be used to "link" files
# that are input to the generated executables.  On Windows, the codes
# using the native Windows system libraries, cannot understand symbolic
# links, and a copy should be used instead of 'ln -s'.
# The result is stored in coin_link_input_cmd

AC_DEFUN([AC_COIN_LINK_INPUT_CMD],
[AC_REQUIRE([AC_PROG_LN_S])
AC_BEFORE([AC_COIN_PROG_CC], [$0])
AC_BEFORE([AC_COIN_ENABLE_DOSCOMPILE], [$0])

AC_MSG_CHECKING([which command should be used to link input files])
coin_link_input_cmd="$LN_S"
if test "$enable_doscompile" = mingw; then
  coin_link_input_cmd=cp
fi
case "$CC" in
  cl* | */cl* | CL* | */CL* | icl* | */icl* | ICL* | */ICL*)
    coin_link_input_cmd=cp ;;
esac
AC_MSG_RESULT($coin_link_input_cmd)
])

###########################################################################
#                              COIN_FINALIZE                              #
###########################################################################

# This macro should be called at the very end of the configure.ac file.
# It creates the output files (by using AC_OUTPUT), and might do some other
# things (such as generating links to data files in a VPATH configuration).
# It also prints the "success" message.
# Note: If the variable coin_skip_ac_output is set to yes, then no output
# files are written.

AC_DEFUN([AC_COIN_FINALIZE],
[
AC_REQUIRE([AC_COIN_LINK_INPUT_CMD])
if test x$coin_skip_ac_output != xyes; then

  FADDLIBS="$ADDLIBS"
  if test x"$coin_need_flibs" = xyes; then
    ADDLIBS="$ADDLIBS $FLIBS"
  fi

  # library extension
  AC_SUBST(LIBEXT)
  case "$CC" in
    cl* | */cl* | CL* | */CL* | icl* | */icl* | ICL* | */ICL*)
         LIBEXT=lib ;;
      *) LIBEXT=a ;;
  esac

  # Define VPATH_DISTCLEANFILES to be everything that needs to be
  # cleaned for distclean in a vpath configuration
  AC_SUBST(VPATH_DISTCLEANFILES)
  VPATH_DISTCLEANFILES="$coin_vpath_link_files"

  # Take out subdirectories if their configuration concluded that they
  # don't need to be compiled
  if test x"$coin_ac_skip_subdirs" != x; then
    new_subdirs=
    for i in $subdirs; do
      skipme=no
      for j in $coin_ac_skip_subdirs; do
        if test $i = $j; then
          skipme=yes;
        fi
      done
      if test $skipme = no; then
        new_subdirs="$new_subdirs $i"
      fi
    done
    subdirs="$new_subdirs"
  fi

  AC_OUTPUT

  if test x"$coin_vpath_link_files" = x; then : ; else
    lnkcmd="$coin_link_input_cmd"
    if test "$lnkcmd" = cp; then
      AC_MSG_NOTICE(Copying data files for VPATH configuration)
    else
      AC_MSG_NOTICE(Creating VPATH links for data files)
    fi
    for file in $coin_vpath_link_files; do
      dir=`AS_DIRNAME(["./$file"])`
      if test -d $dir; then : ; else
        AS_MKDIR_P($dir)
      fi
      rm -f $file
      $lnkcmd $abs_source_dir/$file $file
    done
  fi

  if test x$coin_projectdir = xyes; then
    AC_MSG_NOTICE([Configuration of $PACKAGE_NAME successful])
  else
    AC_MSG_NOTICE([Main configuration of $PACKAGE_NAME successful])
  fi
else
  AC_MSG_NOTICE([No configuration of $PACKAGE_NAME necessary])
fi

]) #AC_COIN_FINALIZE

###########################################################################
#                             COIN_VPATH_LINK                             #
###########################################################################

# This macro makes sure that a symbolic link is created to a file in
# the source code directory tree if we are in a VPATH compilation, and
# if this package is the main package to be installed

AC_DEFUN([AC_COIN_VPATH_LINK],
[AC_REQUIRE([AC_COIN_CHECK_VPATH])
if test $coin_vpath_config = yes; then
  coin_vpath_link_files="$coin_vpath_link_files $1"
fi
]) #AC_COIN_VPATH_LINK

###########################################################################
#                       COIN_ENABLE_GNU_PACKAGES                          #
###########################################################################

# This macro defined the --enable-gnu-packages flag.  This can be used
# to check if a user wants to compile GNU packges (such as readline or
# zlib) into the executable.  By default, GNU packages are disabled.
# This also defines the automake conditional COIN_ENABLE_GNU_PACKAGES

AC_DEFUN([AC_COIN_ENABLE_GNU_PACKAGES],
[AC_ARG_ENABLE([gnu-packages],
               [AC_HELP_STRING([--enable-gnu-packages],
                               [compile with GNU packages (disabled by default)])],
	       [coin_enable_gnu=$enableval],
	       [coin_enable_gnu=no])
]) # AC_COIN_ENABLE_GNU_PACKAGES

###########################################################################
#                           COIN_CHECK_GNU_ZLIB                           #
###########################################################################

# This macro checks for the libz library.

AC_DEFUN([AC_COIN_CHECK_GNU_ZLIB],
[AC_REQUIRE([AC_COIN_ENABLE_GNU_PACKAGES])
AC_BEFORE([AC_COIN_PROG_CXX],[$0])
AC_BEFORE([AC_COIN_PROG_CC],[$0])
AC_BEFORE([AC_COIN_PROG_F77],[$0])
AC_BEFORE([$0],[AC_COIN_FINALIZE])

coin_has_zlib=no
if test $coin_enable_gnu = yes; then
  AC_COIN_CHECK_HEADER([zlib.h],[coin_has_zlib=yes])

  if test $coin_has_zlib = yes; then
    AC_CHECK_LIB([z],[gzopen],
                 [ADDLIBS="-lz $ADDLIBS"],
                 [coin_has_zlib=no])
  fi

  if test $coin_has_zlib = yes; then
    AC_DEFINE([COIN_HAS_ZLIB],[1],[Define to 1 if zlib is available])
  fi
fi

AM_CONDITIONAL(COIN_HAS_ZLIB,test x$coin_has_zlib = xyes)
]) # AC_COIN_CHECK_GNU_ZLIB


###########################################################################
#                          COIN_CHECK_GNU_BZLIB                           #
###########################################################################

# This macro checks for the libbz2 library.

AC_DEFUN([AC_COIN_CHECK_GNU_BZLIB],
[AC_REQUIRE([AC_COIN_ENABLE_GNU_PACKAGES])
AC_BEFORE([AC_COIN_PROG_CXX],[$0])
AC_BEFORE([AC_COIN_PROG_CC],[$0])
AC_BEFORE([AC_COIN_PROG_F77],[$0])
AC_BEFORE([$0],[AC_COIN_FINALIZE])

coin_has_bzlib=no
if test $coin_enable_gnu = yes; then
  AC_COIN_CHECK_HEADER([bzlib.h],[coin_has_bzlib=yes])

  if test $coin_has_bzlib = yes; then
    AC_CHECK_LIB([bz2],[BZ2_bzReadOpen],
                 [ADDLIBS="-lbz2 $ADDLIBS"],
                 [coin_has_bzlib=no])
  fi

  if test $coin_has_bzlib = yes; then
    AC_DEFINE([COIN_HAS_BZLIB],[1],[Define to 1 if bzlib is available])
  fi
fi
]) # AC_COIN_CHECK_GNU_BZLIB


###########################################################################
#                         COIN_CHECK_GNU_READLINE                         #
###########################################################################

# This macro checks for GNU's readline.  It verifies that the header
# readline/readline.h is available, and that the -lreadline library
# contains "readline".  It is assumed that #include <stdio.h> is included
# in the source file before the #include<readline/readline.h>

AC_DEFUN([AC_COIN_CHECK_GNU_READLINE],
[AC_REQUIRE([AC_COIN_ENABLE_GNU_PACKAGES])
AC_BEFORE([AC_COIN_PROG_CXX],[$0])
AC_BEFORE([AC_COIN_PROG_CC],[$0])
AC_BEFORE([AC_COIN_PROG_F77],[$0])
AC_BEFORE([$0],[AC_COIN_FINALIZE])

coin_has_readline=no
if test $coin_enable_gnu = yes; then
  AC_COIN_CHECK_HEADER([readline/readline.h],
                       [coin_has_readline=yes],[],
                       [#include <stdio.h>])

  coin_save_LIBS="$LIBS"
  LIBS=
  # First we check if tputs and friends are available
  if test $coin_has_readline = yes; then
    AC_SEARCH_LIBS([tputs],[ncurses termcap curses],[],
                   [coin_has_readline=no])
  fi

  # Now we check for readline
  if test $coin_has_readline = yes; then
    AC_CHECK_LIB([readline],[readline],
                 [ADDLIBS="-lreadline $LIBS $ADDLIBS"],
                 [coin_has_readline=no])
  fi

  if test $coin_has_readline = yes; then
    AC_DEFINE([COIN_HAS_READLINE],[1],[Define to 1 if readline is available])
  fi

  LIBS="$coin_save_LIBS"
fi
]) # AC_COIN_CHECK_GNU_READLINE

###########################################################################
#                             COIN_DATA_PATH                              #
###########################################################################

# This macro defines a preprocessor macro with the absolute path to a
# subdirectory of Data.  The argument of this macro is the name of the
# subdirectory (in correct case), and the name of the macro is
# COIN_DATA_DIR_PATH, where dir is replaced by the capitalized name of
# the directory.  The path ends with a separator ("/" for linux and
# '\\' for Windows).  The default value for this path can be
# overwritten with the input variable with the same name
# (COIN_DATA_DIR_PATH).  At this point we chech only for the
# $srcdir/../Data subdirectory.

AC_DEFUN([AC_COIN_DATA_PATH],
[AC_MSG_CHECKING([absolute path to data directory $1])

AC_ARG_VAR(m4_toupper(COIN_DATA_$1_PATH),[Set to absolute path to Data/$1 subdirectory])

if test x"$m4_toupper(COIN_DATA_$1_PATH)" = x; then
  m4_toupper(COIN_DATA_$1_PATH)=`cd $srcdir/../Data/$1; pwd`
fi

# Under Cygwin, use Windows path.  Add separator
case $build in
  *-cygwin*)
    m4_toupper(COIN_DATA_$1_PATH)=`cygwin -w $m4_toupper(COIN_DATA_$1_PATH)`\\
    ;;
  *)
    m4_toupper(COIN_DATA_$1_PATH)="$m4_toupper(COIN_DATA_$1_PATH)/"
    ;;
esac

if test -r $m4_toupper(COIN_DATA_$1_PATH); then
  AC_DEFINE_UNQUOTED(m4_toupper(COIN_DATA_$1_PATH),["$m4_toupper(COIN_DATA_$1_PATH)"],
            [Define to absolute path for Data subdirectory $1])
  AC_MSG_RESULT($m4_toupper(COIN_DATA_$1_PATH))
else
  AC_MSG_ERROR(Directory $m4_toupper(COIN_DATA_$1_PATH) does not exist)
fi
]) # AC_COIN_HAS_DATA

###########################################################################
#                       COIN_LINK_FROM_FILELIST                           #
###########################################################################

# This macro creates links (or copies, if necessary) to files listed
# as content in a text file (second argument) into a target directory
# (first argument), which is created if it doesn't exist yet.  If s link
# already exists, nothing happens.

AC_DEFUN([AC_COIN_LINKCOPY_FROM_FILELIST],
[cmd="$3"
if test -e $srcdir/$2 ; then
  my_target_dir="$1"
  my_link_files=`cat $srcdir/$2`
  my_dirname=`AS_DIRNAME($2)`
#  if test -e $my_target_dir; then : ; else
#    AS_MKDIR_P($my_target_dir)
#  fi
  for i in $my_link_files; do
    #rm -rf $my_target_dir/$i
    if test -e $my_target_dir/$i; then : ; else
      dirn2=`AS_DIRNAME($my_target_dir/$i)`
      if test -e $dirn2; then : ; else
        AS_MKDIR_P($dirn2)
      fi
      $cmd $abs_source_dir/$my_dirname/$i $my_target_dir/$i
    fi
  done
else
  AC_MSG_WARN([File list file $2 missing!])
fi
])

AC_DEFUN([AC_COIN_LINK_FROM_FILELIST], 
[
AC_REQUIRE([AC_COIN_LINK_INPUT_CMD])
echo Creating links in $1 ...
AC_COIN_LINKCOPY_FROM_FILELIST($1, $2, $coin_link_input_cmd)
])

###########################################################################
#                       COIN_COPY_FROM_FILELIST                           #
###########################################################################

# Like COIN_LINK_FROM_FILELIST, but copies the files.

AC_DEFUN([AC_COIN_COPY_FROM_FILELIST], 
[
echo Creating copies in $1 ...
AC_COIN_LINKCOPY_FROM_FILELIST($1, $2, [cp])
])

###########################################################################
#                          COIN_EXAMPLE_FILES                             #
###########################################################################

# This macro determines the names of the example files (using the
# argument in an "ls" command) and sets up the variables EXAMPLE_FILES
# and EXAMPLE_CLEAN_FILES.  If this is a VPATH configuration, it also
# creates soft links to the example files.

AC_DEFUN([AC_COIN_EXAMPLE_FILES],
[AC_REQUIRE([AC_COIN_CHECK_GNU_ZLIB])
AC_REQUIRE([AC_COIN_CHECK_VPATH])
files=`cd $srcdir; ls $1`
# We need to do the following loop to make sure that are no newlines
# in the variable
for file in $files; do
  EXAMPLE_FILES="$EXAMPLE_FILES $file"
done
if test $coin_vpath_config = yes; then
  lnkcmd=
  if test "$enable_doscompile" != no; then
    lnkcmd=cp
  fi
  case "$CC" in
    cl* | */cl* | CL* | */CL* | icl* | */icl* | ICL* | */ICL*)
      lnkcmd=cp ;;
  esac
  if test "$lnkcmd" = cp; then
    AC_MSG_NOTICE([Copying example files ($1)])
  else
    AC_PROG_LN_S
    AC_MSG_NOTICE([Creating links to the example files ($1)])
    lnkcmd="$LN_S"
  fi
  for file in $EXAMPLE_FILES; do
    rm -f $file
    $lnkcmd $srcdir/$file $file
  done
  EXAMPLE_CLEAN_FILES="$EXAMPLE_CLEAN_FILES $1"
else
  EXAMPLE_CLEAN_FILES="$EXAMPLE_CLEAN_FILES"
fi

# In case there are compressed files, we create a variable with the
# uncompressed names
EXAMPLE_UNCOMPRESSED_FILES=
for file in $EXAMPLE_FILES; do
  case $file in
    *.gz)
      EXAMPLE_UNCOMPRESSED_FILES="$EXAMPLE_UNCOMPRESSED_FILES `echo $file | sed -e s/.gz//`"
      ;;
  esac
done

AC_SUBST(EXAMPLE_UNCOMPRESSED_FILES)
AC_SUBST(EXAMPLE_FILES)
AC_SUBST(EXAMPLE_CLEAN_FILES)
]) #AC_COIN_EXAMPLE_FILES

###########################################################################
#                            COIN_HAS_PROJECT                             #
###########################################################################

# This macro sets up usage of a Coin package.  It defines the
# PKGSRCDIR and PKGOBJDIR variables, refering to the main source and
# object directory of the package, respectively.  It also defines
# a COIN_HAS_PKG preprocessor macro and makefile conditional.  The
# argument should be the name (Pkg) of the project (in correct lower
# and upper case).
#
# It is also possible to specify a preinstalled version of this Coin
# package.  In this case, only PKGINSTDIR is set to the prefix of that
# installation.  If the user is providing such a preinstalled version,
# this macro checks for a header (second argument) and the existance
# of a library file (third argument - only the stub without extension
# is to be given, e.g., "libCoinUtils").

AC_DEFUN([AC_COIN_HAS_PROJECT],
[AC_MSG_CHECKING([for COIN project $1])
# First check, if the sub-project is actually available (ToDo: allow
# other locations)

m4_tolower(coin_has_$1)=notGiven
if test x"$COIN_SKIP_PROJECTS" != x; then
  for dir in $COIN_SKIP_PROJECTS; do
    if test $dir = $1; then
      m4_tolower(coin_has_$1)=skipping
    fi
  done
fi

AC_SUBST(m4_toupper($1SRCDIR))
AC_SUBST(m4_toupper($1OBJDIR))
AC_SUBST(m4_toupper($1DOCDIR))
AC_SUBST(m4_toupper($1INSTDIR))
m4_toupper($1SRCDIR)=This_Variable_Is_Not_Set
m4_toupper($1OBJDIR)=This_Variable_Is_Not_Set
m4_toupper($1DOCDIR)=This_Variable_Is_Not_Set
m4_toupper($1INSTDIR)=This_Variable_Is_Not_Set

if test $m4_tolower(coin_has_$1) != skipping; then
  if test $PACKAGE_TARNAME = m4_tolower($1); then
    m4_tolower(coin_has_$1)=.
  else
    AC_ARG_WITH([m4_tolower($1)-instdir],
	        AC_HELP_STRING([--with-m4_tolower($1)-instdir],
                               [prefix of installation directory for precompiled $1 package]),
                [m4_tolower(coin_has_$1)=installed
		 if test -d "$withval"; then : ; else
		   AC_MSG_ERROR([argument for --with-m4_tolower($1)-instdir not a directory])
		 fi
                 m4_toupper($1INSTDIR)=`cd $withval; pwd`], [])

    AC_ARG_WITH([m4_tolower($1)-srcdir],
	        AC_HELP_STRING([--with-m4_tolower($1)-srcdir],
                               [source directory for $1 package]),
                [if test "$m4_tolower(coin_has_$1)" = installed; then
		   AC_MSG_ERROR([--with-$m4_tolower($1)-srcdir cannot be given together with --with-m4_tolower($1)-instdir])
		 fi
		 if test -d "$withval"; then : ; else
		   AC_MSG_ERROR([argument for --with-m4_tolower($1)-srcdir not a directory])
		 fi
		 m4_tolower(coin_has_$1)=$withval
                 m4_toupper($1SRCDIR)=`cd $withval; pwd`], [])

    AC_ARG_WITH([m4_tolower($1)-objdir],
	        AC_HELP_STRING([--with-m4_tolower($1)-objdir],
                               [object directory for $1 package (if different from source directory in VPATH compilation)]),
                [if test "$m4_tolower(coin_has_$1)" = installed; then
		   AC_MSG_ERROR([--with-m4_tolower($1)-objdir cannot be given together with --with-m4_tolower($1)-instdir])
		 fi
		 if test "$m4_tolower(coin_has_$1)" = notGiven; then
		   AC_MSG_ERROR([--with-m4_tolower($1)-objdir must be accompanied by --with-m4_tolower($1)-srcdir])
		 fi
		 if test -d "$withval"; then : ; else
		   AC_MSG_ERROR([argument for --with-m4_tolower($1)-objdir not a directory])
		 fi
                 m4_toupper($1OBJDIR)=`cd $withval; pwd`], [m4_toupper($1OBJDIR)="$m4_toupper($1SRCDIR)"])

    # if we still don't have a location, check for generic flag
    if test "$m4_tolower(coin_has_$1)" = notGiven; then
      AC_ARG_WITH([coin-instdir],
	          AC_HELP_STRING([--with-coin-instdir],
                                 [prefix of installation directory for precompiled COIN packages]),
                  [if test -d "$withval"; then : ; else
		     AC_MSG_ERROR([argument for --with-coin-instdir not a directory])
		   fi
		   m4_ifvaln([$2],
                     [if test -r $withval/include/coin/$2; then
                        m4_tolower(coin_has_$1)=installed
			m4_toupper($1INSTDIR)=`cd $withval; pwd`
                      fi],
		     [AC_MSG_WARN([Cannot verify availability of $1 with --with-coin-instdir flag.  Project manager needs to provide name of header to look for])])],
                   [])
    fi

    if test "$m4_tolower(coin_has_$1)" = notGiven; then
      if test -d $srcdir/../$1; then
        m4_tolower(coin_has_$1)=../$1
      fi
    fi
  fi
fi

if test $m4_tolower(coin_has_$1) != notGiven &&
   test $m4_tolower(coin_has_$1) != skipping; then
  # Set the #define if the component is available
  AC_DEFINE(m4_toupper(COIN_HAS_$1),[1],[Define to 1 if the $1 package is used])

  if test $m4_tolower(coin_has_$1) = installed; then
    m4_toupper($1DOCDIR)=$m4_toupper($1INSTDIR)/share/doc/coin/$1
  else
    if test "$m4_toupper($1OBJDIR)" = This_Variable_Is_Not_Set; then
      # Set the variables for source and object code location
      m4_toupper($1SRCDIR)=`cd $srcdir/$m4_tolower(coin_has_$1); pwd`
      m4_toupper($1OBJDIR)=`pwd`/$m4_tolower(coin_has_$1)
      m4_toupper($1DOCDIR)=$abs_lib_dir/../share/doc/coin/$1
    else
      # This is just a guess:
      m4_toupper($1DOCDIR)=$m4_toupper($1OBJDIR)/../share/doc/coin/$1
    fi
  fi
fi

  # Define the Makefile conditional
AM_CONDITIONAL(m4_toupper(COIN_HAS_$1),
               [test $m4_tolower(coin_has_$1) != notGiven &&
                test $m4_tolower(coin_has_$1) != skipping])
AM_CONDITIONAL(m4_toupper(COIN_HAS_$1_PREINSTALLED),
               [test $m4_tolower(coin_has_$1) = installed])

if test $m4_tolower(coin_has_$1) = installed; then
  AC_MSG_RESULT([installed in $m4_toupper($1INSTDIR)])
  # Check for header
  m4_ifvaln([$2],[tmp=$m4_toupper($1INSTDIR)/include/coin
                  AC_MSG_CHECKING(whether header $2 is available in $tmp)
                  if test -r $tmp/$2; then
		    AC_MSG_RESULT([yes])
		  else
		    AC_MSG_RESULT([no])
		    AC_MSG_ERROR([Header $2 not available in $1 installation directory])
		  fi])
  # Check for library (this is not qo reliable)
  m4_ifvaln([$3],[tmp=$m4_toupper($1INSTDIR)/lib
                  AC_MSG_CHECKING([whether library $3.* is available in $tmp])
		  bla=`ls $tmp/$3.* 2>/dev/null | head -n 1`
		  if test x"$bla" = x; then
		    AC_MSG_RESULT([no])
		    AC_MSG_ERROR([Library $3.* not available in $1 installation directory])
		  else
		    AC_MSG_RESULT([found $bla])
		  fi])
else
  AC_MSG_RESULT([$m4_tolower(coin_has_$1)])
fi
]) # AC_COIN_HAS_PROJECT

###########################################################################
#                        COIN_HAS_USER_LIBRARY                            #
###########################################################################
# This macro sets up usage of a user library with header files. The assumption
# is that the header file(s) and library do not reside in standard system
# directories, hence both the include directory and link flags must be
# specified. There are two mandatory arguments and two optional arguments.
#
# The first argument (mandatory) should be a name (LibraryName) for the
# library.  The second argument (mandatory) should be an abbreviation in
# upper case letters (LBRY) for the library. Ultimately, the macro will
# specify two variables, LBRYINCDIR and LBRYLIB, to be substituted in files
# generated during configuration; a preprocessor symbol COIN_HAS_LBRY; and a
# matching automake conditional COIN_HAS_LBRY. LBRYINCDIR should specify the
# directory containing include files for the library. LBRYLIB should specify
# the flags necessary to link to the library. A variable coin_has_lbry will
# be set to true or false, as appropriate. A variable lbry_libcheck will be
# be set to yes or no; no indicates link checks should not be attempted.
#
# The macro defines three configure arguments, --with-libraryname-incdir,
# --with-libraryname-lib, and --disable-libraryname-libcheck, by converting
# LibraryName to lower case.
#
# LBRYINCDIR and LBRYLIB can be specified as environment variables or as
# part of the configure command using --with-libraryname-incdir and
# --with-libraryname-lib, respectively. Command line arguments override
# environment variables.
#
# If a third argument is given, it should specify a file in LBRYINCDIR.  The
# macro will check for the presence of the file. If a fourth argument is given,
# it should specify a function name, `fname'.  The macro will attempt to link a
# trivial program containing a parameterless call to the function, `fname()',
# using the LBRYLIB flags. The link check uses C as the language; this has been
# adequate to date but has limitations. It is possible to disable the link
# check by specifying --disable-libraryname-libcheck. This is a workaround for
# instances where the link check does not work properly, for whatever reason.
# If you're trying to link to a Fortran library, consider using F77_FUNC or
# FC_FUNC to obtain a mangled fname appropriate for use from C code. For a C++
# library, you're on your own unless the library declares some function with
# extern "C" linkage. Otherwise, you'll have to somehow find the mangled C++
# name.

AC_DEFUN([AC_COIN_HAS_USER_LIBRARY],
[ AC_REQUIRE([AC_COIN_PROJECTDIR_INIT])
  AC_MSG_CHECKING(if user provides library for $1)

# Check for header file directory

  AC_ARG_WITH(m4_tolower($1)-incdir,
      AS_HELP_STRING([--with-m4_tolower($1)-incdir],
		     [specify the header file directory for library $1]),
      [$2INCDIR=`cd $withval; pwd`])

# Check for library directory

  AC_ARG_WITH(m4_tolower($1)-lib,
      AS_HELP_STRING([--with-m4_tolower($1)-lib],
		     [specify the flags used to link with the library $1]),
      [$2LIB=$withval])

# Switch to disable library check if requested

  AC_ARG_ENABLE(m4_tolower($1)-libcheck,
      AS_HELP_STRING([--enable-m4_tolower($1)-libcheck],
		     [use disable-m4_tolower($1)-libcheck to skip the link
		      check at configuration time]),
      [m4_tolower($1)_libcheck=$enableval],
      [m4_tolower($1)_libcheck=yes])

# At this point, if we're going to use the library, both LBRYINCDIR and
# LBRYLIB must be defined and not empty.

  if test x"$$2INCDIR" != x || test x"$$2LIB" != x; then
    if test x"$$2INCDIR" = x || test x"$$2LIB" = x; then
      AC_MSG_ERROR([You need to specify both an include directory and link flags to use library $1. Use --with-m4_tolower($1)-incdir of environment variable $$2INCDIR to specify the include directory. Use --with-m4_tolower($1)-lib or environment variable $$2LIB to specify link flags.])
    fi
    m4_tolower(coin_has_$2)=true
    AC_MSG_RESULT(yes)
  else
    m4_tolower(coin_has_$2)=false
    AC_MSG_RESULT(no)
  fi

# If we have instructions for use, consider header and link checks.

  if test $m4_tolower(coin_has_$2) = true; then

# If argument 3 (file) is given, check for the file. Typically this will be a
# header file, but that's not assumed.

    m4_ifval([$3],
        [AC_COIN_CHECK_FILE([$$2INCDIR/$3],[],
	     [AC_MSG_ERROR([Cannot find file $3 in $$2INCDIR])])])

# Now see if we can link the function. There are arguments for and against
# assuming argument 3 is a header file declaring the function. A correct
# function declaration is the main argument in favour. Having to cope with
# possible dependencies or other oddities are the main arguments against.
# Force the use of C as the best single choice amongst C++, C, and Fortran.
# Obviously, this has limits.

    m4_ifvaln([$4],
        [if test x"$m4_tolower($1)_libcheck" != xno; then
	   coin_save_LIBS="$LIBS"
	   LIBS="$$2LIB $ADDLIBS"
	   coin_$2_link=no
	   AC_LANG_PUSH(C)
	   for fnm in $4 ; do
	     AC_MSG_CHECKING([whether symbol $fnm is available with $2])
	     AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],[[$fnm()]])],
		 [AC_MSG_RESULT(yes)
		  coin_$2_link=yes
		  break],
		 [AC_MSG_RESULT(no)])
	   done
	   AC_LANG_POP(C)
	   if test x"$coin_$2_link" = xyes ; then
	     LIBS="$coin_save_LIBS"
	   else
	     AC_MSG_ERROR([Cannot find symbol(s) $4 with $2])
	   fi
	 fi])

# If we make it this far, we've verified the file and linked the function. Add
# the necessary link flags to ADDLIBS and define the preprocessor symbol
# COIN_HAS_LBRY.

    ADDLIBS="$$2LIB $ADDLIBS"
    AC_DEFINE(COIN_HAS_$2,[1],[Define to 1 if the $1 package is available])
  fi

# Arrange for configure to substitute LBRYINCDIR and LBRYLIB and create the
# automake conditional. These actions must occur unconditionally.

  AC_SUBST($2INCDIR)
  AC_SUBST($2LIB)
  AM_CONDITIONAL(COIN_HAS_$2, test $m4_tolower(coin_has_$2) = true)
]) #AC_COIN_HAS_USER_LIBRARY 

###########################################################################
#                               COIN_HAS_ASL                              #
###########################################################################

# This macro checks if the user has provide arguments that say where
# the precompiled ASL files should be found (with the --with-asldir
# flag).  If this is not the case, we check if the ThirdParty/ASL
# directory has been configured, which indicates that the files will
# be in that directory and can be used.

AC_DEFUN([AC_COIN_HAS_ASL],
[coin_aslobjdir=../ThirdParty/ASL
coin_aslsrcdir=$srcdir/$coin_aslobjdir

# Determine the name of the ASL library
case "$CC" in
  cl* | */cl* | CL* | */CL* | icl* | */icl* | ICL* | */ICL*)
    ampllib=amplsolv.lib ;;
  *)
    ampllib=amplsolver.a ;;
esac

AC_ARG_WITH([asldir],
            AC_HELP_STRING([--with-asldir],
                           [specify path to AMPL solver directory (or BUILD for compilation, or "no" for disabling AMPL)]),
            [use_asldir=$withval], [use_asldir=])

if test "$use_asldir" = BUILD; then
  if test "$PACKAGE_NAME" != ThirdPartyASL; then
    # If we are configuring ThirdParty/ASL, don't check
    if test -r $coin_aslobjdir/.MakeOk; then :; else
      AC_MSG_ERROR([option \"BUILD\" specified for asldir, but directory is not configure (sources missing?)])
    fi
  fi
elif test -z "$use_asldir"; then
 # try to find sources - if not given don't compile
  if test "$PACKAGE_NAME" != ThirdPartyASL; then
    if test -r $coin_aslobjdir/.MakeOk; then
      use_asldir=BUILD
    else
      use_asldir=no
    fi
  else
    use_asldir=no
  fi
elif test "$use_asldir" != "no"; then
  AC_COIN_CHECK_FILE([$use_asldir/$ampllib],[],
                     [AC_MSG_ERROR([ASL directory \"$use_asldir\" specified, but library missing])])
  AC_COIN_CHECK_FILE([$use_asldir/asl.h],[],
                     [AC_MSG_ERROR([ASL directory \"$use_asldir\" specified, but header files are missing])])
  use_asldir=`cd $use_asldir; pwd`
  case $build in
    *-cygwin*) use_asldir=`cygpath -w $use_asldir | sed -e sX\\\\\\\\X/Xg` ;;
  esac
fi

# Variable containing ASL library (including full path)
AC_SUBST(ASLLIB)
# Variable containing flags for including ASL header files
AC_SUBST(ASL_CPPFLAGS)

if test "$use_asldir" = BUILD; then
  coin_aslobjdir=`cd $coin_aslobjdir; pwd`
  ASLLIB=`$CYGPATH_W $coin_aslobjdir/$ampllib | sed -e sX\\\\\\\\X/Xg`
  coin_aslsrcdir=`cd $coin_aslsrcdir; pwd`
  ASL_CPPFLAGS="-I"`$CYGPATH_W $coin_aslobjdir | sed -e sX\\\\\\\\X/Xg`" -I"`$CYGPATH_W $coin_aslsrcdir/solvers | sed -e sX\\\\\\\\X/Xg`
elif test "$use_asldir" != no; then
  ASLLIB=`$CYGPATH_W $use_asldir/$ampllib | sed -e sX\\\\\\\\X/Xg`
  ASL_CPPFLAGS="-I"`$CYGPATH_W $use_asldir | sed -e sX\\\\\\\\X/Xg`
fi

if test "$use_asldir" != no; then
  AC_LANG_PUSH(C)
  AC_CHECK_LIB(dl,[dlopen],[ASLLIB="$ASLLIB -ldl"],[])
  coin_has_asl=yes
  AC_DEFINE([COIN_HAS_ASL],[1],
            [If defined, the Ampl Solver Library is available.])
  AC_LANG_POP(C)
else
  coin_has_asl=no
fi
AM_CONDITIONAL(COIN_HAS_ASL, test $coin_has_asl = yes)
]) # AC_COIN_HAS_ASL

###########################################################################
#                            COIN_TRY_FLINK                               #
###########################################################################

# Auxilliary macro to test if a Fortran function name can be linked,
# given the current settings of LIBS.  We determine from the context, what
# the currently active programming language is, and cast the name accordingly.
# The first argument is the name of the function/subroutine, in small letters,
# the second argument are the actions taken when the test works, and the
# third argument are the actions taken if the test fails.

AC_DEFUN([AC_COIN_TRY_FLINK],
[case $ac_ext in
  f)
    AC_TRY_LINK(,[      call $1],[$2],[$3])
    ;;
  c)
    AC_F77_FUNC($1,cfunc$1)
    if test x"$coin_need_flibs" = xyes; then
      flink_try=no;
    else
      AC_TRY_LINK([void $cfunc$1();],[$cfunc$1()],
                  [flink_try=yes],[flink_try=no])
    fi
    if test $flink_try = yes; then
      $2
    else
      if test x"$FLIBS" != x; then
        flink_save_libs="$LIBS"
        LIBS="$LIBS $FLIBS"
        AC_TRY_LINK([void $cfunc$1();],[$cfunc$1()],
                    [LIBS="$flink_save_libs"
                     $2
                     coin_need_flibs=yes],
                    [LIBS="$flink_save_libs"
                     $3])
      else
        $3
      fi
    fi
    ;;
  cc|cpp)
    AC_F77_FUNC($1,cfunc$1)
    if test x"$coin_need_flibs" = xyes; then
      flink_try=no;
    else
      AC_TRY_LINK([extern "C" {void $cfunc$1();}],[$cfunc$1()],
                  [flink_try=yes],[flink_try=no])
    fi
    if test $flink_try = yes; then
      $2
    else
      if test x"$FLIBS" != x; then
        flink_save_libs="$LIBS"
        LIBS="$LIBS $FLIBS"
        AC_TRY_LINK([extern "C" {void $cfunc$1();}],[$cfunc$1()],
                    [LIBS="$flink_save_libs"
                     $2
                     coin_need_flibs=yes],
                    [LIBS="$flink_save_libs"
                     $3])
      else
        $3
      fi
    fi
    ;;
esac
]) # AC_COIN_TRY_FLINK

###########################################################################
#                             COIN_HAS_BLAS                               #
###########################################################################

# This macro checks for a library containing the BLAS library.  It
# tries standard libraries, and if none is found to be working, it
# checks whether the BLAS ThirdParty/Blas directory has been configured.
# It adds to ADDLIBS any flags required to link with an externally provided
# BLAS.  It defines the makefile conditional and preprocessor macro
# COIN_HAS_BLAS, if blas is available, and it defines the makefile conditional
# COIN_BUILD_BLAS, if blas is compiled within COIN.

AC_DEFUN([AC_COIN_HAS_BLAS],
[if test "$PACKAGE_NAME" = ThirdPartyBlas || test "$PACKAGE_NAME" = ThirdPartyLapack || test "$PACKAGE_NAME" = ThirdPartyMumps; then
  coin_blasobjdir=../Blas
else
  coin_blasobjdir=../ThirdParty/Blas
fi
coin_blassrcdir=$srcdir/$coin_blasobjdir

AC_ARG_WITH([blas],
            AC_HELP_STRING([--with-blas],
                           [specify BLAS library (or BUILD for compilation)]),
            [use_blas=$withval], [use_blas=])

MAKEOKFILE=.MakeOk
# Check if user supplied option makes sense
if test x"$use_blas" != x; then
  if test "$use_blas" = "BUILD"; then
    # Don't check for course code if this is executed in ThirdParty/Blas
    if test "$PACKAGE_NAME" != ThirdPartyBlas; then
      if test -r $coin_blasobjdir/.MakeOk; then :; else
        AC_MSG_ERROR([option \"BUILD\" specified for Blas, but $coin_blasobjdir directory is not properly configured])
      fi
    fi
  elif test "$use_blas" != no ; then
    AC_MSG_CHECKING([whether user supplied BLASLIB=\"$use_blas\" works])
    LIBS="$use_blas $LIBS"
    ADDLIBS="$use_blas $ADDLIBS"
    AC_COIN_TRY_FLINK([daxpy],
                      [AC_MSG_RESULT([yes])],
                      [AC_MSG_RESULT([no])
                       AC_MSG_ERROR([user supplied BLAS library \"$use_blas\" does not work])])
  fi
else
# Try to autodetect the library for blas based on build system
  #AC_MSG_CHECKING([default locations for BLAS])
  skip_lblas_check=no
  case $build in
    *-sgi-*) 
      SAVE_LIBS="$LIBS"
      AC_MSG_CHECKING([whether -lcomplib.sgimath has BLAS])
      LIBS="-lcomplib.sgimath $LIBS"
      AC_COIN_TRY_FLINK([daxpy],
                        [AC_MSG_RESULT([yes])
                         use_blas=-lcomplib.sgimath;
                         ADDLIBS="-lcomplib.sgimath $ADDLIBS"],
                        [AC_MSG_RESULT([no])
                         SAVE_LIBS="$LIBS"])
      ;;

# Ideally, we'd use -library=sunperf, but it's an imperfect world. Studio
# cc doesn't recognise -library, it wants -xlic_lib. Studio 12 CC doesn't
# recognise -xlic_lib. Libtool doesn't like -xlic_lib anyway. Sun claims
# that CC and cc will understand -library in Studio 13. The main extra
# function of -xlic_lib and -library is to arrange for the Fortran run-time
# libraries to be linked for C++ and C. We can arrange that explicitly.
    *-*-solaris*)
      SAVE_LIBS="$LIBS"
      AC_MSG_CHECKING([for BLAS in libsunperf])
      LIBS="-lsunperf $FLIBS $LIBS"
      AC_COIN_TRY_FLINK([daxpy],
                        [AC_MSG_RESULT([yes])
                         use_blas='-lsunperf'
                         ADDLIBS="-lsunperf $ADDLIBS"
			 coin_need_flibs=yes],
                        [AC_MSG_RESULT([no])
                         LIBS="$SAVE_LIBS"])
      ;;
    *-cygwin* | *-mingw*)
# On cygwin, consider -lblas only if doscompile is disabled. The prebuilt
# library will want to link with cygwin, hence won't run standalone in DOS.
      if test "$enable_doscompile" = mingw; then
	skip_lblas_check=yes
      fi
      case "$CC" in
        cl* | */cl* | CL* | */CL* | icl* | */icl* | ICL* | */ICL*)
          SAVE_LIBS="$LIBS"
          AC_MSG_CHECKING([for BLAS in MKL])
          LIBS="mkl_intel_c.lib mkl_sequential.lib mkl_core.lib $LIBS"
          AC_COIN_TRY_FLINK([daxpy],
                            [AC_MSG_RESULT([yes])
                             use_blas='mkl_intel_c.lib mkl_sequential.lib mkl_core.lib'
                             ADDLIBS="mkl_intel_c.lib mkl_sequential.lib mkl_core.lib $ADDLIBS"],
                            [AC_MSG_RESULT([no])
                             LIBS="$SAVE_LIBS"])
          ;;
      esac
      ;;
  esac

  if test -z "$use_blas" && test $skip_lblas_check = no; then
    SAVE_LIBS="$LIBS"
    AC_MSG_CHECKING([whether -lblas has BLAS])
    LIBS="-lblas $LIBS"
    AC_COIN_TRY_FLINK([daxpy],
		      [AC_MSG_RESULT([yes])
		       ADDLIBS="-lblas $ADDLIBS"
		       use_blas='-lblas'],
		      [AC_MSG_RESULT([no])
	               LIBS="$SAVE_LIBS"])
  fi

# If we have no other ideas, consider building BLAS.
  if test -z "$use_blas"; then
    if test "$PACKAGE_NAME" != ThirdPartyBlas; then
      if test -r $coin_blasobjdir/.MakeOk; then
        use_blas=BUILD
      fi
    fi
  fi
fi

if test "$use_blas" = BUILD; then
  coin_need_flibs=yes
fi

AM_CONDITIONAL([COIN_HAS_BLAS],[test x"$use_blas" != x])
AM_CONDITIONAL([COIN_BUILD_BLAS],[test "$use_blas" = BUILD])

if test x"$use_blas" = x || test "$use_blas" = no; then
  coin_has_blas=no
else
  coin_has_blas=yes
  AC_DEFINE([COIN_HAS_BLAS],[1],
            [If defined, the BLAS Library is available.])
fi
]) # AC_COIN_HAS_BLAS

###########################################################################
#                            COIN_HAS_LAPACK                              #
###########################################################################

# This macro checks for a library containing the LAPACK library.  It
# tries standard libraries, and if none is found to be working, it
# checks whether the LAPACK ThirdParty/Lapack directory has been
# configured.  It adds to ADDLIBS any flags required to link with an
# externally provided LAPACK.  It defines the makefile conditional and
# preprocessor macro COIN_HAS_LAPACK, if lapack is available, and it
# defines the makefile conditional COIN_BUILD_LAPACK, if lapack is
# compiled within COIN.

AC_DEFUN([AC_COIN_HAS_LAPACK],
[coin_lapackobjdir=../ThirdParty/Lapack
coin_lapacksrcdir=$srcdir/$coin_lapackobjdir

AC_ARG_WITH([lapack],
            AC_HELP_STRING([--with-lapack],
                           [specify LAPACK library (or BUILD for compilation)]),
            [use_lapack=$withval], [use_lapack=])

# Check if user supplied option makes sense
if test x"$use_lapack" != x; then
  if test "$use_lapack" = "BUILD"; then
    # Don't check for course code if this is executed in ThirdParty/Lapack
    if test "$PACKAGE_NAME" != ThirdPartyLapack; then
      if test -r $coin_lapackobjdir/.MakeOk; then :; else
        AC_MSG_ERROR([option \"BUILD\" specified for LAPACK, but $coin_lapackobjdir directory is not configured])
      fi
    fi
  elif test "$use_lapack" != no; then
    AC_MSG_CHECKING([whether user supplied LAPACKLIB=\"$use_lapack\" works])
    LIBS="$use_lapack $LIBS"
    ADDLIBS="$use_lapack $ADDLIBS"
    AC_COIN_TRY_FLINK([dsyev],
                      [AC_MSG_RESULT([yes])],
                      [AC_MSG_RESULT([no])
                       AC_MSG_ERROR([user supplied LAPACK library \"$use_lapack\" does not work])])
  fi
else
  if test x$coin_has_blas = xyes; then
    # First try to see if LAPACK is already available with BLAS library
    AC_MSG_CHECKING([whether LAPACK is already available with BLAS library])
    AC_COIN_TRY_FLINK([dsyev],
                      [AC_MSG_RESULT([yes]); use_lapack=ok],
                      [AC_MSG_RESULT([no])])
  fi
  skip_llapack_check=no
  if test -z "$use_lapack"; then
    # Try to autodetect the library for lapack based on build system
    case $build in
      *-sgi-*) 
        SAVE_LIBS="$LIBS"
        AC_MSG_CHECKING([whether -lcomplib.sgimath has LAPACK])
        LIBS="-lcomplib.sgimath $LIBS"
        AC_COIN_TRY_FLINK([dsyev],
                          [AC_MSG_RESULT([yes])
                           use_lapack=-lcomplib.sgimath;
                           ADDLIBS="-lcomplib.sgimath $ADDLIBS"],
                          [AC_MSG_RESULT([no])
                           SAVE_LIBS="$LIBS"])
        ;;

# See comments in COIN_HAS_BLAS.
      *-*-solaris*)
      SAVE_LIBS="$LIBS"
      AC_MSG_CHECKING([for LAPACK in libsunperf])
      LIBS="-lsunperf $FLIBS $LIBS"
      AC_COIN_TRY_FLINK([dsyev],
                        [AC_MSG_RESULT([yes])
                         use_blas='-lsunperf'
                         ADDLIBS="-lsunperf $ADDLIBS"
			 coin_need_flibs=yes],
                        [AC_MSG_RESULT([no])
                         LIBS="$SAVE_LIBS"])
        ;;
# On cygwin, do this check only if doscompile is disabled. The prebuilt library
# will want to link with cygwin, hence won't run standalone in DOS.
      *-cygwin*)
	if test "$enable_doscompile" = mingw; then
	  skip_llapack_check=yes
	fi
	;;
    esac
  fi

  if test -z "$use_lapack" && test $skip_llapack_check = no; then
    SAVE_LIBS="$LIBS"
    AC_MSG_CHECKING([whether -llapack has LAPACK])
    LIBS="-llapack $LIBS"
    AC_COIN_TRY_FLINK([dsyev],
		      [AC_MSG_RESULT([yes])
		       ADDLIBS="-llapack $ADDLIBS"
		       use_lapack='-llapack'],
		      [AC_MSG_RESULT([no])
		       LIBS="$SAVE_LIBS"])
  fi

# If we have no other ideas, consider building LAPACK.
  if test -z "$use_lapack"; then
    if test "$PACKAGE_NAME" != ThirdPartyLapack; then
      if test -r $coin_lapackobjdir/.MakeOk; then
        use_lapack=BUILD
      fi
    fi
  fi
fi

if test "$use_lapack" = BUILD; then
  coin_need_flibs=yes
fi

AM_CONDITIONAL([COIN_HAS_LAPACK],[test x"$use_lapack" != x])
AM_CONDITIONAL([COIN_BUILD_LAPACK],[test "$use_lapack" = BUILD])

if test x"$use_lapack" = x || test "$use_lapack" = no; then
  coin_has_lapack=no
else
  coin_has_lapack=yes
  AC_DEFINE([COIN_HAS_LAPACK],[1],
            [If defined, the LAPACK Library is available.])
fi
]) # AC_COIN_HAS_LAPACK

###########################################################################
#                            COIN_HAS_MUMPS                               #
###########################################################################

# This macro checks for a library containing the MUMPS library.  It
# checks if the user has provided an argument for the MUMPS library,
# and if not, it checks whether the MUMPS ThirdParty/Mumps directory has
# been configured.  It adds to ADDLIBS any flags required to link with
# an externally provided MUMPS.  It defines the makefile conditional
# and preprocessor macro COIN_HAS_MUMPS, if MUMPS is available, and it
# defines the makefile conditional COIN_BUILD_MUMPS, if MUMPS is
# compiled within COIN.

AC_DEFUN([AC_COIN_HAS_MUMPS],
[
if test "$PACKAGE_NAME" = ThirdPartyMumps; then
  coin_mumpsobjdir=../Mumps
else
  coin_mumpsobjdir=../ThirdParty/Mumps
fi
coin_mumpssrcdir=$abs_source_dir/$coin_mumpsobjdir/MUMPS

MAKEOKFILE=.MakeOk

#check if user provides a MUMPS library (that works)
AC_LANG_PUSH(C)
SAVE_ADDLIBS="$ADDLIBS"
ADDLIBS="$ADDLIBS $FLIBS"
AC_COIN_HAS_USER_LIBRARY(mumps, MUMPS, dmumps_c.h, dmumps_c)
ADDLIBS="$SAVE_ADDLIBS"
AC_LANG_POP(C)

if test "$coin_has_mumps" = "true"; then  # user provided mumps library
  use_mumps=yes
  coin_has_mumps=yes

  MUMPS_INCFLAGS="-I\`\$(CYGPATH_W) $MUMPSINCDIR\`"
  ADDLIBS="$MUMPSLIB $ADDLIBS"

else # no user provided library, so we try to build our own
  use_mumps=BUILD

  # Check if the MUMPS' ThirdParty project has been configured
  if test "$PACKAGE_NAME" != ThirdPartyMumps; then
    if test -r $coin_mumpsobjdir/.MakeOk; then
      use_mumps=BUILD
      # Mumps needs pthreads
      AC_LANG_PUSH(C)
      save_LIBS="$LIBS"
      LIBS="$LIBS $FLIBS"
      AC_CHECK_LIB([pthread],[pthread_create],[LIBS="-lpthread $save_LIBS"; ADDLIBS="-lpthread $ADDLIBS"],[LIBS="save_LIBS"])
      AC_LANG_POP(C)
  
      MUMPS_INCFLAGS="-I\`\$(CYGPATH_W) $coin_mumpssrcdir/libseq\` -I\`\$(CYGPATH_W) $coin_mumpssrcdir/include\`"
    else
      use_mumps=
    fi
  fi

  # if a user provided library is used, then COIN_HAS_USER_LIBRARY takes care of the COIN_HAS_MUMPS conditional and preprocessor symbol
  AM_CONDITIONAL([COIN_HAS_MUMPS],[test x"$use_mumps" != x])
  if test x"$use_mumps" != x; then
    AC_DEFINE([COIN_HAS_MUMPS],[1],[If defined, the MUMPS Library is available.])
    coin_has_mumps=yes
  else
    coin_has_mumps=no
  fi
  AC_MSG_CHECKING([whether MUMPS is available])
  AC_MSG_RESULT([$coin_has_mumps])
fi

if test x"$use_mumps" != x; then
  # we need the Fortran runtime libraries if we want to link with C/C++
  coin_need_flibs=yes

  AC_SUBST(MUMPS_INCFLAGS)
fi

AM_CONDITIONAL([COIN_BUILD_MUMPS],[test "$use_mumps" = BUILD])

]) # AC_COIN_HAS_MUMPS

###########################################################################
#                            COIN_HAS_METIS                               #
###########################################################################

# This macro checks for a library containing the METIS library.  It
# checks if the user has provided an argument for the METIS library,
# and if not, it checks whether the METIS ThirdParty/Metis directory has
# been configured.  It adds to ADDLIBS any flags required to link with
# an externally provided METIS.  It defines the makefile conditional
# and preprocessor macro COIN_HAS_METIS, if METIS is available, and it
# defines the makefile conditional COIN_BUILD_METIS, if METIS is
# compiled within COIN.

AC_DEFUN([AC_COIN_HAS_METIS],
[
case "$PACKAGE_NAME" in
  ThirdParty*)
    coin_metisobjdir=../Metis
    ;;
  *)
    coin_metisobjdir=../ThirdParty/Metis
    ;;
esac

MAKEOKFILE=.MakeOk

#check if user provides a METIS library (that works)
AC_LANG_PUSH(C)
AC_ARG_WITH(metis,
   AS_HELP_STRING([--with-metis], [specify flags to link with METIS library]),
   [METISLIB="$withval"; coin_has_metis=true], [coin_has_metis=no])

if test $coin_has_metis = true; then
  coin_save_LIBS="$LIBS"
  LIBS="$METISLIB $ADDLIBS"
  AC_MSG_CHECKING([whether symbol metis_nodend is available with ])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],[[metis_nodend()]])],
     [AC_MSG_RESULT(yes)],
     [AC_MSG_RESULT(no)
      AC_MSG_ERROR([User-supplied METIS library does not work])])
  LIBS="$coin_save_LIBS"
fi
AC_LANG_POP(C)

if test "$coin_has_metis" = "true"; then  # user provided metis library
  use_metis=yes
  coin_has_metis=yes

  ADDLIBS="$METISLIB $ADDLIBS"

else # no user provided library, so we try to build our own
  use_metis=BUILD

  # Check if the METIS' ThirdParty project has been configured
  if test "$PACKAGE_NAME" != ThirdPartyMetis; then
    if test -r $coin_metisobjdir/.MakeOk; then
      use_metis=BUILD
    else
      use_metis=
    fi
  fi

  if test x"$use_metis" != x; then
    coin_has_metis=yes
  else
    coin_has_metis=no
  fi
  AC_MSG_CHECKING([whether METIS is available])
  AC_MSG_RESULT([$coin_has_metis])
fi

AM_CONDITIONAL([COIN_HAS_METIS],[test $coin_has_metis = yes])
if test $coin_has_metis = yes; then
  AC_DEFINE([COIN_HAS_METIS],[1],[If defined, the METIS library is available.])
fi

]) # AC_COIN_HAS_METIS


###########################################################################
#                             COIN_HAS_GLPK                               #
###########################################################################
#
# This macro checks for the GLPK package. GLPK provides two capabilities,
# an LP and MIP solver (GLPK) and the GNU Mathprog modelling language (GMPL).
# The macro checks for either Glpk or Gmpl, according to the value specified as
# the parameter. Use one of Glpk or Gmpl. Use *exactly* these strings, eh?
#
#
# The macro first uses COIN_HAS_USER_LIBRARY to see if the user has specified
# a preexisting library (this allows the use of any glpk version, if the user
# is fussy). The macro then checks for ThirdParty/Glpk.
#
# This macro will define the following variables for Glpk:
#   coin_has_glpk	true or false
#   GLPKLIB		link flags for GLPK (if user supplied)
#   GLPKINCDIR		location of glpk include files
#   COIN_HAS_GLPK	Preprocessor symbol, defined to 1
#   COIN_HAS_GLPK	Automake conditional
#   COIN_BUILD_GLPK	Automake conditional, defined only if Glpk is to be
#			built in ThirdParty/Glpk
#
# With the exception of COIN_BUILD_GLPK, an identical set of variables is
# defined for Gmpl.

AC_DEFUN([AC_COIN_HAS_GLPK],
[
if test "$PACKAGE_NAME" = ThirdPartyGlpk; then
  coin_glpkobjdir=../Glpk
else
  coin_glpkobjdir=../ThirdParty/Glpk
fi
coin_glpksrcdir=$abs_source_dir/$coin_glpkobjdir

use_thirdpartyglpk=no

# Check for the requested component. If the user specified an external glpk
# library don't force a ThirdParty build, let the error propagate.

m4_if([$1],[Glpk],
[AC_COIN_HAS_USER_LIBRARY([Glpk],[GLPK],[glpk.h],
    [_glp_lpx_simplex glp_lpx_simplex])
 if test x"$coin_has_glpk" = xfalse && test x"$GLPKLIB" = x ; then
   use_thirdpartyglpk=try
 fi])
m4_if([$1],[Gmpl],
[AC_COIN_HAS_USER_LIBRARY([Gmpl],[GMPL],[glpmpl.h],
    [_glp_mpl_initialize glp_mpl_initialize])
 if test x"$coin_has_gmpl" = xfalse && test x"$GMPLLIB" = x ; then
   use_thirdpartyglpk=try
 fi])

# If the user has supplied an external library, use it. Otherwise, consider
# a build in ThirdParty/Glpk. If we build, assume we get both glpk and gmpl.

if test x"$use_thirdpartyglpk" = xtry ; then
  MAKEOKFILE=.MakeOk
  # Check if the Glpk's ThirdParty project has been configured
  if test "$PACKAGE_NAME" != ThirdPartyGlpk; then
    if test -r $coin_glpkobjdir/.MakeOk; then
      use_thirdpartyglpk=build
    else
      use_thirdpartyglpk=no
    fi
  else
    use_thirdpartyglpk=build
  fi

# If we're building, set the library and include directory variables, create a
# preprocessor symbol, define a variable that says we're using glpk/gmpl, and
# another to indicate a link check is a bad idea (hard to do before the library
# exists).

  if test x"$use_thirdpartyglpk" = xbuild ; then
    m4_toupper($1INCDIR)="$coin_glpksrcdir/glpk/include"
    AC_SUBST(m4_toupper($1INCDIR))
    AC_DEFINE(m4_toupper(COIN_HAS_$1),[1],
	      [Define to 1 if $1 package is available])
    m4_tolower(coin_has_$1)=true
    m4_tolower($1_libcheck)=no
    AC_MSG_NOTICE([Using $1 in ThirdParty/Glpk])
  fi
fi

# Define the necessary automake conditionals.

AM_CONDITIONAL(m4_toupper(COIN_HAS_$1),
	       [test x"$m4_tolower(coin_has_$1)" = xtrue])
AM_CONDITIONAL([COIN_BUILD_GLPK],[test x"$use_thirdpartyglpk" = xbuild])

]) # AC_COIN_HAS_GLPK
