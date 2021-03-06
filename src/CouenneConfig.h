/* Copyright (C) 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 *
 * Include file for the configuration of Couenne.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header file is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN-OR packages or third party code) are set
 * by the files config_*default.h. The project maintainer needs to remember
 * to update these files and choose reasonable defines.
 * A user can modify the default setting by editing the config_*default.h files.
 */

#ifndef __COUENNECONFIG_H__
#define __COUENNECONFIG_H__

#ifdef HAVE_CONFIG_H

#ifdef COUENNELIB_BUILD
#include "config.h"
#else
#include "config_couenne.h"
#endif

/* overwrite COUENNELIB_EXPORT from config.h when building Couenne
 * we want it to be __declspec(dllexport) when building a DLL on Windows
 * we want it to be __attribute__((__visibility__("default"))) when building with GCC,
 *   so user can compile with -fvisibility=hidden
 */
#ifdef COUENNELIB_BUILD
#ifdef DLL_EXPORT
#undef COUENNELIB_EXPORT
#define COUENNELIB_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) && __GNUC__ >= 4
#undef COUENNELIB_EXPORT
#define COUENNELIB_EXPORT __attribute__((__visibility__("default")))
#endif
#endif

#else /* HAVE_CONFIG_H */

#ifdef COUENNELIB_BUILD
#include "config_default.h"
#else
#include "config_couenne_default.h"
#endif

#endif /* HAVE_CONFIG_H */

#endif /*__HAVE_COUENNE_CONFIG_H__*/
