
/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_couenne_default.h */
#ifndef COUENNELIB_EXPORT
#if defined(_WIN32) && defined(DLL_EXPORT)
#define COUENNELIB_EXPORT __declspec(dllexport)
#else
#define COUENNELIB_EXPORT
#endif
#endif

/* include the public project specific macros */
#include "config_couenne_default.h"

/***************************************************************************/
/*        HERE DEFINE THE PROJECT SPECIFIC PRIVATE MACROS                  */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/
