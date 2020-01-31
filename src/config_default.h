
/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_ipopt_default.h */
#ifndef COUENNELIB_EXPORT
#ifdef _WIN32
/* assuming we build an Ipopt DLL */
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
