/*!
\file
\brief  Exception handler
\author Martin Boer, Biometris
\date   1998-2007
*/
#ifndef EXCEPTION_HANDLER_HEADER
#define EXCEPTION_HANDLER_HEADER

#include "Args.h"

namespace mbl {

int exception_handler(int (*F)());
int exception_handler(int (*F)(const Args& args), const Args& args);

}

#endif

