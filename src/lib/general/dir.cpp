#include "dir.h"
#include "mblexcept.h"

#ifndef LINUX
#include <direct.h>
#else
#include <sys/stat.h> // for mkdir
#include <sys/types.h> // for mkdir
#include <float.h> // for DBL_MIN
#endif

using namespace mbl;
using namespace std;

//const int _MAX_PATH = 10000;

string mbl::get_current_dir()
{
	string retBuf;
#ifdef LINUX
	retBuf = getenv("PWD");
#else
	char buffer[_MAX_PATH];
	_getcwd(buffer,_MAX_PATH);
	retBuf = string(buffer);
#endif
	return retBuf;
}

void mbl::ChangeDir(const char * dir)
{
#ifdef LINUX
    if (chdir(dir) == 0 )    return;
    if (mkdir(dir,0700) < 0) throw mblib_error("Cannot make directory");  
    if (chdir(dir) == -1)    throw mblib_error("Cannot change to directory");
#else
	if (_chdir(dir) == 0)  return;   
	if (_mkdir(dir) == -1) throw mblib_error("Cannot make directory");  
	if (_chdir(dir) == -1) throw mblib_error("Cannot change to directory");
#endif
}

void mbl::ChangeDir(std::string dir) { ChangeDir(dir.c_str()); }
