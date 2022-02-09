#include "dir.h"
#include "ibdexcept.h"

//#ifndef LINUX
//#include <direct.h>
//#else
#include <sys/stat.h> // for mkdir
#include <sys/types.h> // for mkdir
#include <float.h> // for DBL_MIN
#endif

using namespace ibd;
using namespace std;

string ibd::get_current_dir()
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

void ibd::ChangeDir(const char * dir)
{
#ifdef LINUX
    if (chdir(dir) == 0 )    return;
    if (mkdir(dir,0700) < 0) throw ibd_error("Cannot make directory");  
    if (chdir(dir) == -1)    throw ibd_error("Cannot change to directory");
#else
	if (_chdir(dir) == 0)  return;   
	if (_mkdir(dir) == -1) throw ibd_error("Cannot make directory");  
	if (_chdir(dir) == -1) throw ibd_error("Cannot change to directory");
#endif
}

void ibd::ChangeDir(std::string dir) { ChangeDir(dir.c_str()); }
