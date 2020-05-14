#include <iostream>

#include "mblexcept.h"
#include "excepthnd.h"

using std::cerr;
using std::endl;

int mbl::exception_handler(int (*F)())
{
	// solves problem with comma separated integers in VC8 (Studio 2005)
	std::locale::global(std::locale("C")); 
	try
	{
		F();
	}
	catch (mblib_error& e)
	{
		cerr << "mblib_error: " << e.what() << endl;
	}
	catch (std::exception& e)
	{
		cerr << "std exception: " << e.what() << endl;
	}
	catch (...)
	{
		cerr << "unknown exception: " << endl;
	}
	return 0;
}

int mbl::exception_handler(int (*F)(const Args& args), const Args& args)
{
	// solves problem with comma separated integers in VC8 (Studio 2005)
	std::locale::global(std::locale("C"));
	try
	{
		F(args);
	}
	catch (mblib_error& e)
	{
		cerr << "mblib_error: " << e.what() << endl;
	}
	catch (std::exception& e)
    {
	    cerr << "std exception: " << e.what() << endl;
	}
	catch (...)
	{
		cerr << "unknown exception: " << endl;
	}
	return 0;
}
