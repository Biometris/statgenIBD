#include <iostream>

#include "mblexcept.h"
#include "excepthnd.h"
#include <Rcpp.h>

using std::endl;
using namespace Rcpp;

int mbl::exception_handler(int (*F)())
{
	try
	{
		F();
	}
	catch (mblib_error& e)
	{
		Rcerr << "mblib_error: " << e.what() << endl;
	}
	catch (std::exception& e)
	{
		Rcerr << "std exception: " << e.what() << endl;
	}
	catch (...)
	{
		Rcerr << "unknown exception: " << endl;
	}
	return 0;
}

int mbl::exception_handler(int (*F)(const Args& args),
                           const Args& args)
{
	try
	{
		F(args);
	}
	catch (mblib_error& e)
	{
		Rcerr << "mblib_error: " << e.what() << endl;
	}
	catch (std::exception& e)
    {
	    Rcerr << "std exception: " << e.what() << endl;
	}
	catch (...)
	{
		Rcerr << "unknown exception: " << endl;
	}
	return 0;
}
