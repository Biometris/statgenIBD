/*!
\file
\brief  Exception handling in mbl-library
\author Martin Boer, Biometris
\date   1998-2007
*/
#ifndef MBLIB_EXCEPTION_HEADER
#define MBLIB_EXCEPTION_HEADER

#include <string>
#include <stdexcept>
#include <iostream>

namespace mbl
{

class mblib_error : public std::runtime_error
{
public:
	mblib_error(const std::string& what_arg) : std::runtime_error(what_arg) {}
	virtual ~mblib_error() throw() {;}
};

class mblib_file_error : public mblib_error
{
public:
	mblib_file_error(const std::string& filename,int line_nr, const std::string& what_arg);
	virtual ~mblib_file_error() throw() {;}
};

typedef mblib_error		 MQMlib_error;
typedef mblib_file_error MQMlib_file_error;

}

#endif

