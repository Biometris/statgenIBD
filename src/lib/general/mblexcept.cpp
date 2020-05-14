#include <iostream>

#include "convert.h"
#include "mblexcept.h"

mbl::mblib_file_error::mblib_file_error(const std::string& filename, 
								int line_nr, const std::string& what_arg)
  : mblib_error("file: " + filename + ",line " + stringify(line_nr) + ": " + what_arg) {}      

