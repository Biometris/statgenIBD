/*!
\file
\brief  Read LDF file format: IBD between founders
\author Martin Boer, Biometris
\date   2006-2010
*/
#ifndef READ_LDF_FILE_HEADER
#define READ_LDF_FILE_HEADER

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "matvec.h"
#include "Loc.h"

template<class T1> void read_columns(std::vector<T1>& x, int dim, std::ifstream& inp)
{
	std::string line;
	for (int i=0;i<dim;i++)
	{
		getline(inp,line);
		std::istringstream line_stream(line);
		T1 val_x;
		line_stream >> val_x;
		x.push_back(val_x);
	}
}

template<class T1,class T2> void read_columns(std::vector<T1>& x, std::vector<T2>& y, int dim,
											  std::ifstream& inp)
{
	std::string line;
	for (int i=0;i<dim;i++)
	{
		getline(inp,line);
		std::istringstream line_stream(line);
		T1 val_x;
		T2 val_y;
		line_stream >> val_x >> val_y;
		x.push_back(val_x);
		y.push_back(val_y);
	}
}

template<class T1,class T2,class T3> void read_columns(std::vector<T1>& x, std::vector<T2>& y,
							std::vector<T3>& z, int dim, std::ifstream& inp)
{
	std::string line;
	for (int i=0;i<dim;i++)
	{
		getline(inp,line);
		std::istringstream line_stream(line);
		T1 val_x;
		T2 val_y;
		T3 val_z;
		line_stream >> val_x >> val_y >> val_z;
		x.push_back(val_x);
		y.push_back(val_y);
		z.push_back(val_z);
	}
}

// Jan 21 2009: Check the correct format!!
void read_ldf_file(std::vector<int>& fndname,
                   LinkageMap& markermap,
                   mbl::matrix3D<double>& Q,
                   const std::string filename);

#endif

