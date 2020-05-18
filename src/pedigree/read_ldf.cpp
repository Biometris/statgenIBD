// Martin Boer, Biometris

#include "read_ldf.h"
#include "manipulators.h"
#include "convert.h"

using namespace mbl;
using namespace std;

// Jan 21 2009: Check the correct format!!
void read_ldf_file(vector<int>& fndname,
				   LinkageMap& markermap,
				   matrix3D<double>& Q,
				   const string filename)
{
	ifstream inp;
	OpenFile(inp,filename);

	int nind,npos;
	string line,keyword;
	inp >> keyword >> nind;
	if (keyword != "nind")
		throw mblib_error("Error!");
	skip_rest_of_line(inp);
	skip_lines(inp,1);

	vector<int> nr;
	read_columns(nr,fndname,nind,inp);
	skip_lines(inp,1);

	inp >> keyword >> npos;
	if (keyword != "npos")
		throw mblib_error("Error!");
	skip_rest_of_line(inp);
	skip_lines(inp,1);

	vector<int> posnr,chr;
	vector<double> pos;
	read_columns(posnr,chr,pos,npos,inp);
	skip_lines(inp,2);

	Q = matrix3D<double>(npos,nind,nind);
	for (int m=0;m<npos;m++)
	{
		markermap.push_back(Locus(chr[m],pos[m],EVAL_POS));
		skip_lines(inp,1);
		inp >> Q[m];
	}
}
