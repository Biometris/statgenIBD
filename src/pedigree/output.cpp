// Martin Boer, Biometris
#include <fstream>

#include "output.h"
#include "output_bin.h"
#include "output_tab.h"
#include "output_flx.h"

using namespace std;
using namespace mbl;

void outp_header(ostream& outp, 
					 const LinkageMap& markermap,
					 const vector<IndProp>& pop,
					 const std::set<string>& sel_ind,
					 OutputType type)
{
	if (type == BIN) 
	{	
		outp_bin_header(outp,markermap,pop);
		return;
	}
	if (type == TAB) 
	{
		outp_tab_header(outp,markermap,pop,sel_ind);
		return;
	}
	if (type == FLX) 
	{
		outp_flx_header(outp,markermap,pop,sel_ind);
		return;
	}
}

void outp_fam(ostream& outp,const Pedigree& ped_fam, 
					const matrix3D<double>& Z,
					const LinkageMap& markermap,
					OutputType type)
{
	if (type == BIN) return outp_bin_fam(outp,ped_fam,Z);
	if (type == TAB) return outp_tab_fam(outp,ped_fam,Z,markermap);
}

void outp_ped_ftp(ostream& outp,const vector<int>& ndx, 
					    const IBDped& IBD,const Pedigree& ped,
					    const LinkageMap& markermap,
						OutputType type)
{
	if (type == BIN) return outp_bin_ped_ftp(outp,ndx,IBD,ped);
	if (type == TAB) return outp_tab_ped_ftp(outp,ndx,IBD,ped,markermap);
}

void outp_ped_coa(ostream& outp,const vector<int>& ndx, 
					const IBDped& IBD,
					const LinkageMap& markermap,
					OutputType type)
{
	if (type == BIN) return outp_bin_ped_coa(outp,ndx,IBD);
	if (type == TAB) return outp_tab_ped_coa(outp,ndx,IBD,markermap);
	if (type == FLX) return outp_flx_ped_coa(outp,ndx,IBD,markermap);
}

OutputType open_output(ofstream& outp, const Args& Argu,
					   const string& output, bool coa_parents_fam)
{
	string outputfile;
	OutputType outp_type = TAB;
	if (Argu.GetOption("bin"))
	{
		if (!coa_parents_fam) 
			throw mblib_error("It is not possible to use -bin and -coa together");
		outp_type = BIN;
		outputfile = output + ".bin";
		outp.open(outputfile.c_str(),ios::binary);
	}
	else 
	{
	    int prec=3;
	    Argu.GetOption("prec",prec);
		if (Argu.GetOption("flx"))
		{
			outp_type = FLX;
			outputfile = output + ".ldf";
		}
		else
			outputfile = output + ".txt";
		OpenFile(outp,outputfile);
		outp << setprecision(prec);
	}
	return outp_type;
}
