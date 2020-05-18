// // Martin Boer, Biometris
// #include <fstream>
// #include <sstream>
//
// #include "misc.h"
// #include "read_ped.h"
// #include "main.h"
//
// using namespace mbl;
// using namespace std;
//
// map<string,string> get_family_ID(const string& filename)
// {
// 	string line;
// 	ifstream inp;
// 	OpenFile(inp,filename);
// 	int line_nr = 0;
// 	map<string,string> result;
// 	while (getline(inp,line))
// 	{
// 		line_nr++;
// 		if (line.empty()) continue;
// 		istringstream line_stream(line);
// 		string c1,c2,c3,c4,c5;
// 		line_stream >> c1 >> c2 >> c3 >> c4 >> c5;
// 		result[c2] = c5;
// 	}
// 	return result;
// }
//
// int convert_ped_file()
// {
// 	map<string,string> fam = get_family_ID("FlexQTL.dat");
//
// 	string line;
// 	ofstream outp;
// 	ifstream inp;
// 	OpenFile(inp,"Martin_ped.dat");
// 	OpenFile(outp,"Martin_ped_new.dat");
// 	int line_nr = 0;
// 	while (getline(inp,line))
// 	{
// 		line_nr++;
// 		if (line.empty()) continue;
// 		istringstream line_stream(line);
// 		string ID,type,P1,P2,fam_ID;
// 		line_stream >> ID >> type >> P1 >> P2;
// 		if (type == "RIL" || type == "INBFND")
// 			fam_ID = "*";
// 		else
// 			fam_ID = fam[ID];
// 		outp << ID << '\t' << fam_ID << '\t'
// 			 << type << '\t' << P1 << '\t' << P2 << '\n';
// 	}
// 	return 0;
// }
//
