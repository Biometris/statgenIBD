#include <iostream>

#include "PopName.h"
#include "util_genetics.h"
#include "convert.h"

using namespace ibd;
using namespace std;


/*
pop_base * init_pop(const std::string& poptype)
{
	// biparental crosses:
	if (poptype == "DH")                   // Doubled Haploid
		return new popDH();
	if (poptype.find("F") == 0)            // Fx and FxDH populations
		return init_F(poptype.substr(1));  
	if (poptype.find("BC") == 0)           // Backcross or backcross followed by selfing
	{
		if (poptype.find("BCS") == 0)      // Backcross followed by selfing
			return init_BCS(poptype.substr(3));
		else                               // Backcross
			return init_BC(poptype.substr(2));
	}

	// three-way crosses: 
	if (poptype == "C3")          // three-way cross, no selfing
		return new pop3wSx(0);
	if (poptype == "C3DH")        // three-way cross, no selfing, DH
		return new pop3wSxDH(0);
	if (poptype.find("C3S") == 0) // three-way cross, followed by selfing (and DH)
		return init_C3S(poptype.substr(3));

	// four-way crosses:
	if (poptype == "C4")          // four-way cross, no selfing
		return new pop4wSx(0); 
	if (poptype == "C4DH")        // four-way cross, no selfing, DH
		return new pop4wSxDH(0);
	if (poptype.find("C4S") == 0) // four-way cross, followed by selfing (and DH)
		return init_C4S(poptype.substr(3));

	throw mblib_error("unknown type " + poptype);
	return 0;
}
*/


Pop::Pop(string& name) : x(-1)
{
	if (name == "C3")   name = "C3S0";
	if (name == "C3DH") name = "C3S0DH";
	if (name == "C4")   name = "C4S0";
	if (name == "C4DH") name = "C4S0DH";
	// biparental crosses:
	if (name == "DH")		   { type = DH; return; }
	if (match(x,name,"Fx"))	   { type = Fx; return; }
	if (match(x,name,"FxDH"))  { type = FxDH; return; }
	if (match(x,name,"BCx"))   { type = BCx; return; }

	// if (match(x,name,"C4SxDH")) return new pop4wSxDH(x);
}


int just_test()
{
	string name = "F3DH";
	Pop pop(name);
	switch (pop.GetType())
	{ 
		case Pop::Fx   : cout << "Fx: "   << pop.get_x() << endl; return 0;
		case Pop::FxDH : cout << "FxDH: " << pop.get_x() << endl; return 0;
		default :		 cout << "Wrong! " << endl;
	}
	return 0;
}
