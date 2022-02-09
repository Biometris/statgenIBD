#include "SimPop.h"

using namespace ibd;
using namespace std;

SIMpop_base * init_SIMpop(const std::string& poptype)
{
	int x,y;
	// biparental crosses:
	if (poptype == "DH")				return new SIMpopDH();
	if (match(x,poptype,"Fx"))			return new SIMpopFx(x);
	if (match(x,poptype,"FxDH"))		return new SIMpopFxDH(x);
	if (match(x,poptype,"BCx"))			return new SIMpopBCx(x);
	if (match(x,poptype,"BCxDH"))		return new SIMpopBCxDH(x);

	if (match(x,y,poptype,"BCxSy"))		return new SIMpopBCxSy(x,y);
	if (match(x,y,poptype,"BCxSyDH"))	return new SIMpopBCxSyDH(x,y);

	// three-way crosses:
	if (poptype == "C3") 				return new SIMpopC3Sx(0);   
	if (poptype == "C3DH")				return new SIMpopC3SxDH(0); 
	if (match(x,poptype,"C3Sx"))		return new SIMpopC3Sx(x);
	if (match(x,poptype,"C3SxDH"))		return new SIMpopC3SxDH(x);
	if (match(x,y,poptype,"C3RCxSy"))	return new SIMpopC3RCxSy(x,y);
	if (match(x,y,poptype,"C3RCxSyDH"))	return new SIMpopC3RCxSyDH(x,y);

	// four-way crosses:
	if (poptype == "C4")				return new SIMpopC4Sx(0);   
	if (poptype == "C4DH")				return new SIMpopC4SxDH(0); 
	if (match(x,poptype,"C4Sx"))		return new SIMpopC4Sx(x);
	if (match(x,poptype,"C4SxDH"))		return new SIMpopC4SxDH(x);

	throw ibd_error("unknown type " + poptype);
	return 0;
}

ibd::Genome SIMpopDH::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{ 
	return DoubledHaploid(p1*p2); 
}

ibd::Genome SIMpopFx::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{
	Genome g = Selfing(p1*p2,ngen);
	return g; 
}

ibd::Genome SIMpopFxDH::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{
	Genome g = Selfing(p1*p2,ngen);
	return DoubledHaploid(g); 
}

ibd::Genome SIMpopBCx::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{
	Genome g = p1*p2;
	for (int i=0;i<ngen;i++)
		g = g*p2;
	return g; 
}

ibd::Genome SIMpopBCxDH::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{
	Genome g = p1*p2;
	for (int i=0;i<ngen;i++)
		g = g*p2;
	return DoubledHaploid(g); 
}

ibd::Genome SIMpopBCxSy::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{
	Genome g = p1*p2;
	for (int i=0;i<ngen_BC;i++)
		g = g*p2;
	g = Selfing(g,ngen_S);
	return g; 
}

ibd::Genome SIMpopBCxSyDH::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{
	Genome g = p1*p2;
	for (int i=0;i<ngen_BC;i++)
		g = g*p2;
	g = Selfing(g,ngen_S);
	return DoubledHaploid(g); 
}

ibd::Genome SIMpopC3Sx::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{ 
	Genome g = p1*p2;
	Genome x = Selfing(g,ngen_self);
	return x;
}

ibd::Genome SIMpopC3SxDH::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{ 
	Genome g = p1*p2;
	Genome x = Selfing(g,ngen_self);
	return DoubledHaploid(x);
}

ibd::Genome SIMpopC3RCxSy::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{ 
	Genome g = p1; // p1 = AxB
	for (int i=0;i<ngen_RC;i++)
		g = g*p2;
	Genome x = Selfing(g,ngen_self);
	return x;
}

ibd::Genome SIMpopC3RCxSyDH::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{ 
	Genome g = p1; // p1 = AxB
	for (int i=0;i<ngen_RC;i++)
		g = g*p2;
	Genome x = Selfing(g,ngen_self);
	return DoubledHaploid(x);
}


ibd::Genome SIMpopC4Sx::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{ 
	Genome g = p1*p2;
	Genome x = Selfing(g,ngen_self);
	return x;
}

ibd::Genome SIMpopC4SxDH::sim(const ibd::Genome& p1, const ibd::Genome& p2) const
{ 
	Genome g = p1*p2;
	Genome x = Selfing(g,ngen_self);
	return DoubledHaploid(x);
}

