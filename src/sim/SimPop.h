#ifndef SIMPOP_SIMQTL_HEADER
#define SIMPOP_SIMQTL_HEADER

#include "Genome.h"
#include "SmartPtr.h"
#include <string>

class SIMpop_base
{
public:
	SIMpop_base() {}
	virtual ~SIMpop_base(){}
	virtual ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const = 0;
};

class SIMpopDH : public SIMpop_base
{
public:
	SIMpopDH() {}
	virtual ~SIMpopDH(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
};

class SIMpopFx : public SIMpop_base
{
public:
	SIMpopFx(int x) : ngen(x-1) {}
	virtual ~SIMpopFx(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen;
};

class SIMpopFxDH : public SIMpop_base
{
public:
	SIMpopFxDH(int x) : ngen(x-1) {}
	virtual ~SIMpopFxDH(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen;
};

class SIMpopBCx : public SIMpop_base
{
public:
	SIMpopBCx(int x) : ngen(x) {}
	virtual ~SIMpopBCx(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen;
};

class SIMpopBCxDH : public SIMpop_base
{
public:
	SIMpopBCxDH(int x) : ngen(x) {}
	virtual ~SIMpopBCxDH(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen;
};

class SIMpopBCxSy : public SIMpop_base
{
public:
	SIMpopBCxSy(int x, int y) : ngen_BC(x), ngen_S(y) {}
	virtual ~SIMpopBCxSy(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen_BC;
	int ngen_S;
};

class SIMpopBCxSyDH : public SIMpop_base
{
public:
	SIMpopBCxSyDH(int x, int y) : ngen_BC(x), ngen_S(y) {}
	virtual ~SIMpopBCxSyDH(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen_BC;
	int ngen_S;
};

class SIMpopC3Sx : public SIMpop_base
{
public:
	SIMpopC3Sx(int ngen) : ngen_self(ngen) {}
	virtual ~SIMpopC3Sx(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen_self;
};

class SIMpopC3SxDH : public SIMpop_base
{
public:
	SIMpopC3SxDH(int ngen) : ngen_self(ngen) {}
	virtual ~SIMpopC3SxDH(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen_self;
};

class SIMpopC3RCxSy : public SIMpop_base
{
public:
	SIMpopC3RCxSy(int x, int y) : ngen_RC(x), ngen_self(y) {}
	virtual ~SIMpopC3RCxSy(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen_RC;
	int ngen_self;
};

class SIMpopC3RCxSyDH : public SIMpop_base
{
public:
	SIMpopC3RCxSyDH(int x, int y) : ngen_RC(x), ngen_self(y) {}
	virtual ~SIMpopC3RCxSyDH(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen_RC;
	int ngen_self;
};

class SIMpopC4Sx : public SIMpop_base
{
public:
	SIMpopC4Sx(int ngen) : ngen_self(ngen) {}
	virtual ~SIMpopC4Sx(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen_self;
};

class SIMpopC4SxDH : public SIMpop_base
{
public:
	SIMpopC4SxDH(int ngen) : ngen_self(ngen) {}
	virtual ~SIMpopC4SxDH(){}
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const;
private:
	int ngen_self;
};

SIMpop_base * init_SIMpop(const std::string& poptype);

class PopProp
{
public:
	PopProp(int n,std::string crossname, std::string popt, std::vector<std::string> par)
		: nind(n), npar(par.size()), name(crossname), type(popt), P(par) 
			{ ptr = init_SIMpop(popt); }
	ibd::Genome sim(const ibd::Genome& p1, const ibd::Genome& p2) const
		{ return ptr->sim(p1,p2); }
	int GetNpar() const { return npar; }
	int GetNind() const { return nind; }
	std::string GetType() const { return type; }
	std::string GetName() const { return name; }
	std::string GetParName(int i) const { return P[i]; } 
private:
	int nind,npar;
	std::string name,type;
	std::vector<std::string> P;
	ibd::SmartPtr<SIMpop_base> ptr;
};

#endif