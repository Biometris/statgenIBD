#include <Rcpp.h>
#include <iomanip>

#include "sim.h"
#include "urn.h"
#include "misc.h"

using namespace ibd;
using namespace std;

SimPop sim_population(const PopProp& pop, const map<string,Genome>& genome_par, const string& pop_name)
{
	const int npar = pop.GetNpar();
	const int nhybr = npar - 2;
	const int nheader = npar + nhybr;
	const int ntot = nheader + pop.GetNind();
	vector<Genome> par(npar);
	for (int i=0;i<npar;i++)
	{
		map<string,Genome>::const_iterator it = genome_par.find(pop.GetParName(i));
		if (it!=genome_par.end())
			par[i] = it->second;
		else
			throw ibd_error("Cannot find parent " + pop.GetParName(i));
	}
	Genome P1,P2;
	string P1_name,P2_name;

	SimPop ped(ntot);  
	for (int i=0;i<npar;i++)
		ped[i] = SimInd(IndProp(pop.GetParName(i),"*","INBPAR","0","0"), par[i]);
	if (npar == 4)
	{
		P1 = par[0]*par[1];
		P2 = par[2]*par[3];
		P1_name = "H1";
		P2_name = "H2";
		ped[npar]   = SimInd(IndProp(P1_name,"*","HYBRID",pop.GetParName(0),pop.GetParName(1)),P1);
		ped[npar+1] = SimInd(IndProp(P2_name,"*","HYBRID",pop.GetParName(2),pop.GetParName(3)),P2);
	}
	else if (npar == 3)
	{
		P1 = par[0]*par[1];
		P2 = par[2];
		P1_name = "H";
		P2_name = pop.GetParName(2);
		ped[npar]   = SimInd(IndProp(P1_name,"*","HYBRID",pop.GetParName(0),pop.GetParName(1)),P1);
	}
	else if (npar == 2)
	{
		P1 = par[0];
		P2 = par[1];
		P1_name = pop.GetParName(0);
		P2_name = pop.GetParName(1);
	}
	else
		throw ibd_error("npar: " + npar);

	MakeLabel ID_name(pop_name,4);
	for (int i=0;i<pop.GetNind();i++)
	{
		string ID = ID_name(i);
		Genome x = pop.sim(P1,P2);
		IndProp indprop(ID,"FAM",pop.GetType(),P1_name,P2_name);
		ped[i+nheader] = SimInd(indprop,x);
	}
	return ped;
}


SimPop sim_multiple_populations(const vector<PopProp>& pops, const map<string,Genome>& genome_par)
{
	SimPop sim_pop1,sim_pop2;
	for (vector<PopProp>::const_iterator it=pops.begin();it!=pops.end();++it)
	{
		string pop_name = it->GetName();
		SimPop sim_sub_pop = sim_population(*it,genome_par,pop_name);
		for (SimPop::const_iterator it2=sim_sub_pop.begin();it2!=sim_sub_pop.end();++it2)
		{
			if (it2->first.IsInbredParent())
			{
				bool already_defined = false;
				string ID = it2->first.GetID();
				for (SimPop::const_iterator it3=sim_pop1.begin();it3!=sim_pop1.end();it3++)
				{
					if (it3->first.GetID() == ID)
						already_defined = true;
				}
				if (!already_defined)
					sim_pop1.push_back(*it2);
			}
			else
			{
				sim_pop2.push_back(*it2);
			}
		}
	}		
	for (SimPop::const_iterator it=sim_pop2.begin();it!=sim_pop2.end();it++)
		sim_pop1.push_back(*it);

	return sim_pop1;
}


map<string,Genome> sim_pedigree(const vector<IndProp>& pedigree, const map<string,Genome>& founders)
{
	map<string,Genome> simped;
	for (vector<IndProp>::const_iterator it=pedigree.begin();it!=pedigree.end();it++)
	{
		string id = it->GetID();
		if (it->IsFounder())
		{
			map<string,Genome>::const_iterator it_fnd = founders.find(id);
			if (it_fnd == founders.end())
				throw ibd_error("cannot find founder " + id);
			simped[id] = it_fnd->second;
		}
		else if (it->IsRIL())  
		{
			const int Ngen = 25;
			Genome x = simped[it->GetP1()]*simped[it->GetP2()];
			x = Selfing(x,Ngen);
			x = DoubledHaploid(x);
			simped[id] = x;
		}
		else 
			throw ibd_error("error in sim_pedigree!");
	}	
	return simped;
}


void sim_SS_NS_pedigree(const vector<string>& fnd_names, string filename, int Ngen, int N)
{
	ofstream outp;
	OpenFile(outp,filename);

	for (vector<string>::const_iterator it=fnd_names.begin();it!=fnd_names.end();it++)
	{
		outp << setw(12) << *it << setw(12) << "INBFND" 
			 << setw(12) << "0" << setw(12) << "0" << endl;
	}

	for (int grp=0;grp<2;grp++)
	{
		string grp_name = (grp==0) ? "A" : "B";

    	vector<string> sel_names = fnd_names;
		MakeLabel gen_name(grp_name,1);

		for (int gen=0;gen<Ngen;gen++)
		{
			Urn<string> urn(sel_names);
			string cur_gen_name = gen_name(gen);
			MakeLabel name(cur_gen_name+"_",2);
			vector<string> names;
			for (int i=0;i<N;i++)
			{
				string P1 = urn.random_draw();
				string P2 = urn.random_draw();
				urn.replace_all();

				string ID = name(i);
				names.push_back(ID);
				outp << setw(12) << ID << setw(12) << "RIL" 
					 << setw(12) << P1 << setw(12) << P2 << endl;
			}
			sel_names = names;
		}
	}

}

