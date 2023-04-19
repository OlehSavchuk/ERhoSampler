#ifndef __MSU_ERHO_SAMPLER_H__
#define __MSU_ERHO_SAMPLER_H__
#include <Eigen/Dense>
#include "msu_eos/resonances.h"
#include "msu_sampler/part.h"
#include "msu_sampler/classdefs.h"
#include "msu_sampler/master.h"
#include "msu_sampler/sampler.h"
#include "msu_eos/eos.h"
#include "msu_commonutils/log.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/decay_nbody.h"
using namespace std;

namespace NMSU_ERrhoSampler{
	void FillOutHyperBjorken(Chyper *hyper,double T,double tau,double R,double deleta,double rhoB,double rhoII);
	void DecayParts(Crandy *randyset,CpartList *partlist);
	void IncrementEQtest(CpartList *partlist,Eigen::VectorXd &EQtot,Eigen::VectorXd &EQTarget);
	void Chi4Test(CpartList *partlist,Eigen::MatrixXd &chitest);
}

class CcorrVsEta{
public:
	double DETA;
	virtual void GetCorrVsEta(double eta,Eigen::MatrixXd &corr){
		//This is a dummy
	}
};
class CcorrVsEtaScott : public CcorrVsEta{
public:
	void GetCorrVsEta(double eta,Eigen::MatrixXd &corr);
};

class CcorrVsEtaOleh : public CcorrVsEta{
public:
	void GetCorrVsEta(double eta,Eigen::MatrixXd &corr);
    CcorrVsEtaOleh()
    {
        std::string fol_name= "./CFs/";
        this->init1(fol_name);
    }
private:
    ifstream t_file;
    void init1(std::string);
    static const int Ncharges=7;
    static const int Netas=1002;
    std::string charg_names[Ncharges]={"E","Ux","Uy","Uz","B","Q","S"};
    std::vector<std::vector<double>> corr_input=std::vector<std::vector<double>>(Ncharges*Ncharges,std::vector<double>(Netas,0  ));
    int index_eta(double);
    int index_charge(int,int);
};

class CcorrVsY{
public:
	CparameterMap *parmap;
	double DY;
	int NY;
	Crandy *randy;
	CcorrVsY(CparameterMap *parmap);
	vector<Eigen::MatrixXd> corr;
	double denom;
	void Increment(CpartList *partlista,CpartList *partlistb,CcorrVsEta *corrvseta);
	double GetPairWeight(Cpart *part1,Cpart *part2,Eigen::MatrixXd &CorrMatrix);
	void WriteResults();
	
};




#endif
