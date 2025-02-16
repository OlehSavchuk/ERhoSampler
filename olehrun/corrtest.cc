#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"
#include "msu_erhosampler/erhosampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	Crandy *randy=new Crandy(time(NULL));
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters/parameters.txt");
	int NEVENTS_TOT=parmap.getI("NEVENTS_TOT",10);
	
	CresList *reslist=new CresList(&parmap);
	CresInfo::randy=randy;
	
	Chyper *hyper=new Chyper();
	CcorrVsEtaOleh corrvseta;
	CcorrVsY corrvsy(&parmap);
	
	int ievent;
	double T=0.150,tau=10.0,R=5.0,deleta=0.05;
	double rhoB=0.1,rhoQ;
	rhoQ=0.4*rhoB;
	Csampler *sampler=new Csampler(T,0.093,&parmap,reslist,randy);
	
	CpartList *partlista=new CpartList(&parmap,reslist);
	CpartList *partlistb=new CpartList(&parmap,reslist);
	
	NMSU_ERrhoSampler::FillOutHyperBjorken(hyper,T,tau,R,deleta,rhoB,rhoQ);
	hyper->sampler=sampler;
	
	sampler->CalcDensitiesMu0();
	sampler->GetNHMu0();
	sampler->GetMuNH(hyper);
	sampler->CalcChi4BQS(hyper);
	
	for(ievent=0;ievent<NEVENTS_TOT;ievent++){
		sampler->partlist=partlista;
		sampler->MakeParts(hyper);
		
		partlista->SetEQWeightVec(hyper);
		NMSU_ERrhoSampler::DecayParts(randy,partlista);
		
		
		sampler->partlist=partlistb;
		sampler->MakeParts(hyper);
		
		partlistb->SetEQWeightVec(hyper);
		NMSU_ERrhoSampler::DecayParts(randy,partlistb);
		
		corrvsy.Increment(partlista,partlistb,&corrvseta);
		partlista->Clear();
		partlistb->Clear();
		
		if(10*(ievent+1)%NEVENTS_TOT==0){
			printf("finished %d percent\n",((ievent+1)*100)/NEVENTS_TOT);
		}
	}
	
	corrvsy.WriteResults();
	
	delete partlista;
	delete partlistb;
	
	return 0;
}
