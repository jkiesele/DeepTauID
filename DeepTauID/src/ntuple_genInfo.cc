
#include "../interface/ntuple_genInfo.h"
#include <iostream>
#include <stdexcept>

void ntuple_genInfo::initBranches(TTree* t){

	// TBI FIXME
//	ADDBRANCH(t, genEvtWeight);
//	ADDBRANCH(t, evtWeight);


	ADDBRANCH(t, genTau_pt);
	ADDBRANCH(t, genTau_eta);
	ADDBRANCH(t, genTau_phi);
	ADDBRANCH(t, genTau_M);
	ADDBRANCH(t, genTauDeltaR);

	/*
	ADDBRANCH(t, genVisTau_pt);
	ADDBRANCH(t, genVisTau_eta);
	ADDBRANCH(t, genVisTau_phi);
	ADDBRANCH(t, genVisTau_M);
	ADDBRANCH(t, genVisTauDeltaR);
*/

	ADDBRANCH(t, genTauDecayMode);
	ADDBRANCH(t, genTauMatch);

	ADDBRANCH(t, isTau);
	ADDBRANCH(t, isNoTau);

}

bool ntuple_genInfo::fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau){
	if(genTau){
		isTau=1;
		isNoTau=0;


		// TBI FIXME
		genEvtWeight=1;
		evtWeight=1;

		genTau_pt=genTau->pt();
		genTau_eta=genTau->eta();
		genTau_phi=genTau->phi();
		genTau_M=genTau->mass();
		genTauDeltaR=1;
		if(recJet)
			genTauDeltaR=reco::deltaR(recJet->p4(),genTau->p4());

		genTauMatch=1;
		genTauDecayMode=dec_helper_.getGenTauDecayMode(genTau);

		if(genTau_pt<minpt_)return false;
		if(fabs(genTau_eta)>maxeta_) return false;
	}
	else{
		clear();
		isTau=0;
		isNoTau=1;
	}
	return true;
}


ntuple_genInfo::ntuple_genInfo():ntuple_content(){
	minpt_=1;
	maxeta_=3.1;
	clear();
}

void ntuple_genInfo::clear(){
	genEvtWeight=1;
	evtWeight=1;

	genTau_pt=0;
	genTau_eta=0;
	genTau_phi=0;
	genTau_M=0;
	genTauDeltaR=0;

	genVisTau_pt=0;
	genVisTau_eta=0;
	genVisTau_phi=0;
	genVisTau_M=0;
	genVisTauDeltaR=0;

	genTauDecayMode=0;
	genTauMatch=0;

	isTau=0;
	isNoTau=1;
}
