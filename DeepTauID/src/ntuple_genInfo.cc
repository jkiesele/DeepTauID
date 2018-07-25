
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


	ADDBRANCH(t,  isElectron);
	ADDBRANCH(t,  isMuon);

}

bool ntuple_genInfo::fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau){

	if(genTau){


		bool istau = dec_helper_.isPromptTau(*genTau);
		genDecayHelper::decayModes dec=dec_helper_.getGenTauDecayMode(genTau);
		istau &=  !(dec == genDecayHelper::electron || dec == genDecayHelper::muon);

		//these would be a fake coming from electron -> fill as electron/muon

		isTau=istau;
		isNoTau=!istau;


		isElectron = (dec == genDecayHelper::electron || abs(genTau->pdgId()) == 11);
		isMuon = (dec == genDecayHelper::muon || abs(genTau->pdgId()) == 13);



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

		genTauMatch=isTau;
		genTauDecayMode=dec_helper_.getGenTauDecayMode(genTau);

	}
	else{
		clear();
		isTau=0;
		isNoTau=1;
	}
	return true;
}


ntuple_genInfo::ntuple_genInfo():ntuple_content(){
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
	isElectron=0;
	isMuon=0;
}
