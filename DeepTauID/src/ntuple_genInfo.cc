
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


	ADDBRANCH(t, genVisTau_pt);
	ADDBRANCH(t, genVisTau_eta);
	ADDBRANCH(t, genVisTau_phi);
	ADDBRANCH(t, genVisTau_M);
	ADDBRANCH(t, genVisTauDeltaR);

	ADDBRANCH(t, genVisTauOrJet_pt);
	ADDBRANCH(t, genVisTauOrJet_eta);
	ADDBRANCH(t, genVisTauOrJet_phi);
	ADDBRANCH(t, genVisTauOrJet_M);

	ADDBRANCH(t, genTauDecayMode);
	ADDBRANCH(t, genTauMatch);

	ADDBRANCH(t, isTau);
	ADDBRANCH(t, isNoTau);


	ADDBRANCH(t,  isElectron);
	ADDBRANCH(t,  isMuon);

	ADDBRANCH(t, isTauToElectron);
	ADDBRANCH(t, isTauToMuon);

}

bool ntuple_genInfo::fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau){

	if(genTau){


		bool istau = dec_helper_.isPromptTau(*genTau);
		genDecayHelper::decayModes dec=dec_helper_.getGenTauDecayMode(genTau);
		istau &=  !(dec == genDecayHelper::electron || dec == genDecayHelper::muon);

		//these would be a fake coming from electron -> fill as electron/muon

		isTau=istau;
		isNoTau=!istau;


		isElectron = abs(genTau->pdgId()) == 11;
		isTauToElectron = dec == genDecayHelper::electron ;
		isMuon = abs(genTau->pdgId()) == 13;
		isTauToMuon = dec == genDecayHelper::muon;



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

		reco::Candidate::LorentzVector visP4 = dec_helper_.getVisMomentum(genTau);

		genVisTau_pt=visP4.pt();
		genVisTau_eta=visP4.eta();
		genVisTau_phi=visP4.phi();
		genVisTau_M=visP4.M();
		genVisTauDeltaR=reco::deltaR(visP4,genTau->p4());


	    genVisTauOrJet_pt=visP4.pt();
	    genVisTauOrJet_eta=visP4.eta();
	    genVisTauOrJet_phi=visP4.phi();
	    genVisTauOrJet_M=visP4.M();

	}
	else{
		clear();
		isTau=0;
		isNoTau=1;

		if(recJet && recJet->genJet()){

		    genVisTauOrJet_pt=recJet->genJet()->pt();
		    genVisTauOrJet_eta=recJet->genJet()->eta();
		    genVisTauOrJet_phi=recJet->genJet()->phi();
		    genVisTauOrJet_M=recJet->genJet()->mass();

		}
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

	genVisTauOrJet_pt=0;
	genVisTauOrJet_eta=0;
	genVisTauOrJet_phi=0;
	genVisTauOrJet_M=0;

	genTauDecayMode=-1;
	genTauMatch=0;

	isTau=0;
	isNoTau=1;
	isElectron=0;
	isMuon=0;

    isTauToElectron=0;
    isTauToMuon=0;
}
