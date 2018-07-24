
#include "../interface/ntuple_JetInfo.h"




bool ntuple_JetInfo::fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau){
	if(!recJet){
		clear();
		return true;
	}

	return true;
}
