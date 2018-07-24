

#include "../interface/genDecayHelper.h"


bool genDecayHelper::selectTau(const reco::GenParticle& genParticle)const{

	if ( !genParticle.isPromptDecayed() ) return false;
	if( ! abs(genParticle.pdgId()) == 15)return false;
	return true;
}



genDecayHelper::decayModes genDecayHelper::getGenTauDecayMode(const reco::GenParticle* genTau)const{

	//--- determine generator level tau decay mode
		//
		//    NOTE:
		//        (1) function implements logic defined in PhysicsTools/JetMCUtils/src/JetMCTag::genTauDecayMode
		//            for different type of argument
		//        (2) this implementation should be more robust to handle cases of tau --> tau + gamma radiation
		//
		int numElectrons           = 0;
		int numElecNeutrinos       = 0;
		int numMuons               = 0;
		int numMuNeutrinos         = 0;
		int numChargedHadrons      = 0;
		int numPi0s                = 0;
		int numOtherNeutralHadrons = 0;
		int numPhotons             = 0;

		countDecayProducts<reco::GenParticle>(genTau,
			numElectrons, numElecNeutrinos, numMuons, numMuNeutrinos,
			numChargedHadrons, numPi0s, numOtherNeutralHadrons, numPhotons);

		if      ( numElectrons == 1 && numElecNeutrinos == 1 ) return electron;
		else if ( numMuons     == 1 && numMuNeutrinos   == 1 ) return muon;

		switch ( numChargedHadrons ) {
		case 1 :
			if ( numOtherNeutralHadrons != 0 ) return oneProngOther;
			switch ( numPi0s ) {
			case 0:
				return oneProng0Pi0;
			case 1:
				return oneProng1Pi0;
			case 2:
				return oneProng2Pi0;
			default:
				return oneProngOther;
			}
			case 3 :
				if ( numOtherNeutralHadrons != 0 ) return threeProngOther;
				switch ( numPi0s ) {
				case 0:
					return threeProng0Pi0;
				case 1:
					return threeProng1Pi0;
				default:
					return threeProngOther;
				}
				default:
					return rare;
		}
}
