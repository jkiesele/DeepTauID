

#include "../interface/genDecayHelper.h"
#include "TMath.h"

reco::Candidate::LorentzVector genDecayHelper::getVisMomentum(const reco::GenParticle* genTau)const
{
	std::vector<const reco::GenParticle*> stableDaughters;
	findDaughters(genTau, stableDaughters, 1);
	reco::Candidate::LorentzVector genVisTauP4 = getVisMomentum(stableDaughters, 1);
	return genVisTauP4;
}

bool genDecayHelper::isPromptTau(const reco::GenParticle& genParticle)const{

	if ( !genParticle.isPromptDecayed() ) return false;
	if( ! (abs(genParticle.pdgId()) == 15))return false;
	return true;
}

bool genDecayHelper::isPromptLepton(const reco::GenParticle& genParticle)const{

	if ( !(genParticle.isPromptFinalState() ||  genParticle.isPromptDecayed()) ) return false;
	if(  abs(genParticle.pdgId()) == 15 || abs(genParticle.pdgId()) == 13 || abs(genParticle.pdgId()) == 11)return true;
	return false;
}

bool genDecayHelper::isLooseGenLepton(const reco::GenParticle& genParticle)const{
	if((genParticle.isLastCopy() && genParticle.statusFlags().fromHardProcess())){
		if(  abs(genParticle.pdgId()) == 15 || abs(genParticle.pdgId()) == 13 || abs(genParticle.pdgId()) == 11)return true;
	}
	else if(isPromptLepton(genParticle)) return true;
	return false;
}


bool genDecayHelper::isNeutrino(const reco::GenParticle & daughter)const{
	return ( TMath::Abs(daughter.pdgId()) == 12 || TMath::Abs(daughter.pdgId()) == 14 || TMath::Abs(daughter.pdgId()) == 16 );
}

reco::Candidate::LorentzVector genDecayHelper::getVisMomentum(const std::vector<const reco::GenParticle*>& daughters, int status)const
{
	reco::Candidate::LorentzVector p4Vis(0,0,0,0);
	for ( std::vector<const reco::GenParticle*>::const_iterator daughter = daughters.begin();
			daughter != daughters.end(); ++daughter ) {
		if ( (status == -1 || (*daughter)->status() == status) && !isNeutrino(*(*daughter)) ) {
			p4Vis += (*daughter)->p4();
		}
	}
	return p4Vis;
}

void genDecayHelper::findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters, int status)const
{
	unsigned numDaughters = mother->numberOfDaughters();
	for ( unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
		const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(mother->daughter(iDaughter));
		if ( status == -1 || daughter->status() == status ) daughters.push_back(daughter);
		findDaughters(daughter, daughters, status);
	}
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
		return other;
}
