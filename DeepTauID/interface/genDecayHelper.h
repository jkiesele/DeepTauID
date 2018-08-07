/*
 * genDecayHelper.h
 *
 *  Created on: 24 Jul 2018
 *      Author: jkiesele
 */

#ifndef DEEPTAU_DEEPTAUID_INTERFACE_GENDECAYHELPER_H_
#define DEEPTAU_DEEPTAUID_INTERFACE_GENDECAYHELPER_H_


#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TMath.h"


class genDecayHelper{
	public:
	genDecayHelper(){}

	enum decayModes{
		electron,
		muon,
		oneProngOther,
		oneProng0Pi0,
		oneProng1Pi0,
		oneProng2Pi0,
		threeProngOther,
		threeProng0Pi0,
		threeProng1Pi0,
		rare,
		other
	};

	decayModes getGenTauDecayMode(const reco::GenParticle* genTau)const;

	template<class T>
	void countDecayProducts(const T* genParticle,
			int& numElectrons, int& numElecNeutrinos, int& numMuons, int& numMuNeutrinos,
			int& numChargedHadrons, int& numPi0s, int& numOtherNeutralHadrons, int& numPhotons) const;

	bool isPromptTau(const reco::GenParticle& )const;
	bool isPromptLepton(const reco::GenParticle& )const;

	bool isNeutrino(const reco::GenParticle& )const;

	reco::Candidate::LorentzVector getVisMomentum(const reco::GenParticle* tau)const;


private:
	reco::Candidate::LorentzVector getVisMomentum(const std::vector<const reco::GenParticle*>& daughters, int status)const;
	void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters, int status=1)const;
};

template<class T>
void genDecayHelper::countDecayProducts(const T* genParticle,
		int& numElectrons, int& numElecNeutrinos, int& numMuons, int& numMuNeutrinos,
		int& numChargedHadrons, int& numPi0s, int& numOtherNeutralHadrons, int& numPhotons)const{

	if(!genParticle)return;

	int absPdgId = TMath::Abs(genParticle->pdgId());
	int status   = genParticle->status();
	int charge   = genParticle->charge();

	if ( absPdgId == 111 ) ++numPi0s;
	else if ( status == 1 ) {
		if      ( absPdgId == 11 ) ++numElectrons;
		else if ( absPdgId == 12 ) ++numElecNeutrinos;
		else if ( absPdgId == 13 ) ++numMuons;
		else if ( absPdgId == 14 ) ++numMuNeutrinos;
		else if ( absPdgId == 15 ) {
			edm::LogError ("genDecayHelper::countDecayProducts")
			<< "Found tau lepton with status code 1 !!";
			return;
		}
		else if ( absPdgId == 16 ) return; // no need to count tau neutrinos
		else if ( absPdgId == 22 ) ++numPhotons;
		else if ( charge   !=  0 ) ++numChargedHadrons;
		else                       ++numOtherNeutralHadrons;
	} else {
		unsigned numDaughters = genParticle->numberOfDaughters();
		for ( unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
			const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(genParticle->daughter(iDaughter));

			countDecayProducts<reco::GenParticle>(daughter,
					numElectrons, numElecNeutrinos, numMuons, numMuNeutrinos,
					numChargedHadrons, numPi0s, numOtherNeutralHadrons, numPhotons);
		}
	}
}


#endif /* DEEPTAU_DEEPTAUID_INTERFACE_GENDECAYHELPER_H_ */
