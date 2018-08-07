/*
 * ntuple_genInfo.h
 *
 *  Created on: 20 Jul 2018
 *      Author: jkiesele
 */

#ifndef DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_GENINFO_H_
#define DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_GENINFO_H_

#include "ntuple_content.h"
#include "genDecayHelper.h"
/**
 * Contains gen information for tau / jet
 */
class ntuple_genInfo: public ntuple_content{
public:
	ntuple_genInfo();


    void getInput(const edm::ParameterSet& iConfig){}
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent){}
    void readSetup(const edm::EventSetup& iSetup){}
    //use either of these functions

    bool fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau);


    void clear();

private:

    genDecayHelper dec_helper_;

public:
    //branches
    float genEvtWeight;
    float evtWeight;

    float genTau_pt;
    float genTau_eta;
    float genTau_phi;
    float genTau_M;
    float genTauDeltaR;

    float genVisTau_pt;
    float genVisTau_eta;
    float genVisTau_phi;
    float genVisTau_M;
    float genVisTauDeltaR;

    float genVisTauOrJet_pt;
    float genVisTauOrJet_eta;
    float genVisTauOrJet_phi;
    float genVisTauOrJet_M;

    int genTauDecayMode;
    float genTauMatch;

    int isTau;
    int isNoTau;

    int isElectron;
    int isMuon;

    int isTauToElectron;
    int isTauToMuon;


};



#endif /* DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_GENINFO_H_ */
