/*
 * ntuple_JetInfo.h
 *
 *  Created on: 20 Jul 2018
 *      Author: jkiesele
 */

#ifndef DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_JETINFO_H_
#define DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_JETINFO_H_

#include "ntuple_content.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
/**
 *
 * Contains global jet information
 *
 */
class ntuple_JetInfo: public ntuple_content{
public:
    ntuple_JetInfo();
    ~ntuple_JetInfo();


    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent){}
    void readSetup(const edm::EventSetup& iSetup){}
    //use either of these functions

    bool fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau);


    void clear();
private:

    PFJetIDSelectionFunctor* loosePFJetIdAlgo_;

    float recJet_pt;
    float recJet_eta;
    float recJet_phi;
    float recJet_mass;
    float recJetLooseId;

    int hasGenMatch;


};



#endif /* DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_JETINFO_H_ */
