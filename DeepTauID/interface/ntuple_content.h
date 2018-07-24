/*
 * ntuple_content.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */

#ifndef DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_CONTENT_H_
#define DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_CONTENT_H_


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TTree.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include <iostream>
#include <math.h>
#include <iostream>
#include "TString.h"

/**
 * Base class for modules to inherit from.
 */
class ntuple_content{
public:

    ntuple_content():rhoInfo_(0),read_(false){}
    virtual ~ntuple_content();

    virtual void getInput(const edm::ParameterSet& iConfig)=0;
    virtual void initBranches(TTree* )=0;
    virtual void readEvent(const edm::Event& iEvent)=0;
    virtual void readSetup(const edm::EventSetup& iSetup)=0;
    //use either of these functions

    virtual bool fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau)=0;

    virtual void clear(){};

    void setIsRead(bool isread){read_=isread;}

    void setPrimaryVertices(const reco::VertexCollection* v){
            vertices_=v;
        }

    std::vector<TString> getListOfBranches(){
        if(allbranches_.size())
            return allbranches_;
        else{
            TTree *t=new TTree();
            initBranches(t);
            return allbranches_;
        }
    }

    void setRhoInfo(const double* r){rhoInfo_=r;}
    const double* rhoInfo()const;
    const reco::VertexCollection * vertices()const{return vertices_;}

    static bool useoffsets;
    static bool debug;

protected:



    template <class T>
    void addBranch(TTree* t, const char* name,  T*, const char* leaflist=0);


    static inline const float& catchInfs(const float& in,const float& replace_value){
        if(in==in){
            if(std::isinf(in))
                return replace_value;
            else if(in < -1e32 || in > 1e32)
                return replace_value;
            return in;
        }
        return replace_value;
    }

    static inline float catchInfsAndBound(const float& in,const float& replace_value,
            const float& lowerbound, const float& upperbound,const float offset=0){
        float withoutinfs=catchInfs(in,replace_value);
        if(withoutinfs+offset<lowerbound) return lowerbound;
        if(withoutinfs+offset>upperbound) return upperbound;
        if(useoffsets)
            withoutinfs+=offset;
        return withoutinfs;
    }

private:
    const reco::VertexCollection* vertices_;
    const double* rhoInfo_;
    bool read_;
    std::vector<TString> allbranches_;
};

template <class T>
void ntuple_content::addBranch(TTree* t, const char* name,  T* address, const char* leaflist){

    if(read_ ){
        t->SetBranchAddress(name,address);
    }
    else{
        if(leaflist)
            t->Branch(name  ,address  ,leaflist );
        else
            t->Branch(name  ,address);
    }
    allbranches_.push_back((TString)name);

}

#define ADDBRANCH(t,x) addBranch(t, #x, &x)



#endif /* DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_CONTENT_H_ */
