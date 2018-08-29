/*
 * ntuple_pfCands.h
 *
 *  Created on: 20 Jul 2018
 *      Author: jkiesele
 */

#ifndef DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_PFCANDS_H_
#define DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_PFCANDS_H_



#include "ntuple_content.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include <vector>

/**
 * Contains individual information of pf Candidates around the tau / in the tau/jet
 */
class ntuple_pfCands: public ntuple_content{
public:
	ntuple_pfCands():ntuple_content(),nCpfcan(0),nNpfcand(0){}


	void getInput(const edm::ParameterSet& iConfig);
	void initBranches(TTree* );
	void readEvent(const edm::Event& iEvent){}
	void readSetup(const edm::EventSetup& iSetup);
	//use either of these functions

	void clear();

	bool fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau);

private:

	float catchInfsAndBound_track(const pat::PackedCandidate* cand, const float& in,const float& replace_value,
            const float& lowerbound, const float& upperbound,const float offset=0)const;

	edm::ESHandle<TransientTrackBuilder> builder;

	float nCpfcan;
	std::vector<float>   Cpfcan_pt;
	std::vector<float>   Cpfcan_eta;
	std::vector<float>   Cpfcan_phi;
	std::vector<float>   Cpfcan_ptrel;
	std::vector<float>   Cpfcan_erel;
	std::vector<float>   Cpfcan_phirel;
	std::vector<float>   Cpfcan_etarel;
	std::vector<float>   Cpfcan_deltaR;
	std::vector<float>   Cpfcan_puppiw;
	std::vector<float>   Cpfcan_VTX_ass;

	std::vector<float>   Cpfcan_fromPV;

	std::vector<float>  Cpfcan_vertexChi2;
	std::vector<float>  Cpfcan_vertexNdof;
	std::vector<float>  Cpfcan_vertexNormalizedChi2;
	std::vector<float>  Cpfcan_vertex_rho;
	std::vector<float>  Cpfcan_vertex_phirel;
	std::vector<float>  Cpfcan_vertex_etarel;
	std::vector<float>  Cpfcan_vertexRef_mass;

	// covariance
	std::vector<float>   Cpfcan_dz;
	std::vector<float>   Cpfcan_dxy;

	std::vector<float>   Cpfcan_dxyerrinv;
	std::vector<float>   Cpfcan_dxysig;

	std::vector<float>  Cpfcan_BtagPf_trackMomentum;
	std::vector<float>  Cpfcan_BtagPf_trackEta;
	std::vector<float>  Cpfcan_BtagPf_trackEtaRel;
	std::vector<float>  Cpfcan_BtagPf_trackPtRel;
	std::vector<float>  Cpfcan_BtagPf_trackPPar;
	std::vector<float>  Cpfcan_BtagPf_trackDeltaR;
	std::vector<float>  Cpfcan_BtagPf_trackPtRatio;
	std::vector<float>  Cpfcan_BtagPf_trackPParRatio;
	std::vector<float>  Cpfcan_BtagPf_trackSip3dVal;
	std::vector<float>  Cpfcan_BtagPf_trackSip3dSig;
	std::vector<float>  Cpfcan_BtagPf_trackSip2dVal;
	std::vector<float>  Cpfcan_BtagPf_trackSip2dSig;


	std::vector<float>  Cpfcan_BtagPf_trackJetDistVal;

	std::vector<float>  Cpfcan_isMu;
	std::vector<float>  Cpfcan_isEl;
	std::vector<float>  Cpfcan_pdgID;

	std::vector<float>  Cpfcan_charge;

	// track quality
	std::vector<float>  Cpfcan_lostInnerHits;
	std::vector<float>  Cpfcan_numberOfPixelHits;

	std::vector<float>  Cpfcan_chi2;
	std::vector<float>  Cpfcan_quality;


	float  nNpfcand;
	std::vector<float>   Npfcan_pt;
	std::vector<float>   Npfcan_eta;
	std::vector<float>   Npfcan_phi;
	std::vector<float>   Npfcan_ptrel;
	std::vector<float>   Npfcan_erel;
	std::vector<float>   Npfcan_puppiw;
	std::vector<float>   Npfcan_phirel;
	std::vector<float>   Npfcan_etarel;
	std::vector<float>   Npfcan_deltaR;
	std::vector<float>   Npfcan_isGamma;
	std::vector<float>   Npfcan_HadFrac;
	std::vector<float>   Npfcan_drminsv;
	std::vector<float>   Npfcan_pdgID;



};




#endif /* DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_PFCANDS_H_ */
