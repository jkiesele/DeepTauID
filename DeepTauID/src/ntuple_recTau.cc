/*
 * ntuple_recTau.cc
 *
 *  Created on: 24 Jul 2018
 *      Author: jkiesele
 */


#include "../interface/ntuple_recTau.h"
#include "../interface/ntuple_content.h"

ntuple_recTau::ntuple_recTau():ntuple_content()
{
	clear();
}


void ntuple_recTau::clear(){

	isRecTau=0;
	recTauDecayMode_i=0;
	recTauDecayMode=0;
	recTau_pt=0;
	recTau_eta=0;
	recTau_phi=0;
	recTau_M=0;
	recTauVtxZ=0;
	recImpactParamPCA_x=0;
	recImpactParamPCA_y=0;
	recImpactParamPCA_z=0;
	recImpactParam=0;
	recImpactParamSign=0;

	recImpactParam3D=0;
	recImpactParamSign3D=0;

	recChi2DiffEvtVertex=0;
	//
	hasRecDecayVertex=0;


	recDecayDist_x=0;
	recDecayDist_y=0;
	recDecayDist_z=0;
	//missing recDecayDistCov
	recDecayDistSign=0;

	//missing recEvtVertexCov
	recTauPtWeightedDetaStrip=0;
	recTauPtWeightedDphiStrip=0;
	recTauPtWeightedDrSignal=0;
	recTauPtWeightedDrIsolation=0;
	recTauNphoton=0;
	recTauEratio=0;
	recTauLeadingTrackChi2=0;
	recTauNphotonSignal=0;
	recTauNphotonIso=0;
	recTauNphotonTotal=0;


	recTauLeadPFCand_pt=0;
	recTauLeadPFCand_eta=0;
	recTauLeadPFCand_phi=0;
	recTauLeadPFCand_M=0;

	recTauLeadPFChargedHadrCand_pt=0;
	recTauLeadPFChargedHadrCand_eta=0;
	recTauLeadPFChargedHadrCand_phi=0;
	recTauLeadPFChargedHadrCand_M=0;


    chargedIsoPtSum=0;
    neutralIsoPtSum=0;
    puCorrPtSum=0;
    neutralIsoPtSumWeight=0;
    footprintCorrection=0;
    photonPtSumOutsideSignalCone=0;
    chargedIsoPtSumdR03=0;
    neutralIsoPtSumdR03=0;
    neutralIsoPtSumWeightdR03=0;
    footprintCorrectiondR03=0;
    photonPtSumOutsideSignalConedR03=0;


    demetraIsolation=0;

}

void ntuple_recTau::initBranches(TTree* t){
	ADDBRANCH(t,isRecTau);
	ADDBRANCH(t,recTauDecayMode_i);
	ADDBRANCH(t,recTauDecayMode);
	ADDBRANCH(t,recTau_pt);
	ADDBRANCH(t,recTau_eta);
	ADDBRANCH(t,recTau_phi);
	ADDBRANCH(t,recTau_M);
	ADDBRANCH(t,recTauVtxZ);
	ADDBRANCH(t,recImpactParamPCA_x);
	ADDBRANCH(t,recImpactParamPCA_y);
	ADDBRANCH(t,recImpactParamPCA_z);
	ADDBRANCH(t,recImpactParam);
	ADDBRANCH(t,recImpactParamSign);

	ADDBRANCH(t,recImpactParam3D);
	ADDBRANCH(t,recImpactParamSign3D);

	// none in official producer ADDBRANCH(t,recChi2DiffEvtVertex);

	ADDBRANCH(t,hasRecDecayVertex);


	ADDBRANCH(t,recDecayDist_x);
	ADDBRANCH(t,recDecayDist_y);
	ADDBRANCH(t,recDecayDist_z);

	ADDBRANCH(t,recDecayDistSign);


	ADDBRANCH(t,recTauPtWeightedDetaStrip);
	ADDBRANCH(t,recTauPtWeightedDphiStrip);
	ADDBRANCH(t,recTauPtWeightedDrSignal);
	ADDBRANCH(t,recTauPtWeightedDrIsolation);
	ADDBRANCH(t,recTauNphoton);
	ADDBRANCH(t,recTauEratio);
	ADDBRANCH(t,recTauLeadingTrackChi2);
	ADDBRANCH(t,recTauNphotonSignal);
	ADDBRANCH(t,recTauNphotonIso);
	ADDBRANCH(t,recTauNphotonTotal);

	ADDBRANCH(t,recTauLeadPFCand_pt);
	ADDBRANCH(t,recTauLeadPFCand_eta);
	ADDBRANCH(t,recTauLeadPFCand_phi);
	ADDBRANCH(t,recTauLeadPFCand_M);

	ADDBRANCH(t,recTauLeadPFChargedHadrCand_pt);
	ADDBRANCH(t,recTauLeadPFChargedHadrCand_eta);
	ADDBRANCH(t,recTauLeadPFChargedHadrCand_phi);
	ADDBRANCH(t,recTauLeadPFChargedHadrCand_M);


    ADDBRANCH(t, chargedIsoPtSum);
    ADDBRANCH(t, neutralIsoPtSum);
    ADDBRANCH(t, puCorrPtSum);
    ADDBRANCH(t, neutralIsoPtSumWeight);
    ADDBRANCH(t, footprintCorrection);
    ADDBRANCH(t, photonPtSumOutsideSignalCone);
    ADDBRANCH(t, chargedIsoPtSumdR03);
    ADDBRANCH(t, neutralIsoPtSumdR03);
    ADDBRANCH(t, neutralIsoPtSumWeightdR03);
    ADDBRANCH(t, footprintCorrectiondR03);
    ADDBRANCH(t, photonPtSumOutsideSignalConedR03);

    ADDBRANCH(t, demetraIsolation);
}


bool ntuple_recTau::fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau){

	clear();
	if(!recTau){
		return true;
	}

	isRecTau=1;

	recTauDecayMode_i=recTau->decayMode();
	recTauDecayMode=recTau->decayMode();
	recTau_pt=recTau->pt();
	recTau_eta=recTau->eta();
	recTau_phi=recTau->phi();
	recTau_M=recTau->mass();



	recTauVtxZ = catchInfsAndBound(recTau->vertex().z(),0,-20,20);
	recImpactParamPCA_x=catchInfsAndBound(recTau->dxy_PCA().x(),0,0,1);
	recImpactParamPCA_y=catchInfsAndBound(recTau->dxy_PCA().y(),0,0,1);
	recImpactParamPCA_z=catchInfsAndBound(recTau->dxy_PCA().z(),0,-20,20);
	recImpactParam=catchInfsAndBound(recTau->dxy(),0,-5,5);
	recImpactParam3D=catchInfsAndBound(recTau->ip3d(),0,-100,100);
	recImpactParamSign=catchInfsAndBound(recTau->dxy_Sig(),0,-100,100);
	recImpactParamSign3D=catchInfsAndBound(recTau->ip3d_Sig(),0,-10000,10000);
	//
	recChi2DiffEvtVertex=0;
	//
	hasRecDecayVertex=recTau->hasSecondaryVertex();

	recDecayDist_x=catchInfsAndBound(recTau->flightLength().x(),0,-20,20);
	recDecayDist_y=catchInfsAndBound(recTau->flightLength().y(),0,-20,20);
	recDecayDist_z=catchInfsAndBound(recTau->flightLength().z(),0,-50,50);
	//missing recDecayDistCov
	recDecayDistSign=catchInfsAndBound(recTau->flightLengthSig(),0,-1000,1000);

//recTau->decayMode()
	recTauPtWeightedDetaStrip=catchInfsAndBound(pt_weighted_deta_strip(*recTau, recTau->decayMode()),0,0,1);
	recTauPtWeightedDphiStrip=catchInfsAndBound(pt_weighted_dphi_strip(*recTau, recTau->decayMode()),0,0,1);
	recTauPtWeightedDrSignal=catchInfsAndBound(pt_weighted_dr_signal(*recTau, recTau->decayMode()),0,0,1);
	recTauPtWeightedDrIsolation=catchInfsAndBound(pt_weighted_dr_iso(*recTau, recTau->decayMode()),0,0,1);
	recTauNphoton=n_photons_total(*recTau);
	recTauEratio=catchInfsAndBound(recTau->ecalEnergy()/(recTau->ecalEnergy() + recTau->hcalEnergy()),0,0,1);
	recTauLeadingTrackChi2=catchInfsAndBound(recTau->leadingTrackNormChi2(),0,0,100);
	recTauNphotonSignal=recTau->signalGammaCands().size();
	recTauNphotonIso=recTau->isolationGammaCands().size();
	recTauNphotonTotal=recTau->signalGammaCands().size()+recTau->isolationGammaCands().size();

	if(recTau->leadCand().isNonnull()){
		recTauLeadPFCand_pt=recTau->leadCand()->pt();
		recTauLeadPFCand_eta=recTau->leadCand()->eta();
		recTauLeadPFCand_phi=recTau->leadCand()->phi();
		recTauLeadPFCand_M=recTau->leadCand()->mass();
	}

	if(recTau->leadChargedHadrCand().isNonnull()){
		recTauLeadPFChargedHadrCand_pt=recTau->leadChargedHadrCand()->pt();
		recTauLeadPFChargedHadrCand_eta=recTau->leadChargedHadrCand()->eta();
		recTauLeadPFChargedHadrCand_phi=recTau->leadChargedHadrCand()->phi();
		recTauLeadPFChargedHadrCand_M=recTau->leadChargedHadrCand()->mass();
	}


    chargedIsoPtSum=recTau->tauID("chargedIsoPtSum");
    neutralIsoPtSum=recTau->tauID("neutralIsoPtSum");
    puCorrPtSum=recTau->tauID("puCorrPtSum");
    neutralIsoPtSumWeight=recTau->tauID("neutralIsoPtSumWeight");
    footprintCorrection=recTau->tauID("footprintCorrection");
    photonPtSumOutsideSignalCone=recTau->tauID("photonPtSumOutsideSignalCone");
    chargedIsoPtSumdR03=recTau->tauID("chargedIsoPtSumdR03");
    neutralIsoPtSumdR03=recTau->tauID("neutralIsoPtSumdR03");
    neutralIsoPtSumWeightdR03=recTau->tauID("neutralIsoPtSumWeightdR03");
    footprintCorrectiondR03=recTau->tauID("footprintCorrectiondR03");
    photonPtSumOutsideSignalConedR03=recTau->tauID("photonPtSumOutsideSignalConedR03");




	return true;
}




float ntuple_recTau::pt_weighted_dx(const pat::Tau& tau, int mode , int var , int decaymode )const{

	float sum_pt = 0.;
	float sum_dx_pt = 0.;
	float signalrad = std::max(0.05, std::min(0.1, 3./tau.pt()));
	int is3prong = (decaymode==10);

	auto cands = tau.isolationGammaCands();
	if(mode < 2)
		cands=tau.signalGammaCands();

	if(cands.isNull()) return 0;

	for (const auto& cand : cands) {
		// only look at electrons/photons with pT > 0.5
		if (cand->pt() < 0.5)
			continue;

		float dr = reco::deltaR(*cand, tau);
		float deta = fabs(cand->eta() - tau.eta());
		float dphi = fabs(reco::deltaPhi(cand->phi(), tau.phi()));
		float pt = cand->pt();

		bool flag = dynamic_isInside(pt, deta, dphi);

		if(is3prong==0){
			if (mode == 2 || (mode == 0 && dr < signalrad) || (mode == 1 && dr > signalrad)) {
				sum_pt += pt;
				if (var == 0)
					sum_dx_pt += pt * dr;
				else if (var == 1)
					sum_dx_pt += pt * deta;
				else if (var == 2)
					sum_dx_pt += pt * dphi;
			}
		}else if(is3prong==1){

			if( (mode==2 && flag==false) || (mode==1 && flag==true) || mode==0){
				sum_pt += pt;

				if (var == 0)
					sum_dx_pt += pt * dr;
				else if (var == 1)
					sum_dx_pt += pt * deta;
				else if (var == 2)
					sum_dx_pt += pt * dphi;
			}
		}
	}

	if (sum_pt > 0.)
		return sum_dx_pt/sum_pt;
	return 0.;
}

float ntuple_recTau::pt_weighted_dr_signal(const pat::Tau& tau, int dm)const {
	return pt_weighted_dx(tau, 0, 0, dm);
}

float ntuple_recTau::pt_weighted_deta_strip(const pat::Tau& tau, int dm)const {
	if(dm==10){
		return pt_weighted_dx(tau, 2, 1, dm);
	}else{
		return pt_weighted_dx(tau, 1, 1, dm);
	}
}

float ntuple_recTau::pt_weighted_dphi_strip(const pat::Tau& tau, int dm)const {
	if(dm==10){
		return pt_weighted_dx(tau, 2, 2, dm);
	}else{
		return pt_weighted_dx(tau, 1, 2, dm);
	}
}

float ntuple_recTau::pt_weighted_dr_iso(const pat::Tau& tau, int dm)const {
	return pt_weighted_dx(tau, 2, 0, dm);
}


unsigned int ntuple_recTau::n_photons_total(const pat::Tau& tau)const {
	unsigned int n_photons = 0;
	for (auto& cand : tau.signalGammaCands()) {
		if (cand->pt() > 0.5)
			++n_photons;
	}
	for (auto& cand : tau.isolationGammaCands()) {
		if (cand->pt() > 0.5)
			++n_photons;
	}
	return n_photons;
}


bool ntuple_recTau::dynamic_isInside(float photon_pt, float deta, float dphi)const{

	if(photon_pt==0){return false;}

	if(
			(dphi < TMath::Min(0.3, TMath::Max(0.05, 0.352476*TMath::Power(photon_pt, -0.707716)))) && \
			(deta < TMath::Min(0.15, TMath::Max(0.05, 0.197077*TMath::Power(photon_pt, -0.658701))))
	){
		return true;
	}

	return false;

}



float ntuple_recTau::calculate_demetraIsolation(const pat::Tau& tau)const{
	unsigned int tau_vertex_idxpf=-1;
	pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
	tau_vertex_idxpf = packedLeadTauCand->vertexRef().key();

	float isoDR05pt05dz015=0;
	float isoDR05pt1dz015=0;
	float gamma_DR03sum=0;

	for(const auto& IsoCand: tau.isolationChargedHadrCands()){
		pat::PackedCandidate const* cand = dynamic_cast<pat::PackedCandidate const*>(IsoCand.get());
		if (! cand->charge() )continue;
		//WATCH OUT WHICH VERTICES THESE ARE
		const auto& tau_vertex = (*vertices())[tau_vertex_idxpf];
		if ((cand->pt()<=0.5) || (cand->dxy(tau_vertex.position())>=0.1))continue;
		if (cand->hasTrackDetails()){
			const auto &tt = cand->pseudoTrack();
			if (tt.normalizedChi2()>=100. || cand->numberOfHits()<3)continue;
		}

		if (reco::deltaR2(&tau,cand)<0.5*0.5
				&& fabs(cand->dz(tau_vertex.position()))<0.15){

			isoDR05pt05dz015+=cand->pt();
			if(cand->pt()>1)
				isoDR05pt1dz015+=cand->pt();
		}
	}
	for(const auto&  IsoCand: tau.isolationGammaCands()){
		pat::PackedCandidate const* cand = dynamic_cast<pat::PackedCandidate const*>(IsoCand.get());
		if ( cand->pt() < 0.5 ) continue;
		if (reco::deltaR2(&tau,cand)<0.3*0.3
				&& cand->pt()>1.){
			gamma_DR03sum+=cand->pt();
		}
	}

	return isoDR05pt05dz015 + 0.2 * std::max(0., gamma_DR03sum - 5.);

}





