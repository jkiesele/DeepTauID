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

	ADDBRANCH(t,recChi2DiffEvtVertex);

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
}


bool ntuple_recTau::fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau){
#warning "ntuple_recTau::fillBranches: catchInfsAndBound needs to be implemented"

	clear();
	if(!recTau){
		return true;
	}

	recTauDecayMode_i=recTau->decayMode();
	recTauDecayMode=recTau->decayMode();
	recTau_pt=recTau->pt();
	recTau_eta=recTau->eta();
	recTau_phi=recTau->phi();
	recTau_M=recTau->mass();



	recTauVtxZ = catchInfsAndBound(recTau->vertex().z(),20,-20,20);
	recImpactParamPCA_x=recTau->dxy_PCA().x();
	recImpactParamPCA_y=recTau->dxy_PCA().y();
	recImpactParamPCA_z=recTau->dxy_PCA().z();
	recImpactParam=recTau->dxy();
	recImpactParam3D=recTau->ip3d();
	recImpactParamSign=recTau->dxy_Sig();
	recImpactParamSign3D=recTau->ip3d_Sig();
	//
	recChi2DiffEvtVertex=0;
	//
	hasRecDecayVertex=recTau->hasSecondaryVertex();

	recDecayDist_x=recTau->flightLength().x();
	recDecayDist_y=recTau->flightLength().y();
	recDecayDist_z=recTau->flightLength().z();
	//missing recDecayDistCov
	recDecayDistSign=recTau->flightLengthSig();

//recTau->decayMode()
	recTauPtWeightedDetaStrip=pt_weighted_deta_strip(*recTau, recTau->decayMode());
	recTauPtWeightedDphiStrip=pt_weighted_dphi_strip(*recTau, recTau->decayMode());
	recTauPtWeightedDrSignal=pt_weighted_dr_signal(*recTau, recTau->decayMode());
	recTauPtWeightedDrIsolation=pt_weighted_dr_iso(*recTau, recTau->decayMode());
	recTauNphoton=n_photons_total(*recTau);
	recTauEratio=recTau->ecalEnergy()/(recTau->ecalEnergy() + recTau->hcalEnergy());
	recTauLeadingTrackChi2=recTau->leadingTrackNormChi2();
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


	return true;
}




float ntuple_recTau::pt_weighted_dx(const pat::Tau& tau, int mode = 0, int var = 0, int decaymode = -1)const{

	float sum_pt = 0.;
	float sum_dx_pt = 0.;
	float signalrad = std::max(0.05, std::min(0.1, 3./tau.pt()));
	int is3prong = (decaymode==10);

	const auto cands = tau.isolationGammaCands();
	if(mode < 2)
		cands=tau.signalGammaCands;


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



