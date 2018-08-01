#define signal_pu0_cxx
// The class definition in signal_pu0.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("signal_pu0.C")
// Root > T->Process("signal_pu0.C","some options")
// Root > T->Process("signal_pu0.C+")
//

#include "signal_pu0.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include "TFile.h"
#include "TTree.h"

   ///// Added ////
   static const int maxnum=100;
float pt[maxnum], eta[maxnum], phi[maxnum], dz[maxnum], dxy[maxnum], DR[maxnum];
int q[maxnum], nhits[maxnum], DR_gamma[maxnum], Npart, nmax, n, NumberOfParticles, NumberOfGammas;
float absdz[maxnum], absdxy[maxnum], w[60], weight[maxnum], cand_Iso, cand_Iso_Weight, cand_Iso_cut, cand_Iso1, cand_Iso_Weight1, cand_Iso2, cand_Iso_Weight2, cand_Iso3, cand_Iso_Weight3, cand_Iso_Weight13, gamma_Iso, gamma_Iso1pt, gamma_Iso15pt, gamma_Iso14pt, gamma_Iso13pt, gamma_Iso12pt, gamma_Iso11pt, gamma_iso1, gamma_iso2, gamma_iso3, gamma_iso4, gamma_iso5, gamma_iso45, gamma_DR04sum, gamma_DR03sum, gamma_DR02sum, gamma_DR03sum05, gamma_DR03sum15, gamma_npart, gamma_npart1, gamma_npart2, neutral_iso, NumberOfParticles2, NumberOfCharged, charged_from_PU, isoDR05pt05dz01, isoDR05pt1dz01, isoDR03pt05dz01, isoDR03pt1dz01, isoDR05pt05dz015, isoDR05pt1dz015, isoDR03pt05dz015, isoDR03pt1dz015;

//float bins1[] = {0,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.1,0.2,0.3,0.4};
//int binnum1 = sizeof(bins1)/sizeof(Float_t) - 1;
float bins2[] = {-0.02,0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40};
   int binnum2 = sizeof(bins2)/sizeof(Float_t) - 1;   
   float bins1[] = {0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40};
   int binnum1 = sizeof(bins1)/sizeof(Float_t) - 1;

//Int_t           nbins;
// Int_t           nmax;
//  Float_t         w[60];   //[nmax]
n=0;
 TFile *f1 = TFile::Open("Weight2.root", "READ");
 TTree *fr = (TTree*)f1->Get("weight_tree");
 TBranch *b_nmax=fr->GetBranch("nmax");//!
 TBranch *b_w=fr->GetBranch("w");
   
//TBranch        *b_nmax;   //!
//  TBranch        *b_w;   //!

TFile Histo_output("Histo_signal_pu200_new.root","recreate");
   
void signal_pu0::Begin(TTree * /*tree*/)
  
{
  b_nmax->SetAddress(&nmax);
  b_w->SetAddress(w);
     TTree *tree = new TTree("tree","Example of a ROOT tree");
  //TTree *tree = new TTree("tree","Example of a ROOT tree");
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
  // TFile Histo_output("Histo_signal_pu0.root","recreate");
  
   TString option = GetOption();
   tree->Branch("genTauPt", &genTauPt, "genTauPt/D");
   tree->Branch("tauPt", &tauPt, "tauPt/D");
   tree->Branch("genTauEta", &genTauEta, "genTauEta/D");
   tree->Branch("tauEta", &tauEta, "tauEta/D");
   tree->Branch("tauChargedIsoPtSum ", &tauChargedIsoPtSum, "tauChargedIsoPtSum/F");
   tree->Branch("tauNeutralIsoPtSum", &tauNeutralIsoPtSum, "tauNeutralIsoPtSum/F");
   tree->Branch("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding, "taupfTausDiscriminationByDecayModeFinding/B");
   tree->Branch("genTauMatch", &genTauMatch, "genTauMatch/I");
   tree->Branch("event", &event, "event/I");
   tree->Branch("cand_Iso", &cand_Iso, "cand_Iso/F");
   tree->Branch("cand_Iso_Weight", &cand_Iso_Weight, "cand_Iso_Weight/F");
   tree->Branch("cand_Iso1", &cand_Iso1, "cand_Iso1/F");
   tree->Branch("cand_Iso_Weight1", &cand_Iso_Weight1, "cand_Iso_Weight1/F");
   tree->Branch("cand_Iso2", &cand_Iso2, "cand_Iso2/F");
   tree->Branch("cand_Iso_Weight2", &cand_Iso_Weight2, "cand_Iso_Weight2/F");
   tree->Branch("cand_Iso_Weight3", &cand_Iso_Weight3, "cand_Iso_Weight3/F");
   tree->Branch("cand_Iso3", &cand_Iso3, "cand_Iso3/F");
   tree->Branch("cand_Iso_Weight13", &cand_Iso_Weight13, "cand_Iso_Weight13/F");
   tree->Branch("tauPuCorrPtSum", &tauPuCorrPtSum, "tauPuCorrPtSum/F");
   tree->Branch("cand_Iso_cut", &cand_Iso_cut, "cand_Iso_cut/F");
   tree->Branch("gamma_Iso", &gamma_Iso, "gamma_Iso/F");
   tree->Branch("gamma_Iso1pt", &gamma_Iso1pt, "gamma_Iso1pt/F");
   tree->Branch("gamma_Iso15pt", &gamma_Iso15pt, "gamma_Iso15pt/F");
   tree->Branch("gamma_Iso14pt", &gamma_Iso14pt, "gamma_Iso14pt/F");   
   tree->Branch("gamma_Iso13pt", &gamma_Iso13pt, "gamma_Iso13pt/F");
   tree->Branch("gamma_Iso12pt", &gamma_Iso12pt, "gamma_Iso12pt/F");
   tree->Branch("gamma_Iso11pt", &gamma_Iso11pt, "gamma_Iso11pt/F");    
   tree->Branch("Npart", &Npart, "Npart/I");
   tree->Branch("pt", pt, "pt[Npart]/F");
   tree->Branch("eta", eta, "eta[Npart]/F");
   tree->Branch("phi", phi, "phi[Npart]/F");
   tree->Branch("dz", dz, "dz[Npart]/F");
   tree->Branch("dxy", dxy, "dxy[Npart]/F");
   tree->Branch("DR", DR, "DR[Npart]/F");
   tree->Branch("q", q, "q[Npart]/I");
   tree->Branch("nhits", nhits, "nhits[Npart]/I");
   tree->Branch("gamma_iso1", &gamma_iso1, "gamma_iso1/F");
   tree->Branch("gamma_iso2", &gamma_iso2, "gamma_iso2/F");
   tree->Branch("gamma_iso3", &gamma_iso3, "gamma_iso3/F");
   tree->Branch("gamma_iso4", &gamma_iso4, "gamma_iso4/F");
   tree->Branch("gamma_iso5", &gamma_iso5, "gamma_iso5/F");
   tree->Branch("gamma_iso45", &gamma_iso45, "gamma_iso45/F");
   tree->Branch("gamma_DR04sum", &gamma_DR04sum, "gamma_DR04sum/F");
   tree->Branch("gamma_DR03sum", &gamma_DR03sum, "gamma_DR03sum/F");   
   tree->Branch("gamma_DR02sum", &gamma_DR02sum, "gamma_DR02sum/F");
   tree->Branch("gamma_DR03sum05", &gamma_DR03sum05, "gamma_DR03sum05/F");
   tree->Branch("gamma_DR03sum15", &gamma_DR03sum15, "gamma_DR03sum15/F");    
   tree->Branch("gamma_npart", &gamma_npart, "gamma_npart/F");
   tree->Branch("gamma_npart1", &gamma_npart1, "gamma_npart1/F");
   tree->Branch("gamma_npart2", &gamma_npart2, "gamma_npart2/F");
   tree->Branch("charged_from_PU", &charged_from_PU, "charged_from_PU/F");
   tree->Branch("isoDR05pt05dz01" , &isoDR05pt05dz01, "isoDR05pt05dz01/F");
   tree->Branch("isoDR05pt1dz01" , &isoDR05pt1dz01, "isoDR05pt1dz01/F");
   tree->Branch("isoDR03pt05dz01" , &isoDR03pt05dz01, "isoDR03pt05dz01/F");
   tree->Branch("isoDR03pt1dz01" , &isoDR03pt1dz01, "isoDR03pt1dz01/F");
   tree->Branch("isoDR05pt05dz015" , &isoDR05pt05dz015, "isoDR05pt05dz015/F");
   tree->Branch("isoDR05pt1dz015" , &isoDR05pt1dz015, "isoDR05pt1dz015/F");
   tree->Branch("isoDR03pt05dz015" , &isoDR03pt05dz015, "isoDR03pt05dz015/F");
   tree->Branch("isoDR03pt1dz015" , &isoDR03pt1dz015, "isoDR03pt1dz015/F");
   
   // TBranch        *b_nmax;   //!



   //TBranch        *b_nbins;   //!

   //TBranch        *b_w;   //!
   
  //weight_tree->SetBranchAddress("nbins", &nbins, &b_nbins);

   //cout <<"nmax =  "<<nmax<<endl; 
   //float bins0[] = {0,20,40,60,80,100,120,160,200};
   float bins0[] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200};
   int binnum0 = sizeof(bins0)/sizeof(Float_t) - 1;
   float bins2[] = {-0.02,0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40};
   int binnum2 = sizeof(bins2)/sizeof(Float_t) - 1;
   float bins1[] = {0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40};
   int binnum1 = sizeof(bins1)/sizeof(Float_t) - 1;
   
   TH1F *histo_gen_pt=new TH1F("histo_gen_pt","", binnum0, bins0);
   TH1F *histo_rec_pt=new TH1F("histo_rec_pt","", binnum0, bins0);
   TH1F *histo_gen_eta=new TH1F("histo_gen_eta","gen_eta",16,-4,4);
   TH1F *histo_rec_eta=new TH1F("histo_rec_eta","rec_eta",16,-4,4);
   TH1F *histo_DM=new TH1F("histo_DM","DM",30,-15,15);
   TH1F *histo_iso=new TH1F("histo_iso","iso",10, 0, 5);
   TH1F *histo_iso_norm=new TH1F("histo_iso_norm","",100,0,100);
   TH1F *histo_gamma=new TH1F("histo_gamma","",10,0,5);
   TH1F *histo_isoOverpt1=new TH1F("histo_isoOverpt1","isoOverpt",40,0,2);
   TH1F *histo_isoOverpt2=new TH1F("histo_isoOverpt2","isoOverpt",10,0,0.4);
   THStack *hs1 = new THStack("hs","");
   THStack *hs2 = new THStack("hs","");
   THStack *hs3 = new THStack("hs","");
   THStack *hs4 = new THStack("hs","");
   TH1F *histo_cand_pt=new TH1F("histo_cand_pt","", 20, 0, 10);// binnum0, bins0);
   TH1F *histo_cand_eta=new TH1F("histo_cand_eta","", 18, -4.5, 4.5);
   TH1F *histo_cand_phi=new TH1F("histo_cand_phi","", 14, -3.5, 3.5);
   TH1F *histo_cand_numHits=new TH1F("histo_cand_numHits","", 40, 0, 40);
   TH1F *histo_cand_dz=new TH1F("histo_cand_dz","",binnum1, bins1);//40, -1, 1);// 80, -0.4, 0.4);
   TH1F *histo_new_dz=new TH1F("histo_new_dz","", binnum2, bins2);   
   TH1F *histo_cand_dxy=new TH1F("histo_cand_dxy","", 20, -0.2, 0.2 );
   TH1F *histo_cand_DR=new TH1F("histo_cand_DR","", 30, 0, 0.6);
   TH1F *hiso_cand_Iso=new TH1F("histo_cand_Iso","",10, 0, 5);
   TH1F *hiso_cand_Iso_weight=new TH1F("histo_cand_Iso_weight","",10, 0, 5);
   TH1F *hiso_cand_Iso1=new TH1F("histo_cand_Iso1","",10, 0, 5);
   TH1F *hiso_cand_Iso_weight1=new TH1F("histo_cand_Iso_weight1","",10, 0, 5);
   TH1F *hiso_cand_Iso3=new TH1F("histo_cand_Iso3","",10, 0, 5);   
   TH1F *hiso_cand_Iso_weight3=new TH1F("histo_cand_Iso_weight3","",10, 0, 5);
   TH1F *hiso_cand_Iso_weight13=new TH1F("histo_cand_Iso_weight13","",10, 0, 5);
   TH1F *hiso_cand_Iso_cut=new TH1F("histo_cand_Iso_cut","",10, 0, 5);
   TH1F *histo_gamma_Iso=new TH1F("histo_gamma_Iso","",100,0,100);
   TH1F *histo_gamma_Iso1pt=new TH1F("histo_gamma_Iso1pt","",60,0,60);
   TH1F *histo_gamma_Iso15pt=new TH1F("histo_gamma_Iso15pt","",60,0,60);   
   TH1F *histo_gamma_DR04sum=new TH1F("histo_gamma_DR04sum","",60,0,60);
   TH1F *histo_gamma_DR03sum=new TH1F("histo_gamma_DR03sum","",60,0,60);
   TH1F *histo_gamma_DR02sum=new TH1F("histo_gamma_DR02sum","",60,0,60);
   TH1F *histo_gamma_DR03sum05=new TH1F("histo_gamma_DR03sum05","",60,0,60);
   TH1F *histo_gamma_DR02sum15=new TH1F("histo_gamma_DR03sum15","",60,0,60);   
   TH1F *histo_gamma_iso1=new TH1F("histo_gamma_iso1","",10,0,5);
   TH1F *histo_gamma_iso2=new TH1F("histo_gamma_iso2","",10,0,5);
   TH1F *histo_gamma_iso3=new TH1F("histo_gamma_iso3","",10,0,5);
   TH1F *histo_gamma_iso4=new TH1F("histo_gamma_iso4","",10,0,5);
   TH1F *histo_gamma_iso5=new TH1F("histo_gamma_iso5","",10,0,5);
   TH1F *histo_gamma_iso45=new TH1F("histo_gamma_iso45","",10,0,5);
   TH1F *histo_gamma_pt=new TH1F("histo_gamma_pt","", 20, 0, 10);// binnum0, bins0);
   TH1F *histo_gamma_eta=new TH1F("histo_gamma_eta","", 18, -4.5, 4.5);
   TH1F *histo_gamma_phi=new TH1F("histo_gamma_phi","", 18, -4.5, 4.5);  
   TH1F *histo_gamma_DR=new TH1F("histo_gamma_DR","", 30, 0, 0.6);   
   TH1F *histo_Npart= new TH1F("histo_Npart", "" ,20,0,20);
   TH1F *histo_Npart_new= new TH1F("histo_Npart_new", "" ,20,0,20);
   TH1F *histo_Ngammas= new TH1F("histo_Ngammas", "" ,40,0,40);   
   TH1F *histo_Ngammas_new= new TH1F("histo_Ngammas_new", "" ,40,0,40);
   TH1F *histo_DR_gammas_new=new TH1F("histo_DR_gammas_new","", 30, 0, 0.6);
   TH1F *histo_gammas_pt_new=new TH1F("histo_gammas_pt_new","", 20, 0, 10);
   TH1F *histo_neutral_iso=new TH1F("histo_neutral_iso","",50,0,50);
   TH1F *histo_neutral_iso1=new TH1F("histo_neutral_iso1","",50,0,50);
   TH1F *histo_neutral_iso2=new TH1F("histo_neutral_iso2","",50,0,50);   
   TH1F *histo_numberOfParticles=new TH1F("histo_numberOfParticles","",20,0,20);
   TH1F *histo_Ncharged=new TH1F("histo_Ncharged","",20,0,20);
   TH1F *histo_Ncharged_new=new TH1F("histo_Ncharged_new","",20,0,20);
   TH1F *histo_charged_from_PU=new TH1F("histo_charged_from_PU","",60,0,60);
   TH1F *histo_isoDR05pt05dz01=new TH1F("histo_isoDR05pt05dz01","",10,0,5);
   TH1F *histo_isoDR05pt1dz01=new TH1F("histo_isoDR05pt1dz01","",10,0,5);   
   TH1F *histo_isoDR03pt05dz01=new TH1F("histo_isoDR03pt05dz01","",10,0,5);
   TH1F *histo_isoDR03pt1dz01=new TH1F("histo_isoDR03pt1dz01","",10,0,5);
   TH1F *histo_isoDR05pt05dz015=new TH1F("histo_isoDR05pt05dz015","",10,0,5);
   TH1F *histo_isoDR05pt1dz015=new TH1F("histo_isoDR05pt1dz015","",10,0,5);   
   TH1F *histo_isoDR03pt05dz015=new TH1F("histo_isoDR03pt05dz015","",10,0,5);
   TH1F *histo_isoDR03pt1dz015=new TH1F("histo_isoDR03pt1dz015","",10,0,5);
   
   //TGraphAsymmErrors* Eff_tau = new TGraphAsymmErrors();
   //TGraphAsymmErrors* Eff_eta = new TGraphAsymmErrors();
   
}

void signal_pu0::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   
   TString option = GetOption();


}

   ///////Fin added/////
Bool_t signal_pu0::Process(Long64_t entry)
{
  //fr->SetBranchAddress("nmax", &nmax, &b_nmax);
  //fr->SetBranchAddress("w", w, &b_w[nmax]);
   Int_t nentries_wtn = fChain->GetTree()->GetEntry(entry);
        bool genTau = genTauPt > 20 && genTauEta < 3.0 &&  genTauEta > -3.0 ;//  && (dm!=5&&dm!=6 && dm > -1);
     bool recTau = tauPt > 20 && tauEta < 3.0 && tauEta > -3.0 && taupfTausDiscriminationByDecayModeFinding==1 && genTauMatch == 1 && genTau;// tauChargedIsoPtSum<2.5 && genTau;
     float var_1 = tauChargedIsoPtSum/tauPt;
     float var_2 = tauChargedIsoPtSum/genTauPt;
   
     //cout<<"number of entries =" <<  nentries_wtn << endl;
     //for (Int_t i = 0; i < nentries_wtn; i++){

     /*
     */

     //cout<<"Nparticles ="<<Nparticles<<endl;

     
     /*  histo_tau_iso_pt01->Fill(tau_iso_dz01);
       histo_tau_iso_pt0015->Fill(tau_iso_dz0015);
       histo_tau_iso_dz01->Fill(tau_iso_dz01_dz);
       histo_tau_iso_dz0015->Fill(tau_iso_dz0015_dz);
       histo_tau_iso_dxy01->Fill(tau_iso_dz01_dxy);
       histo_tau_iso_dxy0015->Fill(tau_iso_dz0015_dxy);
     */
	 
     if (genTau){

	 histo_gen_pt->Fill(genTauPt);
	 histo_gen_eta->Fill(genTauEta);
	 histo_isoOverpt2->Fill(var_2);

       
       //cout<<"gettaupt ="<<genTauEta<<endl;
     }
      
     if (recTau){
       histo_rec_pt->Fill(genTauPt);
       histo_rec_eta->Fill(genTauEta);
       histo_iso->Fill(tauChargedIsoPtSum);
       histo_iso_norm->Fill(tauChargedIsoPtSum);       
       histo_gamma->Fill(tauNeutralIsoPtSum);
       histo_isoOverpt1->Fill(var_1);
       histo_DM->Fill(dm);
       histo_Npart->Fill(Nparticles);
       histo_Ngammas->Fill(Ngammas);
       //cout<<"Nparticles =" << Nparticles<<endl;
       Npart=Nparticles;
       //float weght =0;
       //n=0;	      
       int nbins=fr->GetEntries();
       // cout<<"nbins =" << nbins <<endl;
       // cand_Iso[0]=0;
       cand_Iso1=0;
       cand_Iso2=0;
       cand_Iso_Weight1=0;
       cand_Iso_Weight2=0;       
       cand_Iso=0;
       cand_Iso_Weight=0;
       cand_Iso_cut=0;
       cand_Iso3 = 0;
       cand_Iso_Weight3=0;
       cand_Iso_Weight13=0;
       gamma_Iso=0;
       gamma_Iso1pt=0;
       gamma_Iso15pt=0;
       gamma_Iso14pt=0;
       gamma_Iso13pt=0;
       gamma_Iso12pt=0;
       gamma_Iso11pt=0;
       
       gamma_iso1=0;
       gamma_iso2=0;
       gamma_iso3=0;
       gamma_iso4=0;
       gamma_iso5=0;
       gamma_iso45=0;
       gamma_DR04sum=0;
       gamma_DR03sum=0;
       gamma_DR02sum=0;
       gamma_DR03sum05=0;
       gamma_DR03sum15=0;
       gamma_npart2=0;
       gamma_npart=0;
       gamma_npart1=0;
       neutral_iso=0;
       NumberOfParticles=0;
       NumberOfCharged=0;
       charged_from_PU=0;
       isoDR05pt05dz01=0;
       isoDR05pt1dz01=0;
       isoDR03pt05dz01=0;
       isoDR03pt1dz01=0;
       isoDR05pt05dz015=0;
       isoDR05pt1dz015=0;
       isoDR03pt05dz015=0;
       isoDR03pt1dz015=0;
       
      if (Nparticles==0){
	//cout<<"Nparticles =" << Nparticles<<endl;	
	 histo_new_dz->Fill(-0.01);
       }
       for(int i=0; i<Nparticles; i++){
       absdz[i]=cand_dz_new[i];
       if(absdz[i]<0){
	 absdz[i]=-1*absdz[i];
       }
       else {
       	 absdz[i]=absdz[i];
         }

       absdxy[i]=cand_dxy_new[i];
       if (absdxy[i]<0) absdxy[i]=-1*absdxy[i];
       
       //for(int i=0; i<Nparticles; i++){
       // float weight=0.;
       for (int k=0; k<binnum1; k++){
       b_nmax->GetEntry(k);
       b_w->GetEntry(k);
       // cout<<"nbins =   "<<nbins<< "B_W  =   " << w[nbins]<<endl;
       // weight[0]=0.;
        if (absdz[i]>=bins1[k]){
	  
         weight[i] =w[k+1];
	 //cout<< "weight[i] = " <<weight[i]<<endl;
       	 }
	if (absdz[i]>0.5) {
	  //cout<<"Hello"<<endl;
	  weight[i] =0.;
	}
       }

       if (absdz[i]<0.1 && absdxy[i]<0.03){
	cand_Iso_cut+=cand_pt[i];
       }
       
       if (absdz[i]<0.15 && absdxy[i]<0.05){
	 NumberOfParticles=NumberOfParticles+1;  ///////////////////************number of particles passes the definitions //////////////////
	 //histo_numberOfParticles->Fill(Nparticles); 
	 //cand_Iso_cut+=cand_pt[i];
	NumberOfCharged=NumberOfCharged+weight[i];
        cand_Iso+=cand_pt[i];
        cand_Iso_Weight+=cand_pt[i]*weight[i];	
        cand_Iso_Weight1+=cand_pt[i]*weight[i]/0.7755541; ///0.772168;
       }
       
       
       if (absdz[i]<0.15 && absdxy[i]<0.05){
	 //cand_Iso_cut+=cand_pt[i];
        cand_Iso3+=cand_pt[i];
        cand_Iso_Weight3+=cand_pt[i]*weight[i];	
        cand_Iso_Weight13+=cand_pt[i]*weight[i]*1.778122;
       }
       
       if (absdz[i]>0.15 && absdz[i]<0.4 && absdxy[i]<0.05){
	 cand_Iso2+=cand_pt[i];
         //cand_Iso_Weight2+=cand_pt[i]*weight[i];
       }

       if (DR_IsoCandTau[i]<0.5 && absdz[i]>0.15 && absdz[i]<0.4 && cand_pt[i]>1){
	 charged_from_PU+=cand_pt[i];
       }

       if (DR_IsoCandTau[i]<0.5 && absdz[i]<0.1 && cand_pt[i]>0.5){
	 isoDR05pt05dz01+=cand_pt[i];
       }
       if (DR_IsoCandTau[i]<0.5 && absdz[i]<0.1 && cand_pt[i]>1){
	 isoDR05pt1dz01+=cand_pt[i];
       }
       if (DR_IsoCandTau[i]<0.3 && absdz[i]<0.1 && cand_pt[i]>0.5){
	 isoDR03pt05dz01+=cand_pt[i];
       }
       if (DR_IsoCandTau[i]<0.3 && absdz[i]<0.1 && cand_pt[i]>1){
	 isoDR03pt1dz01+=cand_pt[i];
       }
       if (DR_IsoCandTau[i]<0.5 && absdz[i]<0.15 && cand_pt[i]>0.5){
	 isoDR05pt05dz015+=cand_pt[i];
       }
       if (DR_IsoCandTau[i]<0.5 && absdz[i]<0.15 && cand_pt[i]>1){
	 isoDR05pt1dz015+=cand_pt[i];
       }
       if (DR_IsoCandTau[i]<0.3 && absdz[i]<0.15 && cand_pt[i]>0.5){
	 isoDR03pt05dz015+=cand_pt[i];
       }
       if (DR_IsoCandTau[i]<0.3 && absdz[i]<0.15 && cand_pt[i]>1){
	 isoDR03pt1dz015+=cand_pt[i];
       }       
	 
       cand_Iso1=cand_Iso-0.6*cand_Iso2;
       //cout<< "cand_pt[i] = "<< cand_pt[i] <<endl;
       //cout <<  "weight[i] =" << weight[i] <<endl;
       histo_cand_pt->Fill(cand_pt[i],weight[i]);//w[i]);
       histo_cand_eta->Fill(cand_eta[i],weight[i]);//,weight);//,w[i]);
       histo_cand_phi->Fill(cand_phi[i],weight[i]);//,weight);//,w[i]);
       histo_cand_numHits->Fill(cand_numHits[i],weight[i]);//,weight);//,w[i]);
       histo_cand_dz->Fill(absdz[i],weight[i]);//,weight);//cand_dz_new[i]);//,w[i]);
       histo_new_dz->Fill(absdz[i]);//       
       histo_cand_dxy->Fill(cand_dxy_new[i],weight[i]);//,weight);//,w[i]);
       histo_cand_DR->Fill(DR_IsoCandTau[i],weight[i]);//,weight);//,w[i]);
       pt[i]=cand_pt[i];
     eta[i]=cand_eta[i];
     phi[i]=cand_phi[i];
     nhits[i]=cand_numHits[i];
     dz[i]=cand_dz_new[i];
     dxy[i]=cand_dxy_new[i];
     DR[i]=DR_IsoCandTau[i];

     //}
       }
       histo_numberOfParticles->Fill(NumberOfParticles); ///////////////////////////////////////Filling the histo with Number of particles passes the previous cuts /////////////
       histo_Ncharged->Fill(NumberOfCharged);      

       if (cand_Iso1 <= 0){
	 cand_Iso1 = 0.;
	 // cout <<"Hello cand_Iso"<<cand_Iso<<endl;
       }
	 /*
       if (cand_Iso_Weight1 <= 0){
	 cand_Iso_Weight1= 0.;
	 //cout <<"Hello cand_Iso"<<cand_Iso<<endl;
       }       
       
	 */
       histo_cand_Iso->Fill(cand_Iso);
       histo_cand_Iso_weight->Fill(cand_Iso_Weight);
       histo_cand_Iso_cut->Fill(cand_Iso_cut);
       histo_cand_Iso1->Fill(cand_Iso1);
       histo_cand_Iso_weight1->Fill(cand_Iso_Weight1);
       histo_cand_Iso3->Fill(cand_Iso3);       
       histo_cand_Iso_weight3->Fill(cand_Iso_Weight3);
       histo_cand_Iso_weight13->Fill(cand_Iso_Weight13);
       histo_charged_from_PU->Fill(charged_from_PU);
       histo_isoDR05pt05dz01->Fill(isoDR05pt05dz01);
       histo_isoDR05pt1dz01->Fill(isoDR05pt1dz01);
       histo_isoDR03pt05dz01->Fill(isoDR03pt05dz01);
       histo_isoDR03pt1dz01->Fill(isoDR03pt1dz01);       
       histo_isoDR05pt05dz015->Fill(isoDR05pt05dz015);
       histo_isoDR05pt1dz015->Fill(isoDR05pt1dz015);
       histo_isoDR03pt05dz015->Fill(isoDR03pt05dz015);
       histo_isoDR03pt1dz015->Fill(isoDR03pt1dz015);
       
       for(int kk=0; kk<Ngammas; kk++){
	 
	 //if (kk>19)break;
	 if  (DR_gamma[kk]<0.5 && gamma_pt[kk]>0.5){
	   NumberOfGammas=NumberOfGammas+1;
	   gamma_Iso+=gamma_pt[kk];
	 }
       
	 if  (DR_gamma[kk]<0.5 && gamma_pt[kk]>1.0){
	   gamma_Iso1pt+=gamma_pt[kk];
	   //cout<<"pt ="<<gamma_pt[kk]<< "DR ="<< DR_gamma[kk]<<endl;
	 }
	 if  (DR_gamma[kk]<0.5 && gamma_pt[kk]>1.5) gamma_Iso15pt+=gamma_pt[kk];
	 if  (DR_gamma[kk]<0.3 && gamma_pt[kk]>1.1) gamma_Iso11pt+=gamma_pt[kk];
	 if  (DR_gamma[kk]<0.3 && gamma_pt[kk]>1.2) gamma_Iso12pt+=gamma_pt[kk];
	 if  (DR_gamma[kk]<0.3 && gamma_pt[kk]>1.3) gamma_Iso13pt+=gamma_pt[kk];
	 if  (DR_gamma[kk]<0.3 && gamma_pt[kk]>1.4) gamma_Iso14pt+=gamma_pt[kk];	 
	 
	 
	 if  (DR_gamma[kk]<0.4 && gamma_pt[kk]>1.0){
	   gamma_DR04sum+=gamma_pt[kk];
	 }	 
	 if  (DR_gamma[kk]<0.3 && gamma_pt[kk]>1.0){
	   gamma_DR03sum+=gamma_pt[kk];
	 }
	 
	 if  (DR_gamma[kk]<0.2 && gamma_pt[kk]>1.0) gamma_DR02sum+=gamma_pt[kk];
	 if  (DR_gamma[kk]<0.3 && gamma_pt[kk]>0.5) gamma_DR03sum05+=gamma_pt[kk];
	 if  (DR_gamma[kk]<0.3 && gamma_pt[kk]>1.5) gamma_DR03sum15+=gamma_pt[kk];

        

       histo_gamma_pt->Fill(gamma_pt[kk]);
       histo_gamma_eta->Fill(gamma_eta[kk]);
       histo_gamma_phi->Fill(gamma_phi[kk]);
       histo_gamma_DR->Fill(DR_gamma[kk]);
       }
       //histo_Ngammas->Fill(NumberOfGammas);
       gamma_iso1=gamma_Iso-0.1*tauPuCorrPtSum;
       gamma_iso2=gamma_Iso-0.2*tauPuCorrPtSum;
       gamma_iso3=gamma_Iso-0.3*tauPuCorrPtSum;
       gamma_iso4=gamma_Iso-0.4*tauPuCorrPtSum;
       gamma_iso5=gamma_Iso-0.5*tauPuCorrPtSum;
       gamma_iso45=gamma_Iso-0.45*tauPuCorrPtSum;
       
       histo_gamma_Iso->Fill(gamma_Iso);
       histo_gamma_Iso1pt->Fill(gamma_Iso1pt);
       histo_gamma_Iso15pt->Fill(gamma_Iso15pt);      
       histo_gamma_DR04sum->Fill(gamma_DR04sum);
       histo_gamma_DR03sum->Fill(gamma_DR03sum);       
       histo_gamma_DR02sum->Fill(gamma_DR02sum);
       histo_gamma_DR03sum05->Fill(gamma_DR03sum05);
       histo_gamma_DR03sum15->Fill(gamma_DR03sum15);	 
       histo_gamma_iso1->Fill(gamma_iso1);	
       histo_gamma_iso2->Fill(gamma_iso2);
       histo_gamma_iso3->Fill(gamma_iso3);
       histo_gamma_iso4->Fill(gamma_iso4);	
       histo_gamma_iso5->Fill(gamma_iso5);
       histo_gamma_iso45->Fill(gamma_iso45);

/////****************** New part added in the code *********************************

       
       if (cand_Iso_Weight1<4){

	 
	 /*	 for (int mm=0; mm<Nparticles; mm++){
           if (absdz[mm]<0.15 && absdxy[mm]<0.05){
	       NumberOfParticles2=Nparticles;
	       //histo_Npart_new->Fill(Nparticles); 
	   }
	   //histo_Npart_new->Fill(NumberOfParticles2);
	   // histo_Npart_new->Fill(Nparticles); 
	   }*/
	 histo_Npart_new->Fill(NumberOfParticles);
	 histo_Ncharged_new->Fill(NumberOfCharged);
         histo_Ngammas_new->Fill(Ngammas);
	 
	 for (int qq=0; qq<Ngammas; qq++){
           histo_DR_gammas_new->Fill(DR_gamma[qq]);
	   histo_gammas_pt_new->Fill(gamma_pt[qq]);
	   if (DR_gamma[qq]<0.5 && gamma_pt>0.5) neutral_iso+=gamma_pt[qq];
	 }
	 //histo_neutral_iso->Fill(neutral_iso);

       }

//// **************************************************************************************

       if (Ngammas>=NumberOfParticles){
	int Np=NumberOfParticles;
       for(int jj=0; jj<Np; jj++){
	 //if (DR_gamma[jj]<0.5 && gamma_pt[jj]>1.0) gamma_npart+=gamma_pt[jj];
	  if (DR_gamma[jj]<0.5 && gamma_pt[jj]>0.5) gamma_npart1+=gamma_pt[jj];	
       }
       }
       
       if (Ngammas>=(int)NumberOfCharged){
	 int Ng=(int)NumberOfCharged;
       for(int t=0; t<Ng; t++){
	  if (DR_gamma[t]<0.5 && gamma_pt[t]>0.5) gamma_npart+=gamma_pt[t];
	 //if (DR_gamma[jj]<0.5 && gamma_pt[jj]>0.5) gamma_npart1+=gamma_pt[jj];	
       }
       }

       
       if (Ngammas>=5){
	int N=5;
       for(int tt=0; tt<N; tt++){
	 //if (DR_gamma[jj]<0.5 && gamma_pt[jj]>1.0) gamma_npart+=gamma_pt[jj];
	  if (DR_gamma[tt]<0.5 && gamma_pt[tt]>0.5) gamma_npart2+=gamma_pt[tt];	
       }
       }

       histo_neutral_iso2->Fill(gamma_npart2);       
       histo_neutral_iso1->Fill(gamma_npart1);
       histo_neutral_iso->Fill(gamma_npart);
     }

     tree->Fill();
     // tree->Write();
   }

//histo_Npart->Fill(Nparticles);
   return kTRUE;
   //}
//cout<<"tree has %i entries\n"<<tree_test->GetEntries()<<endl;
//printf("tree has %i entries\n",tree_test->GetEntries());
void signal_pu0::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
  //Int_t nentries_wtn = fChain->GetTree()->GetEntry(entry);
  //for (Int_t i = 0; i < nentries_wtn; i++){  


    


}

void signal_pu0::Terminate()
{
  //TFile *rf = new TFile("rootfile.root", "RECREATE");
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
    TCanvas *c1 = new TCanvas("c","c",800,800);
  histo_rec_pt->SetFillColor(kRed);
  histo_gen_pt->SetFillColor(kBlue);
  histo_rec_eta->SetFillColor(kRed);
  histo_gen_eta->SetFillColor(kBlue);
  histo_iso->SetFillColor(kGreen);
  histo_isoOverpt1->SetFillColor(kCyan);
  histo_isoOverpt1->SetFillColor(kCyan-3);
  histo_cand_pt->SetFillColor(kBlue-4);
  histo_cand_eta->SetFillColor(kRed-2);
  histo_cand_phi->SetFillColor(kRed);
  histo_cand_numHits->SetFillColor(kRed);
  histo_cand_dz->SetFillColor(kRed);
  histo_cand_dxy->SetFillColor(kRed);
  histo_cand_DR->SetFillColor(kRed-4);
  histo_cand_Iso->SetFillColor(kBlue-4);
  histo_cand_Iso_weight->SetFillColor(kBlue-2);
  histo_cand_Iso_cut->SetFillColor(kBlue);
  histo_Npart->SetFillColor(kPink);
  histo_new_dz->SetFillColor(kPink);
  histo_cand_Iso1->SetFillColor(kBlue-4);
  histo_cand_Iso_weight1->SetFillColor(kBlue-2);  
  c1->Divide(4,4);

  c1->cd(1);
  histo_rec_pt->Draw("histobar");
  c1->cd(2);
  histo_gen_pt->Draw("histobar");
  c1->cd(3);
  histo_rec_eta->Draw("histobar");
  c1->cd(4);
  histo_gen_eta->Draw("histobar");
  c1->cd(5);
  histo_iso->Draw("histobar");
  c1->cd(6);
  histo_Npart->Draw("histobar");
  c1->cd(7);
  histo_cand_Iso->Draw("histobar");
  c1->cd(8);
  histo_cand_Iso_cut->Draw("histobar");//  histo_DM->Draw("histobar");
  c1->cd(9);
  histo_cand_pt->Draw("histobar");
  c1->cd(10);
  histo_cand_eta->Draw("histobar");
  c1->cd(11);
  histo_cand_phi->Draw("histobar");
  c1->cd(12); 
  histo_cand_numHits->Draw("histobar");
  c1->cd(13);
  histo_cand_dz->Draw("histobar");
  c1->cd(14); 
  histo_cand_dxy->Draw("histobar");  
  c1->cd(15);
  histo_cand_DR->Draw("histobar");
  c1->cd(16);
  histo_cand_Iso_weight->Draw("histobar");
  
 

   histo_rec_pt->Write();
   histo_gen_pt->Write();
   histo_rec_eta->Write();
   histo_gen_eta->Write();
   histo_iso->Write();
   histo_iso_norm->Write();
   histo_isoOverpt1->Write();
   histo_cand_DR->Write();
   histo_DM->Write();
   histo_cand_pt->Write();
   histo_cand_eta->Write();
   histo_cand_phi->Write();
   histo_cand_numHits->Write();
   histo_cand_dz->Write();
   histo_cand_dxy->Write();
   histo_cand_Iso->Write();
   histo_cand_Iso_weight->Write();
   histo_cand_Iso1->Write();
   histo_cand_Iso_weight1->Write();
   histo_cand_Iso3->Write();
   histo_cand_Iso_weight3->Write();
   histo_cand_Iso_weight13->Write();   
   histo_cand_Iso_cut->Write();
   histo_Npart->Write();
   histo_new_dz->Write();
   histo_gamma->Write();
   histo_gamma_Iso->Write();
   histo_gamma_iso1->Write();
   histo_gamma_iso2->Write();
   histo_gamma_iso3->Write();
   histo_gamma_iso4->Write();
   histo_gamma_iso5->Write();
   histo_gamma_iso45->Write();
   histo_gamma_pt->Write();
   histo_gamma_eta->Write();
   histo_gamma_phi->Write();
   histo_gamma_DR->Write();
   histo_gamma_DR04sum->Write();
   histo_gamma_DR03sum->Write();
   histo_gamma_DR02sum->Write();
   histo_gamma_DR03sum05->Write();
   histo_gamma_DR03sum15->Write();
   histo_gamma_Iso1pt->Write();
   histo_gamma_Iso15pt->Write();
   histo_Npart_new->Write();
   histo_Ngammas->Write();
   histo_Ngammas_new->Write();
   histo_DR_gammas_new->Write();
   histo_gammas_pt_new->Write();
   histo_neutral_iso->Write();
   histo_neutral_iso1->Write();
   histo_neutral_iso2->Write();   
   histo_numberOfParticles->Write();
   histo_Ncharged->Write();
   histo_Ncharged_new->Write();
   histo_charged_from_PU->Write();
   histo_isoDR05pt05dz01->Write();
   histo_isoDR05pt1dz01->Write();
   histo_isoDR03pt05dz01->Write();
   histo_isoDR03pt1dz01->Write();
   histo_isoDR05pt05dz015->Write();
   histo_isoDR05pt1dz015->Write();
   histo_isoDR03pt05dz015->Write();
   histo_isoDR03pt1dz015->Write();
   
   tree->Write();

   //TFile rf("rootfile.root", "recreate");

   //   rf.Close();
   
}
   
     Histo_output.Close();
