/*
 * removeNonTauFraction.cc
 *
 *  Created on: 3 Aug 2018
 *      Author: jkiesele
 */


#include "../../../DeepTau/DeepTauID/interface/ntuple_content.h"
#include "../../../DeepTau/DeepTauID/interface/ntuple_JetInfo.h"
#include "../../../DeepTau/DeepTauID/interface/ntuple_genInfo.h"
#include "../../../DeepTau/DeepTauID/interface/ntuple_recTau.h"
#include "../../../DeepTau/DeepTauID/interface/ntuple_pfCands.h"

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include <iostream>

#include "TRandom3.h"

int main(int argc, char *argv[]){

	if(argc<3){
		exit(-1);
	}
	const TString infile=argv[1];
	const TString outfile=argv[2];
	float target_fraction=1;
	if(argc>3){
		target_fraction=atof(argv[3]);
	}
	std::vector<ntuple_content*> branchinfos;

  //  branchinfos.push_back(new ntuple_JetInfo());
    auto geninfo = new ntuple_genInfo();
    branchinfos.push_back(geninfo);
 //   branchinfos.push_back(new ntuple_recTau());
  //  branchinfos.push_back(new ntuple_pfCands());
  //  branchinfos.push_back(new ntuple_global());

    TFile * f = new TFile(infile,"READ");
    TTree * intree=(TTree*) f->Get("tree");

    for(auto& bi:branchinfos){
        bi->setIsRead(true);
        bi->initBranches(intree);
    }



    //count taus
    size_t ntaus=0;
    for(int i=0;i<intree->GetEntries();i++){
    	intree->GetEntry(i);
    	if(geninfo->isTau)
    		ntaus++;
    }
    float fraction = (float)ntaus/(float)intree->GetEntries();

    TRandom3 rand;

    TFile * fout= new TFile(outfile,"RECREATE");
    TTree * outtree = intree->CloneTree(0);

    for(int i=0;i<intree->GetEntries();i++){
    	intree->GetEntry(i);

    	if(geninfo->isTau)
    		outtree->Fill();
    	else if( rand.Uniform(0.,1.) < fraction/target_fraction)
    		outtree->Fill();
    }

    outtree->Write();
    //check file as ok
    fout->Close();

}
