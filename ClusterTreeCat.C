#define ClusterTree_cxx
#include "ClusterTreeCat.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <string>

#include <cstdio>
#include <cmath>

void ClusterTree::Loop(const std::string& outFileName) {
  // -- invalid setup
  if ( fChain == nullptr ) { printf("[ClusterTree::Loop()] ERROR No valid input tree found -> STOP\n"); }

  // -- create a new root file
  TFile *f = new TFile(outFileName.c_str(), "RECREATE");
  if ( f == nullptr ) {
    printf("[ClusterTree::Loop()] ERROR cannot open output file \042%s\042 -> STOP\n",outFileName.c_str()); 
  } else { 
    printf("[ClusterTree::Loop()] INFO  opened output file \042%s\042\n",f->GetName());
  }

  Long64_t nentries(fChain->GetEntries());
  printf("[ClusterTree::Loop()] INFO  processing %lli events\n",nentries);
 
  // -- new output trees 
  TTree *EM_tree  = new TTree("EM_tree" , "The EM tree"      );
  TTree *Had_tree = new TTree("Had_tree", "The Hadronic tree");
  Init_EM (EM_tree ); 
  Init_Had(Had_tree); 

  // -- clone input tree to output
  TTree *Cloned_tree = fChain->CloneTree();
  printf("[ClusterTree::Loop()] INFO  output trees are %s, %s and %s\n",EM_tree->GetName(),Had_tree->GetName(),Cloned_tree->GetName()); 

  // -- loop input
  Int_t ndec((Int_t)std::log10((Double_t)nentries)+1);
  Long64_t kentry(1); 
  for (Long64_t jentry=0; jentry<nentries;jentry++,kentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    Delta_E       = cluster_ENG_CALIB_TOT - CalibratedE;
    Delta_Calib_E = cluster_ENG_CALIB_TOT - clusterECalib;
    // fill original tree 
    Cloned_tree->Fill();
    // check which particle we are talking about
    if ( std::abs(truthPDG) == 211 ) { 
      // this is pi+/- -> Had
      Had_tree->Fill(); 
    } else if ( truthPDG == 111 ) { 
      // this is pi0 -> EM
      EM_tree->Fill(); 
    }
    // print a message
    Long64_t istep = kentry < 10 ? 1 : kentry < 100 ? 10 : kentry < 1000 ? 100 : kentry < 10000 ? 1000 : 10000;    
    if ( kentry % istep == 0 ) { printf("[ClusterTree::Loop()] INFO  entry %*lli/%lli\n",ndec,jentry+1,nentries); }
  }

  Cloned_tree->Write(); 
  EM_tree->Write(); 
  Had_tree->Write(); 

  printf("[ClusterTree::Loop()] INFO  wrote %lli entries for tree %s\n",Cloned_tree->GetEntries(),Cloned_tree->GetName());
  printf("[ClusterTree::Loop()] INFO  wrote %lli entries for tree %s\n",EM_tree->GetEntries()    ,EM_tree->GetName()    );
  printf("[ClusterTree::Loop()] INFO  wrote %lli entries for tree %s\n",Had_tree->GetEntries()   ,Had_tree->GetName()   );

  f->Close();
}
