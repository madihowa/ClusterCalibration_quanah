#include <TStyle.h>
#include <TSystem.h>

#include "convert_csv_ttree.C"
#include "ClusterTreeCat.C"
//#include "PlotHisto.C"

void runT() {
    const char* name_out = "results.root";

    if ( ! gSystem->IsFileInIncludePath(name_out) ) { convert("results.csv", name_out); }

    TFile* inFile = new TFile(name_out,"UPDATE");
    TTree* inTree = (TTree*)inFile->FindObjectAny("ClusterTree"); 
    if ( inTree != nullptr ) { 
      ClusterTree cluster(inTree);
      cluster.Loop();
      //    PlotHisto(gSystem->pwd());
    }
}
