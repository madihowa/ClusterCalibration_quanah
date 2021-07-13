#include "convert_csv_ttree.C"
#include "NewClusterTree.C"
#include <TStyle.h>
#include <TCanvas.h>
#include "PlotHisto.C"

void createROOT(){
    const char* name_out = "results.root";

    if ( ! gSystem->IsFileInIncludePath(name_out) )
    {
        convert("results.csv", name_out);
    }

    NewClusterTree cluster;
    // topo.InitHist()
    cluster.Loop();
}
