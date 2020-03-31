#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "TRint.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TH1.h"
#include "TString.h"

void hist2TKA(TH1F* hist){
   
    int num_bins = hist->GetNbinsX();
    TString name = "data_TKA/";
    name.Append(hist->GetName());
    name.Append(".TKA");
   
    ofstream TKAfile;
    TKAfile.open(name.Data());
    for (int i=0; i<num_bins; i++){
        TKAfile << hist->GetBinContent(i) << endl;
    }
    TKAfile.close();

}
