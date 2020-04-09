////////////////////////////////////////////////////////////////////////////////////////
// deadtime.C
// Calculates the dead time for all MDPP16 channels for the entire run. Additionally does
// this for the cycle time sub-divisions defined in RabVar.h
// To run: root -l "deadtime.C(XXX)" where XXX is run number
// Requires that both mvme2root and process_rabbit have been run. 
///////////////////////////////////////////////////////////////////////////////////////

using namespace std;
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "TRint.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TH1.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TFile.h"

#include "include/processed.h"
#include "include/processed_QDC.h"
#include "RabVar.h"

using namespace RabVar;

struct DT{
    double num_PU= 0;
    double num_events = 0;
    int last_event=0;
};

void deadtime(int run_num){

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //Variables
    double elapsed_time;
    DT SCP[16];

    //in file
    processed rabbit(run_num);
    processed_QDC rabbit_QDC(run_num);

    cout << rabbit.rawfile->Get("start_time")->GetTitle() << endl;
    cout << rabbit.rawfile->Get("stop_time")->GetTitle() << endl;

    //loop over SCP data
    cout << "Looping over SCP data..." << endl;
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }

        if (rabbit.ADC[BCI_chn] >10){

            if ((rabbit.cycle_time>time_irr[0])
              &&(rabbit.cycle_time<time_irr[1])){

            }
        }
    }
    cout << endl;
    nb = rabbit.GetEntry(nentries-1);
    elapsed_time = rabbit.seconds;

    //loop over QDC data
    /*
    cout << "Looping over QDC data..." << endl;
    nentries = rabbit_QDC.fChain->GetEntriesFast();
    nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit_QDC.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }
        if (rabbit_QDC.ADC_long[nmon_chn]>min_nmon_E){

            if ((rabbit_QDC.nmon_PSD>nmon_PSD_cut[0])
                && (rabbit_QDC.nmon_PSD<nmon_PSD_cut[1])){

            }
            if ((rabbit_QDC.cycle_time>time_irr[0])
              &&(rabbit_QDC.cycle_time<time_irr[1])){

            }
        }
    }
    cout << endl;
    */


}
