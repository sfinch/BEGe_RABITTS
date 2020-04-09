////////////////////////////////////////////////////////////////////////////////////////
// beam_corr.C
// Calculates the beam fluctuation correction for an activation run. Does this using the
// (1) BCI, the (2) neutron monitor, and the (3)neutron monitor with a PSD cut. The values
// for these cuts are all defined in RabVar.h
// To run: root -l "beam_corr.C(XXX)" where XXX is run number
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

struct correction{
    double corr = 0;
    double num_events = 0;
    double corr_irr = 0;
    double num_events_irr = 0;
};

void beam_corr(int run_num){

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //Variables
    double half_life = 2.66*3600;
    double lambda = 0.69314718/half_life;

    double elapsed_time;

    double DC_corr = 0;

    correction corr_BCI;
    correction corr_nmon;
    correction corr_nPSD;

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

        if (rabbit.ADC[BCI_chn] >min_BCI){
            corr_BCI.corr += exp(-1*lambda*rabbit.seconds);
            corr_BCI.num_events++;
            if ((rabbit.cycle_time>time_irr[0])
              &&(rabbit.cycle_time<time_irr[1])){
                corr_BCI.corr_irr += exp(-1*lambda*rabbit.seconds);
                corr_BCI.num_events_irr++;
            }
        }
    }
    cout << endl;
    corr_BCI.corr = corr_BCI.corr/corr_BCI.num_events;
    corr_BCI.corr_irr = corr_BCI.corr_irr/corr_BCI.num_events_irr;
    nb = rabbit.GetEntry(nentries-1);
    elapsed_time = rabbit.seconds;
    DC_corr = (1-exp(-1*lambda*elapsed_time));

    //loop over QDC data
    cout << "Looping over QDC data..." << endl;
    nentries = rabbit_QDC.fChain->GetEntriesFast();
    nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit_QDC.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }
        if (rabbit_QDC.ADC_long[nmon_chn]>min_nmon_E){
            corr_nmon.corr += exp(-1*lambda*rabbit_QDC.seconds);
            corr_nmon.num_events++;
            if ((rabbit_QDC.nmon_PSD>nmon_PSD_cut[0])
                && (rabbit_QDC.nmon_PSD<nmon_PSD_cut[1])){
                corr_nPSD.corr += exp(-1*lambda*rabbit_QDC.seconds);
                corr_nPSD.num_events++;
            }
            if ((rabbit_QDC.cycle_time>time_irr[0])
              &&(rabbit_QDC.cycle_time<time_irr[1])){
                corr_nmon.corr_irr += exp(-1*lambda*rabbit_QDC.seconds);
                corr_nmon.num_events_irr++;
                if ((rabbit_QDC.nmon_PSD>nmon_PSD_cut[0])
                    && (rabbit_QDC.nmon_PSD<nmon_PSD_cut[1])){
                    corr_nPSD.corr_irr += exp(-1*lambda*rabbit_QDC.seconds);
                    corr_nPSD.num_events_irr++;
                }
            }
        }
    }
    cout << endl;
    corr_nmon.corr = corr_nmon.corr/corr_nmon.num_events;
    corr_nPSD.corr = corr_nPSD.corr/corr_nPSD.num_events;
    corr_nmon.corr_irr = corr_nmon.corr_irr/corr_nmon.num_events_irr;
    corr_nPSD.corr_irr = corr_nPSD.corr_irr/corr_nPSD.num_events_irr;

    cout << "n-mon:                  "  << corr_nmon.corr << endl;
    cout << "n-mon with PSD cut:     "  << corr_nPSD.corr << endl;
    cout << "BCI:                    "  << corr_BCI.corr << endl;

    cout << "n-mon irradiation only: "  << corr_nmon.corr_irr << endl;
    cout << "n-mon w/ PSD irr only:  "  << corr_nPSD.corr_irr << endl;
    cout << "BCI irradiation only:   "  << corr_BCI.corr_irr << endl;

    cout << "Assuming DC beam:       "  << DC_corr<< endl;

}
