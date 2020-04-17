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
#include "TCanvas.h"
#include "TStyle.h"
#include "TDatime.h"

#include "include/processed.hh"
#include "include/processed_QDC.hh"
#include "include/RabVar.hh"


class correction{
  public:
    double dt = 1; 
    double lambda = 1; 

    double corr = 0;
    double num_events = 0;
    double last_time = 0;
    
    double corr_irr = 0;
    double num_events_irr = 0;
    double last_time_irr = 0;

    vector<int> scalers;
    vector<int> scalers_irr;

    correction(double l, double t);
    void calc();
};

correction::correction(double l, double t){
    dt = t;
    lambda = l;
};

void correction::calc(){

    for (int i=0; i<scalers.size(); i++){
        corr += scalers.at(i)*(1-exp(-1*lambda*dt))
            *exp(-1*lambda*(scalers.size()-i)*dt);
        num_events += scalers.at(i);
    }

    for (int i=0; i<scalers_irr.size(); i++){
        corr_irr += scalers_irr.at(i)*(1-exp(-1*lambda*dt))
            *exp(-1*lambda*(scalers_irr.size()-i)*dt);
        num_events_irr += scalers_irr.at(i);
    }

    corr = corr*scalers.size()/num_events;
    corr_irr = corr_irr*scalers.size()/num_events_irr;
};

void beam_corr(int run_num){

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //Variables
    double half_life = 2.66*3600;
    double lambda = 0.69314718/half_life;
    double dt = 1.00;

    double elapsed_time;
    double SCP_time;
    double QDC_time;

    double DC_corr = 0;

    TDatime tStart;
    TDatime tStop;

    correction corr_BCI(lambda, dt);
    correction corr_nmon(lambda, dt);
    correction corr_nPSD(lambda, dt);

    //in file
    processed rabbit(run_num);
    processed_QDC rabbit_QDC(run_num);
    tStart = TDatime(rabbit.rawfile->Get("start_time")->GetTitle());
    tStop = TDatime(rabbit.rawfile->Get("stop_time")->GetTitle());
    elapsed_time = tStop.Convert() - tStart.Convert();

    corr_BCI.scalers.resize(int(elapsed_time/dt), 0);
    corr_BCI.scalers_irr.resize(int(elapsed_time/dt), 0);
    corr_nmon.scalers.resize(int(elapsed_time/dt), 0);
    corr_nPSD.scalers.resize(int(elapsed_time/dt), 0);
    corr_nmon.scalers_irr.resize(int(elapsed_time/dt), 0);
    corr_nPSD.scalers_irr.resize(int(elapsed_time/dt), 0);

    //loop over SCP data
    cout << "Looping over SCP data..." << endl;
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    nb = rabbit.GetEntry(nentries-1);
    SCP_time = rabbit.seconds;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }
        if (rabbit.seconds > SCP_time){
            SCP_time = rabbit.seconds;
        }

        if (rabbit.ADC[RabVar::BCI_chn] > RabVar::min_BCI){
    
            corr_BCI.scalers.at(int(rabbit.seconds/dt))++;

            if ((rabbit.cycle_time > RabVar::time_irr[0])
              &&(rabbit.cycle_time < RabVar::time_irr[1])){

                corr_BCI.scalers_irr.at(int(rabbit.seconds/dt))++;
            }
        }
    }
    cout << '\r' << nentries << " total SCP events" << endl;
    cout << SCP_time << " s SCP time" << endl;

    //loop over QDC data
    cout << "Looping over QDC data..." << endl;
    nentries = rabbit_QDC.fChain->GetEntriesFast();
    nbytes = 0, nb = 0;
    nb = rabbit_QDC.GetEntry(nentries-1);
    QDC_time = rabbit_QDC.seconds;

    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit_QDC.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }
        if (rabbit_QDC.seconds > QDC_time){
            QDC_time = rabbit_QDC.seconds;
        }
        if (rabbit_QDC.ADC_long[RabVar::nmon_chn] > RabVar::min_nmon_E){

            corr_nmon.scalers.at(int(rabbit_QDC.seconds/dt))++;

            if ((rabbit_QDC.nmon_PSD > RabVar::nmon_PSD_cut[0])
              && (rabbit_QDC.nmon_PSD < RabVar::nmon_PSD_cut[1])){

                corr_nPSD.scalers.at(int(rabbit_QDC.seconds/dt))++;

            }
            if ((rabbit_QDC.cycle_time > RabVar::time_irr[0])
              &&(rabbit_QDC.cycle_time < RabVar::time_irr[1])){

                corr_nmon.scalers_irr.at(int(rabbit_QDC.seconds/dt))++;

                if ((rabbit_QDC.nmon_PSD > RabVar::nmon_PSD_cut[0])
                  && (rabbit_QDC.nmon_PSD < RabVar::nmon_PSD_cut[1])){

                    corr_nPSD.scalers_irr.at(int(rabbit_QDC.seconds/dt))++;

                }
            }
        }
    }
    cout << '\r' << nentries << " total QDC events" << endl;
    cout << QDC_time << " s QDC time" << endl;
    

    //output scalers to file
    FILE *file_ptr = fopen(Form("datafiles/run%i.scalers", run_num), "w");
    for (int i=0; i<corr_nmon.scalers.size(); i++){
        fprintf(file_ptr, "%i\t", i);
        fprintf(file_ptr, "%f\t", i*dt);
        fprintf(file_ptr, "%i\t", corr_BCI.scalers.at(i));
        fprintf(file_ptr, "%i\n", corr_nmon.scalers.at(i));
    }
    fclose(file_ptr);

    //calculate corrections
    DC_corr = (1-exp(-1*lambda*elapsed_time));
    corr_BCI.calc();
    corr_nmon.calc();
    corr_nPSD.calc();

    //print corrections
    cout << "------------------------------------------------------" << endl;
    cout << "Run start:    " << tStart.AsString() << endl;
    cout << "Run stop:     " << tStop.AsString() << endl;
    cout << "Elapsed time: " << elapsed_time/3600. << " h" << endl;
    cout << "------------------------------------------------------" << endl;
    cout.precision(5);
    cout << "\t\t\t" << "Produced \tCorrection %" << endl;
    cout << "Assuming DC beam:       "  << DC_corr
         << "\t\t" << DC_corr/DC_corr << endl;

    cout << "n-mon:                  "  << corr_nmon.corr
         << "\t\t" << corr_nmon.corr/DC_corr << endl;
    cout << "n-mon with PSD cut:     "  << corr_nPSD.corr
         << "\t\t" << corr_nPSD.corr/DC_corr << endl;
    cout << "BCI:                    "  << corr_BCI.corr
         << "\t\t" << corr_BCI.corr/DC_corr << endl;

    cout << "n-mon irradiation only: "  << corr_nmon.corr_irr
         << "\t\t" << corr_nmon.corr_irr/DC_corr << endl;
    cout << "n-mon w/ PSD irr only:  "  << corr_nPSD.corr_irr
         << "\t\t" << corr_nPSD.corr_irr/DC_corr << endl;
    cout << "BCI irradiation only:   "  << corr_BCI.corr_irr
         << "\t\t" << corr_BCI.corr_irr/DC_corr << endl;


}
