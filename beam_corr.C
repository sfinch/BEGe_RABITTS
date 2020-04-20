////////////////////////////////////////////////////////////////////////////////////////
// beam_corr.C
// Calculates the beam fluctuation correction for an activation run. Does this using the
// (1) BCI, the (2) neutron monitor, and the (3)neutron monitor with a PSD cut. The values
// for these cuts are all defined in RabVar.h
// To run: root -l "beam_corr.C(XXX, YYY)" where XXX and YY are run numbers
// Will calculate the beam correction for runs XXX through YYY inclusive 
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

    correction();
    correction(double l, double t);
    void calc();
    void resize(int length);
};

correction::correction(){
}

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

void correction::resize(int length){
    scalers.resize(length, 0);
    scalers_irr.resize(length, 0);
}

void beam_corr(int run_num, int run_num2 = 0){

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //Variables
    double half_life = 2.66*3600;
    double lambda = 0.69314718/half_life;
    double dt = 0.10;

    if (run_num2 < run_num){
        run_num2 = run_num;
    }
    const int num_runs = run_num2 - run_num + 1;

    double elapsed_time[num_runs];
    double SCP_time[num_runs];
    double QDC_time[num_runs];
    double total_time;

    double DC_corr = 0;
    int numTbins = 0;
    int startOffset = 0;

    TDatime tStart[num_runs];
    TDatime tStop[num_runs];

    Long64_t nentries;
    Long64_t nbytes = 0, nb = 0;

    correction corr_BCI(lambda, dt);
    correction corr_nmon(lambda, dt);
    correction corr_nPSD(lambda, dt);

    //in file
    processed *rabbit[num_runs];
    processed_QDC *rabbit_QDC[num_runs];
    for (int i=0; i<num_runs; i++){
        rabbit[i] = new processed(i+run_num);
        rabbit_QDC[i] = new processed_QDC(i+run_num);

        tStart[i] = TDatime(rabbit[i]->rawfile->Get("start_time")->GetTitle());
        tStop[i] = TDatime(rabbit[i]->rawfile->Get("stop_time")->GetTitle());
        elapsed_time[i] = tStop[i].Convert() - tStart[i].Convert();
    }
    total_time = tStop[num_runs-1].Convert() - tStart[0].Convert();

    numTbins =  total_time/dt;
    corr_BCI.resize(numTbins);
    corr_nmon.resize(numTbins);
    corr_nPSD.resize(numTbins);

    //loop over runs
    for (int run=0; run<num_runs; run++){
        cout << "------------------------------------------------------" << endl;
        cout << "Run " << run+run_num << " start:        " << tStart[run].AsString() << endl;
        cout << "Run " << run+run_num << " stop:         " << tStop[run].AsString() << endl;
        cout << "Run " << run+run_num << " elapsed time: " 
             << elapsed_time[run] << " s, " << elapsed_time[run]/3600. << " h" << endl;
        cout << "------------------------------------------------------" << endl;
        
        startOffset = tStart[run].Convert() - tStart[0].Convert();

        //loop over SCP data
        cout << "Looping over SCP data..." << endl;
        nentries = rabbit[run]->fChain->GetEntriesFast();
        nbytes = 0, nb = 0;
        nb = rabbit[run]->GetEntry(nentries-1);
        SCP_time[run] = rabbit[run]->seconds;
        
        for (Long64_t jentry=0; jentry<nentries; jentry++) {
            nb = rabbit[run]->GetEntry(jentry);   nbytes += nb;
            if (jentry%100000==0){
                cout << '\r' << "Processing event " << jentry;
            }
            if (rabbit[run]->seconds > SCP_time[run]){
                SCP_time[run] = rabbit[run]->seconds;
            }

            if (rabbit[run]->ADC[RabVar::BCI_chn] > RabVar::min_BCI){
        
                corr_BCI.scalers.at(int((startOffset + rabbit[run]->seconds)/dt))++;

                if ((rabbit[run]->cycle_time > RabVar::time_irr[0])
                  &&(rabbit[run]->cycle_time < RabVar::time_irr[1])){

                    corr_BCI.scalers_irr.at(int((startOffset + rabbit[run]->seconds)/dt))++;
                }
            }
        } //end loop over SCP
        cout << '\r' << nentries << " total SCP events" << endl;
        cout << SCP_time[run] << " s SCP time" << endl;

        //loop over QDC data
        cout << "Looping over QDC data..." << endl;
        nentries = rabbit_QDC[run]->fChain->GetEntriesFast();
        nbytes = 0, nb = 0;
        nb = rabbit_QDC[run]->GetEntry(nentries-1);
        QDC_time[run] = rabbit_QDC[run]->seconds;

        for (Long64_t jentry=0; jentry<nentries; jentry++) {
            nb = rabbit_QDC[run]->GetEntry(jentry);   nbytes += nb;
            if (jentry%100000==0){
                cout << '\r' << "Processing event " << jentry;
            }
            if (rabbit_QDC[run]->seconds > QDC_time[run]){
                QDC_time[run] = rabbit_QDC[run]->seconds;
            }
            if (rabbit_QDC[run]->ADC_long[RabVar::nmon_chn] > RabVar::min_nmon_E){

                corr_nmon.scalers.at(int((startOffset + rabbit_QDC[run]->seconds)/dt))++;

                if ((rabbit_QDC[run]->nmon_PSD > RabVar::nmon_PSD_cut[0])
                  && (rabbit_QDC[run]->nmon_PSD < RabVar::nmon_PSD_cut[1])){

                    corr_nPSD.scalers.at(int((startOffset + rabbit_QDC[run]->seconds)/dt))++;

                }
                if ((rabbit_QDC[run]->cycle_time > RabVar::time_irr[0])
                  &&(rabbit_QDC[run]->cycle_time < RabVar::time_irr[1])){

                    corr_nmon.scalers_irr.at(int((startOffset + rabbit_QDC[run]->seconds)/dt))++;

                    if ((rabbit_QDC[run]->nmon_PSD > RabVar::nmon_PSD_cut[0])
                      && (rabbit_QDC[run]->nmon_PSD < RabVar::nmon_PSD_cut[1])){

                        corr_nPSD.scalers_irr.at(int((startOffset + rabbit_QDC[run]->seconds)/dt))++;

                    }
                }
            }
        } //end loop over QDC
        cout << '\r' << nentries << " total QDC events" << endl;
        cout << QDC_time[run] << " s QDC time" << endl;
    } //end loop over runs

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
    DC_corr = (1-exp(-1*lambda*total_time));
    corr_BCI.calc();
    corr_nmon.calc();
    corr_nPSD.calc();

    //print corrections
    cout << "------------------------------------------------------" << endl;
    cout << "Total start:    " << tStart[0].AsString() << endl;
    cout << "Total stop:     " << tStop[num_runs-1].AsString() << endl;
    cout << "Elapsed time:   " << total_time/3600. << " h" << endl;
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
