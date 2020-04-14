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
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

#include "include/processed.hh"
#include "include/processed_QDC.hh"
#include "include/RabVar.hh"

struct DT{

  public:
    double rate = 0;
    double num_PU = 0;
    double num_events = 0;
    double num_write = 0;

    double percent_PU = 0;
    double percent_write = 0;
    double percent_tot = 0;

    int last_event = 0;

    double exp_int = 0;
    int rebin = 100;

    TString name;
    TH1F *hTS;
    TF1 *fExpo;
    
    DT(TString);
    void calc();
};

DT::DT(TString n){
    name = n;
    hTS = new TH1F(Form("hTS%s", name.Data()), Form("hTS%s", name.Data()), 
         3000000/rebin, 0, 3000000);
    fExpo = new TF1(Form("fExpo%s", name.Data()), "exp([0]+[1]*x)",
         0, 3000000);
    fExpo->SetParameter(0, 5);
    fExpo->SetParameter(1, -1e-5);
};

void DT::calc(){
    if (num_events>0){
        hTS->GetXaxis()->SetRangeUser(1000, 3000000);
        hTS->Fit(Form("fExpo%s", name.Data()));
        exp_int = fExpo->Eval(0)*(-1.)/(rebin*1.*fExpo->GetParameter(1));

        num_write = exp_int - num_events;

        percent_PU = 100*num_PU/num_events;
        percent_write = 100*num_write/num_events;
        if (percent_write>100){
            percent_write = 0;
            num_write = 0;
        }
        percent_tot = 100*(num_PU+num_write)/num_events;

    }
};

void deadtime(int run_num){

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //Variables
    double elapsed_time;
    double time_win[RabVar::num_win][2];
    for (int i=0; i<RabVar::num_win; i++){
        time_win[i][0] = (RabVar::time_bin*i)+RabVar::time_count[0];
        time_win[i][1] = (RabVar::time_bin*(i+1))+RabVar::time_count[0];
    }

    DT *SCP[16];
    DT *HPGe_count[RabVar::num_det];
    DT *HPGe_win[RabVar::num_win][RabVar::num_det];

    for (int chn=0; chn<16; chn++){
        SCP[chn] = new DT(Form("scp%i", chn));
    }
    for (int chn=0; chn<RabVar::num_det; chn++){
        HPGe_count[chn] = new DT(Form("HPGe%i", chn));
        for (int win=0; win<RabVar::num_win; win++){
            HPGe_win[win][chn] = new DT(Form("HPGe%i%i", win, chn));
        }
    }

    //in file
    processed rabbit(run_num);
    processed_QDC rabbit_QDC(run_num);

    //loop over SCP data
    cout << "Looping over SCP data..." << endl;
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }

        for (int chn=0; chn<16; chn++){
            if (rabbit.TDC[chn] >1){
                SCP[chn]->num_events++;
                SCP[chn]->hTS->Fill(double(rabbit.time_stamp-SCP[chn]->last_event));

                if (rabbit.pileup[chn]){
                    SCP[chn]->num_PU++;
                }

                SCP[chn]->last_event = rabbit.time_stamp;
            }
        }

        for (int chn=0; chn<RabVar::num_det; chn++){
            if (rabbit.TDC[RabVar::det_chn[chn]] >1){
                if ((rabbit.cycle_time>RabVar::time_count[0])
                  &&(rabbit.cycle_time<RabVar::time_count[1])){

                    HPGe_count[chn]->num_events++;
                    HPGe_count[chn]->hTS->Fill(double(rabbit.time_stamp-HPGe_count[chn]->last_event));
                    if (rabbit.pileup[RabVar::det_chn[chn]]){
                        HPGe_count[chn]->num_PU++;
                    }
                    HPGe_count[chn]->last_event = rabbit.time_stamp;

                    for (int window=0; window<RabVar::num_win; window++){
                        if ((rabbit.cycle_time>time_win[window][0]) 
                        && (rabbit.cycle_time<time_win[window][1])){
                            HPGe_win[window][chn]->num_events++;
                            HPGe_win[window][chn]->hTS->Fill(
                              double(rabbit.time_stamp-HPGe_win[window][chn]->last_event));
                            if (rabbit.pileup[RabVar::det_chn[chn]]){
                                HPGe_win[window][chn]->num_PU++;
                            }
                            HPGe_win[window][chn]->last_event = rabbit.time_stamp;
                        }
                    }
                }
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

            if ((rabbit_QDC.nmon_PSD>nmon_PSD_cut[0]) && (rabbit_QDC.nmon_PSD<nmon_PSD_cut[1])){

            }
            if ((rabbit_QDC.cycle_time>RabVar::time_irr[0])
              &&(rabbit_QDC.cycle_time<RabVar::time_irr[1])){

            }
        }
    }
    cout << endl;
    */

    TCanvas *cDT = new TCanvas("cDT", "cDT", 1000, 1000);
    cDT->Divide(4,4);
    for (int chn=0; chn<16; chn++){
        cDT->cd(chn+1);
        SCP[chn]->hTS->Draw();
        SCP[chn]->calc();
    }
    for (int chn=0; chn<RabVar::num_det; chn++){
        cDT->cd(chn+1);
        //HPGe_count[chn]->hTS->Draw();
        HPGe_count[chn]->calc();
        for (int window=0; window<RabVar::num_win; window++){
            HPGe_win[window][chn]->calc();
        }
    }


    for (int chn=0; chn<RabVar::num_det; chn++){
        cout << "-------------------------------------------------------------" << endl;
        cout << "Det " << chn << " time windows" << endl;
        cout << "-------------------------------------------------------------" << endl;
        cout << "Win \tPileup \tWrite \tCombined " << endl;
        for (int window=0; window<RabVar::num_win; window++){
            cout << window << "\t";
            cout.precision(3);
            cout << HPGe_win[window][chn]->percent_PU << "%\t" 
                 << HPGe_win[window][chn]->percent_write << "%\t"
                 << HPGe_win[window][chn]->percent_tot << "%" << endl;
        }
    }

    cout << "-------------------------------------------------------------" << endl;
    cout << "Irradiation only" << endl;
    cout << "-------------------------------------------------------------" << endl;
    cout << "Chn \tPileup \tWrite \tCombined " << endl;
    for (int chn=0; chn<RabVar::num_det; chn++){
        cout << chn << "\t";
        cout.precision(3);
        cout << HPGe_count[chn]->percent_PU << "%\t" 
             << HPGe_count[chn]->percent_write << "%\t"
             << HPGe_count[chn]->percent_tot << "%" << endl;
    }

    cout << "-------------------------------------------------------------" << endl;
    cout << "Total run" << endl;
    cout << "-------------------------------------------------------------" << endl;
    cout << "Chn \tRate\tPileup \tWrite \tCombined " << endl;
    for (int chn=0; chn<16; chn++){
        SCP[chn]->rate = SCP[chn]->num_events/elapsed_time;
        //cout << chn << "\t" << SCP[chn]->num_events<< "\t"
        //     << SCP[chn]->num_PU << "%\t" 
        //     << SCP[chn]->num_write << "%" << endl;
        cout.precision(4);
        cout << chn << "\t" << SCP[chn]->rate << "\t";
        cout.precision(3);
        cout << SCP[chn]->percent_PU << "%\t" 
             << SCP[chn]->percent_write << "%\t"
             << SCP[chn]->percent_tot << "%" << endl;
    }


    cout << "-------------------------------------------------------------" << endl;
    cout << "Run start:    " << rabbit.rawfile->Get("start_time")->GetTitle() << endl;
    cout << "Run stop:     " << rabbit.rawfile->Get("stop_time")->GetTitle() << endl;
    cout << "Elapsed time: " << elapsed_time << " s" << endl;
    cout << "-------------------------------------------------------------" << endl;
    
}
