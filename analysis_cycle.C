///////////////////////////////////////////////////////////////////////////
// analysis_cycle.C
// For analyzing RABITTS runs: Projects BEGe energy spectra for the irradiation time,
// counting time, and sub-divides the counting time into different time divisions. 
// The time divisions may be changed in RabVar.h
// ROOT histograms are saved in data_root. After ananlysis_cycle.C has been run, you may
// use plot_cycle.C to plot the spectra and output the spectra to TKA files
// Requries running of mvme2root, followed by process_rabbit // // To run: root -l "analysis_cycle.C(XXX)" where XXX is run number // 
//
///////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <TStyle.h>
#include <TCanvas.h>

using std::cout;
using std::cerr;
using std::endl;

#include "include/processed.hh"
#include "include/RabVar.hh"


void analysis_cycle(int run_num){

    //varibles
    double time_win[RabVar::num_win][2];
    for (int i=0; i<RabVar::num_win; i++){
        time_win[i][0] = (RabVar::time_bin*i)+RabVar::time_count[0];
        time_win[i][1] = (RabVar::time_bin*(i+1))+RabVar::time_count[0];
    }

    //in file
    processed rabbit(run_num);
    //out file
    TFile *fHist = new TFile(Form("data_hist/RABBITS_%i.root", run_num), "RECREATE");
    
    //Hitsos
    TH1F *hCycle = new TH1F("hCycle", "hCycle", 4500, -10, 440);
    TH1F *hEnIrr[RabVar::num_BEGe];
    TH1F *hEnCount[RabVar::num_BEGe];
    TH1F *hEnWin[RabVar::num_win][RabVar::num_BEGe]; 

    for (int i=0; i<RabVar::num_win; i++){
        for (int j=0; j<RabVar::num_BEGe; j++){
            hEnWin[i][j] = new TH1F(Form("hEn_Time%i_Det%i", i, j), Form("hEn_Time%i_Det%i", i, j), 50000, 0, 5000);
        }
    }
    for (int i=0; i<RabVar::num_BEGe; i++){
        hEnCount[i] = new TH1F(Form("hEn_AllCount_Det%i", i), Form("hEn_AllCount_Det%i", i), 40000, 0, 4000);
        hEnIrr[i] = new TH1F(Form("hEn_Irr_Det%i", i), Form("hEn_Irr_Det%i", i), 40000, 0, 4000);
    }

    //loop over data
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }
        if ((rabbit.cycle_time > 0.00001)||(rabbit.cycle_time<-0.00001)){
            for (int det=0; det<RabVar::num_BEGe; det++){
                if (rabbit.En[det]>10){
                    hCycle->Fill(rabbit.cycle_time);
                }
            }
        }

        if ((rabbit.cycle_time>RabVar::time_irr[0]) && (rabbit.cycle_time<RabVar::time_irr[1])){
            for (int det=0; det<RabVar::num_BEGe; det++){
                if (rabbit.En[det]>10){
                    hEnIrr[det]->Fill(rabbit.En[det]);
                }
            }
        }
        else if ((rabbit.cycle_time>RabVar::time_count[0]) && (rabbit.cycle_time<RabVar::time_count[1])){
            for (int det=0; det<RabVar::num_BEGe; det++){
                if (rabbit.En[det]>10){
                    hEnCount[det]->Fill(rabbit.En[det]);
                }
                for (int window=0; window<RabVar::num_win; window++){
                    if ((rabbit.cycle_time>time_win[window][0]) && (rabbit.cycle_time<time_win[window][1])){
                        if (rabbit.En[det]>10){
                            hEnWin[window][det]->Fill(rabbit.En[det]);
                        }
                    }
                }//end time windows
            }
        }//end counting time
    }//end loop over events
    cout << endl;

    //write histos to file
    fHist->cd();

    hCycle->Write();
    for (int det=0; det<RabVar::num_BEGe; det++){
        hEnIrr[det]->Write();
        hEnCount[det]->Write();
    }
    for (int i=0; i<RabVar::num_win; i++){
        for (int j=0; j<RabVar::num_BEGe; j++){
            hEnWin[i][j]->Write();
        }
    }

    fHist->Write();
    fHist->Close();
    
}
