
#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "include/processed.h"

using std::cout;
using std::cerr;
using std::endl;

void analysis_cycle(int run_num){

    //Variables
    const int num_det = 2;

    const int num_win = 20;
    double irr_time = 4;
    double count_time = 20;
    double transit_time = 1.05;
    
    //in file
    processed rabbit(run_num);
    //out file
    TFile *fHist = new TFile(Form("data_hist/RABBITS_%i.root", run_num), "RECREATE");
    
    //cuts
    double time_bin = count_time/num_win;
    double time_irr[2] = {transit_time, irr_time+transit_time};
    double time_count[2] = {irr_time + 2*transit_time, 
                            irr_time + 2*transit_time + count_time};
    double time_win[num_win][2];
    for (int i=0; i<num_win; i++){
        time_win[i][0] = (time_bin*i)+time_count[0];
        time_win[i][1] = (time_bin*(i+1))+time_count[0];
    }

    //Hitsos
    TH1F *hCycle = new TH1F("hCycle", "hCycle", 4500, -10, 440);
    TH1F *hEnIrr[num_det];
    TH1F *hEnCount[num_det];
    TH1F *hEnWin[num_win][num_det]; 

    for (int i=0; i<num_win; i++){
        for (int j=0; j<num_det; j++){
            hEnWin[i][j] = new TH1F(Form("hEn_Time%i_Det%i", i, j), Form("hEn_Time%i_Det%i", i, j), 50000, 0, 5000);
        }
    }
    for (int i=0; i<num_det; i++){
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
            for (int det=0; det<num_det; det++){
                if (rabbit.En[det]>10){
                    hCycle->Fill(rabbit.cycle_time);
                }
            }
        }

        if ((rabbit.cycle_time>time_irr[0]) && (rabbit.cycle_time<time_irr[1])){
            for (int det=0; det<num_det; det++){
                if (rabbit.En[det]>10){
                    hEnIrr[det]->Fill(rabbit.En[det]);
                }
            }
        }
        else if ((rabbit.cycle_time>time_count[0]) && (rabbit.cycle_time<time_count[1])){
            for (int det=0; det<num_det; det++){
                if (rabbit.En[det]>10){
                    hEnCount[det]->Fill(rabbit.En[det]);
                }
                for (int window=0; window<num_win; window++){
                    if ((rabbit.cycle_time>time_win[window][0]) && (rabbit.cycle_time<time_win[window][1])){
                        if (rabbit.En[det]>10){
                            hEnWin[window][det]->Fill(rabbit.En[det]);
                        }
                    }
                }//end time windows
            }
        }//end counting time

    }//end loop over events

    //write histos to file
    fHist->cd();

    hCycle->Write();
    for (int det=0; det<num_det; det++){
        hEnIrr[det]->Write();
        hEnCount[det]->Write();
    }
    for (int i=0; i<num_win; i++){
        for (int j=0; j<num_det; j++){
            hEnWin[i][j]->Write();
        }
    }

    fHist->Write();
    fHist->Close();
    
}
