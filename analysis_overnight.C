///////////////////////////////////////////////////////////////////////////
// analysis_overnight.C
// For analyzing overnight runs: Produces BEGe energy spectra in 1 hour division for the
// entire run length. The time division (default 1 hour) may be changed in RabVar.h
// ROOT histograms are saved in data_root, and TKA spectra are saved in data_TKA
// Requries running of mvme2root, followed by process_rabbit
//
// To run: root -l "analysis_overnight.C(XXX)" where XXX is run number
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
#include "src/hist2TKA.C"


void analysis_overnight(int run_num){

    //Variables
    TH1F *hEnAll[RabVar::num_BEGe+1];
    TH1F *hEnWin[RabVar::num_overnight_win][RabVar::num_BEGe+1]; //[window 1-3][det 1, 2, summed]

    for (int i=0; i<RabVar::num_overnight_win; i++){
        for (int j=0; j<RabVar::num_BEGe; j++){
            hEnWin[i][j] = new TH1F(Form("run%i_Hour%i_Det%i", run_num, (i+1), j), Form("hEnWin%i%i", i, j), 30000, 0, 3000);
        }
        hEnWin[i][RabVar::num_BEGe] = new TH1F(Form("run%i_Hour%i_BothDet", run_num, i+1), Form("hEnWin%iBoth", i), 30000, 0, 3000);
    }

    for (int i=0; i<RabVar::num_BEGe; i++){
        hEnAll[i] = new TH1F(Form("run%i_All_Det%i", run_num, i), Form("hEnAll%i", i), 30000, 0, 3000);
    }
    hEnAll[RabVar::num_BEGe] = new TH1F(Form("run%i_All_BothDet", run_num), Form("hEnAll%i", RabVar::num_BEGe), 30000, 0, 3000);

    //cuts
    double time_win[RabVar::num_overnight_win][2];
    for (int i=0; i<RabVar::num_overnight_win; i++){
        time_win[i][0] = i*RabVar::overnight_win_time;
        time_win[i][1] = (i+1)*RabVar::overnight_win_time;
    }

    //in file
    processed rabbit(run_num);

    //out file
    TFile *fHist = new TFile(Form("data_hist/RABBIT_%i.root", run_num), "RECREATE");

    //loop over data
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }

        for (int det=0; det<RabVar::num_BEGe; det++){
            if (rabbit.En[det]>10){
                hEnAll[det]->Fill(rabbit.En[det]);
                hEnAll[RabVar::num_BEGe]->Fill(rabbit.En[det]);
            }
        }

        for (int window=0; window<RabVar::num_overnight_win; window++){
            if ((rabbit.seconds>time_win[window][0]) && (rabbit.seconds<time_win[window][1])){
                //clover 1
                if (rabbit.En[0]>10){
                    hEnWin[window][0]->Fill(rabbit.En[0]);
                    hEnWin[window][RabVar::num_BEGe]->Fill(rabbit.En[0]);
                }
                //clover 2
                if (rabbit.En[1]>10){
                    hEnWin[window][1]->Fill(rabbit.En[1]);
                    hEnWin[window][RabVar::num_BEGe]->Fill(rabbit.En[1]);
                }
            }
        }//end time windows

    }//end loop over events
    cout << endl;
    int cycles_complete = rabbit.seconds/RabVar::overnight_win_time;
    cout << rabbit.seconds << " s elapsed" <<  endl;
    cout << rabbit.seconds/3600. << " h elapsed" <<  endl;
    cout << rabbit.seconds/RabVar::overnight_win_time << " total time cycles" <<  endl;
    cout << cycles_complete << " complete cycles" <<  endl;
    if (cycles_complete>RabVar::num_overnight_win){
        cycles_complete=RabVar::num_overnight_win;
    }

    TCanvas *cDet1 = new TCanvas("cDet1","Det 1",1000, 400);
    hEnWin[0][0]->Draw();
    for (int i=1; i<cycles_complete; i++){
        hEnWin[i][0]->SetLineColor(i+2);
        hEnWin[i][0]->Draw("same");
     }

    TCanvas *cDet2 = new TCanvas("cDet2","Det 2",1000, 400);
    hEnWin[0][1]->Draw();
    for (int i=1; i<cycles_complete; i++){
        hEnWin[i][1]->SetLineColor(i+2);
        hEnWin[i][1]->Draw("same");
    }

    TCanvas *cBothDet= new TCanvas("cDet","Both dets",1000, 400);
    hEnWin[0][2]->Draw();
    for (int i=1; i<cycles_complete; i++){
        hEnWin[i][2]->SetLineColor(i+1);
        hEnWin[i][2]->Draw("same");
    }

    //write histos to file
    fHist->cd();

    for (int det=0; det<RabVar::num_BEGe+1; det++){
        hEnAll[det]->Write();
        hist2TKA(hEnAll[det]);
    }
    for (int i=0; i<cycles_complete; i++){
        for (int j=0; j<RabVar::num_BEGe+1; j++){
            hEnWin[i][j]->Write();
            hist2TKA(hEnWin[i][j]);
        }
    }
    for (int j=0; j<RabVar::num_BEGe+1; j++){
        hEnWin[cycles_complete][j]->SetName(Form("run%i_Part%i_Det%i", run_num, (cycles_complete+1), j+1));
        hEnWin[cycles_complete][j]->Write();
        hist2TKA(hEnWin[cycles_complete][j]);
    }

    fHist->Write();
    fHist->Close();
    
}
