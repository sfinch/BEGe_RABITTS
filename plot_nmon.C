
#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "include/processed.h"

using std::cout;
using std::cerr;
using std::endl;

void plot_nmon(int run_num){

    //Variables
    TH1F *hCycNmon  = new TH1F("hCycNmon", "hCycNmon", 9000, -1, 89);
    TH1F *hNmonT;
    TH1F *hCycFC[2];
    TH1F *hHPGe[2];
    for (int i=0; i<2; i++){
        hCycFC[i] = new TH1F(Form("hCycFC%i", i), Form("hCycFC%i", i), 9000, -1, 89);
        hHPGe[i] = new TH1F(Form("hHPGe%i", i), Form("hHPGe%i", i), 9000, -1, 89);
    }

    //cuts
    int min_nmon_E = 100;
    int min_FC_E[2] = {6000, 6000};
    int min_HPGe_E[2] = {100, 100};

    //in file
    processed rabbit(run_num);

    //loop over data
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;


    rabbit.GetEntry(nentries-1);
    //Double_t total_time = rabbit.seconds;
    //hNmonT  = new TH1F("hNmonT", "hNmonT", 1*total_time, 0, total_time);
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }
        if (rabbit.ADC[14]>min_nmon_E){
            hCycNmon->Fill(rabbit.cycle_time);
            //hNmonT->Fill(rabbit.seconds);
        }
        if (rabbit.ADC[10]>min_FC_E[0]){
            hCycFC[0]->Fill(rabbit.cycle_time);
        }
        if (rabbit.ADC[11]>min_FC_E[1]){
            hCycFC[1]->Fill(rabbit.cycle_time);
        }
        if (rabbit.En[0]>min_HPGe_E[0]){
            hHPGe[0]->Fill(rabbit.cycle_time);
        }
        if (rabbit.En[1]>min_HPGe_E[1]){
            hHPGe[1]->Fill(rabbit.cycle_time);
        }

    }//end loop over events

    //plot
    TCanvas *cCycle = new TCanvas("cCycle", "cCycle", 500, 1000);
    cCycle->Divide(1,5);
    cCycle->cd(1);
    hCycNmon->Draw();
    gPad->SetLogy();
    cCycle->cd(2);
    hCycFC[0]->Draw();
    cCycle->cd(3);
    hCycFC[1]->Draw();
    cCycle->cd(4);
    hHPGe[0]->Draw();
    cCycle->cd(5);
    hHPGe[1]->Draw();

    //save to file
    cCycle->SaveAs(Form("time_spec/run%i.png", run_num));
    cCycle->SaveAs(Form("time_spec/run%i.C", run_num));
    
}
