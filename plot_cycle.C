
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
#include "TRandom3.h"
#include "TFile.h"
#include "TCut.h"

#include "src/hist2TKA.C"

void plot_cycle(int run_num, int run_num2 = 0){

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    if (run_num2 < run_num){
        run_num2 = run_num;
    }

    //Histograms
    TH1F *hCycle = new TH1F("hCycle", "hCycle", 4500, -10, 440);
    TH1F *hEnIrr[num_BEGe+1]; //last index is detectors summed together
    TH1F *hEnCount[num_BEGe+1];
    TH1F *hEnWin[num_win][num_BEGe+1]; 

    for (int i=0; i<num_win; i++){
        for (int j=0; j<num_BEGe; j++){

            hEnWin[i][j] = new TH1F(Form("runs%i-%i_Time%i_Det%i", run_num, run_num2, i, j), Form("hEn_Time%i_Det%i", i, j), 50000, 0, 5000);
            if (run_num2 == run_num){
                hEnWin[i][j]->SetName(Form("run%i_Time%i_Det%i", run_num, i, j));
            }
        }
        hEnWin[i][num_BEGe] = new TH1F(Form("runs%i-%i_Time%i_BothDet", run_num, run_num2, i), Form("hEn_Time%i_BothDet", i), 50000, 0, 5000);
        if (run_num2 == run_num){
            hEnWin[i][num_BEGe]->SetName(Form("run%i_Time%i_BothDet", run_num, i));
        }
    }
    for (int i=0; i<num_BEGe; i++){
        hEnCount[i] = new TH1F(Form("runs%i-%i_AllCount_Det%i", run_num, run_num2, i), Form("hEn_AllCount_Det%i", i), 40000, 0, 4000);
        hEnIrr[i] = new TH1F(Form("runs%i-%i_Irr_Det%i", run_num, run_num2, i), Form("hEn_Irr_Det%i", i), 40000, 0, 4000);
        if (run_num2 == run_num){
            hEnCount[i]->SetName(Form("run%i_AllCount_Det%i", run_num, i));
            hEnIrr[i]->SetName(Form("run%i_Irr_Det%i", run_num, i));
        }
    }
    hEnCount[num_BEGe] = new TH1F(Form("runs%i-%i_AllCount_BothDet", run_num, run_num2), Form("hEn_AllCount_BothDet"), 40000, 0, 4000);
    hEnIrr[num_BEGe] = new TH1F(Form("runs%i-%i_Irr_BothDet", run_num, run_num2), Form("hEn_Irr_BothDet"), 40000, 0, 4000);
    if (run_num2 == run_num){
        hEnCount[num_BEGe]->SetName(Form("run%i_AllCount_BothDet", run_num));
        hEnIrr[num_BEGe]->SetName(Form("run%i_Irr_BothDet", run_num));
    }

    //get histos
    for (int k=run_num; k<=run_num2; k++){
        cout << "Run number:  " << k << endl;
        TFile *fHist = new TFile(Form("data_hist/RABBITS_%i.root",k));

        hCycle->Add((TH1F*) fHist->Get("hCycle"));
        for (int i=0; i<num_win; i++){
            for (int j=0; j<num_BEGe; j++){
                hEnWin[i][j]->Add((TH1F*) fHist->Get(Form("hEn_Time%i_Det%i", i, j)));
                hEnWin[i][num_BEGe]->Add((TH1F*) fHist->Get(Form("hEn_Time%i_Det%i", i, j)));
            }
        }
        for (int i=0; i<num_BEGe; i++){
            hEnCount[i]->Add((TH1F*) fHist->Get(Form("hEn_AllCount_Det%i", i)));
            hEnIrr[i]->Add((TH1F*) fHist->Get(Form("hEn_Irr_Det%i", i)));
            hEnCount[num_BEGe]->Add((TH1F*) fHist->Get(Form("hEn_AllCount_Det%i", i)));
            hEnIrr[num_BEGe]->Add((TH1F*) fHist->Get(Form("hEn_Irr_Det%i", i)));
        }

    //get histos

        fHist->Close();
        delete fHist;
    }

    TCanvas *cDet1 = new TCanvas("cDet1","Det 1",1000, 400);
    hEnCount[0]->Rebin(energy_rebin);
    hEnCount[0]->Draw();
    for (int i=0; i<num_win; i++){
        hEnWin[i][0]->Rebin(energy_rebin);
        hEnWin[i][0]->SetLineColor(i+2);
        hEnWin[i][0]->Draw("same");
    }

    TCanvas *cDet2 = new TCanvas("cDet2","Det 2",1000, 400);
    hEnCount[1]->Rebin(energy_rebin);
    hEnCount[1]->Draw();
    for (int i=0; i<num_win; i++){
        hEnWin[i][1]->Rebin(energy_rebin);
        hEnWin[i][1]->SetLineColor(i+2);
        hEnWin[i][1]->Draw("same");
    }

    TCanvas *cBothDet2 = new TCanvas("cBothDet2","Both Det",1000, 400);
    hEnCount[2]->Rebin(energy_rebin);
    hEnCount[2]->Draw();
    for (int i=0; i<num_win; i++){
        hEnWin[i][num_BEGe]->Rebin(energy_rebin);
        hEnWin[i][num_BEGe]->SetLineColor(i+2);
        hEnWin[i][num_BEGe]->Draw("same");
    }

    hist2TKA(hEnCount[num_BEGe]);
    for (int j=0; j<num_BEGe+1; j++){
        for (int i=0; i<num_win; i++){
            hist2TKA(hEnWin[i][j]);
            hist2TKA(hEnWin[i][num_BEGe]);
        }
    }

}
