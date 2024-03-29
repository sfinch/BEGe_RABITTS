///////////////////////////////////////////////////////////////////////////
// plot_FC.C
// Plots the histograms for the two fission chambers, and integrats the number of fission
// events above threshold. The values for these cuts are defined in RabVar.h
// Only requires raw data has been converted to ROOT using mvme2root
//
// To run: root -l "beam_corr.C(XXX)" where XXX is run number
//
//////////////////////////////////////////////////////////////////////////


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

#include "include/RabVar.hh"


void plot_FC(int run_num){

    gStyle->SetOptStat(0);

    //Variables
    TH1F *hFC[RabVar::num_FC];
    TLine *l[RabVar::num_FC];

    //get histos
    TFile *fHist;
    if (run_num<10){
        fHist = new TFile(Form("data_root/RABITTS_00%i.root", run_num));
    }
    else if (run_num<100){
        fHist = new TFile(Form("data_root/RABITTS_0%i.root", run_num));
    }
    else if (run_num<1000){
        fHist = new TFile(Form("data_root/RABITTS_%i.root", run_num));
    }
    else if (run_num<1010){
        run_num = run_num-1000;
        fHist = new TFile(Form("data_root/RABITTS2_00%i.root", run_num));
    }
    else if (run_num<1100){
        run_num = run_num-1000;
        fHist = new TFile(Form("data_root/RABITTS2_0%i.root", run_num));
    }
    else{
        run_num = run_num-1000;
        fHist = new TFile(Form("data_root/RABITTS2_%i.root", run_num));
    }

    TCanvas *cFC = new TCanvas("cFC","Fisson chamber energy spectra", 800, 800);
    cFC->Divide(1, RabVar::num_FC);

    for (int j=0; j<RabVar::num_FC; j++){
        cFC->cd(j+1);

        hFC[j] = (TH1F*) (fHist->Get(Form("histos_SCP/hADC%i", RabVar::FC_chn[j])))->Clone();
        hFC[j]->SetTitle(Form("Run %i, FC%i", run_num, j+1));
        cout << "FC" << j+1 << " counts: " << hFC[j]->Integral(RabVar::FC_threshold[j], 65535) << endl;

        hFC[j]->Rebin(RabVar::FC_rebin);
        hFC[j]->GetXaxis()->SetRangeUser(2000, 65535);
        hFC[j]->GetYaxis()->SetTitle("Counts");
        hFC[j]->GetXaxis()->SetTitle("ADC channel");
        hFC[j]->Draw();

        l[j] = new TLine(RabVar::FC_threshold[j], 0, RabVar::FC_threshold[j], hFC[j]->GetMaximum());
        l[j]->SetLineColor(2);
        l[j]->SetLineWidth(2);
        l[j]->SetLineStyle(2);
        l[j]->Draw("same");

    }




}
