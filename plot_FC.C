
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

void plot_FC(int run_num){


    //Variables
    const int num_det = 2;
    const int rebin = 1;
    int dets[num_det] = {10, 11};
    int threshold[num_det] = {7000, 7000};
    TLine *l[num_det];
    
    //Histograms
    TH1F *hFC[num_det];

    //get histos
    TFile *fHist;
    if (run_num<10){
        fHist = new TFile(Form("data_root/RABITTS_00%i.root", run_num));
    }
    else if (run_num<100){
        fHist = new TFile(Form("data_root/RABITTS_0%i.root", run_num));
    }
    else{
        fHist = new TFile(Form("data_root/RABITTS_%i.root", run_num));
    }

    TCanvas *cFC = new TCanvas("cFC","Summed segments",1000, 400);
    cFC->Divide(num_det);

    for (int j=0; j<num_det; j++){
        cFC->cd(j+1);

        hFC[j] = (TH1F*) (fHist->Get(Form("histos_SCP/hADC%i", dets[j])))->Clone();
        hFC[j]->SetName(Form("run%i_Det%i", run_num, j+1));
        cout << hFC[j]->Integral(threshold[j], 65535) << endl;
        hFC[j]->Rebin(16);
        hFC[j]->GetXaxis()->SetRangeUser(5000, 65535);
        hFC[j]->Draw();

        l[j] = new TLine(threshold[j], 0, threshold[j], hFC[j]->GetMaximum());
        l[j]->SetLineColor(2);
        l[j]->SetLineWidth(2);
        l[j]->SetLineStyle(2);
        l[j]->Draw("same");
    }




}
