
#include <iostream>
#include <fstream>
#include <vector>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "include/MDPP16_SCP.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::vector;

void check_PU(int run_num){

    TH1F* hEn[2];
    TH1F* hEnPU[2];
    TFile *fRABITTS = new TFile(Form("data_root/RABITTS_0%i.root", run_num), "READ");
    TTree *tRABITTS = (TTree*)fRABITTS->Get("MDPP16_SCP");
    int det[2] = {0, 2};

    for (int i=0; i<2; i++){
        hEn[i] = new TH1F(Form("hEn%i", i), Form("hEn%i", i), 40000, 1, 40001);
        hEnPU[i] = new TH1F(Form("hEnPU%i", i), Form("hEnPU%i", i), 40000, 1, 40001);
    }


    TCanvas *c1 = new TCanvas();
    c1->Divide(2,1);

    for (int i=0; i<2; i++){
        tRABITTS->Project(Form("hEn%i", i), Form("ADC[%i]", i));
        tRABITTS->Project(Form("hEnPU%i", i), Form("ADC[%i]", det[i]), 
                          Form("pileup[%i]==0", det[i]));

        c1->cd(i+1);
        hEn[i]->Draw();
        hEnPU[i]->SetLineColor(2);
        hEnPU[i]->Draw("same");
    }


}

