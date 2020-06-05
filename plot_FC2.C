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

#include "include/MDPP16_SCP.hh"
#include "include/RabVar.hh"


void plot_FC2(int run_num){

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //Variables
    double percent_PU[RabVar::num_FC];
    double counts[RabVar::num_FC];
    double counts_PU[RabVar::num_FC];

    TLine *lThresh[RabVar::num_FC]; 

    TH1F *hFC[RabVar::num_FC];
    TH1F *hFC_pu[RabVar::num_FC];
    TH1F *hFC_rej[RabVar::num_FC];

    for (int fc=0; fc<RabVar::num_FC; fc++){
        hFC[fc] = new TH1F(Form("hFC%i", fc), Form("hFC%i", fc), 16*4096, 0, 16*4096);
        hFC_pu[fc] = new TH1F(Form("hFC_pu%i", fc), Form("hFC_pu%i", fc), 16*4096, 0, 16*4096);
        hFC_rej[fc] = new TH1F(Form("hFC_rej%i", fc), Form("hFC_rej%i", fc), 16*4096, 0, 16*4096);
    }

    //in file
    MDPP16_SCP rabbit(run_num);
    
    //loop over SCP data
    cout << "Looping over SCP data..." << endl;
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }

        for (int fc=0; fc<RabVar::num_FC; fc++){
            if (rabbit.TDC[RabVar::FC_chn[fc]] >1){
                hFC[fc]->Fill(rabbit.ADC[RabVar::FC_chn[fc]]);

                if (!(rabbit.pileup[RabVar::FC_chn[fc]])){
                    hFC_pu[fc]->Fill(rabbit.ADC[RabVar::FC_chn[fc]]);
                }
                else if ((rabbit.pileup[RabVar::FC_chn[fc]])){
                    hFC_rej[fc]->Fill(rabbit.ADC[RabVar::FC_chn[fc]]);
                }

            }
        }
    }
    cout << endl;

    TCanvas *cFC = new TCanvas("cFC","Summed segments", 800, 1200);
    cFC->Divide(1, RabVar::num_FC);

    cout << "----- Run " << run_num << "-----" << endl;
    cout << "      All counts\tPU reject\tPercent PU" << endl;
    for (int j=0; j<RabVar::num_FC; j++){

        counts[j] = hFC[j]->Integral(RabVar::FC_threshold[j], 65500);
        counts_PU[j] = hFC_pu[j]->Integral(RabVar::FC_threshold[j], 65500);
        percent_PU[j] = 100*(counts[j]-counts_PU[j])/counts[j];
        
        cout << "FC" << RabVar::FC_chn[j] << ":\t" 
             << counts[j] << "\t"
             << counts_PU[j] << "\t"
             << percent_PU[j] << "%" << endl;

        hFC[j]->Rebin(4*RabVar::FC_rebin);
        hFC_pu[j]->Rebin(4*RabVar::FC_rebin);
        hFC_rej[j]->Rebin(4*RabVar::FC_rebin);

        hFC[j]->GetXaxis()->SetRangeUser(2000, 65500);
        hFC[j]->GetXaxis()->SetRangeUser(100, 65500);

        cFC->cd(j+1);
        gPad->SetLogy();
        hFC[j]->Draw();
        hFC_pu[j]->SetLineColor(2);
        hFC_pu[j]->Draw("same");
        hFC_rej[j]->SetLineColor(4);
        hFC_rej[j]->Draw("same");

        lThresh[j] = new TLine(RabVar::FC_threshold[j], 0, RabVar::FC_threshold[j], hFC[j]->GetMaximum());
        lThresh[j]->SetLineColor(2);
        lThresh[j]->SetLineWidth(2);
        lThresh[j]->Draw("same");

    }


}
