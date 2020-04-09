// plot.C
// Plots histograms for all detectors
// To run: root -l "plot.C(XXX)" where XXX is run number
// or: root -l "plot.C(XXX, YYY)" sums all runs from XXX to YYY

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

#include "include/processed.h"
#include "include/processed_QDC.h"
#include "RabVar.h"

using namespace RabVar;

void run_info(int run_num){

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //Variables
    double elapsed_time;
    double FC_int[2];        //integral value
    double FC_irr[2];
    double FC_ratio;
    double FC_ratio_err;
    double FC_ratio_irr;
    double FC_ratio_irr_err;

    double BCI_int;
    double NE_int;
    double NPSD_int;
    double BCI_irr;
    double NE_irr;
    double NPSD_irr;

    //Histos
    TH1F *hFC[num_FC];
    TH1F *hFC_irr[num_FC];
    TH1F *hBCI = new TH1F("hBCI", "hBCI", 16*4096, 0, 16*4096);
    TH1F *hNmon_E = new TH1F("hNmon_E", "hNmon_E", 16*4096, 0, 16*4096);
    TH1F *hNmon_PSD = new TH1F("hNmon_PSD", "hNmon_PSD", 16*4096, 0, 16*4096);
    TH1F *hBCI_irr = new TH1F("hBCI_irr", "hBCI_irr", 16*4096, 0, 16*4096);
    TH1F *hNmon_E_irr = new TH1F("hNmon_E_irr", "hNmon_E_irr", 16*4096, 0, 16*4096);
    TH1F *hNmon_PSD_irr = new TH1F("hNmon_PSD_irr", "hNmon_PSD_irr", 16*4096, 0, 16*4096);
    for (int i=0; i<num_FC; i++){
        hFC[i] = new TH1F(Form("hFC%i",i), Form("FC %i",i+1), 
            16*4096, 0, 16*4096);
        hFC_irr[i] = new TH1F(Form("hFC_irr%i",i), Form("FC_irr %i",i+1), 
            16*4096, 0, 16*4096);
    }
    
    //in file
    processed rabbit(run_num);
    processed_QDC rabbit_QDC(run_num);

    //out file
    FILE *file_ptr;
    file_ptr = fopen("datafiles/run_info.dat","a");
    
    //lines
    TLine *line[num_FC];

    cout << rabbit.rawfile->Get("start_time")->GetTitle() << endl;
    cout << rabbit.rawfile->Get("stop_time")->GetTitle() << endl;

    //loop over SCP data
    cout << "Looping over SCP data..." << endl;
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }

        if (rabbit.ADC[BCI_chn] >10){
            hBCI->Fill(rabbit.ADC[BCI_chn]);
            if ((rabbit.cycle_time>time_irr[0])
              &&(rabbit.cycle_time<time_irr[1])){
                hBCI_irr->Fill(rabbit.ADC[BCI_chn]);
            }
        }
        for (int i=0; i<num_FC; i++){
            if (rabbit.ADC[FC_chn[i]] >10){
                hFC[i]->Fill(rabbit.ADC[FC_chn[i]]);
                if ((rabbit.cycle_time>time_irr[0])
                  &&(rabbit.cycle_time<time_irr[1])){
                    hFC_irr[i]->Fill(rabbit.ADC[FC_chn[i]]);
                }
            } 
        }
    }
    cout << endl;
    nb = rabbit.GetEntry(nentries-1);
    elapsed_time = rabbit.seconds/3600.;

    //loop over QDC data
    cout << "Looping over QDC data..." << endl;
    nentries = rabbit_QDC.fChain->GetEntriesFast();
    nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit_QDC.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }
        if (rabbit_QDC.ADC_long[nmon_chn]>min_nmon_E){
            hNmon_E->Fill(rabbit_QDC.ADC_long[nmon_chn]);
            if ((rabbit_QDC.nmon_PSD>nmon_PSD_cut[0])
                && (rabbit_QDC.nmon_PSD<nmon_PSD_cut[1])){
                hNmon_PSD->Fill(rabbit_QDC.ADC_long[nmon_chn]);
            }
            if ((rabbit_QDC.cycle_time>time_irr[0])
              &&(rabbit_QDC.cycle_time<time_irr[1])){
                hNmon_E_irr->Fill(rabbit_QDC.ADC_long[nmon_chn]);
                if ((rabbit_QDC.nmon_PSD>nmon_PSD_cut[0])
                    && (rabbit_QDC.nmon_PSD<nmon_PSD_cut[1])){
                    hNmon_PSD_irr->Fill(rabbit_QDC.ADC_long[nmon_chn]);
                }
            }
        }
    }
    cout << endl;

    //plot histograms
    TCanvas *cEnergy = new TCanvas("cFC","FC", 800, 800);
    cEnergy->Divide(1,2);

    for (int i=0; i<num_FC; i++){

        cEnergy->cd(i+1);
        //gPad->SetLogy();
        hFC[i]->Draw();
        hFC_irr[i]->SetLineColor(4);
        hFC_irr[i]->Draw("same");
        hFC[i]->GetXaxis()->SetRange(3000,16*4096);
        hFC[i]->GetXaxis()->SetTitle("Chn");
        hFC[i]->GetYaxis()->SetTitle("Counts");

        int max = hFC[i]->GetMaximumBin();
        max = hFC[i]->GetBinContent(max);
        
        line[i] = new TLine(FC_threshold[i], 0., FC_threshold[i], 1.*max);
        line[i]->SetLineColor(2);
        line[i]->SetLineStyle(2);
        line[i]->Draw("same");

        FC_int[i] = hFC[i]->Integral(FC_threshold[i], 65536);
        FC_irr[i] = hFC_irr[i]->Integral(FC_threshold[i], 65536);
        cout << "FC" << i+1 << ":       " << FC_int[i] << endl;
        cout << "FC" << i+1 << " irr:   " << FC_irr[i] << endl;
        
    }

    BCI_int = hBCI->Integral(10, 65000);
    NE_int = hNmon_E->Integral(10, 65000);
    NPSD_int = hNmon_PSD->Integral(10, 65000);
    BCI_irr = hBCI_irr->Integral(10, 65000);
    NE_irr = hNmon_E_irr->Integral(10, 65000);
    NPSD_irr = hNmon_PSD_irr->Integral(10, 65000);

    FC_ratio = FC_int[1]/FC_int[0];
    FC_ratio_err = FC_ratio*sqrt(pow(FC_int[0],-1)+pow(FC_int[1], -1));

    FC_ratio_irr = FC_irr[1]/FC_irr[0];
    FC_ratio_irr_err = FC_ratio_irr*sqrt(pow(FC_irr[0],-1)+pow(FC_irr[1], -1));

    cout << "N-mon E:       "  << NE_int << endl;
    cout << "N-mon E irr:   "  << NE_irr << endl;
    cout << "N-mon PSD:     "  << NPSD_int << endl;
    cout << "N-mon PSD irr: "  << NPSD_irr << endl;
    cout << "BCI:           "  << BCI_int << endl;
    cout << "BCI irr:       "  << BCI_irr << endl;
    
    //output to file
    fprintf(file_ptr, "%i\t", run_num);
    fprintf(file_ptr, "\t"); // target
    fprintf(file_ptr, "\t"); // FC 
    fprintf(file_ptr, "\t"); // cycle time 
    fprintf(file_ptr, "%s\t", rabbit.rawfile->Get("start_time")->GetTitle());
    fprintf(file_ptr, "%s\t", rabbit.rawfile->Get("stop_time")->GetTitle());
    fprintf(file_ptr, "%f\t", elapsed_time);
    //entire run
    for (int i=0; i<num_FC; i++){
        fprintf(file_ptr,"%f\t", FC_int[i]); 
    }
    fprintf(file_ptr,"%f\t", FC_ratio); 
    fprintf(file_ptr,"%f\t", FC_ratio_err); 

    fprintf(file_ptr,"%f\t", NE_int); 
    fprintf(file_ptr,"%f\t", NPSD_int); 
    fprintf(file_ptr,"%f\t", BCI_int); 
    //beam on only
    for (int i=0; i<num_FC; i++){
        fprintf(file_ptr,"%f\t", FC_irr[i]); 
    }
    fprintf(file_ptr,"%f\t", FC_ratio_irr); 
    fprintf(file_ptr,"%f\t", FC_ratio_irr_err); 
    fprintf(file_ptr,"%f\t", NE_irr); 
    fprintf(file_ptr,"%f\t", NPSD_irr); 
    fprintf(file_ptr,"%f\n", BCI_irr); 
    fclose(file_ptr);
    cout << "Run " << run_num << " Output complete!" << endl;

}
