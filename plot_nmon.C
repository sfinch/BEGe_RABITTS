///////////////////////////////////////////////////////////////////////////
// plot_nmon.C
// Plots the detector count rate as a function of cycle time for the BCI, neutron monitor
// with and without PSD cut, the two fission chambers, and the two BEGe detectors. 
// The spectra are saved in time_spectra for quick reference
// The channel numbers, thresholds, and PSD cuts are all stores in RabVar.C
// Requries running of mvme2root, followed by process_rabbit
//
// To run: root -l "plot_nmon.C(XXX)" where XXX is run number
//
//////////////////////////////////////////////////////////////////////////


#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "include/processed.h"
#include "include/processed_QDC.h"
#include "plot_FC.C"
#include "RabVar.h"

using namespace RabVar;

using std::cout;
using std::cerr;
using std::endl;

void plot_nmon(int run_num){

    //Histograms
    TH1F *hCycFC[num_FC];
    TH1F *hHPGe[num_BEGe];
    TH1F *hCycNmon  = new TH1F("hCycNmon", "Neutron monitor", 
                               cycle_time_bins, min_cycle_time, max_cycle_time);
    TH1F *hCycNPSD  = new TH1F("hCycNPSD", "Neutron monitor with PSD", 
                               cycle_time_bins, min_cycle_time, max_cycle_time);
    TH1F *hCycBCI   = new TH1F("hCycBCI", "BCI",
                               cycle_time_bins, min_cycle_time, max_cycle_time);
    for (int i=0; i<num_FC; i++){
        hCycFC[i] = new TH1F(Form("hCycFC%i", i), Form("FC %i", i), 
                             cycle_time_bins, min_cycle_time, max_cycle_time);
    }
    for (int i=0; i<num_BEGe; i++){
        hHPGe[i] = new TH1F(Form("hHPGe%i", i), Form("BEGe %i", i),
                            cycle_time_bins, min_cycle_time, max_cycle_time);
    }

    //TLines
    TLine *lIrr[2][5];
    TLine *lCount[2][5];

    //in file
    processed rabbit(run_num);
    processed_QDC rabbit_QDC(run_num);

    //loop over SCP data
    cout << "Looping over SCP data..." << endl;
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }
        if (rabbit.ADC[BCI_chn]>min_BCI){
            hCycBCI->Fill(rabbit.cycle_time);
        }
        for (int i=0; i<num_FC; i++){
            if (rabbit.ADC[FC_chn[i]]>FC_threshold[i]){
                hCycFC[i]->Fill(rabbit.cycle_time);
            }
        }
        for (int i=0; i<num_FC; i++){
            if (rabbit.En[i]>min_BEGe_E[i]){
                hHPGe[i]->Fill(rabbit.cycle_time);
            }
        }

    }//end loop over events
    cout << endl;

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
            hCycNmon->Fill(rabbit_QDC.cycle_time);
            if ((rabbit_QDC.nmon_PSD>nmon_PSD_cut[0])
                && (rabbit_QDC.nmon_PSD<nmon_PSD_cut[1])){
                hCycNPSD->Fill(rabbit_QDC.cycle_time);
            }
        }

    }//end loop over events
    cout << endl;

    for (int i=0; i<2; i++){
        lIrr[i][0] = new TLine(time_irr[i], 0, time_irr[i], hCycBCI->GetMaximum());
        lIrr[i][1] = new TLine(time_irr[i], 0, time_irr[i], hCycNmon->GetMaximum());
        lIrr[i][2] = new TLine(time_irr[i], 0, time_irr[i], hCycFC[0]->GetMaximum());
        lIrr[i][3] = new TLine(time_irr[i], 0, time_irr[i], hCycFC[1]->GetMaximum());
        lIrr[i][4] = new TLine(time_irr[i], 0, time_irr[i], hHPGe[0]->GetMaximum());

        lCount[i][0] = new TLine(time_count[i], 0, time_count[i], hCycBCI->GetMaximum());
        lCount[i][1] = new TLine(time_count[i], 0, time_count[i], hCycNmon->GetMaximum());
        lCount[i][2] = new TLine(time_count[i], 0, time_count[i], hCycFC[0]->GetMaximum());
        lCount[i][3] = new TLine(time_count[i], 0, time_count[i], hCycFC[1]->GetMaximum());
        lCount[i][4] = new TLine(time_count[i], 0, time_count[i], hHPGe[0]->GetMaximum());

        for (int j=0; j<5; j++){
            lIrr[i][j]->SetLineColor(2);
            lIrr[i][j]->SetLineWidth(2);
            lIrr[i][j]->SetLineStyle(2);

            lCount[i][j]->SetLineColor(4);
            lCount[i][j]->SetLineWidth(2);
            lCount[i][j]->SetLineStyle(2);
        }
    }

    //plot
    TCanvas *cCycle = new TCanvas("cCycle", "Cycle time counts", 800, 1000);
    cCycle->Divide(1,5);

    cCycle->cd(1);
    hCycBCI->Draw();
    gPad->SetLogy();
    for (int i=0; i<2; i++){
        lIrr[i][0]->Draw("same");
        lCount[i][0]->Draw("same");
    }

    cCycle->cd(2);
    hCycNmon->Draw();
    hCycNPSD->Draw("same");
    gPad->SetLogy();
    for (int i=0; i<2; i++){
        lIrr[i][1]->Draw("same");
        lCount[i][1]->Draw("same");
    }

    cCycle->cd(3);
    hCycFC[0]->Draw();
    for (int i=0; i<2; i++){
        lIrr[i][2]->Draw("same");
        lCount[i][2]->Draw("same");
    }

    cCycle->cd(4);
    hCycFC[1]->Draw();
    for (int i=0; i<2; i++){
        lIrr[i][3]->Draw("same");
        lCount[i][3]->Draw("same");
    }

    cCycle->cd(5);
    for (int i=0; i<num_BEGe; i++){
        hHPGe[i]->SetLineColor(1+i);
        hHPGe[i]->Draw("same");
    }
    for (int i=0; i<2; i++){
        lIrr[i][4]->Draw("same");
        lCount[i][4]->Draw("same");
    }

    plot_FC(run_num);

    //save to file
    cCycle->SaveAs(Form("time_spec/run%i.png", run_num));
    cCycle->SaveAs(Form("time_spec/run%i.C", run_num));
    
}
