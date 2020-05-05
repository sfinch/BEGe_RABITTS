////////////////////////////////////////////////////////////////////////////////////////
// beam_corr.C
// Calculates the beam fluctuation correction for an activation run. Does this using the
// (1) BCI, the (2) neutron monitor, and the (3)neutron monitor with a PSD cut. The values
// for these cuts are all defined in RabVar.h
// To run: root -l "beam_corr.C(XXX, YYY)" where XXX and YY are run numbers
// Will calculate the beam correction for runs XXX through YYY inclusive 
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
#include "TDatime.h"

#include "include/processed.hh"
#include "include/processed_QDC.hh"
#include "include/RabVar.hh"


class correction{
  public:
    double dt = 1; 
    double lambda = 1; 

    double activity = 0;
    double beam_avg = 0;
    double counts = 0;
    
    double activity_irr = 0;
    double beam_avg_irr = 0;
    double counts_irr = 0;

    double DC_corr = 0;
    int total_cycles = 0;

    vector<int> scalers;
    vector<int> scalers_irr;
    vector<bool> counting;
    vector<bool> irradiation;

    vector<float> vec_act;
    vector<float> vec_counts;

    correction();
    correction(double t);

    void resize(int length);
    void calc(double l);
    void DCcalc(double l);
    void write_scalers(char* filename);
    void plot_hist();
    void print();
};

correction::correction(){
}

correction::correction(double t){
    dt = t;
};

void correction::DCcalc(double l){
    lambda = l;

    double cycleT = RabVar::irr_time + RabVar::count_time + 2*RabVar::transit_time;
    DC_corr = (1/lambda)
             *(1-exp(-1.*lambda*RabVar::irr_time))
             *exp(-1.*lambda*RabVar::transit_time)
             *(1-exp(-1.*lambda*RabVar::count_time))
             *((total_cycles/(1-exp(-1.*lambda*cycleT))) 
              - (exp(-1.*lambda*cycleT)*(1-exp(-1.*total_cycles*lambda*cycleT))
                 /pow((1-exp(-1.*lambda*cycleT)),2)));

};

void correction::calc(double l){
    
    lambda = l;
    DCcalc(lambda);
    int num_beam_events = 0;

    num_beam_events = 0;
    activity = 0;
    beam_avg = 0;
    counts = 0;
    vec_act.resize(0);
    vec_counts.resize(0);
    for (int i=0; i<scalers.size(); i++){
        activity = activity*(exp(-1.*lambda*dt))
                 + scalers.at(i)*(1-exp(-1*lambda*dt));
        beam_avg += scalers.at(i);
        counts += counting.at(i)*activity*(1-exp(-1*lambda*dt))/lambda;

        if (irradiation.at(i)){
            num_beam_events++;
        }

        vec_act.push_back(activity);
        vec_counts.push_back(counts);
    }
    beam_avg = beam_avg/num_beam_events;

    num_beam_events = 0;
    activity_irr = 0;
    beam_avg_irr = 0;
    counts_irr = 0;
    for (int i=0; i<scalers_irr.size(); i++){
        activity_irr = activity_irr*(exp(-1.*lambda*dt))
                 + scalers_irr.at(i)*(1-exp(-1*lambda*dt));
        beam_avg_irr += scalers_irr.at(i);
        counts_irr += counting.at(i)*activity_irr*(1-exp(-1*lambda*dt))/lambda;
        
        if (irradiation.at(i)){
            num_beam_events++;
        }
    }
    beam_avg_irr = beam_avg_irr/num_beam_events;

    counts = counts/beam_avg;
    counts_irr = counts_irr/beam_avg_irr;

};

void correction::resize(int length){
    scalers.resize(length, 0);
    scalers_irr.resize(length, 0);
    counting.resize(length, 0);
    irradiation.resize(length, 0);
}

void correction::plot_hist()
{
    TH1F *hAct = new TH1F("hAct", "Target activity", scalers.size(), 0, scalers.size()*dt);
    TH1F *hCounts = new TH1F("hCounts", "Decays in detector", scalers.size(), 0, scalers.size()*dt);
    TH1F *hCounting = new TH1F("hCounting", "hCounting", scalers.size(), 0, scalers.size()*dt);
    TH1F *hIrr = new TH1F("hIrr", "hIrr", scalers.size(), 0, scalers.size()*dt);

    for (int i=0; i<scalers.size(); i++){
        hAct->SetBinContent(i, vec_act.at(i)/beam_avg);
        hCounts->SetBinContent(i, vec_counts.at(i)/beam_avg);
        hCounting->SetBinContent(i, counting.at(i)*vec_act.at(i)/beam_avg);
        hIrr->SetBinContent(i, irradiation.at(i)*vec_act.at(i)/beam_avg);
    }
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->Divide(1,2);

    c1->cd(1);
    hAct->GetXaxis()->SetTitle("Time (s)");
    hAct->SetLineWidth(2);
    hAct->Draw();

    hCounting->SetFillColor(4);
    hCounting->Draw("same");

    hIrr->SetFillColor(2);
    hIrr->Draw("same");

    c1->cd(2);
    hCounts->SetLineWidth(2);
    hCounts->Draw();
    hCounts->GetXaxis()->SetTitle("Time (s)");

}

void correction::write_scalers(char* filename){
    //output scalers to file
    FILE *file_ptr = fopen(filename, "w");
    for (int i=0; i<scalers.size(); i++){
        fprintf(file_ptr, "%i\t", i);
        fprintf(file_ptr, "%f\t", i*dt);
        fprintf(file_ptr, "%i\t", scalers.at(i));
    }
    fclose(file_ptr);
}

void correction::print(){
    cout << "Full rate:              "  << counts
         << "\t\t" << counts/DC_corr << endl;
    cout << "Irradiation only:       "  << counts_irr
         << "\t\t" << counts_irr/DC_corr << endl;
}

void cyclic_corr(int run_num, int run_num2 = 0){

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //Variables
    double half_life = 10;
    //double half_life = 2.66*3600;
    double lambda = 0.69314718/half_life;
    double dt = 0.01;

    if (run_num2 < run_num){
        run_num2 = run_num;
    }
    const int num_runs = run_num2 - run_num + 1;

    double elapsed_time[num_runs];
    double SCP_time[num_runs];
    double QDC_time[num_runs];
    double total_time;

    int num_cycles[num_runs];
    int total_cycles = 0;

    int numTbins = 0;
    int startOffset = 0;

    int irr_length = RabVar::irr_time/dt;
    int count_length = RabVar::count_time/dt;
    int tbin = 0;
    double maxT = 0;

    TDatime tStart[num_runs];
    TDatime tStop[num_runs];
    vector<double> irr_start;

    Long64_t nentries;
    Long64_t nbytes = 0, nb = 0;

    correction corr_BCI(dt);
    correction corr_nmon(dt);
    correction corr_nPSD(dt);

    //in file
    processed *rabbit[num_runs];
    processed_QDC *rabbit_QDC[num_runs];
    for (int i=0; i<num_runs; i++){
        rabbit[i] = new processed(i+run_num);
        rabbit_QDC[i] = new processed_QDC(i+run_num);

        tStart[i] = TDatime(rabbit[i]->rawfile->Get("start_time")->GetTitle());
        tStop[i] = TDatime(rabbit[i]->rawfile->Get("stop_time")->GetTitle());
        elapsed_time[i] = tStop[i].Convert() - tStart[i].Convert();

        num_cycles[i] = ((TH1F*)rabbit[i]->file->Get("histos/hIrrTime"))->Integral();
        total_cycles += num_cycles[i];
    }
    total_time = tStop[num_runs-1].Convert() - tStart[0].Convert();

    numTbins =  total_time/dt;
    corr_BCI.resize(numTbins);
    corr_nmon.resize(numTbins);
    corr_nPSD.resize(numTbins);
    corr_BCI.total_cycles = total_cycles;
    corr_nmon.total_cycles = total_cycles;
    corr_nPSD.total_cycles = total_cycles;

    //loop over runs
    for (int run=0; run<num_runs; run++){
        cout << "------------------------------------------------------" << endl;
        cout << "Run " << run+run_num << " start:        " << tStart[run].AsString() << endl;
        cout << "Run " << run+run_num << " stop:         " << tStop[run].AsString() << endl;
        cout << "Run " << run+run_num << " elapsed time: " 
             << elapsed_time[run] << " s, " << elapsed_time[run]/3600. << " h" << endl;
        cout << "Run " << run+run_num << " num cycles:   " << num_cycles[run] << endl;
        cout << "------------------------------------------------------" << endl;
        
        startOffset = tStart[run].Convert() - tStart[0].Convert();

        //loop over SCP data
        cout << "Looping over SCP data..." << endl;
        nentries = rabbit[run]->fChain->GetEntriesFast();
        nbytes = 0, nb = 0;
        nb = rabbit[run]->GetEntry(nentries-1);
        SCP_time[run] = rabbit[run]->seconds;
        
        for (Long64_t jentry=0; jentry<nentries; jentry++) {
            nb = rabbit[run]->GetEntry(jentry);   nbytes += nb;
            if (jentry%100000==0){
                cout << '\r' << "Processing event " << jentry;
            }
            if (rabbit[run]->seconds > SCP_time[run]){
                SCP_time[run] = rabbit[run]->seconds;
            }

            if (rabbit[run]->ADC[RabVar::BCI_chn] > RabVar::min_BCI){
        
                corr_BCI.scalers.at(int((startOffset + rabbit[run]->seconds)/dt))++;

                if ((rabbit[run]->cycle_time > RabVar::time_irr[0])
                  &&(rabbit[run]->cycle_time < RabVar::time_irr[1])){

                    corr_BCI.scalers_irr.at(int((startOffset + rabbit[run]->seconds)/dt))++;
                }
            }
        } //end loop over SCP
        cout << '\r' << nentries << " total SCP events" << endl;
        cout << SCP_time[run] << " s SCP time" << endl;

        //loop over QDC data
        cout << "Looping over QDC data..." << endl;
        nentries = rabbit_QDC[run]->fChain->GetEntriesFast();
        nbytes = 0, nb = 0;
        nb = rabbit_QDC[run]->GetEntry(nentries-1);
        QDC_time[run] = rabbit_QDC[run]->seconds;

        for (Long64_t jentry=0; jentry<nentries; jentry++) {
            nb = rabbit_QDC[run]->GetEntry(jentry);   nbytes += nb;
            if (jentry%100000==0){
                cout << '\r' << "Processing event " << jentry;
            }
            if (rabbit_QDC[run]->seconds > QDC_time[run]){
                QDC_time[run] = rabbit_QDC[run]->seconds;
            }
            if (rabbit_QDC[run]->ADC_long[RabVar::nmon_chn] > RabVar::min_nmon_E){

                corr_nmon.scalers.at(int((startOffset + rabbit_QDC[run]->seconds)/dt))++;

                if ((rabbit_QDC[run]->nmon_PSD > RabVar::nmon_PSD_cut[0])
                  && (rabbit_QDC[run]->nmon_PSD < RabVar::nmon_PSD_cut[1])){

                    corr_nPSD.scalers.at(int((startOffset + rabbit_QDC[run]->seconds)/dt))++;

                }
                if ((rabbit_QDC[run]->cycle_time > RabVar::time_irr[0])
                  &&(rabbit_QDC[run]->cycle_time < RabVar::time_irr[1])){

                    corr_nmon.scalers_irr.at(int((startOffset + rabbit_QDC[run]->seconds)/dt))++;

                    if ((rabbit_QDC[run]->nmon_PSD > RabVar::nmon_PSD_cut[0])
                      && (rabbit_QDC[run]->nmon_PSD < RabVar::nmon_PSD_cut[1])){

                        corr_nPSD.scalers_irr.at(int((startOffset + rabbit_QDC[run]->seconds)/dt))++;

                    }
                }
            }

        } //end loop over QDC
        cout << '\r' << nentries << " total QDC events" << endl;
        cout << QDC_time[run] << " s QDC time" << endl;
        
        maxT = SCP_time[run];
        if (QDC_time[run] > maxT){
            maxT = QDC_time[run];
        }

        irr_start.resize(0);
        for (int ncycle=0; ncycle<rabbit[run]->irr_start_times->GetNoElements(); ncycle++){
            irr_start.push_back( ( (*rabbit[run]->irr_start_times) )[ncycle] );
        }

        for (int ncycle=0; ncycle<irr_start.size(); ncycle++){
            tbin = (irr_start.at(ncycle) + RabVar::time_irr[0])/dt;
            for (int t=0; t<irr_length; t++){
                if (tbin+t < maxT/dt){
                    corr_BCI.irradiation.at((startOffset/dt) + tbin + t) = 1;
                    corr_nmon.irradiation.at((startOffset/dt) + tbin + t) = 1;
                    corr_nPSD.irradiation.at((startOffset/dt) + tbin + t) = 1;
                }
            }

            tbin = (irr_start.at(ncycle) + RabVar::time_count[0])/dt;
            for (int t=0; t<count_length; t++){
                if (tbin+t < maxT/dt){
                    corr_BCI.counting.at((startOffset/dt) + tbin + t) = 1;
                    corr_nmon.counting.at((startOffset/dt) + tbin + t) = 1;
                    corr_nPSD.counting.at((startOffset/dt) + tbin + t) = 1;
                }
            }
        }

    } //end loop over runs

    if (num_runs == 1){
        corr_BCI.write_scalers(Form("datafiles/run%i.scalers", run_num));
        corr_nmon.write_scalers(Form("datafiles/run%i.scalers", run_num));
        corr_nPSD.write_scalers(Form("datafiles/run%i.scalers", run_num));
    }
    else{
        corr_BCI.write_scalers(Form("datafiles/runs%i-%i.scalers", run_num, run_num2));
        corr_nmon.write_scalers(Form("datafiles/runs%i-%i.scalers", run_num, run_num2));
        corr_nPSD.write_scalers(Form("datafiles/runs%i-%i.scalers", run_num, run_num2));
    }

    //calculate corrections
    corr_BCI.calc(lambda);
    corr_nmon.calc(lambda);
    corr_nPSD.calc(lambda);

    corr_BCI.plot_hist();

    //print corrections
    cout << "------------------------------------------------------" << endl;
    cout << "Total start:    " << tStart[0].AsString() << endl;
    cout << "Total stop:     " << tStop[num_runs-1].AsString() << endl;
    cout << "Total cycles:   " << total_cycles << endl;
    cout << "Elapsed time:   " << total_time/3600. << " h" << endl;
    cout << "------------------------------------------------------" << endl;
    cout.precision(5);
    cout << "\t\t\t" << "Produced \tCorrection %" << endl;
    cout << "Assuming DC beam:       "  << corr_BCI.DC_corr
         << "\t\t" << corr_BCI.DC_corr/corr_BCI.DC_corr << endl;

    cout << "--- BCI ---"  << endl;
    corr_BCI.print();
    cout << "--- n-mon ---"  << endl;
    corr_nmon.print();
    cout << "--- n-mon with PSD cut ---"  << endl;
    corr_nPSD.print();


}
