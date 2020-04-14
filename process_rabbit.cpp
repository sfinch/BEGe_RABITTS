///////////////////////////////////////////////////////////////////////////
// process_rabbit
// Reads in the MDPP16_SCP and MDPP_QDC trees from mvme2root, then outputs a processed
// tree with branches for energy calibrated HPGe detecors, cycle time, and neutron PSD. 
// output ROOT files saved in data_processed
// Experimental configuration set in RabVar.h
// Requries running of mvme2root. 
//
// To run: first make, then
// ./process_rabbit [-v] <run numbers>
//
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <vector>

#include <TStyle.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "RabVar.hh"
#include "MDPP16_SCP.hh"
#include "MDPP16_QDC.hh"

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::vector;


struct calibration{     //struct storing calibration data for both dets
    Float_t m[RabVar::num_det];
    Float_t b[RabVar::num_det];
};
    
calibration read_in_cal(int run_num); // function to read in calibration from file
int process_rabbit(int run_num, bool opt_verbose);

int main(int argc, char *argv[]){

    bool optverbose = 0;
    int startindex = 1;
    int num_errors = 0;

    if (argc<2)
    {
        cerr << "Invalid number of arguments" << endl;
        cerr << "Usage: " << argv[0] << " [-v] <run numbers>" << endl;
        return 1;
    }
    
    if (!strcmp(argv[1], "-v")) //verbose option
    {
        optverbose = 1;
        startindex++;
    }

    //loop over all given files
    for (int file=startindex; file<argc; file++){
        cout << "----- Processing run number " << argv[file] << " -----" << endl;

        //process mvmelst file
        try
        {
            process_rabbit(atoi(argv[file]), optverbose);
        }
        catch (const std::exception &e)
        {
            cerr << "Error processing run: " << e.what() << endl;
            num_errors++;
        }

        cout << "----- Run " << argv[file] << " complete -----" << endl;
    }
    cout << "----- " << argc-startindex-num_errors << "/" << argc-startindex << " files converted -----" << endl;

    return 0;
}

int process_rabbit(int run_num, bool opt_verbose){

    //cycle time variables
    bool source_run = 1;
    int pos = 0;            //-1 = counting, 0 not def, 1 = irradiation

    double last_move = 0;
    double last_irr_length = 0;
    double last_count_length = 0;
    double new_move = 0;
    vector<double> irr_start;
    vector<double> count_start;

    double start_offset = 0;
    double end_offset = 0;
    double start_event= 0;
    double end_event= 0;

    int num_irr = 0;
    int num_count = 0;
    int num_missed_irr = 0;
    int num_missed_count = 0;

    //in file
    MDPP16_SCP rabbit(run_num);
    MDPP16_QDC rabbit_QDC(run_num);
    //check for any additional calibration corrections
    calibration cal = read_in_cal(run_num);

    //output tree variables
    Float_t cycle_time = 0;
    Float_t nmon_PSD;
    Float_t En[RabVar::num_det];
    TRandom3 r;

    //out file
    TFile *fOut = new TFile(Form("data_processed/RABITTS_%i.root", run_num), "RECREATE");
    TTree *tProcessed = new TTree("processed", "Processed rabbit data");
    TTree *tProcessed_QDC = new TTree("processed_QDC", "Processed rabbit QDC data");

    tProcessed->Branch(Form("En[%i]", RabVar::num_det), &En, Form("En[%i]/F", RabVar::num_det));
    tProcessed->Branch("cycle_time", &cycle_time, "cycle_time/F");

    tProcessed_QDC->Branch("cycle_time", &cycle_time, "cycle_time/F");
    tProcessed_QDC->Branch(Form("nmon_PSD"), &nmon_PSD, Form("nmon_PSD/F"));

    //Histrograms
    TH1F *hIrrTime = new TH1F("hIrrTime", "hIrrTime", 100*100, 0, 100);
    TH1F *hCountTime = new TH1F("hCountTime", "hCountTime", 100*100, 0, 100);
    TH1F *hCycleTime = new TH1F("hCycleTime", "hCycleTime", 10*110, -10, 100);
    TH1F *hNmon_PSD = new TH1F("hNmon_PSD", "hNmon_PSD", 2500, -5, 5);
    TH1F *hEn[RabVar::num_det];
    for (int det=0; det<RabVar::num_det; det++){
        hEn[det] = new TH1F(Form("hEn%i", det), Form("hEn%i", det), 10*3000, 0, 3000);
    }

    //check if Rabbit used during run
    TH1F *hRabbit = (TH1F*)rabbit.file->Get(Form("histos_SCP/hADC%i", RabVar::rabbit_chn));
    int counts = hRabbit->Integral(10,65000);
    if (counts>1){
        source_run = 0;
    }

    //loop over data
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nentries_QDC = rabbit_QDC.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    //find start offset
    if (!(source_run)){
        for (Long64_t jentry=0; jentry<nentries;jentry++) { rabbit.GetEntry(jentry);
            if (rabbit.ADC[RabVar::rabbit_chn]>10){
                start_offset = rabbit.seconds; 
                start_event = jentry-1; 
                jentry = nentries;
            }
        }
    }
    cout << start_offset << " seconds data before RABITT turned on." << endl;

    //Find ending offset
    end_event = nentries;
    if (!(source_run)){
        for (Long64_t jentry=nentries-1; jentry>0; jentry--) {
            rabbit.GetEntry(jentry);
            if (rabbit.ADC[RabVar::rabbit_chn]>10){
                end_offset = rabbit.seconds; 
                end_event = jentry+1; 
                jentry = 0;
            }
        }
    }
    rabbit.GetEntry(nentries-1);
    cout << rabbit.seconds-end_offset << " seconds data after RABITT turned off." << endl;

    //Find all movements
    cout << "Finding all motor movements..." << endl;
    if (source_run){
        cout << "...No motor movements: this was as source or overnight run!" << endl;
    }
    else if (!(source_run)){
        for (Long64_t jentry=start_event; jentry<end_event; jentry++) {
            nb = rabbit.GetEntry(jentry);   nbytes += nb;
            if (jentry%100000==0){
                if (!opt_verbose){
                    cout << '\r' << "Processing event " << jentry;
                }
            }

            //calculate cycle time
            if (pos == 0){//Rabbit has yet to be moved
                if ((rabbit.ADC[RabVar::rabbit_chn]>10)||(rabbit.TDC[RabVar::rabbit_chn]>10)){

                    pos = 1;
                    irr_start.push_back(rabbit.seconds);
                    //pos = -1;
                    //count_start.push_back(rabbit.seconds);
                    
                    last_move = rabbit.seconds;
                }
            }
            else{//rabbit has moved
                if ((rabbit.ADC[RabVar::rabbit_chn]>10)||(rabbit.TDC[RabVar::rabbit_chn]>10)){
                    //rabbit currently moving
                    if ((rabbit.seconds-last_move)>RabVar::min_time){ //min time filter

                        if (pos == 1){          //irradiation
                            if ((rabbit.seconds-irr_start.back())>(last_irr_length*RabVar::min_var)){//min var filter
                                pos = -1;
                                last_irr_length = rabbit.seconds - irr_start.back();
                                count_start.push_back(rabbit.seconds);
                                last_move = rabbit.seconds;

                                hIrrTime->Fill(last_irr_length);
                                num_irr++;
                            }
                        }
                        else if (pos == -1){    //counting
                            if ((rabbit.seconds-count_start.back())>(last_count_length*RabVar::min_var)){//min var filter
                                pos = 1;
                                last_count_length = rabbit.seconds - count_start.back();
                                irr_start.push_back(rabbit.seconds);
                                last_move = rabbit.seconds;

                                hCountTime->Fill(last_count_length);
                                num_count++;
                            }
                        }

                        if (opt_verbose){
                            cout << pos << "\t" << rabbit.seconds<< "\t"
                                 << rabbit.seconds-irr_start.back() << "\t" 
                                 << rabbit.seconds-count_start.back() << endl;
                            cout << "\t" << cycle_time << "\t"
                                 << last_irr_length << "\t" 
                                 << last_count_length << endl;
                        }
                    }

                }
                
                //calculate cycle time
                cycle_time = rabbit.seconds - irr_start.back();

                //Check for missed transition
                if ((last_irr_length>RabVar::min_time) && (last_count_length>RabVar::min_time)){
                    if (rabbit.seconds<end_offset){
                        //not at the end of the run, after rabbit turned off
                        if ((pos==1) && 
                           ((rabbit.seconds-irr_start.back())>(last_irr_length*RabVar::max_var))){
                            //missed irradiation pulse
                            if (opt_verbose){
                                cout << "== MISSED COUNTING PULSE ========" << endl;
                                cout << jentry << endl;
                                cout << "== BEFORE =======================" << endl;
                                cout << pos << "\t" << rabbit.seconds<< "\t"
                                     << rabbit.seconds-irr_start.back() << "\t" 
                                     << rabbit.seconds-count_start.back() << endl;
                                cout << "\t" << cycle_time << "\t"
                                     << last_irr_length << "\t" 
                                     << last_count_length << endl;
                            }

                            //rewind
                            do{
                                jentry--;
                                rabbit.GetEntry(jentry);
                            }while(rabbit.seconds>(irr_start.back()+last_irr_length));
                            //insert a movement equal to the last cycle length
                            pos = -1;
                            new_move = (count_start.back() + last_count_length + last_irr_length);
                            count_start.push_back(new_move);
                            last_move = count_start.back();
                            num_missed_irr++;

                            if (opt_verbose){
                                cout << "== AFTER ========================" << endl;
                                cout << pos << "\t" << rabbit.seconds<< "\t"
                                     << rabbit.seconds-irr_start.back() << "\t" 
                                     << rabbit.seconds-count_start.back() << endl;
                                cout << "\t" << cycle_time << "\t"
                                     << last_irr_length << "\t" 
                                     << last_count_length << endl;
                                cout << "---------------------------------" << endl;
                            }
                        }
                        else if ((pos==-1) && 
                                ((rabbit.seconds-count_start.back())>(last_count_length*RabVar::max_var))){
                            //missed counting pulse
                            if (opt_verbose){
                                cout << "== MISSED COUNTING PULSE ========" << endl;
                                cout << jentry << endl;
                                cout << "== BEFORE =======================" << endl;
                                cout << pos << "\t" << rabbit.seconds<< "\t"
                                     << rabbit.seconds-irr_start.back() << "\t" 
                                     << rabbit.seconds-count_start.back() << endl;
                                cout << "\t" << cycle_time << "\t"
                                     << last_irr_length << "\t" 
                                     << last_count_length << endl;
                            }
                                 
                            //rewind
                            do{
                                jentry--;
                                rabbit.GetEntry(jentry);
                            }while(rabbit.seconds>(count_start.back()+last_count_length));
                            //insert a movement equal to the last cycle length
                            pos = 1;
                            new_move = (irr_start.back() + last_irr_length + last_count_length);
                            irr_start.push_back(new_move);
                            last_move = irr_start.back();
                            num_missed_count++;

                            if (opt_verbose){
                                cout << "== AFTER ========================" << endl;
                                cout << pos << "\t" << rabbit.seconds<< "\t"
                                     << rabbit.seconds-irr_start.back() << "\t" 
                                     << rabbit.seconds-count_start.back() << endl;
                                cout << "\t" << cycle_time << "\t"
                                     << last_irr_length << "\t" 
                                     << last_count_length << endl;
                                cout << "=================================" << endl;
                            }
                        }
                    }
                }

            }//end else

        }//end loop over events
        cout << endl;
    }//end !source run

    //loop over all events and write output tree
    cout << "Processing MDPP16 SCP..." << endl;
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            if (!opt_verbose){
                cout << '\r' << "Processing event " << jentry;
            }
        }

        cycle_time = 0;

        if (!source_run){
            if (jentry<start_event){
                //cycle time should be negative
                cycle_time = rabbit.seconds - start_offset;
            }
            else{
                for (auto & cycle_start : irr_start){
                    if (rabbit.seconds > cycle_start){
                        cycle_time = rabbit.seconds - cycle_start;
                    }
                }
            }
        }

        //energy calibrate detectors
        for (int det=0; det<RabVar::num_det; det++){ //loop over detectors
            En[det] = 0;
            if ((rabbit.ADC[RabVar::det_chn[det]]>10)&&(!(rabbit.overflow[RabVar::det_chn[det]]))){
                //Calibrate
                En[det] = (rabbit.ADC[RabVar::det_chn[det]]+r.Rndm()-0.5)*( (*rabbit.m) )[det] + ( (*rabbit.b) )[det];
                En[det] = int(100*(En[det]*cal.m[det] + cal.b[det]));
                En[det] = En[det]/100.;
                //fill histos
                hEn[det]->Fill(En[det]);
                hCycleTime->Fill(cycle_time);
            }
        }

        //fill tree
        tProcessed->Fill();
    }
    cout << endl;
    cout << nentries << " SCP events processed" << endl;

    //loop over all QDC events and write output tree
    cout << "Processing MDPP16 QDC..." << endl;
    for (Long64_t jentry=0; jentry<nentries_QDC; jentry++) {
        nb = rabbit_QDC.GetEntry(jentry);
        if (jentry%100000==0){
            if (!opt_verbose){
                cout << '\r' << "Processing event " << jentry;
            }
        }

        cycle_time = 0;

        if (!source_run){
            if (jentry<start_event){
                //cycle time should be negative
                cycle_time = rabbit_QDC.seconds - start_offset;
            }
            else{
                for (auto & cycle_start : irr_start){
                    if (rabbit_QDC.seconds > cycle_start){
                        cycle_time = rabbit_QDC.seconds - cycle_start;
                    }
                }
            }
        }

        //Calculated values (PSD)
        nmon_PSD = 0;
        if ((rabbit_QDC.ADC_long[RabVar::nmon_chn]>1)&&(!(rabbit_QDC.overflow[RabVar::nmon_chn]))){
            nmon_PSD = int((r.Rndm() + rabbit_QDC.ADC_long[RabVar::nmon_chn] - rabbit_QDC.ADC_short[RabVar::nmon_chn])
                       *100./(r.Rndm() + rabbit_QDC.ADC_long[RabVar::nmon_chn]));
            nmon_PSD = nmon_PSD/100.;
            hNmon_PSD->Fill(nmon_PSD);
        }

        //fill tree
        tProcessed_QDC->Fill();
    }
    cout << endl;
    cout << nentries_QDC << " QDC events processed" << endl;

    cout << "---------------------------------" << endl;
    cout << "Counting cycles:           " << num_count << endl;
    cout << "Irradiation cycles:        " << num_irr << endl;
    cout << "Counting cycles missed:    " << num_missed_count << endl;
    cout << "Irradiation cycles missed: " << num_missed_irr  << endl;
    cout << "---------------------------------" << endl;
    
    //write data to file
    fOut->cd();
    tProcessed->Write();
    tProcessed_QDC->Write();

    fOut->mkdir("histos");
    fOut->cd("histos");
    hIrrTime->Write();
    hCountTime->Write();
    hCycleTime->Write();
    hNmon_PSD->Write();
    for (int det=0; det<RabVar::num_det; det++){
        hEn[det]->Write();
    }

    fOut->Close();
    return 0;
    
}

// Reads in the calibration file and finds the calibration for the given run
calibration read_in_cal(int run_num){

    //variables to read in from file
    Float_t m_file[RabVar::num_det], b_file[RabVar::num_det];
    int run, run2;
    
    //calibration data
    ifstream infile;
    infile.open("datafiles/det_cal.dat");

    //return values
    calibration cal;
    for (int i = 0; i<RabVar::num_det; i++){
        cal.m[i] = 1.;
        cal.b[i] = 0;
    }

    //read in data from file
    do{
        infile>>run;
        if (run<0 || run>10000){
            break;
        }
        for (int j=0; j<RabVar::num_det; j++){
            infile>>m_file[j];
        }
        infile>>run2;
        if ((run != run2) && (run2 != 0)){
            cout << "ERROR IN DATA FILE! Entry:  " << run << endl;
        }
        for (int j=0; j<RabVar::num_det; j++){
            infile>>b_file[j];
        }
        if (run == run_num){
            break;
        }
    }while(!infile.eof());
    infile.close();
    
    //There exists calibration points for the given run
    if (run == run_num){
        //print calibration values
        cout << "Found run " << run << ", using calibration:" << endl;
        cout << "m:  ";
        for (int i = 0; i<RabVar::num_det; i++){
            cout << m_file[i] << "  ";
        }
        cout << endl << "b:  ";
        for (int i = 0; i<RabVar::num_det; i++){
            cout << b_file[i] << "  ";
        }
        cout << endl;
        //check with user
        //cout << "Is this OK? (y or n)  ";
        char ans = 'y';
        //cin >> ans;
        if  (ans == 'y'){
            for (int i = 0; i<RabVar::num_det; i++){
                cal.m[i] = m_file[i];
                cal.b[i] = b_file[i];
            }
        }
    }
    else{
        cout << "No additional calibration found in datafiles/det_cal.dat. Using only MVME calibration." 
             << endl;
    }
    return cal; //return calibration
}

