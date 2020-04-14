//////////////////////////////////////////////////////////////////////////
//  bkgd_cal.C                                                          //
//  Finds the chn number corresponding to the centriod of n peaks,      //
//  then finds a linear calibration from channel to energy              //
//  which is saved in a text file                                       //
//      datafiles/det_cal.dat                                           //
//  Does this for all detectors.                                        //
//  To run: root -l "bkgd_cal.C(XXX)" where  XXX is the run number      //
//  Requires that data has been converted to ROOT using mvme2root       //
//////////////////////////////////////////////////////////////////////////


#include <vector>
#include "RabVar.h"

using namespace RabVar; //bad habit, need to remove

void bkgd_cal(int run_num = 0){

    //get run number
    if (run_num < 1){
        cout << "Run Number?   ";
        cin >> run_num;
    }
    cout << "Run Number " << run_num << endl;

    std::vector<double> energy;
    std::vector<double> energy_error;
    std::vector<double> startrange[num_det];
    std::vector<double> endrange[num_det];

    // Set initial ranges for each fit. Will fit the largest ADC peak in the region and
    // calibrate it to the corresponding energy
    //185
    energy.push_back(185.72);
    energy_error.push_back(0.02);
    for (int i=0; i<num_det; i++){
        startrange[i].push_back(175);
        endrange[i].push_back(195);
    }

    //351
    energy.push_back(351.92);
    energy_error.push_back(0.02);
    for (int i=0; i<num_det; i++){
        startrange[i].push_back(347);
        endrange[i].push_back(355);
    }

    //511
    /*
    energy.push_back(511.0);
    energy_error.push_back(0.02);
    for (int i=0; i<num_det; i++){
        startrange[i].push_back(500);
        endrange[i].push_back(520);
    }
    */

    //609
    energy.push_back(609.320);
    energy_error.push_back(0.005);
    for (int i=0; i<num_det; i++){
        startrange[i].push_back(605);
        endrange[i].push_back(619);
    }

    //1001
    energy.push_back(1001.03);
    energy_error.push_back(0.10);
    for (int i=0; i<num_det; i++){
        startrange[i].push_back(995);
        endrange[i].push_back(1006);
    }

    //1120
    energy.push_back(1120.4);
    energy_error.push_back(0.10);
    for (int i=0; i<num_det; i++){
        startrange[i].push_back(1105);
        endrange[i].push_back(1130);
    }

    //1461
    energy.push_back(1460.822);
    energy_error.push_back(0.06);
    for (int i=0; i<num_det; i++){
        startrange[i].push_back(1440);
        endrange[i].push_back(1490);
    }

    //1764
    energy.push_back(1764.491);
    energy_error.push_back(0.01);
    for (int i=0; i<num_det; i++){
        startrange[i].push_back(1750);
        endrange[i].push_back(1775);
    }

    //2204
    energy.push_back(2204.059);
    energy_error.push_back(0.022);
    for (int i=0; i<num_det; i++){
        startrange[i].push_back(2190);
        endrange[i].push_back(2220);
    }

    //Pb
    energy.push_back(2614.511);
    energy_error.push_back(0.01);
    for (int i=0; i<num_det; i++){
        startrange[i].push_back(2590);
        endrange[i].push_back(2630);
    }


    //variables
    const int num_points = energy.size();
    double centroid[num_det][num_points];
    double centroid_error[num_det][num_points];
    double FWHM[num_det][num_points];
    double FWHM_error[num_det][num_points];
    double m[num_det];
    double b[num_det];
    double residual[num_det][num_points];
    double residual_error[num_det][num_points];
    for (int i=0; i<num_det; i++){
        m[i]=0;
        b[i]=0;
    }

    //root histos
    TCanvas *cPeak[num_points];
    TH1F *hADC[num_det][num_points];
    TF1 *peak[num_det][num_points];

    // file with root data
    TFile *fRun; 
    if (run_num<10){
        fRun = new TFile(Form("data_root/RABITTS_00%i.root",run_num));
    }
    else if (run_num<100){
        fRun = new TFile(Form("data_root/RABITTS_0%i.root",run_num));
    }
    else{
        fRun = new TFile(Form("data_root/RABITTS_%i.root",run_num));
    }

    //make pretty plots
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    
    //loop
    for (int j=0; j<num_points; j++){ // loop over number of points
        int kev = energy[j];
        cPeak[j] = new TCanvas(Form("cPeak%i",j),Form("Calibration Peak %d", kev),1200,600);
        cPeak[j]->Divide(4,2);
        for (int i=0; i<num_det; i++){      // loop over all detectors
            cPeak[j]->cd(i+1);

            if (j==0){     //project ADC histos from file on first loop
                hADC[i][j] = (TH1F*) fRun->Get(Form("histos_SCP/hEn%i", det_chn[i]));
            }
            else{           //copy histos on following loops
                hADC[i][j] = (TH1F*)hADC[i][0]->Clone(Form("hADC%i%i",i,j));
                hADC[i][j]->SetTitle(Form("Det %i peak %i", i, kev));
            }
            hADC[i][j]->Rebin(2);
            hADC[i][j]->GetXaxis()->SetRangeUser(startrange[i][j],endrange[i][j]);
            if (hADC[i][j]->GetMaximum() > 2){ // check to make sure histogram isn't empty

                // fit for peak = linear background + gaussian
                TF1 *fit = new TF1("fit","[0] + [1]*x + [2]*exp(-(x-[3])^2/(2*[4]^2))",startrange[i][j],endrange[i][j]);

                // find estimates for fit
                int bkgd1 = hADC[i][j]->GetBinContent(hADC[i][j]->FindBin(startrange[i][j]));
                int bkgd2 = hADC[i][j]->GetBinContent(hADC[i][j]->FindBin(endrange[i][j]));
                double slope = (bkgd2-bkgd1)/(endrange[i][j]-startrange[i][j]);
                double y_int = bkgd1 - slope*startrange[i][j];
                double gauss_height =  hADC[i][j]->GetMaximum();
                int maxbin = hADC[i][j]->GetMaximumBin();
                double findpeak = hADC[i][j]->GetXaxis()->GetBinCenter(maxbin);

                
                // set estimates for fit
                fit->SetParameter(0,y_int);
                fit->SetParameter(1,slope);
                fit->SetParameter(2,gauss_height);
                fit->SetParameter(3,findpeak);
                fit->SetParameter(4,1);
                
                fit->SetParLimits(2,1,10*gauss_height);
                fit->SetParLimits(3,startrange[i][j],endrange[i][j]);
                fit->SetParLimits(4,0.1,10);

                // fit
                hADC[i][j]->Fit("fit");

                //get centroid of peak
                peak[i][j] = hADC[i][j]->GetFunction("fit");
                centroid[i][j] = peak[i][j]->GetParameter(3);
                centroid_error[i][j] = peak[i][j]->GetParError(3);
                FWHM[i][j] = 2.355*(peak[i][j]->GetParameter(4));
                FWHM_error[i][j] = 2.355*(peak[i][j]->GetParError(4));

            }
        }
    }
    
    // graphs of calibrations
    TCanvas *cCal = new TCanvas("cCal","Calibration Peaks and Fit",1200,600);
    cCal->Divide(4,2);
    TGraphErrors *tCal[num_det];
    
    // graphs for calibration residuals
    TCanvas *cRes = new TCanvas("cRes","Calibration residuals",1200,600);
    cRes->Divide(4,2);
    TGraphErrors *tRes[num_det];
    TF1 *fCal[num_det];

    // make graphs
    for (int i=0; i<num_det; i++){
        // calibration
        cCal->cd(i+1);
        tCal[i] = new TGraphErrors(num_points, centroid[i], energy.data(), 
                                   centroid_error[i], energy_error.data());
        tCal[i]->SetTitle(Form("ADC linear calibration Det %i", i));
        tCal[i]->SetMarkerStyle(20);
        tCal[i]->SetMarkerColor(2);
        tCal[i]->Draw("AP");

        tCal[i]->Fit("pol1");   //fit to a polynomial of order 1 (y = m*x + b)
        b[i] = tCal[i]->GetFunction("pol1")->GetParameter(0); //get fit values
        m[i] = tCal[i]->GetFunction("pol1")->GetParameter(1);

        fCal[i] = new TF1(Form("fCal%i", i), "[0]+[1]*x", 0, 30000);
        //fCal[i] = new TF1(Form("fCal%i", i), "[0]+[1]*x+[2]*x^2", 0, 30000);
        //fCal[i] = new TF1(Form("fCal%i", i), "[0]+[1]*x+[2]*x^2+[3]*x^3", 0, 30000);
        fCal[i]->SetParameter(0, 0);
        fCal[i]->SetParameter(1, 1);
        fCal[i]->SetParameter(2, 0);
        fCal[i]->SetParameter(3, 0);
        //tCal[i]->Fit(Form("fCal%i", i));

        // residuals
        cRes->cd(i+1);
        for (int j=0; j<num_points; j++){
            //residual[i][j] = energy[j]-b[i]-m[i]*centroid[i][j];
            residual[i][j] = energy[j]-tCal[i]->GetFunction("pol1")->Eval(centroid[i][j]);
            //residual[i][j] = energy[j]-fCal[i]->Eval(centroid[i][j]);
            residual_error[i][j] = centroid_error[i][j];
        }
        tRes[i] = new TGraphErrors(num_points, energy.data(), residual[i], 
                                   energy_error.data(), residual_error[i]);
        tRes[i]->SetTitle(Form("Calibration residual Det %i", i));
        tRes[i]->SetMarkerStyle(20);
        tRes[i]->SetMarkerColor(4);
        tRes[i]->Draw("AP");

    }

    for (int j=0; j<num_det; j++){
        cout << "------------------------------------------------------------" << endl;
        cout << "                    Detector " << j+1 << endl;
        cout << "------------------------------------------------------------" << endl;
        cout << "En [kev] \tCentroid \tFWHM\t\tFWHM uncert" << endl;
        for (int i=0; i<num_points; i++){
            cout << energy[i] << "\t\t" << centroid[j][i] << "\t\t" 
                 << FWHM[j][i] << "\t\t" << FWHM_error[j][i] << endl;
        }
    }
    
    gSystem->ProcessEvents();   //update plots for user

    // check calibration with user
    char ans = 'n';
    cout << endl;
    cout << "Accept this calibration (y/n)?  ";
    cin >> ans;
    if (ans == 'y'){    //write calibration to file
        FILE *file_ptr;
        file_ptr = fopen("datafiles/det_cal.dat","a");
        fprintf(file_ptr, "%i", run_num);
        for (int i=0; i<num_det; i++){
            fprintf(file_ptr,"\t%f", m[i]); 
        }
        fprintf(file_ptr, "\n%i", run_num);
        for (int i=0; i<num_det; i++){
            fprintf(file_ptr,"\t%f", b[i]); 
        }
        fprintf(file_ptr, "\n");
        fclose(file_ptr);
        cout << "Run " << run_num << " Output complete!" << endl;
    }
}   

