///////////////////////////////////////////////////////////////////////////
// RabVar.h
// All variables for the runs are set here. This includes all the number of detectors,
// detector channel numbers, timing cuts, PSD cuts, and FC thresholds.
// If you make changes here, you must run 'make' before running process_rabbit. 
// All other scripts are not compiled, and do not require re-making. 
//
//////////////////////////////////////////////////////////////////////////

#ifndef RabVar_h
#define RabVar_h 1

namespace RabVar{
    
    //cycle time variables
    double irr_time = 10;
    double count_time = 60;
    double transit_time = 1.05;
    
    //analysis_cycle
    const int num_win = 10;
    const int energy_rebin = 1;
    
    //analysis_overnight
    const int num_overnight_win = 16;
    const double overnight_win_time =3600; //in sec

    //Hardware configuration
    const int num_det = 6;
    int det_chn[num_det] = {0, 2, 4, 5, 6, 7};
    const int rabbit_chn = 12;

    //FC Variables
    const int num_FC = 2;
    int FC_chn[num_FC] = {10, 11};
    int FC_threshold[num_FC] = {6000, 6000};
    const int FC_rebin = 16;

    //BEge
    const int num_BEGe = 2;
    int min_BEGe_E[2] = {100, 100};

    //neutron monitor
    int nmon_chn = 2;
    int min_nmon_E = 10;
    double nmon_PSD_cut[2] = {-0.60, 0};

    //BCI
    int BCI_chn = 14;
    int min_BCI = 10;

    //Process rabbit cycle time variables
    double min_time = 0.2;  //filter for duplicate in sec
    double max_var = 1.1;   //filter for missed signals
    double min_var = 0.9;   //filter for extra signals

    //calculated values
    double time_bin = count_time/num_win;
    double time_irr[2] = {transit_time, irr_time+transit_time};
    double time_count[2] = {irr_time + 2*transit_time, 
                            irr_time + 2*transit_time + count_time};
}

#endif
