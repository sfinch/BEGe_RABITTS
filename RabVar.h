#ifndef RabVar_h
#define RabVar_h 1

namespace RabVar{

    
    /*
    //analysis_cycle variables
    const int num_det = 2;

    const int num_win = 10;
    double irr_time = 10;
    double count_time = 20;
    double transit_time = 1.05;
    
    const int rebin = 1;
    
    //cuts
    double time_bin = count_time/num_win;
    double time_irr[2] = {transit_time, irr_time+transit_time};
    double time_count[2] = {irr_time + 2*transit_time, 
                            irr_time + 2*transit_time + count_time};
    double time_win[num_win][2];
    for (int i=0; i<num_win; i++){
        time_win[i][0] = (time_bin*i)+time_count[0];
        time_win[i][1] = (time_bin*(i+1))+time_count[0];
    }
    
    //analysis_overnight
    const int num_det = 2;

    const int num_win = 16;
    const double win_time =3600; //in sec
    
    //run_info
    const int num_det = 2;
    int FC_chn[num_det] = {10, 11};
    int nE_chn= 8;
    int nPSD_chn = 9;
    int BCI_chn = 14;

    double time_irr[2] = {0.0, 10.0};

    const int num_lines = 2;
    double lx[num_lines] = {6000, 66000}; //the range for the integration
    */
    //
    

    //Process rabbit Variables
    const int num_det = 6;
    int det_chn[num_det] = {0, 2, 4, 5, 6, 7};
    const int rabbit_chn = 12;
    bool opt_verbose = 0;

    //cycle time variables
    double min_time = 0.2;  //filter for duplicate in sec
    double max_var = 1.1;   //filter for missed signals
    double min_var = 0.9;   //filter for extra signals

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

    //Cycle time ranges (for plots)
    double min_cycle_time = -1;
    double max_cycle_time = 89;
    double cycle_bin_per_sec = 100;
    int cycle_time_bins = cycle_bin_per_sec*(max_cycle_time-min_cycle_time);

}

#endif
