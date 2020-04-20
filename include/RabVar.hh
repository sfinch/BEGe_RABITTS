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
    
    ///////////////////////////////////////////////////////////////
    // Cycle time variables
    // Values in seconds
    ///////////////////////////////////////////////////////////////
    const double irr_time = 10;
    const double count_time = 60;
    const double transit_time = 1.05;
    
    ///////////////////////////////////////////////////////////////
    // Analysis variables for producing time-binned spectra
    ///////////////////////////////////////////////////////////////

    // Default energy binning is 0.1 keV. 
    // Set to 2 for 0.2 keV, 5 for 0.5 keV, and 10 for 1 keV binning
    const int energy_rebin = 1;

    // analysis_cycle.C
    // Divides the counting time into num_win equal cycles and saves spectra
    const int num_win = 10;
    
    // analysis_overnight.C
    // Saves histograms every 3600 s (1 hour)
    const double overnight_win_time = 3600; //in sec

    ///////////////////////////////////////////////////////////////
    // Hardware configuration
    // chn is the channel number of the signal in the MDPP16
    ///////////////////////////////////////////////////////////////

    // Total number of HPGe detectors (in SCP)
    const int num_det = 6;
    const int det_chn[num_det] = {0, 2, 4, 5, 6, 7};

    // BEGe detectors (in SCP)
    // BEGes are always the first two HPGe detectors in the above array.
    const int num_BEGe = 2;
    const int min_BEGe_E[2] = {100, 100}; // in keV

    // Fission chambers (in SCP)
    const int num_FC = 2;
    const int FC_chn[num_FC] = {10, 11};
    const int FC_threshold[num_FC] = {6000, 6000}; //in ADC channels
    const int FC_rebin = 16;
    
    // Rabbit motor signal (in SCP) 
    const int rabbit_chn = 12;

    //BCI (in SCP)
    const int BCI_chn = 14;
    const int min_BCI = 10; //in ADC channel

    // Neutron monitor (in QDC) 
    const int nmon_chn = 2;
    const int min_nmon_E = 10;
    const double nmon_PSD_cut[2] = {-0.60, 0};

    ///////////////////////////////////////////////////////////////
    // Filter settings for processing rabbit motor signals
    // Used in process_rabbit
    ///////////////////////////////////////////////////////////////
    const double min_time = 0.2;  // filter for removing duplicate in sec
    // The cycle time can vary by at most 10%
    const double max_var = 1.1;   // filter for missed signals in %
    const double min_var = 0.9;   // filter for extra signals in %

    ///////////////////////////////////////////////////////////////
    // Calculated values
    // You shouldn't have to change these values.
    // If the motor signal cables are switched, you can substract 
    // transit_time from all time_irr and time_count to compensate
    ///////////////////////////////////////////////////////////////
    double time_bin = count_time/num_win;
    double time_irr[2] = {transit_time, irr_time+transit_time};
    double time_count[2] = {irr_time + 2*transit_time, 
                            irr_time + 2*transit_time + count_time};

}

#endif
