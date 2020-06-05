
///////////////////////////////////////////////////////////////////////////
// HistVar.h
// If you make changes here, you must run 'make' before running process_rabbit. 
// All other scripts are not compiled, and do not require re-making. 
//
//////////////////////////////////////////////////////////////////////////


#ifndef HistVar_h
#define HistVar_h 1

namespace HistVar{
    
    //Energy spectra ranges (for plots)
    const int start_E = 0; //in keV
    const int end_E = 4000; //in keV
    const double binsperkev = 0.10; 
    const int num_E_bins = (end_E-start_E)/(binsperkev);
    
    //Cycle time ranges (for plots)
    const double min_cycle_time = -1;
    const double max_cycle_time = 89;
    const double cycle_bin_per_sec = 100;
    const int cycle_time_bins = cycle_bin_per_sec*(max_cycle_time-min_cycle_time);

}

#endif
