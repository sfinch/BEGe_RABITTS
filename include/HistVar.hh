
///////////////////////////////////////////////////////////////////////////
// HistVar.h
// If you make changes here, you must run 'make' before running process_rabbit. 
// All other scripts are not compiled, and do not require re-making. 
//
//////////////////////////////////////////////////////////////////////////


#ifndef HistVar_h
#define HistVar_h 1

namespace HistVar{
    
    //Cycle time ranges (for plots)
    double min_cycle_time = -1;
    double max_cycle_time = 89;
    double cycle_bin_per_sec = 100;
    int cycle_time_bins = cycle_bin_per_sec*(max_cycle_time-min_cycle_time);

}

#endif
