#ifndef SINGLE_RUN_RESULTS_H
#define SINGLE_RUN_RESULTS_H

#include "GlobalParameters.h"
#include "GraphicOutputManager.h"
#include "SignalOperations.h"
#include "SingleRunData.h"

class AnalysisManager;
class SingleRunData;
class AllRunsResults;

class SingleRunResults
{
public:
	enum Status { Ok, NotLoaded, NoPMTsignal, Empty, NoGEMsignal };
protected:
	Status _current_status;
	ParameterPile::experiment_area curr_area;
	DVECTOR xs_GEM;
	DVECTOR ys_GEM;
	STD_CONT<STD_CONT<peak>> mppc_peaks; //only for mppc channels
	STD_CONT<DVECTOR> mppc_baseline_xs;
	STD_CONT<DVECTOR> mppc_baseline_ys;
	DVECTOR mppc_S2_peaks_area;
	DVECTOR mppc_S2_start_t;
	DVECTOR mppc_S2_finish_t;
	//DVECTOR mppc_sum_peaks_area;
	DVECTOR mppc_double_I;
	STD_CONT<int> mppc_channels;

	double PMT3_summed_peaks_area;
	int PMT3_n_peaks;
	STD_CONT<peak> PMT3_peaks;
	STD_CONT<peak> PMT1_peaks;

	bool is_valid;
public:
	SingleRunResults(SingleRunData *of_what);
	SingleRunResults::Status getStatus(void) const;
	bool isValid(void) const;
	void setValid(bool val);

friend SingleRunData;
friend AnalysisManager;
friend AllRunsResults;
};

#endif
