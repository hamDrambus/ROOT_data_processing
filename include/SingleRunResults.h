#ifndef SINGLE_RUN_RESULTS_H
#define SINGLE_RUN_RESULTS_H

#include "GlobalDefinitions.h"
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
	std::vector<double> xs_GEM;
	std::vector<double> ys_GEM;
	std::vector<std::vector<peak>> mppc_peaks; //only for mppc channels
	std::vector<DVECTOR> mppc_baseline_xs;
	std::vector<DVECTOR> mppc_baseline_ys;
	DVECTOR mppc_S2_peaks_area;
	DVECTOR mppc_S2_start_t;
	DVECTOR mppc_S2_finish_t;
	DVECTOR mppc_sum_peaks_area;

	double PMT3_summed_peaks_area;
	int PMT3_n_peaks;

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
