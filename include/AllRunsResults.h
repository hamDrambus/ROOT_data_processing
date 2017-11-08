#ifndef ALL_RUNS_RESULTS_H
#define ALL_RUNS_RESULTS_H

#include "TH1D.h"
#include "TH1I.h"
#include "TLine.h"
#include "GlobalDefinitions.h"
#include "SingleRunResults.h"

class AnalysisManager;
class SingleRunData;
class SingleRunResults;

class AllRunsResults
{
protected:
	DVECTOR _ns;
	DVECTOR _Ss;
	DVECTOR _xs_GEM_sum;
	DVECTOR _ys_GEM_sum;
	int N_of_runs;
	int Iteration_N;
	ParameterPile::experiment_area _exp;
	GraphicOutputManager graph_manager;
	double S_peaks_cutoff;
	double N_peaks_cutoff;
	double S_peaks_max_cutoff;

	std::vector<DVECTOR> mppc_peaks_in_S2_area; //size == mppc channels (depends on experiment area)
	std::vector<DVECTOR> mppc_S2_start_time; //size == mppc channels (depends on experiment area)
	std::vector<DVECTOR> mppc_S2_finish_time; //size == mppc channels (depends on experiment area)
	std::vector<DVECTOR> mppc_all_peaks_Ss; //size == mppc channels (depends on experiment area)

	void find_GEM_start_time(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start, int N_trust, GraphicOutputManager &man);
	void find_S_cutoff(void); //in: _Ss, out: S_peaks_cutoff
	//void find_S_cutoff_v2(void);
public:
	AllRunsResults(ParameterPile::experiment_area experiment);//only experiment and channells are important here
	void processAllRuns(std::vector<SingleRunResults> &single_results);
	//For multithreading:
	void Merge(AllRunsResults* with);
	void Merged(void);
	void Clear(void);
	int Iteration(void) const;

	friend AnalysisManager;
	friend SingleRunData;
};


#endif
