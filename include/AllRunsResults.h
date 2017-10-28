#ifndef ALL_RUNS_RESULTS_H
#define ALL_RUNS_RESULTS_H

#include "TH1D.h"
#include "TH1I.h"
#include "GlobalDefinitions.h"
#include "SingleRunResults.h"

class AnalysisManager;
class AllRunsResults
{
protected:
	DVECTOR _ns;
	DVECTOR _Ss;
	DVECTOR _xs_GEM_sum;
	DVECTOR _ys_GEM_sum;
	ParameterPile::experiment_area _exp;
	GraphicOutputManager graph_manager;
	double S_peaks_cutoff;
	void find_GEM_start_time(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start, int N_trust, GraphicOutputManager &man);
public:
	AllRunsResults(ParameterPile::experiment_area experiment);//only experiment and channells are important here
	void processAllRuns(std::vector<SingleRunResults> &single_results);

	friend AnalysisManager;
};


#endif
