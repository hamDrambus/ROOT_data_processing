#ifndef SINGLE_RUN_DATA_H
#define SINGLE_RUN_DATA_H

#include "GlobalParameters.h"
#include "GraphCollection.h"
#include "SignalOperations.h"
#include "AllRunsResults.h"

class AnalysisManager;
class AllRunsResults;
class TVirtualFFT;

class SingleRunData
{
protected:
	double trigger_offset;

	STD_CONT<DVECTOR>xs_channels;
	STD_CONT<DVECTOR>ys_channels;
	
	DVECTOR found_base_lines;

	double PMT3_summed_peaks_area;
	int PMT3_n_peaks;

	ParameterPile::experiment_area curr_area; //sets channel area

	void file_to_vector(std::string fname, DVECTOR &xs, DVECTOR &ys, int index);
	double find_spreaded_peaks_threshold(DVECTOR &x_peaks_spreaded, DVECTOR &y_peaks_spreaded, double &apr_threshold);

	void readOneRun(AllRunsResults *results, int channel);
	void clearOneRun(int channel);
	void calculate_MPPC_threshold_and_baseline(DVECTOR &xs, DVECTOR &ys, double &threshold, double &baseline,STD_CONT<peak> &peaks_before_S1);
	void calculate_PMT_threshold_and_baseline(DVECTOR &xs, DVECTOR &ys, double &threshold, double &threshold_2, double &baseline, STD_CONT<peak> &peaks_before_S1, int channel);
	void push_event (AllRunsResults *all_runs_results);
	void push_average (int ch, bool is_first_call, AllRunsResults *all_runs_results);
	void push_dispersion (int ch, AllRunsResults *all_runs_results);
public:
	SingleRunData(ParameterPile::experiment_area area);
	void processSingleRun_Iter_0(AllRunsResults *all_runs_results);
	void processSingleRun_Iter_1(AllRunsResults *all_runs_results);
	void processSingleRun_Iter_2(AllRunsResults *all_runs_results);
	void processSingleRun(AllRunsResults *all_runs_results);
	//^ must be called only if the previous SingleRunResults was valid
	void runProcessedProc(AllRunsResults *all_runs_results);
	void clear_memory(void); //clears only 'input' data, preserves processing results
	ParameterPile::experiment_area getArea(void) const;

	std::size_t real_size(void);
	friend AnalysisManager;
	friend AllRunsResults;
};

#endif
