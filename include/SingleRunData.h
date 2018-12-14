#ifndef SINGLE_RUN_DATA_H
#define SINGLE_RUN_DATA_H

#include "GlobalParameters.h"
#include "GraphicOutputManager.h"
#include "SignalOperations.h"
#include "AllRunsResults.h"

class AnalysisManager;
class AllRunsResults;
class TVirtualFFT;

class SingleRunData
{
public:
	enum Status { Ok, NotLoaded, NoPMTsignal, Empty, NoGEMsignal, PMT_mismatch, MPPC_mismatch, AVG_mismatch };
protected:
	bool _valid;
	Status _status;

	STD_CONT<DVECTOR>xs_channels;
	STD_CONT<DVECTOR>ys_channels;
	
	DVECTOR found_base_lines;

	double PMT3_summed_peaks_area;
	int PMT3_n_peaks;

	ParameterPile::experiment_area curr_area; //sets channel area
	GraphicOutputManager graph_manager;

	//int get_channel_index(int ch);
	void extract_one_subrun(STD_CONT<DVECTOR> &xxs, STD_CONT<DVECTOR> &yys, int sub_index);
	void add_draw_data(std::string prefix, GraphicOutputManager& graph_manager,
		ParameterPile::DrawEngine de = ParameterPile::DrawEngine::Gnuplot, int ch = -1); //only draws specified runs
	void add_draw_baselines(std::string prefix, GraphicOutputManager& graph_manager,
		ParameterPile::DrawEngine de = ParameterPile::DrawEngine::Gnuplot);

	void file_to_vector(std::string fname, DVECTOR &xs, DVECTOR &ys, int index);
	bool test_PMT_signal(int _N_threshold, double _S_threshold, double _S_max_threshold, AllRunsResults *results); //returns false if the PMT signal is empty
	double find_spreaded_peaks_threshold(DVECTOR &x_peaks_spreaded, DVECTOR &y_peaks_spreaded, double &apr_threshold);

	void readOneRun(AllRunsResults *results);//TODO: rework this function using readOneRunMPPCs (or remove altogether)
#ifdef _HOTFIX_DECREASE_MPPC_MEMORY_USAGE
	void readOneRunMPPCs(AllRunsResults *results, int channel);//TODO: rename function and get rid of _HOTFIX_DECREASE_MEMORY at all
	void clearOneRunMPPCs(int channel);
#endif
	void calculate_MPPC_threshold_and_baseline(DVECTOR &xs, DVECTOR &ys, double &threshold, double &baseline,STD_CONT<peak> &peaks_before_S1);
	void calculate_PMT_threshold_and_baseline(DVECTOR &xs, DVECTOR &ys, double &threshold, double &baseline, STD_CONT<peak> &peaks_before_S1, int channel);
	void push_average (int ch, bool is_first_call, AllRunsResults *all_runs_results);
	void push_dispersion (int ch, bool is_first_call, AllRunsResults *all_runs_results);
public:
	SingleRunData(ParameterPile::experiment_area area);
	void processSingleRun_Iter_0(AllRunsResults *all_runs_results);
	void processSingleRun_Iter_1(AllRunsResults *all_runs_results);
	void processSingleRun_Iter_2(AllRunsResults *all_runs_results);
	void processSingleRun(AllRunsResults *all_runs_results);
	//^ must be called only if the previous SingleRunResults was valid
	void runProcessedProc(void);
	void clear_memory(void); //clears only 'input' data, preserves processing results
	ParameterPile::experiment_area getArea(void) const;

	std::size_t real_size(void);
};

#endif
