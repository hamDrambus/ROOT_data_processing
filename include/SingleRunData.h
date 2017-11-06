#ifndef SINGLE_RUN_DATA_H
#define SINGLE_RUN_DATA_H

#include "GlobalDefinitions.h"
#include "GraphicOutputManager.h"
#include "SignalOperations.h"
#include "SingleRunResults.h"
#include "AllRunsResults.h"

class AnalysisManager;
class AllRunsResults;
class SingleRunResults;

class SingleRunData
{
//public:
	//enum Status {Ok, NotLoaded, NoPMTsignal, Empty, NoGEMsignal};
protected:
	//Status _current_status;

	std::vector<std::vector<double>>xs_channels;
	std::vector<std::vector<double>>ys_channels;
	
	std::vector<double> found_base_lines;//TODO: modify to shaped baselines
	std::vector<std::vector<double>> xs_integral_ch;
	std::vector<std::vector<double>> ys_integral_ch;
	//std::vector<double> xs_GEM;
	//std::vector<double> ys_GEM;

	double PMT3_summed_peaks_area;
	int PMT3_n_peaks;
	bool PMT_peaks_analysed;

	ParameterPile::experiment_area curr_area; //sets channel area

	//int get_channel_index(int ch);
	void extract_one_subrun(std::vector<std::vector<double>> &xxs, std::vector<std::vector<double>> &yys, int sub_index);
	void add_draw_data(std::string prefix, GraphicOutputManager& graph_manager,
		ParameterPile::DrawEngine de = ParameterPile::DrawEngine::Gnuplot); //only draws specified runs
	void SingleRunData::add_draw_baselines(std::string prefix, GraphicOutputManager& graph_manager,
		ParameterPile::DrawEngine de = ParameterPile::DrawEngine::Gnuplot);

	void file_to_vector(std::string fname, std::vector<double> &xs, std::vector<double> &ys, int index);
	bool test_PMT_signal(int _N_threshold, double _S_threshold, double _S_max_threshold, SingleRunResults &results); //returns false if the PMT signal is empty
	void find_time_limits(void);//TODO: this is much more complicated operation, the most difficult part of this analysis actually

	//bool is_valid;
	//SingleRunResults * _results;
	GraphicOutputManager graph_manager;
	void readOneRun(SingleRunResults &results);
public:
	SingleRunData(ParameterPile::experiment_area area);
	SingleRunResults processSingleRun_Iter_0(const AllRunsResults *all_runs_results);
	SingleRunResults processSingleRun_Iter_1(const AllRunsResults *all_runs_results);
	SingleRunResults processSingleRun(const AllRunsResults *all_runs_results);
	//^ must be called only if the previous SingleRunResults was valid
	void runProcessedProc(void);
	void clear_memory(void); //clears only 'input' data, preserves processing results
	ParameterPile::experiment_area getArea(void) const;
};

#endif