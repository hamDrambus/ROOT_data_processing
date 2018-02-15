#ifndef ALL_RUNS_RESULTS_H
#define ALL_RUNS_RESULTS_H

#include "TH1D.h"
#include "TH1I.h"
#include "TLine.h"
#include "TStyle.h"
#include "GlobalParameters.h"
#include "SingleRunResults.h"

class AnalysisManager;
class SingleRunData;
class SingleRunResults;

#ifdef _USE_TIME_STATISTICS
struct time_results {
	long long t_PMT_proc;//0st iteration
	int n_PMT_proc;
	long long t_PMT_baseline;//0st iteration
	int n_PMT_baseline;
	long long t_PMT_peaks;//0st iteration
	int n_PMT_peaks;
	long long t_PMT_file_reading;//0st
	int n_PMT_file_reading;
	long long t_PMT_filtering;//0st
	int n_PMT_filtering;
	//====================================
	long long t_MPPC_proc; //all till ==...= is the 1st iteration
	int n_MPPC_proc;
	long long t_MPPC_file_reading;
	int n_MPPC_file_reading;
	long long t_MPPC_filtering;
	int n_MPPC_filtering;
	long long t_MPPC_threshold_and_first_baseline;
	int n_MPPC_threshold_and_first_baseline;
	//long long t_MPPC_threshold_and_first_baseline_peaks;
	//int n_MPPC_threshold_and_first_baseline_peaks;
	long long t_MPPC_curved_baseline;
	int n_MPPC_curved_baseline;
	long long t_MPPC_curved_baseline_v2;
	int n_MPPC_curved_baseline_v2;
	long long t_MPPC_curved_baseline_v3;
	int n_MPPC_curved_baseline_v3;
	long long t_MPPC_curved_baseline_v4;
	int n_MPPC_curved_baseline_v4;
	long long t_MPPC_curved_baseline_v5;
	int n_MPPC_curved_baseline_v5;
	long long t_MPPC_curved_baseline_v6;
	int n_MPPC_curved_baseline_v6;
	long long t_MPPC_curved_baseline_v7;
	int n_MPPC_curved_baseline_v7;
	long long t_MPPC_curved_baseline_v8;
	int n_MPPC_curved_baseline_v8;

	long long t_MPPC_curved_baseline_baseline;
	int n_MPPC_curved_baseline_baseline;
	long long t_MPPC_baseline_substraction;
	int n_MPPC_baseline_substraction;
	long long t_MPPC_peaks_finding;
	int n_MPPC_peaks_finding;
	long long t_MPPC_peaks_processing;
	int n_MPPC_peaks_processing;
	long long t_MPPC_double_I;
	int n_MPPC_double_I;
	//====================================
	long long t_RUN_proc;  //both iterations, must be set at merged() proc
	int n_RUN_proc;
	long long t_RUN_proc_single_iter;//both iterations, used at merge() proc, set in analysis Manager, must be cleared
	int n_RUN_proc_single_iter;
};
#endif
class AllRunsResults
{
protected:
	DVECTOR _ns;
	DVECTOR _Ss;
	DVECTOR _xs_GEM_sum;
	DVECTOR _ys_GEM_sum;
	int N_of_runs;
	int N_of_valid_runs;//i.e. accepted by PMT cut
	int Iteration_N;
	ParameterPile::experiment_area _exp;
	GraphicOutputManager graph_manager;
	double S_peaks_cutoff;
	double N_peaks_cutoff;
	double S_peaks_max_cutoff;

	STD_CONT<DVECTOR> mppc_peaks_in_S2_area; //size == mppc channels (depends on experiment area)
	STD_CONT<DVECTOR> mppc_S2_start_time; //size == mppc channels (depends on experiment area)
	STD_CONT<DVECTOR> mppc_S2_finish_time; //size == mppc channels (depends on experiment area)
	//STD_CONT<DVECTOR> mppc_all_peaks_Ss; //size == mppc channels (depends on experiment area) //depr: now using mppc_peaks
	STD_CONT<DVECTOR> mppc_double_Is; //size == mppc channels (depends on experiment area)
	STD_CONT<int> mppc_channels;
	STD_CONT<STD_CONT<STD_CONT<peak>>> mppc_peaks; //[channel][run#][peaks]. The number of runs must be equal to the size of DVECTOR above.
	STD_CONT<STD_CONT<peak>> PMT3_peaks; //[run#][peaks]
	STD_CONT<STD_CONT<peak>> PMT1_peaks; //[run#][peaks]

	void find_GEM_start_time(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start, int N_trust, GraphicOutputManager &man);
	void find_S_cutoff(void); //in: _Ss, out: S_peaks_cutoff
	//void find_S_cutoff_v2(void);
	TH1D* createMPPCHist(DVECTOR &what, std::string name, double left_cutoff, double right_cutoff_from_RMS, int N_bins = 0);
	TH1D* createMPPCHist_peaks_S(STD_CONT<STD_CONT<peak>> &what, std::string name, double left_cutoff, double right_cutoff_from_RMS, int N_bins = 0);
	void vector_to_file(DVECTOR &what, std::string fname, std::string title="MPPC");
	void vector_to_file(STD_CONT<STD_CONT<peak>> &pks, std::string fname, std::string title = "MPPC_peaks");
	TF1* createMPPCFitFunc(TH1D* hist, std::string name);

#ifdef _USE_TIME_STATISTICS
	time_results time_stat;
	void report_time_statistics();
#endif

public:
	AllRunsResults(ParameterPile::experiment_area experiment);//only experiment and channells are important here
	void processAllRuns(STD_CONT<SingleRunResults> &single_results);
	//For multithreading:
	void Merge(AllRunsResults* with);
	void Merged(void);
	void Clear(void);
	int Iteration(void) const;

	friend AnalysisManager;
	friend SingleRunData;
};


#endif
