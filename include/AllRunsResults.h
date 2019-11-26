#ifndef ALL_RUNS_RESULTS_H
#define ALL_RUNS_RESULTS_H

#include "TH1D.h"
#include "TH1I.h"
#include "TLine.h"
#include "TStyle.h"
#include "GlobalParameters.h"
#include "SingleRunData.h"

class AnalysisManager;
class SingleRunData;

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
public:
	enum Status { Ok, Empty, NoPMTsignal, PMT_mismatch, MPPC_mismatch, AVG_mismatch};
protected:
	int N_of_runs;
	int N_of_valid_runs;//i.e. accepted by PMT cut
	int Iteration_N;
	ParameterPile::experiment_area _exp;
	GraphCollection graph_manager;
	//[#run]
	std::vector<bool> _valid;
	std::vector<Status> _status;
	//AVERAGES:
	//TODO: ensure that only events which is valid in final are used. (N_of_valid_runs can be different for each iteration)
	STD_CONT<DVECTOR> _xs_sum;  //[channel]
	STD_CONT<DVECTOR> _ys_sum;  //[channel]
	STD_CONT<DVECTOR> _ys_disp; //[channel]
	STD_CONT<int> avr_channels;
	//DONE: The number of runs in MPPC and PMT blocks is always equal
	//MPPC:
	//STD_CONT<STD_CONT<double> > mppc_peaks_in_S2_area; //[run#][channel], size of mppc channels (depends on experiment area)
	//STD_CONT<STD_CONT<double> > mppc_S2_start_time;	 //[run#][channel]
	//STD_CONT<STD_CONT<double> > mppc_S2_finish_time;	 //[run#][channel]
	STD_CONT<STD_CONT<double> > mppc_double_Is;		 //[run#][channel]
	STD_CONT<STD_CONT<STD_CONT<peak> > > mppc_peaks; //[run#][channel][peaks].
	STD_CONT<int> mppc_channels;
	//PMT:
	DVECTOR _ns; //[run#]
	DVECTOR _Ss; //[run#]
	STD_CONT<STD_CONT<STD_CONT<peak> > > pmt_peaks;	//[run#][channel][peaks]
	STD_CONT<int> pmt_channels; //[channel]

	STD_CONT<STD_CONT<double> > pmt_S2_integral; //[run#][channel]
	STD_CONT<int> pmt_integrated_channels; //[channel]

	void find_GEM_start_time(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start, int N_trust, GraphCollection &man);
	TH1D* createMPPCHist(STD_CONT<STD_CONT<double> > &what, int ch_ind, std::string name, double left_cutoff, double right_cutoff_from_RMS, int N_bins = 0);
	TH1D* createMPPCHist_peaks_S(STD_CONT<STD_CONT<STD_CONT<peak> > > &what, int ch_ind, std::string name, double left_cutoff, double right_cutoff_from_RMS, int N_bins = 0);
	void vector_to_file(STD_CONT<STD_CONT<double> > &what, int ch_ind, std::string fname, std::string title="MPPC");
	void vector_to_file(STD_CONT<STD_CONT<STD_CONT<peak> > > &pks, int ch_ind, std::string fname, std::string title = "MPPC_peaks");
	void vector_to_file(DVECTOR &xs, DVECTOR &ys, DVECTOR &ydisps, std::string fname, std::string title = "Average");
	TF1* createMPPCFitFunc(TH1D* hist, std::string name);

	double Mean(STD_CONT<STD_CONT<STD_CONT<peak> > > &peaks, int ch_ind, std::function<double (peak& pk)> &value_picker);
	double RMS(STD_CONT<STD_CONT<STD_CONT<peak> > > &peaks, int ch_ind, std::function<double (peak& pk)> &value_picker);
	double Mean(STD_CONT<STD_CONT<double> > &vals, int ch_ind);
	double RMS(STD_CONT<STD_CONT<double> > &vals, int ch_ind);

#ifdef _USE_TIME_STATISTICS
	time_results time_stat;
	void report_time_statistics();
#endif
public:
	AllRunsResults(ParameterPile::experiment_area experiment);//only experiment and channels are important here
	//For multithreading:
	void Merge(AllRunsResults* with);
	void Merged(void);
	void Clear(void);//called at the end of Merged() and in Merge(with) for 'with'
	void ClearMerged(void); //called for Merged results in multithread case
	int Iteration(void) const;

	AllRunsResults& operator=(const AllRunsResults& right);

	friend AnalysisManager;
	friend SingleRunData;
};


#endif
