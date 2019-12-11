#ifndef ALL_EVENTS_RESULTS_H
#define ALL_EVENTS_RESULTS_H

#include "GlobalParameters.h"
#include "SingleEventData.h"

class AnalysisManager;
class SingleEventData;

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
class AllEventsResults
{
public:
	class AverageData {
	public:
		DVECTOR xs_sum;
		DVECTOR ys_sum;
		DVECTOR ys_disp;
		std::size_t average_event_n;
		std::size_t dispersion_event_n;
		AverageData() : average_event_n(0), dispersion_event_n(0) {}
	};
protected:
	STD_CONT<SingleEventData> events_data;
	int Iteration_N;
	const ParameterPile::experiment_manifest* processing_manifest;
	GraphCollection graph_manager;
	indexed_info<GraphCollection> pictures; //per channel, all plots for the same channel (for all runs/events) are stored in single GraphCollection
	indexed_info<AverageData> averages;

#ifdef _USE_TIME_STATISTICS
	time_results time_stat;
	void report_time_statistics();
#endif
public:
	AllEventsResults(const ParameterPile::experiment_manifest* to_process);//only experiment and channels are important here
	//For multithreading:
	void Merge(AllEventsResults* with);
	void Merged(void);
	void Clear(void);//called at the end of Merged() and in Merge(with) for 'with'
	void ClearMerged(void); //called for Merged results in multithread case
	int Iteration(void) const;

	AllEventsResults& operator=(const AllEventsResults& right);

	friend AnalysisManager;
	friend SingleEventData;
};


#endif
