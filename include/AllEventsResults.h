#ifndef ALL_EVENTS_RESULTS_H
#define ALL_EVENTS_RESULTS_H

#include "GlobalParameters.h"
#include "SingleEventData.h"

class AnalysisManager;
class SingleEventData;

#ifdef _USE_TIME_STATISTICS
struct time_results {
	//proc - processing
	long long t_total_proc;
	std::size_t n_total_proc;
	long long t_file_reading;
	std::size_t n_file_reading;
	long long t_filtering;
	std::size_t n_filtering;
	long long t_simple_baseline;
	std::size_t n_simple_baseline;
	long long t_curved_baseline;
	std::size_t n_curved_baseline;
	long long t_peaks;
	std::size_t n_peaks;
	long long t_integrals;
	std::size_t n_integrals;
	long long t_double_integrals;
	std::size_t n_double_integrals;
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
