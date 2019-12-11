#ifndef SINGLE_EVENT_DATA_H
#define SINGLE_EVENT_DATA_H

#include "GlobalParameters.h"
#include "GraphCollection.h"
#include "SignalOperations.h"
#include "AllEventsResults.h"

class AnalysisManager;
class AllEventsResults;
class TVirtualFFT;

class SingleEventData
{
public:
	enum Status {Ok, Empty, NoPMTsignal, PMT_mismatch, SiPM_mismatch, AVG_mismatch, ExternalRejected };
	struct ChannelData {
	public:
		DVECTOR xs;
		DVECTOR ys;
		double found_baseline;
		DVECTOR curved_bl_xs;
		DVECTOR curved_bl_ys;
		double curved_bl_baseline;
		double integral;
		double double_integral;
		STD_CONT<peak> peaks;
	};
	struct EventIndex {
		int run;
		int subrun;
	};
	Status status;
	bool isValid(void) const {
		return status == Status::Ok;
	}
	std::string Status(void) const {
		switch (status) {
		case Ok:
			return "Ok";
		case Empty:
			return "Empty";
		case NoPMTsignal:
			return "No PMT signal";
		case PMT_mismatch:
			return "PMT mismatch";
		case SiPM_mismatch:
			return "SiPM mismatch";
		case AVG_mismatch:
			return "Averages mismatch";
		case ExternalRejected:
			return "Externally rejected";
		}
	}

protected:

	indexed_info<ChannelData> channel_data;

	double trigger_offset;

	const ParameterPile::experiment_manifest *manifest;
	EventIndex index;
	
	void file_to_vector(std::string fname, DVECTOR &xs, DVECTOR &ys, int subrun);

	void readOneRun(AllEventsResults *results, int channel);
	void clearOneRun(int channel);
	void calculate_threshold_and_baseline(DVECTOR &xs, DVECTOR &ys, double &threshold, double &threshold_edges, double &baseline, STD_CONT<peak> &peaks_before_S1, int channel);
	void push_event (AllEventsResults *all_runs_results);
	void push_average (int ch, AllEventsResults *all_runs_results);
	void push_dispersion (int ch, AllEventsResults *all_runs_results);
public:
	SingleEventData(const ParameterPile::experiment_manifest *to_process);
	void processSingleEvent_Iter_0(AllEventsResults *all_runs_results);
	void processSingleEvent_Iter_1(AllEventsResults *all_runs_results);
	void processSingleEvent(AllEventsResults *all_runs_results);
	void eventProcessedProc(AllEventsResults *all_runs_results);
	void clear_memory(void); //clears only 'input' data, preserves processing results

	std::size_t real_size(void);
	friend AnalysisManager;
	friend AllEventsResults;
};

#endif
