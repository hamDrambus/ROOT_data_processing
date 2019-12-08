#ifndef ANALYSIS_MANAGER_H
#define ANALYSIS_MANAGER_H

#include "TH1D.h"
#include "TH1I.h"
#include "GlobalParameters.h"
#include "SingleEventData.h"
#include "AllEventsResults.h"

class AnalysisManager
{
public:
	AnalysisManager(ParameterPile::analysis_manifest to_process);
protected:
	ParameterPile::analysis_manifest manifest_all;
	ParameterPile::experiment_manifest manifest_under_processing; //corresponds to single experiment/folder
	ParameterPile::experiment_manifest manifest_single_event; //Contains only single run and subrun
	std::size_t index_manifest_under_processing;

	//STD_CONT<SingleRunData> one_run_data;
	STD_CONT<AllEventsResults> all_events_results;

	TCondition* _cond;
	TMutex* _thread_mutex;

	enum NextEventIs {NewSubRun, NewRun, NewExperiment, Null} curr_run;
	virtual void processOneEvent_first_iteration(AllEventsResults *_all_results);
	virtual void nextEvent(void);
	virtual void nextExperiment(void);
	virtual void processAllEvents(void);
	virtual void loopAllEvents_first_iteration(AllEventsResults *_all_results);
	virtual void loopAllEvents(AllEventsResults *_all_results);
	ParameterPile::experiment_manifest refine_exp_area(ParameterPile::experiment_manifest area);//looks up the existing runs in data directory 
	//and intersects them with input area (from ParameterPile::experiment_manifest). This is required in order to split runs between threads equally
public:
	virtual void processAllExperiments(void);
	virtual void proceessAllEventsOneThread(void);//Process only single experiment (folder)
	STD_CONT<AllEventsResults>* getAllEventsResults(void);
	void setAllEventsResults(AllEventsResults* to_what);
	void setCondition(TCondition* cond);
	TCondition* getCondition(void);
	void setThreadMutex(TMutex* mutex);
	TMutex* getThreadMutex(void);
};

#endif
