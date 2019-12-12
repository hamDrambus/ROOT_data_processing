#ifndef MT_ANALYSIS_MANAGER_H
#define MT_ANALYSIS_MANAGER_H

#include "TMutex.h"
#include "TThread.h"
#include "AnalysisManager.h"

void process_runs_in_thread(void* manager);

class MTAnalysisManager : public AnalysisManager
{
protected:
	
	virtual void processOneEvent(void);
	virtual void nextEvent(void);
	virtual void processAllEvents(void);
	virtual void loopAllEvents(void);
	virtual void loopAllEvents(AllEventsResults *_all_results);

	STD_CONT<ParameterPile::experiment_manifest> split_exp_area(ParameterPile::experiment_manifest area_to_split, int N);
public:
	MTAnalysisManager(ParameterPile::analysis_manifest area);

	virtual void processAllExperiments(void);
};

#endif
