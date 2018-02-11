#ifndef MT_ANALYSIS_MANAGER_H
#define MT_ANALYSIS_MANAGER_H

#include "TThread.h"
#include "AnalysisManager.h"

void process_runs_in_thread(void* manager);

class MTAnalysisManager : public AnalysisManager
{
protected:
	
	virtual void processOneRun(void);
	virtual void nextRun(void);
	virtual void processAllRuns(void);
	virtual void loopAllRuns(void);
	virtual void loopAllRuns(AllRunsResults *_all_results);

	STD_CONT<ParameterPile::experiment_area> split_exp_area(ParameterPile::experiment_area area_to_split, int N);
public:
	MTAnalysisManager(ParameterPile::experiment_area area);

	virtual void processAllExperiments(void);
	//void* processAllExperiments(void*);

	//virtual void save_all_runs(void);
	//virtual void save_all_exps(void);

	//virtual bool check_run_processed(void);
	//virtual void save_one_run_results(void);
	//virtual void load_one_run_results(void);
};

#endif