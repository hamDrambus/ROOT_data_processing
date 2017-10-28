#ifndef ANALYSIS_MANAGER_H
#define ANALYSIS_MANAGER_H

#include "TH1D.h"
#include "TH1I.h"
#include "GlobalDefinitions.h"
#include "SingleRunData.h"
#include "AllRunsResults.h"

class AnalysisManager
{
public:
	AnalysisManager(ParameterPile::experiment_area area);
protected:
	ParameterPile::experiment_area _exp_area;
	//ParameterPile::experiment_area first_level_processed;  //per run processed - e.g. filters
	//ParameterPile::experiment_area second_level_processed; //per experiment processed - averaged per runs.
	//ParameterPile::experiment_area third_level_processed;  //everything
	ParameterPile::experiment_area current_under_processing;

	std::vector<SingleRunData> one_run_data;
	std::vector<SingleRunResults> one_run_results;
	std::vector<AllRunsResults> all_runs_results;
#ifdef _TEMP_CODE
	void* all_exps_data;

	//void* all_exps_result;
#endif

	enum NextRunIs { FirstRun, NewSubRun, NewRun, NewExperiment, Null} curr_run;
	virtual void processOneRun(void);
	virtual void nextRun(void);
	virtual void processAllRuns(void);
	virtual void loopAllRuns(void);

public:
	//virtual void nextExperiment(void);

	virtual void processAllExperiments(void);
	
	virtual void save_all_runs(void);
	virtual void save_all_exps(void);

	virtual bool check_run_processed(void);
	virtual void save_one_run_results(void);
	virtual void load_one_run_results(void);

	//double minimize_function(double value);
};

#endif