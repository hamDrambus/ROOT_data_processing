#include "MTAnalysisManager.h"

void process_runs_in_thread(void* manager)
{
	((AnalysisManager*)manager)->proceessAllEventsOneThread();
	//std::cout << "Exiting thread" << std::endl;
}

void MTAnalysisManager::processOneEvent(void)
{
	//not called
}

void MTAnalysisManager::processAllEvents(void)
{
	std::vector<std::thread> pThreads;
	std::vector<AnalysisManager*> _submanagers;
	
	ParameterPile::experiment_manifest actual_area = refine_exp_area(manifest_under_processing);
	STD_CONT<ParameterPile::experiment_manifest> areas = split_exp_area(actual_area, ParameterPile::threads_number);
	all_events_results.push_back(AllEventsResults(&actual_area));
	for (int n = 0; n<areas.size(); ++n) {
		_submanagers.push_back(new AnalysisManager(ParameterPile::analysis_manifest(areas[n]) ) );
		pThreads.push_back(std::thread());
		//pThreads.push_back (new TThread(("AnManager_" + manifest_under_processing.name + std::to_string(n)).c_str(), \
			&process_runs_in_thread, _submanagers[n]));
	}

	while (all_events_results.back().Iteration() <= ParameterPile::Max_iteration_N) {
		for (int n = 0; n < areas.size(); ++n) {
			_submanagers[n]->setAllEventsResults(&all_events_results.back()); //Creates AllEventsResults for each thread at first iteration.
			//Used for passing parameters determined from all runs to the next iteration. Assignment is overloaded, so only some of the data is copied!
			pThreads[n] = std::thread(process_runs_in_thread, _submanagers[n]);//if it is the last iteration, submanager (AnalysisManager) clears its data
		}
		for (int n = 0; n<areas.size(); ++n) {
			pThreads[n].join();
			STD_CONT<AllEventsResults> *res = _submanagers[n]->getAllEventsResults();
			if (!res->empty())
				all_events_results.back().Merge(&res->back());
		}
		all_events_results.back().Merged();
		all_events_results.back().ClearMerged();
	}
	
	for (int n = 0; n<areas.size(); ++n) {
		delete _submanagers[n];
	}
	nextExperiment();
}

void MTAnalysisManager::loopAllEvents(void)
{
	//not called
}
void MTAnalysisManager::loopAllEvents(AllEventsResults *_all_results)
{
	//not called
}

MTAnalysisManager::MTAnalysisManager(ParameterPile::analysis_manifest area) : AnalysisManager(area)
{
	className = "MTAnalysisManager";
}

void MTAnalysisManager::processAllExperiments(void)
{
	AnalysisManager::processAllExperiments();
}

STD_CONT<ParameterPile::experiment_manifest> MTAnalysisManager::split_exp_area(ParameterPile::experiment_manifest area_to_split, int N)
{
	STD_CONT <ParameterPile::experiment_manifest> out_;
	STD_CONT<ParameterPile::area_vector> Runs = area_to_split.runs.split_area(N);
	for (int h = 0; h < Runs.size(); h++) {
		out_.push_back(area_to_split);
		out_[h].runs = Runs[h];
	}
	return out_;
}
