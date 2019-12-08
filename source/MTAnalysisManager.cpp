#include "MTAnalysisManager.h"

void process_runs_in_thread(void* manager)
{
	((AnalysisManager*)manager)->proceessAllEventsOneThread();
	TCondition* cond = ((AnalysisManager*)manager)->getCondition();
	TMutex* mutex = ((AnalysisManager*)manager)->getThreadMutex();
	if (0 != mutex->TryLock()){//means that the main thread is waiting for the signal
		cond->Signal();
		//std::cout << "Signal()" << std::endl;
	}
	//std::cout << "Exiting thread" << std::endl;
}

void MTAnalysisManager::processOneRun(void)
{
	//not called
}

void MTAnalysisManager::nextRun(void)
{
	if (NextEventIs::Null == curr_run) {
		index_manifest_under_processing = 0;
		if (index_manifest_under_processing >= manifest_all.manifests.size()) {
			curr_run = NextEventIs::Null;
			return;
		}
		manifest_under_processing = manifest_all.manifests[index_manifest_under_processing];
		curr_run = NextEventIs::NewExperiment;
		return;
	}
	if (++index_manifest_under_processing < manifest_all.manifests.size()) {
		manifest_under_processing = manifest_all.manifests[index_manifest_under_processing];
		curr_run = NextEventIs::NewExperiment;
		return;
	} else {
		manifest_under_processing = ParameterPile::experiment_manifest();
		curr_run = NextEventIs::Null;
		return;
	}
	//No runs or subruns switch, the MultithreadAnalysisManager is only responsible for experiments switching and runs splitting
}

void MTAnalysisManager::processAllRuns(void)
{
	std::vector<TThread*> pThreads;
	std::vector<AnalysisManager*> _submanagers;
	std::vector<TMutex*> mutexes;
	std::vector<TMutex*> thread_mutexes;
	std::vector<TCondition*> conditions;
	
	ParameterPile::experiment_manifest actual_area = refine_exp_area(manifest_under_processing);
	STD_CONT<ParameterPile::experiment_manifest> areas = split_exp_area(actual_area, ParameterPile::threads_number);
	all_events_results.push_back(AllEventsResults(&actual_area));
	for (int n = 0; n<areas.size(); ++n) {
		mutexes.push_back(new TMutex());
		conditions.push_back(new TCondition(mutexes[n]));
		thread_mutexes.push_back (new TMutex());
		_submanagers.push_back(new AnalysisManager(ParameterPile::analysis_manifest(areas[n]) ) );
		_submanagers[n]->setCondition(conditions[n]);
		_submanagers[n]->setThreadMutex(thread_mutexes[n]);
		pThreads.push_back (new TThread(("AnManager_" + manifest_under_processing.name + std::to_string(n)).c_str(),
			&process_runs_in_thread, _submanagers[n]));
	}

	while (all_events_results.back().Iteration() <= ParameterPile::Max_iteration_N) {
		for (int n = 0; n < areas.size(); ++n) {
			_submanagers[n]->setAllEventsResults(&all_events_results.back()); //Creates AllEventsResults for each thread at first iteration.
			//Used for passing parameters determined from all runs to the next iteration. Assignment is overloaded, so only some of the data is copied!
			pThreads[n]->Run(); //if it is the last iteration, submanager (AnalysisManager) clears its data
		}
		//TThread::Ps();
		for (int n = 0; n<areas.size(); ++n) {
			if (0 != thread_mutexes[n]->TryLock()) { //thread is already executed, so no wait required
			} else {
				conditions[n]->Wait();
			}
			thread_mutexes[n]->UnLock();
			STD_CONT<AllEventsResults> *res = _submanagers[n]->getAllEventsResults();
			if (!res->empty())
				all_events_results.back().Merge(&res->back());
		}
		all_events_results.back().Merged();
		all_events_results.back().ClearMerged();
	}
	
	for (int n = 0; n<areas.size(); ++n) {
		pThreads[n]->Delete();
		delete _submanagers[n];
		conditions[n]->Delete();
		delete mutexes[n];
		delete thread_mutexes[n];
	}
	nextRun();
}

void MTAnalysisManager::loopAllRuns(void)
{
	//not called
}
void MTAnalysisManager::loopAllRuns(AllEventsResults *_all_results)
{
	//not called
}

MTAnalysisManager::MTAnalysisManager(ParameterPile::analysis_manifest area) : AnalysisManager(area)
{}

void MTAnalysisManager::processAllExperiments(void)
{
	AnalysisManager::processAllExperiments();
}

STD_CONT<ParameterPile::experiment_manifest> MTAnalysisManager::split_exp_area(ParameterPile::experiment_manifest area_to_split, int N)
{
	STD_CONT <ParameterPile::experiment_manifest> out_;
	STD_CONT<ParameterPile::area_vector> Runs = area_to_split.runs.split_area(N);
	for (int h = 0; h < Runs.size(); h++){
		out_.push_back(area_to_split);
		out_[h].runs = Runs[h];
	}
	return out_;
}
