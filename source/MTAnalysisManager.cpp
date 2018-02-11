#include "MTAnalysisManager.h"

void process_runs_in_thread(void* manager)
{
	((AnalysisManager*)manager)->proceessAllRunsOneThread();
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
	if (NextRunIs::Null == curr_run){
		current_under_processing.experiments.push_back(_exp_area.experiments[0]);
		current_under_processing.runs=_exp_area.runs;
		current_under_processing.sub_runs = _exp_area.sub_runs;
		current_under_processing.channels = _exp_area.channels;
		curr_run = NextRunIs::FirstRun;//change to NewSubRun ?
		return;
	}
	for (auto i = _exp_area.experiments.begin(); i != _exp_area.experiments.end(); i++)
		if (*i == current_under_processing.experiments.back())
			if (++i != _exp_area.experiments.end()){
				current_under_processing.experiments.back() = *i;
				curr_run = NextRunIs::NewExperiment;
				return;
			} else {
				current_under_processing.experiments.clear();
				current_under_processing.runs.erase();
				current_under_processing.sub_runs.erase();
				current_under_processing.channels.erase();
				curr_run = NextRunIs::Null;
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
	
	ParameterPile::experiment_area actual_area = refine_exp_area(current_under_processing);
	STD_CONT<ParameterPile::experiment_area> areas = split_exp_area(actual_area, ParameterPile::threads_number);
	all_runs_results.push_back(AllRunsResults(actual_area));
	for (int n =0;n<ParameterPile::threads_number; ++n) {
		mutexes.push_back(new TMutex());
		conditions.push_back(new TCondition(mutexes[n]));
		thread_mutexes.push_back (new TMutex());
		_submanagers.push_back(new AnalysisManager(areas[n]));
		_submanagers[n]->setCondition(conditions[n]);
		_submanagers[n]->setThreadMutex(thread_mutexes[n]);
		pThreads.push_back (new TThread(("AnManager_" + current_under_processing.experiments.back() + std::to_string(n)).c_str(),
			&process_runs_in_thread, _submanagers[n]));
	}

	while (all_runs_results.back().Iteration() <= ParameterPile::Max_iteration_N) {
		for (int n = 0; n < ParameterPile::threads_number; ++n) {
			_submanagers[n]->setAllRunsResults(&all_runs_results.back());
			pThreads[n]->Run(); //if it is the last iteration, submanager (AnalysisManager) clears its data
		}
		//TThread::Ps();
		all_runs_results.back().Clear();
		for (int n = 0; n<ParameterPile::threads_number; ++n) {
			if (0 != thread_mutexes[n]->TryLock()){ //thread is already executed, so no wait required
			} else {
				conditions[n]->Wait();
			}
			thread_mutexes[n]->UnLock();
			STD_CONT<AllRunsResults> *res = _submanagers[n]->getAllRunsResults();
			if (!res->empty())
				all_runs_results.back().Merge(&res->back());
		}
		all_runs_results.back().Merged();
	}
	
	for (int n = 0; n<ParameterPile::threads_number; ++n){
		pThreads[n]->Delete();
		delete _submanagers[n];
		conditions[n]->Delete();
		mutexes[n]->Delete();
		thread_mutexes[n]->Delete();
	}
	nextRun();
}

void MTAnalysisManager::loopAllRuns(void)
{
	//not called
}
void MTAnalysisManager::loopAllRuns(AllRunsResults *_all_results)
{
	//not called
}

MTAnalysisManager::MTAnalysisManager(ParameterPile::experiment_area area) : AnalysisManager(area)
{
}

void MTAnalysisManager::processAllExperiments(void)
{
	AnalysisManager::processAllExperiments();
}

STD_CONT<ParameterPile::experiment_area> MTAnalysisManager::split_exp_area(ParameterPile::experiment_area area_to_split, int N)
{
	STD_CONT <ParameterPile::experiment_area> out_;
	STD_CONT<ParameterPile::area_vector> Runs = area_to_split.runs.split_area(N);
	for (int h = 0; h < N; h++){
		out_.push_back(area_to_split);
		out_[h].runs = Runs[h];
	}
	return out_;
}
