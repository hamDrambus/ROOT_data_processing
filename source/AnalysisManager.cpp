#include "AnalysisManager.h"

AnalysisManager::AnalysisManager(ParameterPile::experiment_area exp_area)
{
	_exp_area = exp_area;
	curr_run = NextRunIs::Null;
	_cond = NULL;
	_thread_mutex = NULL;
}

void AnalysisManager::nextRun(void)
{
	if (NextRunIs::Null==curr_run){
		_exp_area.channels.reset();//clears flags, but not the elements, for the valid get_next_index
		_exp_area.runs.reset();
		_exp_area.sub_runs.reset();
		current_under_processing.experiments.push_back(_exp_area.experiments[0]);
		current_under_processing.runs.push_back(_exp_area.runs.get_next_index());
		current_under_processing.sub_runs.push_back(_exp_area.sub_runs.get_next_index());
		current_under_processing.channels= _exp_area.channels;
		curr_run = NextRunIs::FirstRun;//change to NewSubRun ?
		return;
	}
	current_under_processing.sub_runs.back() = _exp_area.sub_runs.get_next_index();
	if (current_under_processing.sub_runs.back() < 0){
		current_under_processing.sub_runs.back() = _exp_area.sub_runs.get_next_index();
		current_under_processing.runs.back() = _exp_area.runs.get_next_index();
		if (current_under_processing.runs.back() < 0){
			current_under_processing.runs.back() = _exp_area.runs.get_next_index();
			for (auto i = _exp_area.experiments.begin(); i != _exp_area.experiments.end(); i++)
				if (*i == current_under_processing.experiments.back())
					if (++i != _exp_area.experiments.end()){
						current_under_processing.experiments.back() = *i;
						curr_run = NextRunIs::NewExperiment;
						return;
					} else {
						current_under_processing.experiments.pop_back();
						current_under_processing.runs.erase();
						current_under_processing.sub_runs.erase();
						current_under_processing.channels.erase();
						curr_run = NextRunIs::Null;
						return;
					}
		} else {
			curr_run = NextRunIs::NewRun;
			return;
		}
	} else {
		curr_run = NextRunIs::NewSubRun;
		return;
	}
}

void AnalysisManager::processOneRun_first_iteration(AllRunsResults *_all_results) //current_under_processing must be preset
{
	one_run_data.push_back(SingleRunData(current_under_processing));
	one_run_results.push_back(one_run_data.back().processSingleRun(_all_results));
	if (!one_run_results.back().isValid()){
		/*std::cout << "invalid: " << current_under_processing.experiments.back() << "_run_" << current_under_processing.runs.back() << "_sub_"
			<< current_under_processing.sub_runs.back() << "_processed" << std::endl;
		std::cout << "reason: " << one_run_results.back().getStatus()<<std::endl;*/
		//std::cout << "S=" << one_run_results.back().PMT3_summed_peaks_area << "| N = " << one_run_results.back().PMT3_n_peaks << std::endl;
		one_run_results.pop_back();
		one_run_data.pop_back();
	} else {
		std::cout<<"processed: "<< current_under_processing.experiments.back() << "_run_" << current_under_processing.runs.back() << "_sub_"
					<< current_under_processing.sub_runs.back() << std::endl;
	}
}

void AnalysisManager::loopAllRuns_first_iteration(AllRunsResults *_all_results)
{
	while ((curr_run != NextRunIs::Null) && (curr_run != NextRunIs::NewExperiment)) {
		processOneRun_first_iteration(_all_results);
		nextRun();
	}
	//processAllRuns();
}

void AnalysisManager::loopAllRuns(AllRunsResults *_all_results)
{
	if (0 == _all_results->Iteration())
		return loopAllRuns_first_iteration(_all_results);
	auto i = one_run_results.begin();
	auto j = one_run_data.begin();
	for (; ((i != one_run_results.end()) && (j != one_run_data.end())); ++i, ++j){
		*i = j->processSingleRun(_all_results);
		std::cout << "processed: " << j->getArea().experiments.back() << "_run" << j->getArea().runs.back() << "_sub"
			<< j->getArea().sub_runs.back() << std::endl;
	}
}

void AnalysisManager::processAllRuns(void)
{
	all_runs_results.push_back(AllRunsResults(current_under_processing));
	while (all_runs_results.back().Iteration() <= ParameterPile::Max_iteration_N){
#ifdef _USE_TIME_STATISTICS
		time_t runs_start_timer;
		time_t runs_end_timer;
		time(&runs_start_timer);
#endif
		loopAllRuns(&all_runs_results.back());
		all_runs_results.back().Clear();
		all_runs_results.back().processAllRuns(one_run_results);
#ifdef _USE_TIME_STATISTICS
		time(&runs_end_timer);
		all_runs_results.back().time_stat.n_RUN_proc = all_runs_results.back().N_of_runs;
		all_runs_results.back().time_stat.t_RUN_proc += difftime(runs_end_timer, runs_start_timer);
#endif
		all_runs_results.back().Merged();//all in one thread, so it is already joined. Merge == iteration++
	}
#ifdef _HOTFIX_CLEAR_MEMORY
	STD_CONT<SingleRunResults>().swap(one_run_results);
	STD_CONT<SingleRunData>().swap(one_run_data);
#else
	one_run_data.clear();
	one_run_results.clear();
#endif
}

void AnalysisManager::processAllExperiments(void)
{
	nextRun();
	while (curr_run != NextRunIs::Null){
		curr_run = NextRunIs::NewSubRun;
		processAllRuns();
	}
	//TODO: average, sumarize over all of the experiments
}

void AnalysisManager::proceessAllRunsOneThread(void)
{
#ifdef _USE_TIME_STATISTICS
	time_t runs_start_timer;
	time_t runs_end_timer;
	time(&runs_start_timer);
#endif
	if (all_runs_results.empty()){
		nextRun();
		all_runs_results.push_back(AllRunsResults(current_under_processing)); //actually MultiThreadManager must set it before call of this method
	} else {
		if (0 == all_runs_results.back().Iteration())
			nextRun();
	}
	
	loopAllRuns(&all_runs_results.back());
	all_runs_results.back().processAllRuns(one_run_results);
	//std::cout << "one_run_results.size()==" << one_run_results.size()<<std::endl;

	if (0 == all_runs_results.back().Iteration())
		while (curr_run != NextRunIs::Null)
			nextRun();//skip the rest of experiments if any

	if (ParameterPile::Max_iteration_N == all_runs_results.back().Iteration()) {
#ifdef _HOTFIX_CLEAR_MEMORY
		STD_CONT<SingleRunData>().swap(one_run_data);
		STD_CONT<SingleRunResults>().swap(one_run_results);
#else
		one_run_data.clear();
		one_run_results.clear();
#endif
	}
#ifdef _USE_TIME_STATISTICS
	time(&runs_end_timer);
	all_runs_results.back().time_stat.n_RUN_proc_single_iter = all_runs_results.back().N_of_runs;
	all_runs_results.back().time_stat.t_RUN_proc_single_iter = difftime(runs_end_timer,runs_start_timer);
#endif
	//std::cout << "finish_proceessAllRunsOneThread" << std::endl;
}

STD_CONT<AllRunsResults>* AnalysisManager::getAllRunsResults(void)
{
	return &all_runs_results;
}

void AnalysisManager::setAllRunsResults(AllRunsResults* to_what)
{
	if (all_runs_results.empty())
		all_runs_results.push_back(*to_what);
	else
		all_runs_results.back() = *to_what;
}

void AnalysisManager::setCondition(TCondition* cond)
{	_cond = cond;}
TCondition* AnalysisManager::getCondition(void)
{	return _cond;}
void AnalysisManager::setThreadMutex(TMutex* mutex)
{	_thread_mutex = mutex;}
TMutex* AnalysisManager::getThreadMutex(void)
{	return _thread_mutex;}
