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
		current_under_processing.experiments.push_back(_exp_area.experiments[0]);
		current_under_processing.runs.push_back(_exp_area.runs[0]);
		current_under_processing.sub_runs.push_back(_exp_area.sub_runs[0]);
		current_under_processing.channels= _exp_area.channels;
		curr_run = NextRunIs::FirstRun;//change to NewSubRun ?
		return;
	}
	current_under_processing.sub_runs.back() = get_next_index(_exp_area.sub_runs, current_under_processing.sub_runs.back());
	if (current_under_processing.sub_runs.back() < 0){
		current_under_processing.sub_runs.back() = _exp_area.sub_runs[0];
		current_under_processing.runs.back() = get_next_index(_exp_area.runs, current_under_processing.runs.back());
		if (current_under_processing.runs.back() < 0){
			current_under_processing.runs.back() = _exp_area.runs[0];
			for (auto i = _exp_area.experiments.begin(); i != _exp_area.experiments.end(); i++)
				if (*i == current_under_processing.experiments.back())
					if (++i != _exp_area.experiments.end()){
						current_under_processing.experiments.back() = *i;
						curr_run = NextRunIs::NewExperiment;
						return;
					} else {
						current_under_processing.experiments.pop_back();
						current_under_processing.runs.pop_back();
						current_under_processing.sub_runs.pop_back();
						current_under_processing.channels.clear();
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

void AnalysisManager::processOneRun(void) //current_under_processing must be preset
{
	one_run_data.push_back(SingleRunData(current_under_processing));
	one_run_results.push_back(one_run_data.back().processSingleRun());
	if (!one_run_results.back().isValid()){
		/*std::cout << "invalid: " << current_under_processing.experiments.back() << "_run_" << current_under_processing.runs.back() << "_sub_"
			<< current_under_processing.sub_runs.back() << "_processed" << std::endl;*/
		//std::cout << "reason: " << one_run_results.back().getStatus()<<std::endl;
		//std::cout << "S=" << one_run_results.back().PMT3_summed_peaks_area << "| N = " << one_run_results.back().PMT3_n_peaks << std::endl;
		one_run_results.pop_back();
		one_run_data.pop_back();
	} else {
		/*std::cout<<"processed: "<< current_under_processing.experiments.back() << "_run_" << current_under_processing.runs.back() << "_sub_"
					<< current_under_processing.sub_runs.back() << "_processed" << std::endl;*/
	}
}

void AnalysisManager::loopAllRuns(void)
{
	while ((curr_run != NextRunIs::Null) && (curr_run != NextRunIs::NewExperiment)) {
		processOneRun();
		nextRun();
	}
	//processAllRuns();
}

void AnalysisManager::loopAllRuns(AllRunsResults *_all_results)
{
	auto i = one_run_results.begin();
	auto j = one_run_data.begin();
	for (; ((i != one_run_results.end())&&(j!=one_run_data.end())); ++i,++j){
		*i = j->processSingleRun(_all_results);
		//actually one_run_results and one_run_data must be the same size, so the use of "of_what" may be unnecessary
		//UPD: of_what is invalid actually
	}
}

void AnalysisManager::processAllRuns(void)
{
	loopAllRuns();
	all_runs_results.push_back(AllRunsResults(one_run_results.back().curr_area));
	all_runs_results.back().processAllRuns(one_run_results);
	all_runs_results.back().Joined();//all in one thread, so it is already joined
	loopAllRuns(&all_runs_results.back());
	all_runs_results.back().processAllRuns(one_run_results);
	all_runs_results.back().Joined();//all in one thread, so it is already joined

	one_run_data.clear();
	one_run_results.clear();
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
	nextRun();
	ParameterPile::experiment_area area = current_under_processing;
	if (all_runs_results.empty()) { //first iteration
		loopAllRuns();
		while (curr_run != NextRunIs::Null)
			nextRun();//skip the rest of experiments if any
		//std::cout << "one_run_results.size()==" << one_run_results.size()<<std::endl;
		all_runs_results.push_back(AllRunsResults(area));
		all_runs_results.back().processAllRuns(one_run_results); //one_run_results may be empty
	} else {
		loopAllRuns(&all_runs_results.back());
		all_runs_results.back().processAllRuns(one_run_results);
		one_run_data.clear();
		one_run_results.clear();
	}
	//std::cout << "finish_proceessAllRunsOneThread" << std::endl;
}

std::vector<AllRunsResults>* AnalysisManager::getAllRunsResults(void)
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