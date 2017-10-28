#include "AnalysisManager.h"

AnalysisManager::AnalysisManager(ParameterPile::experiment_area exp_area)
{
	_exp_area = exp_area;
	curr_run = NextRunIs::Null;
#ifdef _TEMP_CODE //TODO: move this to the all_runs_data
	std::ofstream file;
	open_output_file(std::string(OUTPUT_DIR) + OUTPUT_GEMS, file);//creates folder AND trucates the file
	file << "Experiment\tIntegral[V*us]\tt_from\tt_to" << std::endl;
	file.close();
#endif
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
		one_run_results.pop_back();
		one_run_data.pop_back();
	} else {
		std::cout<<"processed: "<< current_under_processing.experiments.back() << "_run_" << current_under_processing.runs.back() << "_sub_"
					<< current_under_processing.sub_runs.back() << "_processed" << std::endl;
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
	loopAllRuns(&all_runs_results.back());
	all_runs_results.back().processAllRuns(one_run_results);

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

void AnalysisManager::save_all_runs(void)
{

}
void AnalysisManager::save_all_exps(void)
{

}

bool AnalysisManager::check_run_processed(void)
{
	return false;
}
void AnalysisManager::save_one_run_results(void)
{

}
void AnalysisManager::load_one_run_results(void)
{

}
