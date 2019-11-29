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
		available_runs.reset();
		current_under_processing.experiments.push_back(_exp_area.experiments[0]);
		current_under_processing.runs = _exp_area.runs;
		current_under_processing.sub_runs.push_back(_exp_area.sub_runs.get_next_index());
		current_under_processing.channels = _exp_area.channels;
		current_under_processing = refine_exp_area(current_under_processing);
		available_runs = current_under_processing.runs;
		current_under_processing.runs.erase();
		current_under_processing.runs.push_back(available_runs.get_next_index());
		curr_run = NextRunIs::FirstRun;//change to NewSubRun ?
		return;
	}
	current_under_processing.sub_runs.back() = _exp_area.sub_runs.get_next_index();
	if (current_under_processing.sub_runs.back() < 0){
		current_under_processing.sub_runs.back() = _exp_area.sub_runs.get_next_index();
		current_under_processing.runs.back() = available_runs.get_next_index();
		if (current_under_processing.runs.back() < 0){
			for (auto i = _exp_area.experiments.begin(); i != _exp_area.experiments.end(); i++)
				if (*i == current_under_processing.experiments.back())
					if (++i != _exp_area.experiments.end()){
						current_under_processing.experiments.back() = *i;
						curr_run = NextRunIs::NewExperiment;
						available_runs.erase();
						break;
					} else {
						current_under_processing.experiments.pop_back();
						current_under_processing.runs.erase();
						current_under_processing.sub_runs.erase();
						current_under_processing.channels.erase();
						available_runs.erase();
						curr_run = NextRunIs::Null;
						return;
					}
			current_under_processing.runs = _exp_area.runs;
			current_under_processing = refine_exp_area(current_under_processing);
			available_runs = current_under_processing.runs;
			current_under_processing.runs.push_back(available_runs.get_next_index());
			return;
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
	one_run_data.back().processSingleRun(_all_results);
	if (!_all_results->_valid[_all_results->N_of_runs-1]){
		if (_all_results->_status[_all_results->N_of_runs-1] != AllRunsResults::ExternalRejected) {
			std::cout << "invalid: " << current_under_processing.experiments.back() << "/run_" << current_under_processing.runs.back() << "_sub_"
				<< current_under_processing.sub_runs.back() << " processed" << std::endl;
			std::cout << "reason: " << _all_results->_status[_all_results->N_of_runs-1]<<std::endl;
		}
		//one_run_data.pop_back();
	} else {
		//std::size_t data_size =sizeof(one_run_data);
		//for (std::size_t g=0,_end_=one_run_data.size();g!=_end_;++g)
		//	data_size += one_run_data[g].real_size();

		std::cout<<"processed"<<_all_results->Iteration()<<": "<< current_under_processing.experiments.back() << "_run_" << current_under_processing.runs.back() << "_sub_"
			<< current_under_processing.sub_runs.back() << std::endl;
			//<<"\tSingleRunData.size: "<<data_size <<"\tN: "<<one_run_data.size()<< std::endl;

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
	//std::size_t data_size =sizeof(one_run_data);
	for (auto j = one_run_data.begin(); (j != one_run_data.end()); ++j) {
		j->processSingleRun(_all_results);
		if (!_all_results->_valid[_all_results->N_of_runs-1]) {
			if (_all_results->_status[_all_results->N_of_runs-1] != AllRunsResults::ExternalRejected) {
				std::cout << "invalid: " << j->getArea().experiments.back() << "/run_" << j->getArea().runs.back() << "_sub_"
					<< j->getArea().sub_runs.back() << " processed" << std::endl;
				std::cout << "reason: " << _all_results->_status[_all_results->N_of_runs-1]<<std::endl;
			}
		}
		//data_size += j->real_size();
		std::cout << "processed"<<_all_results->Iteration()<<": "<< j->getArea().experiments.back() << "_run" << j->getArea().runs.back() << "_sub"
			<< j->getArea().sub_runs.back() << std::endl; //<< "\tSingleRunData.size: "<<data_size<<std::endl;
	}
}

void AnalysisManager::processAllRuns(void)
{
	all_runs_results.push_back(AllRunsResults(current_under_processing));
	while (all_runs_results.back().Iteration() <= ParameterPile::Max_iteration_N) {
#ifdef _USE_TIME_STATISTICS
		time_t runs_start_timer;
		time_t runs_end_timer;
		time(&runs_start_timer);
#endif
		loopAllRuns(&all_runs_results.back());
#ifdef _USE_TIME_STATISTICS
		time(&runs_end_timer);
		all_runs_results.back().time_stat.n_RUN_proc = all_runs_results.back().N_of_runs;
		all_runs_results.back().time_stat.t_RUN_proc += difftime(runs_end_timer, runs_start_timer);
#endif
		all_runs_results.back().Merged();//all in one thread, so it is already joined. Merge == iteration++
	}
#ifdef _HOTFIX_CLEAR_MEMORY
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

	if (0 == all_runs_results.back().Iteration())
		while (curr_run != NextRunIs::Null)
			nextRun();//skip the rest of experiments if any

	if (ParameterPile::Max_iteration_N == all_runs_results.back().Iteration()) {
		STD_CONT<SingleRunData>().swap(one_run_data);
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


ParameterPile::experiment_area AnalysisManager::refine_exp_area(ParameterPile::experiment_area area)//looks up the existing runs in data directory 
//and intersects them with input area (from ParameterPile::exp_area). This is required in order to split runs between threads equally
//depr: TODO: maybe also move to the AnalysisManager
{
	ParameterPile::experiment_area out_area = area;
	out_area.runs.erase();
	std::vector<int> runs;
	int from = -1, to = -1;
#if defined(__WIN32__)
	HANDLE dir;
	WIN32_FIND_DATA file_data;
	std::string path = DATA_PREFIX;
	path += area.experiments.back();
	if ((dir = FindFirstFile((path + "/*").c_str(), &file_data)) == INVALID_HANDLE_VALUE)
		return out_area;
	do {
		std::string file_name = file_data.cFileName;
		if ((file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0)
			continue;
		if (file_name.size() < 4)
			continue;
		if (file_name[0] == '.')
			continue;
		file_name.erase(file_name.begin(), file_name.begin() + 4); //erase "run_"
		int n_underscore = file_name.find("_");
		if (n_underscore == std::string::npos)
			continue;
		file_name.erase(file_name.begin() + n_underscore, file_name.end());
		if (file_name.empty())
			continue;
		int run = std::stoi(file_name);

		if (from < 0){
			from = run;
			to = run;
			continue;
		}
		if (to == run)
			continue;
		if (to == (run - 1))
			to++;
		else {
			out_area.runs.push_pair(from, to);
			from = run;
			to = from;
		}
	} while (FindNextFile(dir, &file_data));
	FindClose(dir);
#else //__WIN32__

	DIR *dp;
    struct dirent *dirp;
	std::string path = DATA_PREFIX;
	path += area.experiments.back();
    if((dp  = opendir(path.c_str())) == NULL) {
        std::cout << "Error(" << errno << ") opening " << path << std::endl;
        return out_area;
    }
    std::vector<int> run_vector_sorted;
    while ((dirp = readdir(dp)) != NULL) {
		std::string file_name = dirp->d_name;
		if (file_name.size() < 4)
			continue;
		if (file_name[0] == '.')
			continue;
		file_name.erase(file_name.begin(), file_name.begin() + 4); //erase "run_"
		int n_underscore = file_name.find("_");
		if (n_underscore == std::string::npos)
			continue;
		file_name.erase(file_name.begin() + n_underscore, file_name.end());
		if (file_name.empty())
			continue;
		run_vector_sorted.push_back(std::stoi(file_name));
    }
    closedir(dp);
    std::sort(run_vector_sorted.begin(), run_vector_sorted.end());
	for (std::size_t r=0, r_end_= run_vector_sorted.size(); r!=r_end_; ++r) {
        int run = run_vector_sorted[r];
		if (from < 0){
			from = run;
			to = run;
			continue;
     	}
      	if (to == run)
       		continue;
       	if (to == (run - 1))
        	to++;
       	else {
       		out_area.runs.push_pair(from, to);
       		from = run;
       		to = from;
       	}
	}

#endif //__WIN32__
    if ((from >= 0) && (to >= 0))
    		out_area.runs.push_pair(from, to);
    out_area.runs = out_area.runs.intersect(area.runs);
	return out_area;
}
