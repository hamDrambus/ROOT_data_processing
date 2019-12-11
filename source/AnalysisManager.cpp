#include "AnalysisManager.h"

AnalysisManager::AnalysisManager(ParameterPile::analysis_manifest exp_area) : manifest_all(exp_area), 
	index_manifest_under_processing(0), curr_run(NextEventIs::Null), _cond(NULL), _thread_mutex(NULL)
{}

void AnalysisManager::nextExperiment(void)
{
	++index_manifest_under_processing;
	while (index_manifest_under_processing < manifest_all.manifests.size()) {
		manifest_under_processing = manifest_all.manifests[index_manifest_under_processing];
		manifest_under_processing = refine_exp_area(manifest_under_processing);
		manifest_under_processing.runs.reset(); //clears internal flags for correct get_next_index()
		manifest_under_processing.sub_runs.reset();
		manifest_single_event = manifest_under_processing;
		if (manifest_under_processing.runs.empty() || manifest_under_processing.sub_runs.empty()) {
			++index_manifest_under_processing;
			std::cout << "AnalysisManager::nextExperiment: Warnining: Folder \"" << manifest_under_processing.in_folder << "\" does not contain data to process!" << std::endl;
			continue;
		}
		manifest_single_event.runs.erase();
		manifest_single_event.runs.push_back(manifest_under_processing.runs.get_next_index());
		manifest_single_event.sub_runs.erase();
		manifest_single_event.sub_runs.push_back(manifest_under_processing.sub_runs.get_next_index());
		curr_run = NextEventIs::NewExperiment;
		return;
	}
	index_manifest_under_processing = 0;
	curr_run = NextEventIs::Null;
	return;
}

void AnalysisManager::nextEvent(void)
{
	if (NextEventIs::Null==curr_run) {
		index_manifest_under_processing = 0;
		--index_manifest_under_processing;
		nextExperiment();
		curr_run = NextEventIs::NewSubRun;
		return;
	}
	manifest_single_event.sub_runs.back() = manifest_under_processing.sub_runs.get_next_index();
	if (manifest_single_event.sub_runs.back() < 0) {
		manifest_single_event.sub_runs.back() = manifest_under_processing.sub_runs.get_next_index();
		manifest_single_event.runs.back() = manifest_under_processing.runs.get_next_index();
		if (manifest_single_event.runs.back() < 0) {
			nextExperiment();
		} else {
			curr_run = NextEventIs::NewRun;
			return;
		}
	} else {
		curr_run = NextEventIs::NewSubRun;
		return;
	}
}

void AnalysisManager::processOneEvent_first_iteration(AllEventsResults *_all_results) //manifest_single_event must be preset
{
	_all_results->events_data.push_back(SingleEventData(&manifest_single_event));
	_all_results->events_data.back().processSingleEvent(_all_results);
	if (!_all_results->events_data.back().isValid()) {
		if (_all_results->events_data.back().status != SingleEventData::ExternalRejected) {
			std::cout << "invalid: " << manifest_single_event.name << "/run_" << manifest_single_event.runs.back() << "_sub_"
				<< manifest_single_event.sub_runs.back() << " processed" << std::endl;
			std::cout << "reason: " << _all_results->events_data.back().Status() <<std::endl;
		}
		//_all_results->events_data.pop_back();
	} else {
		std::cout<<"processed"<<_all_results->Iteration()<<": "<< manifest_single_event.name << "_run_" << manifest_single_event.runs.back() << "_sub_"
			<< manifest_single_event.sub_runs.back() << std::endl;
	}
}

void AnalysisManager::loopAllEvents_first_iteration(AllEventsResults *_all_results)
{
	while ((curr_run != NextEventIs::Null) && (curr_run != NextEventIs::NewExperiment)) {
		processOneEvent_first_iteration(_all_results);
		nextEvent();
	}
}

void AnalysisManager::loopAllEvents(AllEventsResults *_all_results)
{
	if (0 == _all_results->Iteration())
		return loopAllEvents_first_iteration(_all_results);
	for (auto j = _all_results->events_data.begin(); (j != _all_results->events_data.end()); ++j) {
		j->processSingleEvent(_all_results);
		if (!j->isValid()) {
			if ( j->status != SingleEventData::ExternalRejected) {
				std::cout << "invalid: " << j->manifest->name << "/run_" << j->index.run << "_sub_"	<< j->index.subrun << " processed" << std::endl;
				std::cout << "reason: " << j->Status() << std::endl;
			}
		}
		std::cout << "processed"<<_all_results->Iteration()<<": "<< j->manifest->name << "_run" << j->index.run << "_sub"
			<< j->index.subrun << std::endl;
	}
}

void AnalysisManager::processAllEvents(void)
{
	all_events_results.push_back(AllEventsResults(&manifest_under_processing));
	while (all_events_results.back().Iteration() <= ParameterPile::Max_iteration_N) {
#ifdef _USE_TIME_STATISTICS
		time_t runs_start_timer;
		time_t runs_end_timer;
		time(&runs_start_timer);
#endif
		loopAllEvents(&all_events_results.back());
#ifdef _USE_TIME_STATISTICS
		time(&runs_end_timer);
		all_runs_results.back().time_stat.n_RUN_proc = all_runs_results.back().N_of_runs;
		all_runs_results.back().time_stat.t_RUN_proc += difftime(runs_end_timer, runs_start_timer);
#endif
		all_events_results.back().Merged();//all in one thread, so it is already joined. Merge == iteration++
	}
}

void AnalysisManager::processAllExperiments(void)
{
	nextEvent();
	while (curr_run != NextEventIs::Null){
		curr_run = NextEventIs::NewSubRun;
		processAllEvents();
	}
}

void AnalysisManager::proceessAllEventsOneThread(void)
{
#ifdef _USE_TIME_STATISTICS
	time_t runs_start_timer;
	time_t runs_end_timer;
	time(&runs_start_timer);
#endif
	if (all_events_results.empty()) {
		nextEvent();
		all_events_results.push_back(AllEventsResults(&manifest_under_processing)); //actually MultiThreadManager must set it before call of this method
	} else {
		if (0 == all_events_results.back().Iteration())
			nextEvent();
	}
	
	loopAllEvents(&all_events_results.back());

	if (0 == all_events_results.back().Iteration())
		while (curr_run != NextEventIs::Null)
			nextEvent();//skip the rest of experiments if any

#ifdef _USE_TIME_STATISTICS
	time(&runs_end_timer);
	all_runs_results.back().time_stat.n_RUN_proc_single_iter = all_runs_results.back().N_of_runs;
	all_runs_results.back().time_stat.t_RUN_proc_single_iter = difftime(runs_end_timer,runs_start_timer);
#endif
	//std::cout << "finish_proceessAllRunsOneThread" << std::endl;
}

STD_CONT<AllEventsResults>* AnalysisManager::getAllEventsResults(void)
{
	return &all_events_results;
}

void AnalysisManager::setAllEventsResults(AllEventsResults* to_what)
{
	if (all_events_results.empty())
		all_events_results.push_back(*to_what);
	else
		all_events_results.back() = *to_what;
}

void AnalysisManager::setCondition(TCondition* cond)
{	_cond = cond;}
TCondition* AnalysisManager::getCondition(void)
{	return _cond;}
void AnalysisManager::setThreadMutex(TMutex* mutex)
{	_thread_mutex = mutex;}
TMutex* AnalysisManager::getThreadMutex(void)
{	return _thread_mutex;}


ParameterPile::experiment_manifest AnalysisManager::refine_exp_area(ParameterPile::experiment_manifest area)//looks up the existing runs in data directory 
//and intersects them with input area (from ParameterPile::exp_area). This is required in order to split runs between threads equally
//depr: TODO: maybe also move to the AnalysisManager
{
	ParameterPile::experiment_manifest out_area = area;
	out_area.runs.erase();
	std::vector<int> runs;
	int from = -1, to = -1;
#if defined(__WIN32__)
	HANDLE dir;
	WIN32_FIND_DATA file_data;
	std::string path = area.in_folder;
	if ((dir = FindFirstFile((path + "*").c_str(), &file_data)) == INVALID_HANDLE_VALUE)
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
	std::string path = area.in_folder;
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
