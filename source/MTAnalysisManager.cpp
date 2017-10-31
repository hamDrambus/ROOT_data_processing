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
				current_under_processing.runs.clear();
				current_under_processing.sub_runs.clear();
				current_under_processing.channels.clear();
				curr_run = NextRunIs::Null;
				return;
			}
	//No runs or subruns switch, the MultithreadAnalysisManager is only responsible for experiments switching and runs splitting
}
void MTAnalysisManager::processAllRuns(void)
{
	std::vector<TThread*> pThreads;
	std::vector<AnalysisManager*> _submanagers;
	pThreads.resize(ParameterPile::threads_number, NULL);
	_submanagers.resize(ParameterPile::threads_number, NULL);
	std::vector<TMutex*> mutexes;
	std::vector<TMutex*> thread_mutexes;
	std::vector<TCondition*> conditions;
	mutexes.resize(ParameterPile::threads_number, NULL);
	thread_mutexes.resize(ParameterPile::threads_number, NULL);
	conditions.resize(ParameterPile::threads_number, NULL);
	for (int n =0;n<ParameterPile::threads_number; ++n) {
		mutexes[n] = new TMutex();
		conditions[n] = new TCondition(mutexes[n]);
		thread_mutexes[n] = new TMutex();
	}

	ParameterPile::experiment_area actual_area = refine_exp_area(current_under_processing);
	std::vector<ParameterPile::experiment_area> areas = split_exp_area(actual_area, ParameterPile::threads_number);
	for (int n=0; n<ParameterPile::threads_number;++n) {
		//TODO: mind that analysis manager clears the processed data, fix this
		_submanagers[n] = new AnalysisManager(areas[n]);
		_submanagers[n]->setCondition(conditions[n]);
		_submanagers[n]->setThreadMutex(thread_mutexes[n]);
		pThreads[n] = new TThread(("AnManager_" + current_under_processing.experiments.back() + std::to_string(n)).c_str(),
			&process_runs_in_thread, _submanagers[n]);
		pThreads[n]->Run();
	}
	//TThread::Ps();
	all_runs_results.push_back(AllRunsResults(actual_area));
	for (int n = 0;n<ParameterPile::threads_number; ++n){
		//pThreads[n]->Join(); doesn't work
		//std::cout << "Start of waiting for " << n << std::endl;
		if (0 != thread_mutexes[n]->TryLock()){ //thread is already executed, so no wait required
		} else {
			conditions[n]->Wait();
		}
		thread_mutexes[n]->UnLock();
		//std::cout << "End of waiting for " << n << std::endl;
		std::vector<AllRunsResults> *res = _submanagers[n]->getAllRunsResults();
		if (!res->empty())
			all_runs_results.back().Join(&res->back());
	}
	all_runs_results.back().Joined();
	//start second iteration
	for (int n = 0; n<ParameterPile::threads_number; ++n){
		_submanagers[n]->setAllRunsResults(&all_runs_results.back());
		pThreads[n]->Run();
	}
	all_runs_results.pop_back();
	all_runs_results.push_back(AllRunsResults(actual_area));//clear
	for (int n = 0; n<ParameterPile::threads_number; ++n){
		//pThreads[n]->Join(); doesn't work
		//std::cout << "Start of waiting for " << n << std::endl;
		if (0 != thread_mutexes[n]->TryLock()){ //thread is already executed, so no wait required
		} else {
			conditions[n]->Wait();
		}
		thread_mutexes[n]->UnLock();
		//std::cout << "End of waiting for " << n << std::endl;
		std::vector<AllRunsResults> *res = _submanagers[n]->getAllRunsResults();
		if (!res->empty())
			all_runs_results.back().Join(&res->back());
	}
	all_runs_results.back().Joined();
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

ParameterPile::experiment_area MTAnalysisManager::refine_exp_area(ParameterPile::experiment_area area)//looks up the existing runs in data directory 
//and intersects them with input area (from ParameterPile::exp_area). This is required in order to split runs between threads equally
//TODO: maybe move to the AnalysisManager
{
	ParameterPile::experiment_area out_area = area;
	out_area.runs.clear();
	std::vector<int> runs;
	int from = -1, to = -1;
	HANDLE dir;
	WIN32_FIND_DATA file_data;
	std::string path =  DATA_PREFIX;
	path += "event_x-ray_" + area.experiments.back();
	if ((dir = FindFirstFile((path + "/*").c_str(), &file_data)) == INVALID_HANDLE_VALUE)
		return out_area;
	do {
		std::string file_name = file_data.cFileName;
		if ((file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0)
			continue;
		if (file_name.size() < 4)
			continue;
		if (file_name[0]=='.')
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
			out_area.runs.push_back(from);
			out_area.runs.push_back(to);
			from = run;
			to = from;
		}
	} while (FindNextFile(dir, &file_data));
	if ((from >= 0) && (to >= 0)){
		out_area.runs.push_back(from);
		out_area.runs.push_back(to);
	}
	FindClose(dir);

	//intersect files with input runs:
	bool even_out = true;
	int left_out, right_out;
	bool finished = false;
	for (auto g = out_area.runs.begin(); g != out_area.runs.end(); ++g, even_out = !even_out){
		if (even_out){
			left_out = *g;
		} else {
			right_out = *g;
			bool even = true;
			int left, right;
			for (auto h = area.runs.begin(); h != area.runs.end(); ++h, even = !even){
				if (even){
					left = *h;
				} else {
					right = *h;
					if ((right_out > right) && (left_out <= right)){
						*g = right;
					}
					if ((left_out < left) && (right_out >= left)){
						*(g - 1) = left;
					}
					left_out = *(g - 1);
					right_out = *(g);
				}
			}
		}
	}
	return out_area;
}

std::vector<ParameterPile::experiment_area> MTAnalysisManager::split_exp_area(ParameterPile::experiment_area area_to_split, int N)
{
	std::vector <ParameterPile::experiment_area> out_;
	for (int h = 0; h < N; h++){
		out_.push_back(area_to_split);
		out_[h].runs.erase();
	}
	int N_runs=0;
	bool even = true;
	int l, r;
	for (auto h = area_to_split.runs.begin(); h != area_to_split.runs.end(); ++h,even=!even){
		if (even){
			l = *h;
		} else {
			r = *h;
			N_runs += r - l + 1;
		}
	}
	int N_accum = 0;
	int curr_index = 0;
	even = true;
	for (auto h = area_to_split.runs.begin(); h != area_to_split.runs.end(); ++h, even = !even){
		if (even){
			l = *h;
		} else {
			r = *h;
			N_accum += r - l + 1;
			int new_l = l , new_r = r;
			double delta = N_runs / (double)N;
			while (N_accum >(int)((curr_index + 1)*delta)) {
				if ((N_accum - (int)((curr_index + 1)*delta)) > (int)(delta))
					new_r = new_l + (int)(delta) - 1;
				else
					new_r = N_accum - (int)((curr_index + 1)*delta) + new_l - 1;
				out_[curr_index].runs.push_back(new_l);
				out_[curr_index].runs.push_back(new_r);
				new_l = new_r + 1;
				curr_index++;
			}
			new_r = r;
			if (N_accum <= (int)((curr_index + 1)*delta) && (new_l <= new_r))
			{
				out_[curr_index].runs.push_back(new_l);
				out_[curr_index].runs.push_back(new_r);
			}
		}
	}
	return out_;
}
