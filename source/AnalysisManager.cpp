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
//void AnalysisManager::nextExperiment(void)
//{
//
//}
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
	processAllRuns();
}

void AnalysisManager::processAllRuns(void)
{
	all_runs_results.push_back(AllRunsResults(one_run_results.back().curr_area));
	all_runs_results.back().processAllRuns(one_run_results);

//	int runs_processed = 0;
//	int runs_no_PMT = 0; //false triggers
//	TH1D *hist_S = new TH1D("PMT_S_peaks", "PMT_S_peaks", 60, 0, 2);
//	TH1I *hist_n = new TH1I("PMT_N_peaks", "PMT_N_peaks", 30, 0, 30);
//	DVECTOR ns, Ss;
//	std::vector<double> xs_GEM, ys_GEM;
//	ParameterPile::experiment_area area_;
//	while ((curr_run != NextRunIs::Null) && (curr_run != NextRunIs::NewExperiment)) {
//		processOneRun();
//		if ((one_run_data.back().isValid()) || (SingleRunData::Status::NoPMTsignal == one_run_data.back().getStatus())){		
//			hist_S->Fill(one_run_data.back().PMT3_summed_peaks_area);
//			hist_n->Fill(one_run_data.back().PMT3_n_peaks);
//			Ss.push_back(one_run_data.back().PMT3_summed_peaks_area);
//			ns.push_back(one_run_data.back().PMT3_n_peaks);
//		}
//		if (one_run_data.back().isValid()){
//			std::cout << current_under_processing.experiments.back() << "_run_" << current_under_processing.runs.back() << "_sub_"
//				<< current_under_processing.sub_runs.back() << "_processed" << std::endl;
//			if (xs_GEM.empty()){
//				xs_GEM = one_run_data.back().xs_GEM;
//				ys_GEM = one_run_data.back().ys_GEM;
//			} else {
//				for (auto i = ys_GEM.begin(), j = one_run_data.back().ys_GEM.begin();
//					(i != ys_GEM.end() && j != one_run_data.back().ys_GEM.end()); ++i, ++j)
//					*i += *j;
//			}
//			runs_processed++;
//		} else {
//			if (SingleRunData::Status::NoPMTsignal == one_run_data.back().getStatus())
//				runs_no_PMT++;
//			one_run_data.pop_back();
//		}
//		area_.experiments = current_under_processing.experiments;
//		nextRun();
//	}
//
//	TCanvas *c1 = new TCanvas(("S_peaks_distribution " + area_.experiments.back()).c_str(),
//		("S_peaks_distribution " + area_.experiments.back()).c_str());
//	c1->cd();
//	hist_S->Draw();
//	c1->Update();
//	TCanvas *c2 = new TCanvas(("n_peaks_distribution " + area_.experiments.back()).c_str(),
//		("n_peaks_distribution " + area_.experiments.back()).c_str());
//	c2->cd();
//	hist_n->Draw();
//	c2->Update();
//
//	TCanvas *c3 = new TCanvas(("N-S Scatter " + area_.experiments.back()).c_str(), ("N-S Scatter " + area_.experiments.back()).c_str());
//	double *NNs = new double[ns.size()];
//	double *SSs = new double[Ss.size()];
//	for (int _n = 0, _s = 0; (_n < ns.size()) && (_s < Ss.size()); ++_s, ++_n){
//		NNs[_n] = ns[_n];
//		SSs[_s] = Ss[_s];
//	}
//	TGraph *gr = new TGraph(std::min(ns.size(),Ss.size()),SSs,NNs);
//	c3->cd();
//	gr->Draw("ap");
//	c3->Update();
//	delete[] NNs;
//	delete[] SSs;
//
//	for (auto i = ys_GEM.begin(); i != ys_GEM.end(); ++i)
//		*i /= runs_processed;
//	area_.runs.push_back(ParameterPile::areas_to_draw.back().runs.back());
//	area_.sub_runs.push_back(ParameterPile::areas_to_draw.back().sub_runs.back());
//	area_.channels.push_back(ParameterPile::areas_to_draw.back().channels.back());
//	if (ParameterPile::draw_required(area_)&&runs_processed){
//		//xs_GEM.insert(xs_GEM.begin(), 0);
//		//ys_GEM.insert(ys_GEM.begin(), 0);
//		GraphicOutputManager man;
//		std::string exp = area_.experiments.back();
//		bool invalid = false;
//		for (auto i = exp.begin(); i != exp.end(); invalid ? i =exp.begin():++i){
//			invalid = false;
//			if (*i == '_'){
//				if (i != exp.begin())
//					if (*(i - 1) == '\\')
//						continue;
//				exp.insert(i, '\\');
//				invalid = true;
//			}
//		}
//		Drawing *dr= man.GetDrawing("GEM_" + area_.experiments.back(), 0,ParameterPile::DrawEngine::Gnuplot);
//		dr->AddToDraw(xs_GEM, ys_GEM, "GEM\\\_" + exp,"", 0);
//		std::vector<double> GEM_int;
//		SignalOperations::integrate(xs_GEM, ys_GEM, GEM_int);
//		dr->AddToDraw(xs_GEM, GEM_int, "GEM\\\_I\\\_" + exp,"axes x1y2", 0);
//
//		//find start GEM time 
//		std::vector<double>::iterator x_start = xs_GEM.begin();
//		find_GEM_start_time(xs_GEM, ys_GEM, x_start, ParameterPile::GEM_N_of_averaging, man);
//		if (x_start != xs_GEM.end())
//			dr->AddToDraw_vertical(*x_start, -0.04, 0.04, "lc rgb \"#FF0000\"", 0);
//		//find finish GEM time
//		std::vector<double>::iterator x_finish = xs_GEM.begin();
//		double temp_y_max;
//		DVECTOR temp_xs = xs_GEM, temp_ys_I=GEM_int;
//		SignalOperations::apply_time_limits(temp_xs, temp_ys_I, *x_start, xs_GEM.back());
//		SignalOperations::get_max(temp_xs, temp_ys_I, x_finish, temp_y_max, ParameterPile::GEM_N_of_averaging);
//		if (x_finish != temp_xs.end())
//			dr->AddToDraw_vertical(*x_finish, -0.04, 0.04, "lc rgb \"#FF0000\"", 0);
//		std::cout << "Experiment " << area_.experiments.back() << " processed"<<std::endl;
//		std::cout << "# of runs " << runs_processed << std::endl;
//		std::cout << "# of runs with empty PMT signal" << runs_no_PMT << std::endl;
//		if (x_finish != temp_xs.end() && x_start != xs_GEM.end()){
//			x_finish = SignalOperations::find_x_iterator_by_value(xs_GEM.begin(), xs_GEM.end() - 1, *x_finish);
//			double Integ = *(GEM_int.begin() + (x_finish - xs_GEM.begin())) - *(GEM_int.begin() + (x_start - xs_GEM.begin()));
//			std::cout << "GEM integral is " << Integ << std::endl;
//#ifdef _TEMP_CODE
//			std::ofstream GEM_results_file;
//			GEM_results_file.open(std::string(OUTPUT_DIR) + OUTPUT_GEMS, std::ios_base::app);
//			GEM_results_file << area_.experiments.back() << "\t" << Integ << "\t" << *x_start << "\t" << *x_finish << std::endl;
//			GEM_results_file.close();
//#endif
//			std::cout << "GEM timelimits are: [" << *x_start << ";" << *x_finish << "]" << std::endl;
//		} else {
//			std::cout << "GEM time boundaries are invalid" << std::endl;
//		}
//		dr->DrawData();
//	}
//	//TODO: Add averaging over the runs etc.
//	//TODO: Probably will have to use average value for per run processing again.
//	//te best way I guess is to reverse nextRun (add prevRun, use reverse iterators) and 
//	//call processOneRun (averaging results);
//	//the thing is the same may be requred for experiment processing
//	//The other way is do not clean SingleRunData, and simply rerun processOneRun for every member of one_run_data with external pars
//	one_run_data.clear();
}
void AnalysisManager::processAllExperiments(void)
{
	nextRun();
	while (curr_run != NextRunIs::Null){
		curr_run = NextRunIs::NewSubRun;
		loopAllRuns();
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

//double AnalysisManager::minimize_function(double value)
//{
//	if (one_run_data.back().isValid())
//		return one_run_data.back().minimize_function(value);
//	else
//		return 0;
//}