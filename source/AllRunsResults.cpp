#include "AllRunsResults.h"

AllRunsResults::AllRunsResults(ParameterPile::experiment_area experiment)
{
	_exp = experiment;
	S_peaks_cutoff = ParameterPile::PMT_SArea_peaks_acceptance;
	S_peaks_max_cutoff = S_peaks_cutoff - 1;
	N_peaks_cutoff = ParameterPile::PMT_N_peaks_acceptance;
	N_of_runs = 0;
	N_of_valid_runs = 0;
	Iteration_N = 0;
#ifdef _USE_TIME_STATISTICS
	time_stat.t_PMT_proc=0;//0st iteration
	time_stat.n_PMT_proc=0;
	time_stat.t_PMT_baseline=0;//0st iteration
	time_stat.n_PMT_baseline=0;
	time_stat.t_PMT_peaks=0;//0st iteration
	time_stat.n_PMT_peaks=0;
	time_stat.t_PMT_file_reading=0;//0st
	time_stat.n_PMT_file_reading=0;
	time_stat.t_PMT_filtering = 0;//0st
	time_stat.n_PMT_filtering = 0;
	//====================================
	time_stat.t_MPPC_proc=0; //all till ==...= is the 1st iteration
	time_stat.n_MPPC_proc=0;
	time_stat.t_MPPC_file_reading=0;
	time_stat.n_MPPC_file_reading=0;
	time_stat.t_MPPC_filtering = 0;
	time_stat.n_MPPC_filtering = 0;
	time_stat.t_MPPC_threshold_and_first_baseline=0;
	time_stat.n_MPPC_threshold_and_first_baseline=0;
	//time_stat.t_MPPC_threshold_and_first_baseline_peaks=0;
	//time_stat.n_MPPC_threshold_and_first_baseline_peaks=0;
	time_stat.t_MPPC_curved_baseline=0;
	time_stat.n_MPPC_curved_baseline=0;
	time_stat.t_MPPC_curved_baseline_v2 = 0;
	time_stat.n_MPPC_curved_baseline_v2 = 0;
	time_stat.t_MPPC_curved_baseline_v3 = 0;
	time_stat.n_MPPC_curved_baseline_v3 = 0;
	time_stat.t_MPPC_curved_baseline_v4 = 0;
	time_stat.n_MPPC_curved_baseline_v4 = 0;
	time_stat.t_MPPC_curved_baseline_v5 = 0;
	time_stat.n_MPPC_curved_baseline_v5 = 0;
	time_stat.t_MPPC_curved_baseline_v6 = 0;
	time_stat.n_MPPC_curved_baseline_v6 = 0;
	time_stat.t_MPPC_curved_baseline_v7 = 0;
	time_stat.n_MPPC_curved_baseline_v7 = 0;
	time_stat.t_MPPC_curved_baseline_v8 = 0;
	time_stat.n_MPPC_curved_baseline_v8 = 0;
	time_stat.t_MPPC_curved_baseline_baseline=0;
	time_stat.n_MPPC_curved_baseline_baseline=0;
	time_stat.t_MPPC_baseline_substraction=0;
	time_stat.n_MPPC_baseline_substraction=0;
	time_stat.t_MPPC_peaks_finding=0;
	time_stat.n_MPPC_peaks_finding=0;
	time_stat.t_MPPC_peaks_processing=0;
	time_stat.n_MPPC_peaks_processing=0;
	time_stat.t_MPPC_double_I = 0;
	time_stat.n_MPPC_double_I = 0;
	//====================================
	time_stat.t_RUN_proc=0;  //both iterations, must be set at merge() proc
	time_stat.n_RUN_proc=0;
	time_stat.t_RUN_proc_single_iter = 0;//both iterations, used at merge() proc, set in analysis Manager, must be cleared
	time_stat.n_RUN_proc_single_iter = 0;
#endif
}

void AllRunsResults::processAllRuns(STD_CONT<SingleRunResults> &single_results)
{
	//this->Clear() called before in the AnalysisManager. 
	int n_valid_runs = 0;
	int __N_of_runs = single_results.size();
	N_of_runs = single_results.size();
	N_of_valid_runs = 0;
	for (auto i = single_results.begin(); i != single_results.end(); ++i){
		if (i->isValid()) { //may be false only at the secondary function call (when test_PMT cut is applied)
			++n_valid_runs;
			++N_of_valid_runs;
			_Ss.push_back(i->PMT3_summed_peaks_area);
			_ns.push_back(i->PMT3_n_peaks);
			PMT3_peaks.push_back(i->PMT3_peaks);
			PMT1_peaks.push_back(i->PMT1_peaks);
			if (_xs_PMT3_sum.empty()){
				_xs_PMT3_sum = i->xs_PMT3;
				_ys_PMT3_sum = i->ys_PMT3;
			} else {
				for (auto ys1 = _ys_PMT3_sum.begin(), ys2 = i->ys_PMT3.begin();
					(ys1 != _ys_PMT3_sum.end() && ys2 != i->ys_PMT3.end()); ++ys1, ++ys2)
						*ys1 += *ys2;
			}
			if (_xs_PMT1_sum.empty()){
				_xs_PMT1_sum = i->xs_PMT1;
				_ys_PMT1_sum = i->ys_PMT1;
			} else {
				for (auto ys1 = _ys_PMT1_sum.begin(), ys2 = i->ys_PMT1.begin();
					(ys1 != _ys_PMT1_sum.end() && ys2 != i->ys_PMT1.end()); ++ys1, ++ys2)
						*ys1 += *ys2;
			}
#ifdef _PROCESS_GEMS
			if (_xs_GEM_sum.empty()){
				_xs_GEM_sum = (*i).xs_GEM;
				_ys_GEM_sum = (*i).ys_GEM;
			} else {
				for (auto ys1 = _ys_GEM_sum.begin(), ys2 = (*i).ys_GEM.begin();
					(ys1 != _ys_GEM_sum.end() && ys2 != (*i).ys_GEM.end()); ++ys1, ++ys2)
					*ys1 += *ys2;
			}
#endif
			bool first = false, valid = true;
			if (mppc_peaks_in_S2_area.empty())
				first = true;
			if ((mppc_peaks_in_S2_area.size() != mppc_S2_start_time.size()) || (mppc_peaks_in_S2_area.size() != mppc_S2_finish_time.size())
				/*|| (mppc_peaks_in_S2_area.size() != mppc_all_peaks_Ss.size())*/ || (mppc_peaks_in_S2_area.size()!=mppc_channels.size())
				|| (mppc_peaks.size() != mppc_peaks_in_S2_area.size()))
				valid = false;
			if ((i->mppc_S2_peaks_area.size() != i->mppc_S2_start_t.size()) || (i->mppc_S2_peaks_area.size() != i->mppc_S2_finish_t.size())
				|| (i->mppc_S2_peaks_area.size() != i->mppc_peaks.size()) || (i->mppc_S2_peaks_area.size() != i->mppc_channels.size())
				|| (i->mppc_S2_peaks_area.size() != i->mppc_double_I.size()))
				valid = false;
			if (!valid) {
				std::cout << "SingleRunResulst MPPC channels size mismatch" << std::endl;
				continue; 
			}
			int sz = i->mppc_channels.size();
			if (first){
				mppc_peaks_in_S2_area.resize(sz, DVECTOR());
				mppc_S2_start_time.resize(sz, DVECTOR());
				mppc_S2_finish_time.resize(sz, DVECTOR());
				//mppc_all_peaks_Ss.resize(sz, DVECTOR());
				mppc_double_Is.resize(sz, DVECTOR());
				mppc_peaks.resize(sz, STD_CONT<STD_CONT<peak>>());
#ifndef _USE_DEQUE
				for (int ch = 0; ch != sz; ++ch) {
					mppc_peaks_in_S2_area[ch].reserve(__N_of_runs);
					mppc_S2_start_time[ch].reserve(__N_of_runs);
					mppc_S2_finish_time[ch].reserve(__N_of_runs);
					mppc_double_Is[ch].reserve(__N_of_runs);
				}		
#endif
				mppc_channels = i->mppc_channels;
			}
			for (int ch = 0; ch != sz; ++ch) {
				mppc_peaks_in_S2_area[ch].push_back(i->mppc_S2_peaks_area[ch]);
				mppc_S2_start_time[ch].push_back(i->mppc_S2_start_t[ch]);
				mppc_S2_finish_time[ch].push_back(i->mppc_S2_finish_t[ch]);
				mppc_double_Is[ch].push_back(i->mppc_double_I[ch]);
				mppc_peaks[ch].push_back(i->mppc_peaks[ch]);
				/*for (auto peaks = i->mppc_peaks[ch].begin(); peaks != i->mppc_peaks[ch].end(); ++peaks){
					mppc_all_peaks_Ss[ch].push_back(peaks->S);
				}*/
			}
		}
#ifdef _HOTFIX_CLEAR_MEMORY
		STD_CONT<DVECTOR>().swap(i->mppc_baseline_xs);
		STD_CONT<DVECTOR>().swap(i->mppc_baseline_ys);
		STD_CONT<int>().swap(i->mppc_channels);
		STD_CONT<STD_CONT<peak>>().swap(i->mppc_peaks);
		STD_CONT<peak>().swap(i->PMT3_peaks);
		STD_CONT<peak>().swap(i->PMT1_peaks);
		DVECTOR().swap(i->mppc_S2_finish_t);
		DVECTOR().swap(i->mppc_S2_peaks_area);
		DVECTOR().swap(i->mppc_S2_start_t);
		//DVECTOR().swap(i->mppc_sum_peaks_area);
		DVECTOR().swap(i->mppc_double_I);
		DVECTOR().swap(i->xs_GEM);
		DVECTOR().swap(i->ys_GEM);
#else
		i->mppc_baseline_xs.clear();
		i->mppc_baseline_ys.clear();
		i->mppc_channels.clear();
		i->mppc_peaks.clear();
		i->mppc_S2_finish_t.clear();
		i->mppc_S2_peaks_area.clear();
		i->mppc_S2_start_t.clear();
		//i->mppc_sum_peaks_area.clear();
		i->mppc_double_I.clear();
		i->xs_GEM.clear();
		i->ys_GEM.clear();
#endif
	}
	
}

void AllRunsResults::find_GEM_start_time(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start, int N_trust, GraphicOutputManager &man)
{
	DVECTOR xs_before_S1 = xs, ys_before_S1 = ys;
	SignalOperations::apply_time_limits(xs_before_S1, ys_before_S1, *(xs.begin()), ParameterPile::S1_start_time);
	DITERATOR x_befS1_max;
	double noize_amp;
	SignalOperations::get_max(xs_before_S1, ys_before_S1, x_befS1_max, noize_amp, N_trust);
	Drawing *dr = man.GetDrawing(0);
	if (dr)
		dr->AddToDraw_vertical(*x_befS1_max, -0.4, 0.4, "lc rgb \"#000000\"", 0);
	noize_amp *= ParameterPile::GEM_threshold_to_noise;
	if (dr)
		dr->AddToDraw_baseline(noize_amp, "threshold", "lc rgb \"#000000\"", 0);
	DITERATOR x_S1_left = xs.begin();
	DITERATOR x_S1_right = xs.begin();
	SignalOperations::find_next_peak(xs, ys, x_S1_left, x_S1_right, noize_amp, N_trust);
	if (x_S1_left == xs.end()){
		std::cout << "1st peak not found, threshold " << noize_amp << std::endl;
		x_start = xs.end();
		return;
	}
	if (dr)
		dr->AddToDraw_vertical(*x_S1_right, -0.4, 0.4, "lc rgb \"#00AA00\"", 0);
	//x_S1_right += N_trust;
	SignalOperations::find_next_extremum(xs, ys, x_S1_right, N_trust);
	if (dr)
		dr->AddToDraw_vertical(*x_S1_right, -0.4, 0.4, "lc rgb \"#0000AA\"", 0);
	x_S1_right += N_trust;
	SignalOperations::find_next_extremum(xs, ys, x_S1_right, N_trust);
	if (dr)
		dr->AddToDraw_vertical(*x_S1_right, -0.4, 0.4, "lc rgb \"#00AAAA\"", 0);
	x_S1_right += N_trust;
	SignalOperations::find_next_extremum(xs, ys, x_S1_right, N_trust);
	if (x_S1_right == xs.end()){
		x_start = xs.end();
		return;
	}
	if (*(ys.begin() + (x_S1_right - xs.begin())) > 0){
		x_start = x_S1_right;

		return;
	}
	x_S1_left = x_S1_right;
	SignalOperations::find_next_peak(xs, ys, x_S1_left, x_S1_right, 0, N_trust); //effectively finds intersection with 0
	x_start = x_S1_left;
}

void AllRunsResults::find_S_cutoff(void)//TODO: explain the algorithm and mb test it
//"mb" because the actual test is S histograms pictures with cutoffs. If pictures are ok, there is no need to
//test the algorithm excessively
{
	if (_Ss.size() < ParameterPile::PMT_min_statistics) //not statistically significant. TODO: figure out the number,
		//and ParameterPile, however this is not important
		return;
	DVECTOR areas = _Ss;
	std::sort(areas.begin(), areas.end());
	double G_mean = TMath::Mean(areas.begin(), areas.end());
	double G_rms = TMath::RMS(areas.begin(), areas.end());
	double G_min = *(areas.begin());
	double G_max = areas.back();
	//max_cutoff:
	double max_cutoff = G_mean + ParameterPile::PMT_right_cutoff_from_RMS * G_rms;
	S_peaks_max_cutoff = max_cutoff;
	DITERATOR erase_from;
	for (erase_from = areas.begin(); erase_from != areas.end(); ++erase_from)
		if ((*erase_from) > max_cutoff)
			break;
	areas.erase(erase_from, areas.end());
	G_max = areas.back();
	G_mean = TMath::Mean(areas.begin(), areas.end());
	G_rms = TMath::RMS(areas.begin(), areas.end());
	//applied max cutoff

	int N_min_above_cutoff = areas.size()*ParameterPile::PMT_min_fraction_above_cutoff;
	int N_above_global_mean_acceptable = areas.size()*ParameterPile::PMT_mean_above_cutoff_acceptance; //see usage below for meaning
	DVECTOR areas_full = areas;
	
	double max_ds = 0;
	DITERATOR i_max_ds = areas.begin();
	bool found_max_delta_S = false;
	while (!found_max_delta_S){
		for (auto i = areas.begin(); i != (areas.end() - 1); ++i) {
			if ((*(i + 1) - *i) > max_ds){
				max_ds = (*(i + 1) - *i);
				i_max_ds = i + 1;
			}
		}
		if ((areas.end() - i_max_ds) < N_min_above_cutoff) {//it is not the max of interest
			areas.erase(i_max_ds, areas.end());
			max_ds = 0;
			i_max_ds = areas.begin();
		} else {
			found_max_delta_S = true;
		}
		if (areas.size() < N_min_above_cutoff){
			i_max_ds = areas.begin();
			break;
		}
	}
	double mean = TMath::Mean(areas_full.begin() + (i_max_ds - areas.begin()), areas_full.end());
	double rms = TMath::RMS(areas_full.begin() + (i_max_ds - areas.begin()), areas_full.end());
	double cutoff = mean - rms * ParameterPile::PMT_left_cutoff_from_RMS;
	DITERATOR _last = areas.end()-1;
	DITERATOR _begin = areas.begin();
	DITERATOR i_mean = SignalOperations::find_x_iterator_by_value(_begin, _last, G_mean);
	DITERATOR i_cutoff = SignalOperations::find_x_iterator_by_value(_begin, _last, cutoff);
	if (G_mean < cutoff)
		if ((i_cutoff - i_mean) >= N_above_global_mean_acceptable)
			cutoff = G_mean;//or set to S_peaks_cutoff ?
	S_peaks_cutoff = std::max(cutoff, S_peaks_cutoff);
#ifdef _NO_PMT_SELECTION
	S_peaks_cutoff = ParameterPile::PMT_SArea_peaks_acceptance;
	S_peaks_max_cutoff = S_peaks_cutoff - 1;
#endif
}

//void AllRunsResults::find_S_cutoff_v2(void)
//{
//	if (_Ss.size() < 10) //not statistically significant. TODO: figure out the number, and ParameterPile, however this is not important
//		return;
//	DVECTOR areas = _Ss;
//	std::sort(areas.begin(), areas.end());
//	DVECTOR ys;
//	ys.resize(areas.size(), 1);
//	DVECTOR x_spreaded, y_spreaded;
//	SignalOperations::spread_peaks(areas, ys, x_spreaded, y_spreaded);
//	Drawing *dr = graph_manager.GetDrawing(_exp.experiments.back() + " S_spreaded", 1, ParameterPile::DrawEngine::Gnuplot);
//	dr->AddToDraw(x_spreaded, y_spreaded, _exp.experiments.back() + " S_spreaded", "w lines", 0);
//	//dr->DrawData();
//}

void AllRunsResults::Merge(AllRunsResults* with)
{
	if (0 == Iteration_N) {
		_Ss.insert(_Ss.end(), with->_Ss.begin(), with->_Ss.end());
		_ns.insert(_ns.end(), with->_ns.begin(), with->_ns.end());
		PMT1_peaks.insert(PMT1_peaks.end(), with->PMT1_peaks.begin(), with->PMT1_peaks.end());
		PMT3_peaks.insert(PMT3_peaks.end(), with->PMT3_peaks.begin(), with->PMT3_peaks.end());
	}
	N_peaks_cutoff = with->N_peaks_cutoff;
	S_peaks_cutoff = with->S_peaks_cutoff;
	N_of_runs += with->N_of_runs;
	N_of_valid_runs += with->N_of_valid_runs;
	//Iteration_N = with->Iteration_N;
	if (1 == Iteration_N) {
		if (_xs_GEM_sum.empty()) {
			_xs_GEM_sum = with->_xs_GEM_sum;
			_ys_GEM_sum = with->_ys_GEM_sum;
		} else {
			for (auto i = _ys_GEM_sum.begin(), j = with->_ys_GEM_sum.begin(); (i != _ys_GEM_sum.end()) && (j != with->_ys_GEM_sum.end()); ++j, ++i)
				*i += *j;
		}
		if (_xs_PMT3_sum.empty()) {
			_xs_PMT3_sum = with->_xs_PMT3_sum;
			_ys_PMT3_sum = with->_ys_PMT3_sum;
		} else {
			for (auto i = _ys_PMT3_sum.begin(), j = with->_ys_PMT3_sum.begin(); (i != _ys_PMT3_sum.end()) && (j != with->_ys_PMT3_sum.end()); ++j, ++i)
				*i += *j;
		}
		if (_xs_PMT1_sum.empty()) {
			_xs_PMT1_sum = with->_xs_PMT1_sum;
			_ys_PMT1_sum = with->_ys_PMT1_sum;
		} else {
			for (auto i = _ys_PMT1_sum.begin(), j = with->_ys_PMT1_sum.begin(); (i != _ys_PMT1_sum.end()) && (j != with->_ys_PMT1_sum.end()); ++j, ++i)
				*i += *j;
		}


		bool empty = false, valid = true;
		if (mppc_peaks_in_S2_area.empty())
			empty = true;
		if ((mppc_peaks_in_S2_area.size() != mppc_S2_start_time.size()) || (mppc_peaks_in_S2_area.size() != mppc_S2_finish_time.size())
			/*|| (mppc_peaks_in_S2_area.size() != mppc_all_peaks_Ss.size())*/ || (mppc_peaks_in_S2_area.size() != mppc_channels.size())
			|| (mppc_peaks_in_S2_area.size() != mppc_double_Is.size()) || (mppc_peaks_in_S2_area.size() != mppc_peaks.size()))
			valid = false;
		if ((with->mppc_peaks_in_S2_area.size() != with->mppc_S2_start_time.size()) || (with->mppc_peaks_in_S2_area.size() != with->mppc_S2_finish_time.size())
			/*|| (with->mppc_peaks_in_S2_area.size() != with->mppc_all_peaks_Ss.size())*/
			|| (with->mppc_peaks_in_S2_area.size() != with->mppc_channels.size())
			|| (with->mppc_peaks_in_S2_area.size() != with->mppc_double_Is.size()) || (with->mppc_peaks_in_S2_area.size() != with->mppc_peaks.size()))
			valid = false;
		if (!valid) {
			std::cout << "AllRunsResults contains invalid MPPC data: channels' size mismatches" << std::endl;
			return;
		}
		if (!empty && (mppc_peaks_in_S2_area.size() != with->mppc_peaks_in_S2_area.size())){
			std::cout << "WARNING Two AllRunsResults have MPPC channels' size mismatches: not Merging" << std::endl;
			return;
		}
		if (empty) {
			mppc_peaks_in_S2_area = with->mppc_peaks_in_S2_area;
			mppc_S2_start_time = with->mppc_S2_start_time;
			mppc_S2_finish_time = with->mppc_S2_finish_time;
			//mppc_all_peaks_Ss = with->mppc_all_peaks_Ss;
			mppc_channels = with->mppc_channels;
			mppc_double_Is = with->mppc_double_Is;
			mppc_peaks = with->mppc_peaks;
		} else {
			for (int ch = 0; ch < mppc_peaks_in_S2_area.size(); ++ch) {
				mppc_peaks_in_S2_area[ch].insert(mppc_peaks_in_S2_area[ch].end(), with->mppc_peaks_in_S2_area[ch].begin(), with->mppc_peaks_in_S2_area[ch].end());
				mppc_S2_start_time[ch].insert(mppc_S2_start_time[ch].end(), with->mppc_S2_start_time[ch].begin(), with->mppc_S2_start_time[ch].end());
				mppc_S2_finish_time[ch].insert(mppc_S2_finish_time[ch].end(), with->mppc_S2_finish_time[ch].begin(), with->mppc_S2_finish_time[ch].end());
				//mppc_all_peaks_Ss[ch].insert(mppc_all_peaks_Ss[ch].end(), with->mppc_all_peaks_Ss[ch].begin(), with->mppc_all_peaks_Ss[ch].end());
				mppc_peaks[ch].insert(mppc_peaks[ch].end(), with->mppc_peaks[ch].begin(), with->mppc_peaks[ch].end());
				mppc_double_Is[ch].insert(mppc_double_Is[ch].end(), with->mppc_double_Is[ch].begin(), with->mppc_double_Is[ch].end());
			}
		}
	}
#ifdef _USE_TIME_STATISTICS
	time_stat.t_RUN_proc_single_iter += with->time_stat.t_RUN_proc_single_iter;
	time_stat.n_RUN_proc_single_iter += with->time_stat.n_RUN_proc_single_iter;
	
	if (0 == Iteration_N) {
		time_stat.t_PMT_proc += with->time_stat.t_PMT_proc;//0st iteration
		time_stat.n_PMT_proc += with->time_stat.n_PMT_proc;
		time_stat.t_PMT_baseline += with->time_stat.t_PMT_baseline;//0st iteration
		time_stat.n_PMT_baseline += with->time_stat.n_PMT_baseline;
		time_stat.t_PMT_peaks += with->time_stat.t_PMT_peaks;//0st iteration
		time_stat.n_PMT_peaks += with->time_stat.n_PMT_peaks;
		time_stat.t_PMT_file_reading += with->time_stat.t_PMT_file_reading;//0st
		time_stat.n_PMT_file_reading += with->time_stat.n_PMT_file_reading;
		time_stat.t_PMT_filtering += with->time_stat.t_PMT_filtering;//0st
		time_stat.n_PMT_filtering += with->time_stat.n_PMT_filtering;
	}
	//====================================
	if (1 == Iteration_N) {
		time_stat.t_MPPC_proc += with->time_stat.t_MPPC_proc; //all till ==...= is the 1st iteration
		time_stat.n_MPPC_proc += with->time_stat.n_MPPC_proc;
		time_stat.t_MPPC_file_reading += with->time_stat.t_MPPC_file_reading;
		time_stat.n_MPPC_file_reading += with->time_stat.n_MPPC_file_reading;
		time_stat.t_MPPC_filtering += with->time_stat.t_MPPC_filtering;
		time_stat.n_MPPC_filtering += with->time_stat.n_MPPC_filtering;
		time_stat.t_MPPC_threshold_and_first_baseline += with->time_stat.t_MPPC_threshold_and_first_baseline;
		time_stat.n_MPPC_threshold_and_first_baseline += with->time_stat.n_MPPC_threshold_and_first_baseline;
		//time_stat.t_MPPC_threshold_and_first_baseline_peaks += with->time_stat.t_MPPC_threshold_and_first_baseline_peaks;
		//time_stat.n_MPPC_threshold_and_first_baseline_peaks += with->time_stat.n_MPPC_threshold_and_first_baseline_peaks;
		time_stat.t_MPPC_curved_baseline += with->time_stat.t_MPPC_curved_baseline;
		time_stat.n_MPPC_curved_baseline += with->time_stat.n_MPPC_curved_baseline;
		time_stat.t_MPPC_curved_baseline_v2 += with->time_stat.t_MPPC_curved_baseline_v2;
		time_stat.n_MPPC_curved_baseline_v2 += with->time_stat.n_MPPC_curved_baseline_v2;
		time_stat.t_MPPC_curved_baseline_v3 += with->time_stat.t_MPPC_curved_baseline_v3;
		time_stat.n_MPPC_curved_baseline_v3 += with->time_stat.n_MPPC_curved_baseline_v3;
		time_stat.t_MPPC_curved_baseline_v4 += with->time_stat.t_MPPC_curved_baseline_v4;
		time_stat.n_MPPC_curved_baseline_v4 += with->time_stat.n_MPPC_curved_baseline_v4;
		time_stat.t_MPPC_curved_baseline_v5 += with->time_stat.t_MPPC_curved_baseline_v5;
		time_stat.n_MPPC_curved_baseline_v5 += with->time_stat.n_MPPC_curved_baseline_v5;
		time_stat.t_MPPC_curved_baseline_v6 += with->time_stat.t_MPPC_curved_baseline_v6;
		time_stat.n_MPPC_curved_baseline_v6 += with->time_stat.n_MPPC_curved_baseline_v6;
		time_stat.t_MPPC_curved_baseline_v7 += with->time_stat.t_MPPC_curved_baseline_v7;
		time_stat.n_MPPC_curved_baseline_v7 += with->time_stat.n_MPPC_curved_baseline_v7;
		time_stat.t_MPPC_curved_baseline_v8 += with->time_stat.t_MPPC_curved_baseline_v8;
		time_stat.n_MPPC_curved_baseline_v8 += with->time_stat.n_MPPC_curved_baseline_v8;
		time_stat.t_MPPC_curved_baseline_baseline += with->time_stat.t_MPPC_curved_baseline_baseline;
		time_stat.n_MPPC_curved_baseline_baseline += with->time_stat.n_MPPC_curved_baseline_baseline;
		time_stat.t_MPPC_baseline_substraction += with->time_stat.t_MPPC_baseline_substraction;
		time_stat.n_MPPC_baseline_substraction += with->time_stat.n_MPPC_baseline_substraction;
		time_stat.t_MPPC_peaks_finding += with->time_stat.t_MPPC_peaks_finding;
		time_stat.n_MPPC_peaks_finding += with->time_stat.n_MPPC_peaks_finding;
		time_stat.t_MPPC_peaks_processing += with->time_stat.t_MPPC_peaks_processing;
		time_stat.n_MPPC_peaks_processing += with->time_stat.n_MPPC_peaks_processing;
		time_stat.t_MPPC_double_I += with->time_stat.t_MPPC_double_I;
		time_stat.n_MPPC_double_I += with->time_stat.n_MPPC_double_I;
	}
#endif
	with->Clear();
}

void AllRunsResults::Merged(void)
{
	std::cout << "Iteration " << Iteration_N << std::endl;
	std::cout << "N of runs " << N_of_runs << std::endl;
	std::cout << "N of valid runs " << N_of_valid_runs << std::endl;
	if (!_Ss.empty() && (0==Iteration_N)) {//_xs_GEM is empty at the first run
		find_S_cutoff();

		std::string PMT_output_prefix = std::string(ParameterPile::this_path) + "/" + std::string(OUTPUT_DIR) + OUTPUT_PMTS + _exp.experiments.back()
			+ "/" + "PMT_0" + "/" + "PMT_0_";
		std::string PMT1_output_prefix = std::string(ParameterPile::this_path) + "/" + std::string(OUTPUT_DIR) + OUTPUT_PMTS + _exp.experiments.back()
			+ "/" + "PMT_1" + "/" + "PMT_1_";
		vector_to_file(_Ss, PMT_output_prefix + "S2s.dat");
		if (ParameterPile::exp_area.channels.contains(0))
			vector_to_file(PMT3_peaks, PMT_output_prefix + "peaks.dat");
		if (ParameterPile::exp_area.channels.contains(1))
			vector_to_file(PMT1_peaks, PMT1_output_prefix + "peaks.dat");

		if (ParameterPile::exp_area.channels.contains(0)){
			DITERATOR S2_max = std::max_element(_Ss.begin(), _Ss.end());
			double hist_S2_max = (S_peaks_max_cutoff > S_peaks_cutoff) ? std::min(S_peaks_max_cutoff, *S2_max) : *S2_max;
			TH1D *hist_S2 = new TH1D("3PMT_S2_peaks", "3PMT_S2_peaks", 60, 0, hist_S2_max);

			for (int _n = 0, _s = 0; (_n < _ns.size()) && (_s < _Ss.size()); ++_s, ++_n){
				if ((S_peaks_max_cutoff > S_peaks_cutoff) && (_Ss[_s] > S_peaks_max_cutoff))
					continue; //do not draw if maximum limit is imposed && element exceeds it
				hist_S2->Fill(_Ss[_s]);
			}

			TCanvas *c1 = new TCanvas(("3PMT_S2_peaks_distribution " + _exp.experiments.back()).c_str(),
				("3PMT_S2_peaks_distribution " + _exp.experiments.back()).c_str());
			c1->cd();
			hist_S2->Draw();
			TLine *cutoff = new TLine(S_peaks_cutoff, c1->GetUymin(), S_peaks_cutoff, c1->GetUymax());
			cutoff->SetLineColor(kRed);
			cutoff->Draw();
			c1->Update();
			c1->SaveAs((PMT_output_prefix+"S2s.png").c_str(),"png");

			double S_tot_max =0;
			for (auto i= PMT3_peaks.begin(),_end_=PMT3_peaks.end();i!=_end_;++i){
				double val=0;
				for (auto j= i->begin(),_end1_=i->end();j!=_end1_;++j)
					val+= j->S > 0 ? j->S : 0;
				S_tot_max = std::max(S_tot_max,val);
			}
			TH1D *hist_S_tot = new TH1D("3PMT_S_tot_peaks", "3PMT_S_tot_peaks", 60, 0, S_tot_max);

			for (auto i= PMT3_peaks.begin(),_end_=PMT3_peaks.end();i!=_end_;++i){
				double val=0;
				for (auto j= i->begin(),_end1_=i->end();j!=_end1_;++j)
					val+= j->S > 0 ? j->S : 0;
				hist_S_tot->Fill(val);
			}
			TCanvas *c2 = new TCanvas(("3PMT_S_total_peaks_distribution " + _exp.experiments.back()).c_str(),
				("3PMT_S_total_peaks_distribution " + _exp.experiments.back()).c_str());
			c2->cd();
			hist_S_tot->Draw();
			c2->Update();
			c2->SaveAs((PMT_output_prefix+"S_tot.png").c_str(),"png");
		}
		if (ParameterPile::exp_area.channels.contains(1)){
			double S2_max =0;
			for (auto i= PMT1_peaks.begin(),_end_=PMT1_peaks.end();i!=_end_;++i){
				double val=0;
				for (auto j= i->begin(),_end1_=i->end();j!=_end1_;++j){
					if ((j->left>ParameterPile::S2_start_time.find(_exp.experiments.back())->second)&&(j->right<ParameterPile::S2_finish_time.find(_exp.experiments.back())->second))
						val+= j->S > 0 ? j->S : 0;
				}
				S2_max = std::max(S2_max,val);
			}
			TH1D *hist_S2 = new TH1D("PMT#1_S_tot_peaks", "PMT#1_S2_tot_peaks", 60, 0, S2_max);
			for (auto i= PMT1_peaks.begin(),_end_=PMT1_peaks.end();i!=_end_;++i){
				double val=0;
				for (auto j= i->begin(),_end1_=i->end();j!=_end1_;++j){
					if ((j->left>ParameterPile::S2_start_time.find(_exp.experiments.back())->second)&&(j->right<ParameterPile::S2_finish_time.find(_exp.experiments.back())->second))
						val+= j->S > 0 ? j->S : 0;
				}
				hist_S2->Fill(val);
			}
			TCanvas *c1 = new TCanvas(("PMT#1_S2_peaks_distribution " + _exp.experiments.back()).c_str(),
				("PMT#1_S2_peaks_distribution " + _exp.experiments.back()).c_str());
			c1->cd();
			hist_S2->Draw();
			c1->Update();
			c1->SaveAs((PMT_output_prefix+"S2.png").c_str(),"png");

			double S_tot_max =0;
			for (auto i= PMT1_peaks.begin(),_end_=PMT1_peaks.end();i!=_end_;++i){
				double val=0;
				for (auto j= i->begin(),_end1_=i->end();j!=_end1_;++j)
					val+= j->S > 0 ? j->S : 0;
				S_tot_max = std::max(S_tot_max,val);
			}
			TH1D *hist_S_tot = new TH1D("PMT#1_S_tot_peaks", "PMT#1_S_tot_peaks", 60, 0, S_tot_max);

			for (auto i= PMT1_peaks.begin(),_end_=PMT1_peaks.end();i!=_end_;++i){
				double val=0;
				for (auto j= i->begin(),_end1_=i->end();j!=_end1_;++j)
					val+= j->S > 0 ? j->S : 0;
				hist_S_tot->Fill(val);
			}
			TCanvas *c2 = new TCanvas(("PMT#1_S_total_peaks_distribution " + _exp.experiments.back()).c_str(),
				("PMT#1_S_total_peaks_distribution " + _exp.experiments.back()).c_str());
			c2->cd();
			hist_S_tot->Draw();
			c2->Update();
			c2->SaveAs((PMT_output_prefix+"S_tot.png").c_str(),"png");
		}
	}
#ifdef _PROCESS_GEMS
	if (!_xs_GEM_sum.empty() && (0 != N_of_valid_runs)&&(1==Iteration_N)){
		//x_GEM is not empty only at the second processing (when proper cutoffs were applied)
		for (auto i = _ys_GEM_sum.begin(), _end_ = _ys_GEM_sum.end(); i != _end_; ++i)
			*i /= N_of_valid_runs;
#ifdef GEM_V2_
		STD_CONT<peak> no_peaks;
		no_peaks.push_back(peak());
		no_peaks.back().left = ParameterPile::S1_start_time;
		no_peaks.back().right = _xs_GEM_sum.back();
		double gem_baseline =
			SignalOperations::find_baseline_by_median(0, _xs_GEM_sum, _ys_GEM_sum, no_peaks);
		SignalOperations::substract_baseline(_ys_GEM_sum, gem_baseline);
#endif
		ParameterPile::experiment_area area_ = ParameterPile::areas_to_draw.back().to_point();
		area_.experiments = _exp.experiments;
		if (ParameterPile::draw_required(area_)){ //TODO: remove this condition?
			//xs_GEM.insert(xs_GEM.begin(), 0);
			//ys_GEM.insert(ys_GEM.begin(), 0);
			std::string exp = _exp.experiments.back();
			bool invalid = false;
			for (auto i = exp.begin(); i != exp.end(); invalid ? i = exp.begin() : ++i){
				invalid = false;
				if (*i == '_'){
					if (i != exp.begin())
						if (*(i - 1) == '\\')
							continue;
					exp.insert(i, '\\');
					invalid = true;
				}
			}
			Drawing *dr = graph_manager.GetDrawing("GEM_" + area_.experiments.back(), 2, ParameterPile::DrawEngine::Gnuplot);
			dr->AddToDraw(_xs_GEM_sum, _ys_GEM_sum, "GEM\\\_" + exp, "", 0);
			DVECTOR GEM_int;
			SignalOperations::integrate(_xs_GEM_sum, _ys_GEM_sum, GEM_int,0/*baseline*/);
			dr->AddToDraw(_xs_GEM_sum, GEM_int, "GEM\\\_I\\\_" + exp, "axes x1y2", 0);

			//find start GEM time 
			DITERATOR x_start = _xs_GEM_sum.begin();
			find_GEM_start_time(_xs_GEM_sum, _ys_GEM_sum, x_start, ParameterPile::GEM_N_of_averaging, graph_manager);
			if (x_start != _xs_GEM_sum.end()) {
				dr->AddToDraw_vertical(*x_start, -1, 1, "lc rgb \"#FF0000\"", 0);
				//find finish GEM time
				DITERATOR x_finish = _xs_GEM_sum.begin();
				double temp_y_max;
				DVECTOR temp_xs = _xs_GEM_sum, temp_ys_I = GEM_int;
				SignalOperations::apply_time_limits(temp_xs, temp_ys_I, *x_start, _xs_GEM_sum.back());
				SignalOperations::get_max(temp_xs, temp_ys_I, x_finish, temp_y_max, ParameterPile::GEM_N_of_averaging);
				if (x_finish != temp_xs.end())
					dr->AddToDraw_vertical(*x_finish, -1, 1, "lc rgb \"#FF0000\"", 0);
				std::cout << "Experiment " << area_.experiments.back() << " GEM processed" << std::endl;
				std::cout << "# of runs " << N_of_runs << std::endl;
				std::cout << "# of valid runs " << N_of_valid_runs << std::endl;
				//std::cout << "# of runs with empty PMT signal" << runs_no_PMT << std::endl;
				if (x_finish != temp_xs.end() && x_start != _xs_GEM_sum.end()){
					DITERATOR _begin_ =_xs_GEM_sum.begin();
					DITERATOR _end_ = _xs_GEM_sum.end() - 1;
					x_finish = SignalOperations::find_x_iterator_by_value(_begin_,_end_, *x_finish);
					double Integ = *(GEM_int.begin() + (x_finish - _xs_GEM_sum.begin())) - *(GEM_int.begin() + (x_start - _xs_GEM_sum.begin()));
					std::cout << "GEM integral is " << Integ << std::endl;

					std::ofstream GEM_results_file;
					GEM_results_file.open(std::string(OUTPUT_DIR) + OUTPUT_GEMS + "/" + OUTPUT_GEMS+".txt", std::ios_base::app);
					GEM_results_file << area_.experiments.back() << "\t" << Integ << "\t" << *x_start << "\t" << *x_finish << "\t" << N_of_valid_runs << "\t"
						<< N_of_runs << "\t" << S_peaks_cutoff << std::endl;
					GEM_results_file.close();

					std::cout << "GEM timelimits are: [" << *x_start << ";" << *x_finish << "]" << std::endl;
				} else {
					std::cout << "GEM time boundaries are invalid" << std::endl;
				}
				std::ofstream GEM_results_file;
				open_output_file(std::string(OUTPUT_DIR) + OUTPUT_GEMS + "/" + OUTPUT_GEMS+"_" + area_.experiments.back() + ".dat",GEM_results_file, std::ios_base::trunc);
				for (int i = 0, size_ = _xs_GEM_sum.size(); i != size_; ++i)
					GEM_results_file << _xs_GEM_sum[i] << "\t" << _ys_GEM_sum[i] << std::endl;
				GEM_results_file.close();
			} else {
				std::cout << "GEM x-y size mismatch" << std::endl;
			}
		}
	}
#endif
	bool valid = true;
	if ((mppc_peaks_in_S2_area.size() != mppc_S2_start_time.size()) || (mppc_peaks_in_S2_area.size() != mppc_S2_finish_time.size())
		/*|| (mppc_peaks_in_S2_area.size() != mppc_all_peaks_Ss.size())*/||(mppc_peaks_in_S2_area.size()!=mppc_channels.size())
		|| (mppc_peaks_in_S2_area.size() != mppc_double_Is.size()) || (mppc_peaks_in_S2_area.size() != mppc_peaks.size()))
		valid = false;
	if(1==Iteration_N){
		if (!_xs_PMT3_sum.empty()){
			for (auto i= _ys_PMT3_sum.begin(),_end_=_ys_PMT3_sum.end();i!=_end_;++i)
				*i/=N_of_valid_runs;
			DVECTOR xs_before_S1 = _xs_PMT3_sum, ys_before_S1=_ys_PMT3_sum;
			SignalOperations::apply_time_limits(xs_before_S1,ys_before_S1,*xs_before_S1.begin(),ParameterPile::S1_start_time);
			double baseline = SignalOperations::find_baseline_by_integral(0,xs_before_S1,ys_before_S1);
			SignalOperations::substract_baseline(_ys_PMT3_sum, baseline);
			SignalOperations::integrate(_xs_PMT3_sum,_ys_PMT3_sum,ys_before_S1, 0);
			Drawing* dr = graph_manager.GetDrawing("3PMT_"+_exp.experiments.back()+"\\_AVR\\_",0,ParameterPile::DrawEngine::Gnuplot);
			dr->AddToDraw(_xs_PMT3_sum,_ys_PMT3_sum,"3PMT average signal "+_exp.experiments.back());
			dr->AddToDraw(_xs_PMT3_sum,ys_before_S1,"3PMT I of average signal "+_exp.experiments.back(),"axes x1y2");
		}
		if (!_xs_PMT1_sum.empty()){
			for (auto i= _ys_PMT1_sum.begin(),_end_=_ys_PMT1_sum.end();i!=_end_;++i)
				*i/=N_of_valid_runs;
			DVECTOR xs_before_S1 = _xs_PMT1_sum, ys_before_S1=_ys_PMT1_sum;
			SignalOperations::apply_time_limits(xs_before_S1,ys_before_S1,*xs_before_S1.begin(),ParameterPile::S1_start_time);
			double baseline = SignalOperations::find_baseline_by_integral(0,xs_before_S1,ys_before_S1);
			SignalOperations::substract_baseline(_ys_PMT1_sum, baseline);
			SignalOperations::integrate(_xs_PMT1_sum,_ys_PMT1_sum,ys_before_S1, 0);
			Drawing* dr = graph_manager.GetDrawing("PMT#1_"+_exp.experiments.back()+"\\\_AVR\\\_",1,ParameterPile::DrawEngine::Gnuplot);
			dr->AddToDraw(_xs_PMT1_sum,_ys_PMT1_sum,"PMT#1 average signal "+_exp.experiments.back());
			dr->AddToDraw(_xs_PMT1_sum,ys_before_S1,"PMT#1 I of average signal "+_exp.experiments.back(),"axes x1y2");
		}
	}
	if (valid && !(mppc_peaks_in_S2_area.empty())) {
		ParameterPile::experiment_area area_(ParameterPile::experiment_area::Type::Point);
		area_.experiments.push_back(_exp.experiments.back()); //done: TODO: account for MPPC channel
		if (ParameterPile::draw_required(area_)) {
			std::string fname = std::string(OUTPUT_DIR) + OUTPUT_MPPCS + area_.experiments.back()+".dat";
			std::ofstream output;
			open_output_file(fname, output);
			output << "MPPC#\tS_1pe_avr\tS_2pe_avr\tS2_S_avr\tdouble_I_avr\tS_1pe_sigma\tS_2pe_sigma\tS2_S_sigma\tdouble_I_sigma\ttime_left_avr\ttime_right_avr" << std::endl;
			for (int ch = 0; ch < mppc_peaks_in_S2_area.size(); ++ch) {
				std::string Ss_name = area_.experiments.back()+"_MPPC#" + std::to_string(mppc_channels[ch]) + "_peaks_S";
				std::string S2_S_name = area_.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_S";
				std::string S2_start_t_name = area_.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_start_t";
				std::string S2_finish_t_name = area_.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_finish_t";
				std::string double_I_name = area_.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_double_I";
				
				//done: TODO: Figure out MPPC channel. And add historam creating function
				//TODO: ParameterPile
				//TH1D *hist_S = createMPPCHist(mppc_all_peaks_Ss[ch], Ss_name, 0, 4.0, 60);
				TH1D *hist_S = createMPPCHist_peaks_S(mppc_peaks[ch], Ss_name, 0, 4.0, 60);
				TH1D *hist_S2_S = createMPPCHist(mppc_peaks_in_S2_area[ch], S2_S_name, 0, 4.0, 60);
				TH1D *hist_S2_start_t = createMPPCHist(mppc_S2_start_time[ch], S2_start_t_name, ParameterPile::S2_start_time.find(_exp.experiments.back())->second, 4.0, 60);
				TH1D *hist_S2_finish_t = createMPPCHist(mppc_S2_finish_time[ch], S2_finish_t_name, ParameterPile::S2_start_time.find(_exp.experiments.back())->second, 6.0, 60);
				TH1D *hist_double_I = createMPPCHist(mppc_double_Is[ch], double_I_name, -1, 6.0, 60);

				/*TF1 *g1 = new TF1("m1", "gaus", hist_S->GetMinimum(), hist_S->GetMaximum());
				TF1 *g2 = new TF1("m2", "gaus", hist_S->GetMinimum(), hist_S->GetMaximum());*/
				TF1 *_S_fit = new TF1(Ss_name.c_str(), "gaus(0)+gaus(3)", hist_S ? hist_S->GetXaxis()->GetXmin() : -1,
					hist_S ? hist_S->GetXaxis()->GetXmax() : 1);
				if (NULL != hist_S){
					_S_fit->SetParLimits(0, 0, hist_S->GetMaximum()*1.2);//TODO: ParameterPile
					_S_fit->SetParLimits(1, hist_S->GetXaxis()->GetXmin(), hist_S->GetXaxis()->GetXmax());
					_S_fit->SetParLimits(3, 0, hist_S->GetMaximum() * 1);//TODO: ? ParameterPile
					_S_fit->SetParLimits(4, hist_S->GetBinCenter(hist_S->GetMaximumBin()), hist_S->GetXaxis()->GetXmax());

					_S_fit->SetParameter(1, hist_S->GetBinCenter(hist_S->GetMaximumBin()));
					_S_fit->SetParameter(2, hist_S->GetBinWidth(0));
					_S_fit->SetParameter(4, 2 * hist_S->GetBinCenter(hist_S->GetMaximumBin())); //'2*' because 2 photoelectron peak
					_S_fit->SetParameter(5, hist_S->GetBinWidth(0));
				}
				TF1 *_S2_S_fit = createMPPCFitFunc(hist_S2_S, S2_S_name);
				TF1 *_S2_start_t_fit = createMPPCFitFunc(hist_S2_start_t, S2_start_t_name);
				TF1 *_S2_finish_t_fit = createMPPCFitFunc(hist_S2_finish_t, S2_finish_t_name);
				TF1 *_double_I_fit = createMPPCFitFunc(hist_double_I, double_I_name);

#ifdef _TEMP_CODE
				if (hist_S2_start_t)
					hist_S2_start_t->Fit(_S2_start_t_fit, "Q"); //affects the name of the previous to this code fit in canvas, So not rendered fits
				//should be placed before all of the canvaces
				if (hist_S2_finish_t)
					hist_S2_finish_t->Fit(_S2_finish_t_fit, "Q");
#endif

				gStyle->SetOptFit(102);
				TCanvas *c1 = new TCanvas(Ss_name.c_str(), Ss_name.c_str());
				c1->cd();
				c1->SetTitle(Ss_name.c_str());
				if (hist_S){
					hist_S->Fit(_S_fit, "Q");
					hist_S->Draw();
				}
				std::cout << "GaussSum pars: " << _S_fit->GetParameter(0) << "; "<<_S_fit->GetParameter(1) << "; "<<_S_fit->GetParameter(2) << std::endl;
				//_S_fit->Draw();
				c1->Update();

				TCanvas *c2 = new TCanvas(S2_S_name.c_str(), S2_S_name.c_str());
				c2->cd();
				c2->SetTitle(S2_S_name.c_str());
				if (hist_S2_S) {
					hist_S2_S->Fit(_S2_S_fit, "Q");
					hist_S2_S->Draw();
				}
				//_S2_S_fit->Draw();
				c2->Update();

#ifndef _TEMP_CODE
				TCanvas *c3 = new TCanvas(S2_start_t_name.c_str(), S2_start_t_name.c_str());
				c3->cd();
				c3->SetTitle(S2_start_t_name.c_str());
				if (hist_S2_start_t) {
					hist_S2_start_t->Fit(_S2_start_t_fit,"Q");
					hist_S2_start_t->Draw();
				}
				//_S2_start_t_fit->Draw();
				c3->Update();

				TCanvas *c4 = new TCanvas(S2_finish_t_name.c_str(), S2_finish_t_name.c_str());
				c4->cd();
				c4->SetTitle(S2_finish_t_name.c_str());
				if (hist_S2_finish_t) {
					hist_S2_finish_t->Fit(_S2_finish_t_fit,"Q");
					hist_S2_finish_t->Draw();
				}
				//_S2_finish_t_fit->Draw();
				c4->Update();		
#endif
				TCanvas *c5 = new TCanvas(double_I_name.c_str(), double_I_name.c_str());
				c5->cd();
				c5->SetTitle(double_I_name.c_str());
				if (hist_double_I) {
					hist_double_I->Fit(_double_I_fit,"Q");
					hist_double_I->Draw();
				}
				//_S2_finish_t_fit->Draw();
				c5->Update();
#ifdef OUTPUT_MPPCS_PICS
				std::string output_prefix =std::string(ParameterPile::this_path)+"\\"+ std::string(OUTPUT_DIR) + OUTPUT_MPPCS_PICS + area_.experiments.back() 
					+ "\\" + OUTPUT_MPPCS + std::to_string(mppc_channels[ch]) + "\\" + OUTPUT_MPPCS + std::to_string(mppc_channels[ch]) + "_";
				//vector_to_file(mppc_all_peaks_Ss[ch], output_prefix + "Ss.dat");
				vector_to_file(mppc_peaks_in_S2_area[ch], output_prefix + "S2_S.dat","MPPC_S2");
				vector_to_file(mppc_S2_start_time[ch], output_prefix + "S2_start_t.dat", "MPPC_st");
				vector_to_file(mppc_S2_finish_time[ch], output_prefix + "S2_finish_t.dat", "MPPC_fin");
				vector_to_file(mppc_double_Is[ch], output_prefix + "double_I.dat", "MPPC_II");
				vector_to_file(mppc_peaks[ch], output_prefix + "peaks.dat","MPPC_peaks");
				c1->SaveAs((output_prefix + "Ss.png").c_str(),"png");
				c2->SaveAs((output_prefix + "S2_S.png").c_str(),"png");
#ifndef _TEMP_CODE
				c3->SaveAs((output_prefix + "S2_stat_t.png").c_str(),"png");
				c4->SaveAs((output_prefix + "S2_finish_t.png").c_str(),"png");
#endif
				c5->SaveAs((output_prefix + "double_I.png").c_str(),"png");
#endif
				//output << "MPPC#\tS_1pe_avr\tS_2pe_avr\tS2_S_avr\tdouble_I_avr\tS_1pe_sigma\tS_2pe_sigma\tS2_S_sigma\tdouble_I_sigma\ttime_left_avr\ttime_right_avr" << std::endl;
				output << mppc_channels[ch] << "\t" << _S_fit->GetParameter(1) << "\t" << _S_fit->GetParameter(4)
					<< "\t" << _S2_S_fit->GetParameter(1) << "\t" << _double_I_fit->GetParameter(1)
					<< "\t" << _S_fit->GetParameter(2) << "\t" << _S_fit->GetParameter(5) << "\t" << _S2_S_fit->GetParameter(2)
					<< "\t" << _double_I_fit->GetParameter(2) << "\t" << _S2_start_t_fit->GetParameter(1)
					<< "\t" << _S2_finish_t_fit->GetParameter(1) << std::endl;
				std::cout << "MPPC " << area_.experiments.back() << " processed" << std::endl;
				std::cout << "# of runs " << N_of_runs << std::endl;
				std::cout << "# of valid runs " << N_of_valid_runs << std::endl;
				
				c1->Close();
				c2->Close();
#ifndef _TEMP_CODE
				c3->Close();
				c4->Close();
#endif //_TEMP_CODE			
				c5->Close();
				}

			output.close();
		}
	}
#ifdef _USE_TIME_STATISTICS
	time_stat.t_RUN_proc += time_stat.t_RUN_proc_single_iter;
	if (0==Iteration_N)
		time_stat.n_RUN_proc = time_stat.n_RUN_proc_single_iter;
	if (ParameterPile::Max_iteration_N == Iteration_N)
		report_time_statistics();
#endif
	++Iteration_N;
	N_of_runs = 0;
	N_of_valid_runs = 0;
	graph_manager.Draw();
	graph_manager.Clear();
}

void AllRunsResults::vector_to_file(STD_CONT<STD_CONT<peak>> &pks, std::string fname, std::string title)
{
	std::ofstream str;
	open_output_file(fname, str, std::ios_base::trunc | std::ios_base::binary);
	std::size_t sz = pks.size();
	str.write((char*)&sz,sizeof(std::size_t));
	auto _end_ = pks.end();
	for (auto run = pks.begin(); run != _end_; ++run){
		auto run_end_ = run->end();
		sz = run->size();
		str.write((char*)&sz,sizeof(std::size_t));
		for (auto pp = run->begin(); pp != run_end_; ++pp){
			str.write((char*)&pp->left, sizeof(double));
			str.write((char*)&pp->right, sizeof(double));
			str.write((char*)&pp->S, sizeof(double));
			str.write((char*)&pp->A,sizeof(double));
		}
	}
	str.close();
}

void AllRunsResults::vector_to_file(DVECTOR &what, std::string fname, std::string title)
{
	//ROOT's TTree fails for unknown reason (not dictionary-related I guess). Crashes at Fill() (libRIO.dll) but not 100% of the time,
	//so it is somehow has undefined behaviour. gDebug results in no usefull info. Plus the files that were written were significantly
	//larger than mine, so I opted againts TTree usage.

	//TFile *file = new TFile(fname.c_str(), "RECREATE");
	//TTree *t = new TTree(title.c_str(), title.c_str());
	//t->Branch(title.c_str(), &what);
	//gDebug = 2;
	//t->Fill();
	////t->GetBranch(title.c_str())->SetFile(fname.c_str());
	//t->Write(/*fname.c_str()*/);
	////t->Delete();
	//file->Close();
	//file->Delete();
	std::ofstream str;
	open_output_file(fname, str, std::ios_base::trunc | std::ios_base::binary);
	std::size_t sz = what.size();
	str.write((char*)&sz, sizeof(std::size_t));
	auto _end_ = what.end();
	for (auto i = what.begin(); i != _end_; ++i)
		str.write((char*)&(*i),sizeof(double));
	str.close();
}

TF1* AllRunsResults::createMPPCFitFunc(TH1D* hist, std::string name)
{
	TF1 *out = new TF1(name.c_str(), "gaus",hist? hist->GetXaxis()->GetXmin():-1, hist? hist->GetXaxis()->GetXmax():1);
	out->SetTitle(name.c_str());
	if (NULL!=hist) {
		out->SetParLimits(0, 0, hist->GetMaximum()*1.1);//TODO: ParameterPile
		out->SetParLimits(1, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
		out->SetParLimits(2, hist->GetBinWidth(0), hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin());
		out->SetParameter(2, 2 * hist->GetBinWidth(0));
	}
	return out;
}

TH1D* AllRunsResults::createMPPCHist_peaks_S(STD_CONT<STD_CONT<peak>> &what, std::string name, double left_cutoff, double right_cutoff_from_RMS, int N_bins)
{
	if (what.empty())
		return NULL;
	double Val_max, Val_min;
	int first_run = 0;
	for (int runs = 0; runs != what.size(); ++runs){
		auto _end_ = what[runs].end();
		STD_CONT<peak>::iterator V_max = std::max_element(what[runs].begin(), what[runs].end(), [](peak a, peak b) {
			return a.S < b.S;
		});
		STD_CONT<peak>::iterator V_min = std::min_element(what[runs].begin(), what[runs].end(), [left_cutoff](peak a, peak b) {
			if (a.S <= left_cutoff)
				return false;
			if (b.S <= left_cutoff)
				return true;
			return a.S < b.S;
		});
		if (first_run == runs) { //initializing Val_max/min. run may contain no peaks so it must be accounted for
			if (V_max == _end_ || V_min == _end_)
				++first_run;
			else {
				Val_max = V_max->S;
				Val_min = std::max(V_min->S, left_cutoff);
			}
		} else {
			Val_max = (V_max == _end_) ? Val_max : std::max(V_max->S, Val_max);
			Val_min = (V_min == _end_) ? Val_min : std::min(std::max(V_min->S, left_cutoff), Val_min);
		}
	}
	if (first_run == what.size()){
		Val_max = 10;
		Val_min = 0;
	}
	//TODO: Figure out MPPC channel. And add historam creating function
	
	std::function<double(peak&)> picker = [](peak &pk){return pk.S; };
	double V_mean = SignalOperations::Mean(what.begin(), what.end(), picker);
	double V_RMS = SignalOperations::RMS(what.begin(), what.end(), picker);
	
	double V_left_cutoff = std::max(left_cutoff, Val_min);
	double V_right_cutoff = std::min(V_mean + V_RMS * right_cutoff_from_RMS, Val_max);

	if (N_bins <= 0) {
		N_bins = 0;
		for (int _run = 0; _run < what.size(); ++_run)
			for (int _n = 0; _n < what[_run].size();++_n)
				if ((what[_run][_n].S >= V_left_cutoff) && (V_right_cutoff > what[_run][_n].S))
					N_bins++;
		if (N_bins <= 0)
			N_bins = 1;
		N_bins = std::sqrt(N_bins);
	}
	TH1D *hist = new TH1D(name.c_str(), name.c_str(), N_bins, V_left_cutoff, V_right_cutoff);
	hist->SetTitle(name.c_str());

	for (int _run = 0; _run < what.size(); ++_run)
		for (int _n = 0; _n < what[_run].size(); ++_n)
			if ((what[_run][_n].S >= V_left_cutoff) && (V_right_cutoff > what[_run][_n].S))
				hist->Fill(what[_run][_n].S);

	//TCanvas *c1 = new TCanvas(name.c_str(), name.c_str());
	//c1->cd();
	//hist->Draw();
	//c1->Update();
	return hist;
}

TH1D* AllRunsResults::createMPPCHist(DVECTOR &what, std::string name, double left_cutoff, double right_cutoff_from_RMS, int N_bins)
{
	if (what.empty())
		return NULL;
	DITERATOR V_max = std::max_element(what.begin(), what.end());
	DITERATOR V_min = std::min_element(what.begin(), what.end(), [left_cutoff](double a, double b) {
		if (a <= left_cutoff)
			return false;
		if (b <= left_cutoff)
			return true;
		return a < b;
	});
	//TODO: Figure out MPPC channel. And add historam creating function
	double V_mean = TMath::Mean(what.begin(), what.end());
	double V_RMS = TMath::RMS(what.begin(), what.end());
	double V_left_cutoff = std::max(left_cutoff, *V_min);
	double V_right_cutoff = std::min(V_mean + V_RMS * right_cutoff_from_RMS,*V_max);
	if (0==V_RMS)
		V_right_cutoff = 2*V_mean;
	if (N_bins <= 0) {
		N_bins = 0;
		for (int _n = 0; _n < what.size(); ++_n)
			if ((what[_n]>=V_left_cutoff) && (V_right_cutoff > what[_n]))
				N_bins++;
		if (N_bins <= 0)
			N_bins = 1;
		N_bins = std::sqrt(N_bins);
	}
	TH1D *hist = new TH1D(name.c_str(), name.c_str(), N_bins, V_left_cutoff, V_right_cutoff);
	hist->SetTitle(name.c_str());

	for (int _n = 0; _n < what.size(); ++_n)
		if ((what[_n]>=V_left_cutoff)&&(V_right_cutoff>what[_n]))
			hist->Fill(what[_n]);

	//TCanvas *c1 = new TCanvas(name.c_str(), name.c_str());
	//c1->cd();
	//hist->Draw();
	//c1->Update();
	return hist;
}

int AllRunsResults::Iteration(void) const
{	return Iteration_N;}

void AllRunsResults::Clear(void)
{
#ifdef _HOTFIX_CLEAR_MEMORY
	DVECTOR().swap(_Ss);
	DVECTOR().swap(_ns);
	DVECTOR().swap(_xs_GEM_sum);
	DVECTOR().swap(_ys_GEM_sum);
	STD_CONT<DVECTOR>().swap(mppc_peaks_in_S2_area);
	STD_CONT<DVECTOR>().swap(mppc_S2_start_time);
	STD_CONT<DVECTOR>().swap(mppc_S2_finish_time);
	//STD_CONT<DVECTOR>().swap(mppc_all_peaks_Ss);
	STD_CONT<STD_CONT<STD_CONT<peak>>>().swap(mppc_peaks);
	STD_CONT<DVECTOR>().swap(mppc_double_Is);
	STD_CONT<int>().swap(mppc_channels);
#else
	_Ss.clear();
	_ns.clear();
	_xs_GEM_sum.clear();
	_ys_GEM_sum.clear();
	graph_manager.Clear();
	mppc_peaks_in_S2_area.clear();
	mppc_S2_start_time.clear();
	mppc_S2_finish_time.clear();
	//mppc_all_peaks_Ss.clear();
	mppc_peaks.clear();
	mppc_double_Is.clear();
	mppc_channels.clear();
#endif
	//N_peaks_cutoff preserve
	//S_peaks_cutoff preserve
	N_of_runs = 0;
	//Iteration_N preserve;
#ifdef _USE_TIME_STATISTICS
	time_stat.t_RUN_proc_single_iter=0;
	time_stat.n_RUN_proc_single_iter=0;
#endif
}

#ifdef _USE_TIME_STATISTICS
void AllRunsResults::report_time_statistics()
{
	double coeff = 1e-6;//to milliseconds
	std::cout << "________________________________________________" << std::endl;
	std::cout << "RUN TIMES:" << std::endl;
	std::cout << "Experiment: " << _exp.experiments.back()<<std::endl;
	std::cout << "All times are in milliseconds" << std::endl;
	std::cout << "Number of runs: " << N_of_runs<<std::endl;
	std::cout << "T:" << (double)time_stat.t_RUN_proc*coeff << "\tN: "
		<< time_stat.n_RUN_proc << "\tAvr: " << ((double)time_stat.t_RUN_proc*coeff) / time_stat.n_RUN_proc << std::endl;
	std::cout << "|" << std::endl;
	std::cout << "|____" << "PMT processing:" << std::endl;
	std::cout << "   | " << "T: " << (double)time_stat.t_PMT_proc*coeff << "\tN: "
		<< time_stat.n_PMT_proc << "\tAvr: " << ((double)time_stat.t_PMT_proc*coeff) / time_stat.n_PMT_proc << std::endl;
	std::cout << "   | " << "  |" << std::endl;
	std::cout << "   | " << "  |____" << "Reading:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_PMT_file_reading*coeff << "\tN: "
		<< time_stat.n_PMT_file_reading << "\tAvr: " << ((double)time_stat.t_PMT_file_reading*coeff) / time_stat.n_PMT_file_reading << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Filtering:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_PMT_filtering*coeff << "\tN: "
		<< time_stat.n_PMT_filtering << "\tAvr: " << ((double)time_stat.t_PMT_filtering*coeff) / time_stat.n_PMT_filtering << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Baseline:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_PMT_baseline*coeff << "\tN: "
		<< time_stat.n_PMT_baseline << "\tAvr: " << ((double)time_stat.t_PMT_baseline*coeff) / time_stat.n_PMT_baseline << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Peaks:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_PMT_peaks*coeff << "\tN: "
		<< time_stat.n_PMT_peaks << "\tAvr: " << ((double)time_stat.t_PMT_peaks*coeff) / time_stat.n_PMT_peaks << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Residual time:" << std::endl;
	std::cout << "   | " << "     | " << coeff*(double)(time_stat.t_PMT_proc
		- time_stat.t_PMT_file_reading - time_stat.t_PMT_filtering - time_stat.t_PMT_baseline - time_stat.t_PMT_peaks) << std::endl;
	std::cout << "   | " << std::endl;
	std::cout << "   | " << std::endl;
	std::cout << "   |_" << "MPPC processing:" << std::endl;
	std::cout << "   | " << "T: " << (double)time_stat.t_MPPC_proc*coeff << "\tN: "
		<< time_stat.n_MPPC_proc << "\tAvr: " << ((double)time_stat.t_MPPC_proc*coeff) / time_stat.n_MPPC_proc << std::endl;
	std::cout << "   | " << "  |" << std::endl;
	std::cout << "   | " << "  |____" << "Reading:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_file_reading*coeff << "\tN: "
		<< time_stat.n_MPPC_file_reading << "\tAvr: " << ((double)time_stat.t_MPPC_file_reading*coeff) / time_stat.n_MPPC_file_reading << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Filtering:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_filtering*coeff << "\tN: "
		<< time_stat.n_MPPC_filtering << "\tAvr: " << ((double)time_stat.t_MPPC_filtering*coeff) / time_stat.n_MPPC_filtering << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Threshold and first baseline:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_threshold_and_first_baseline*coeff << "\tN: "
		<< time_stat.n_MPPC_threshold_and_first_baseline << "\tAvr: " <<
		((double)time_stat.t_MPPC_threshold_and_first_baseline*coeff) / time_stat.n_MPPC_threshold_and_first_baseline << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Curved baseline:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_curved_baseline*coeff << "\tN: "
		<< time_stat.n_MPPC_curved_baseline << "\tAvr: " <<
		((double)time_stat.t_MPPC_curved_baseline*coeff) / time_stat.n_MPPC_curved_baseline << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Curved baseline v2:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_curved_baseline_v2*coeff << "\tN: "
		<< time_stat.n_MPPC_curved_baseline_v2 << "\tAvr: " <<
		((double)time_stat.t_MPPC_curved_baseline_v2*coeff) / time_stat.n_MPPC_curved_baseline_v2 << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Curved baseline v3:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_curved_baseline_v3*coeff << "\tN: "
		<< time_stat.n_MPPC_curved_baseline_v3 << "\tAvr: " <<
		((double)time_stat.t_MPPC_curved_baseline_v3*coeff) / time_stat.n_MPPC_curved_baseline_v3 << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Curved baseline v4:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_curved_baseline_v4*coeff << "\tN: "
		<< time_stat.n_MPPC_curved_baseline_v4 << "\tAvr: " <<
		((double)time_stat.t_MPPC_curved_baseline_v4*coeff) / time_stat.n_MPPC_curved_baseline_v4 << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	/*std::cout << "   | " << "     |_" << "Curved baseline v5:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_curved_baseline_v5*coeff << "\tN: "
		<< time_stat.n_MPPC_curved_baseline_v5 << "\tAvr: " <<
		((double)time_stat.t_MPPC_curved_baseline_v5*coeff) / time_stat.n_MPPC_curved_baseline_v5 << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Curved baseline v6:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_curved_baseline_v6*coeff << "\tN: "
		<< time_stat.n_MPPC_curved_baseline_v6 << "\tAvr: " <<
		((double)time_stat.t_MPPC_curved_baseline_v6*coeff) / time_stat.n_MPPC_curved_baseline_v6 << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Curved baseline v7:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_curved_baseline_v7*coeff << "\tN: "
		<< time_stat.n_MPPC_curved_baseline_v7 << "\tAvr: " <<
		((double)time_stat.t_MPPC_curved_baseline_v7*coeff) / time_stat.n_MPPC_curved_baseline_v7 << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Curved baseline v8:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_curved_baseline_v8*coeff << "\tN: "
		<< time_stat.n_MPPC_curved_baseline_v8 << "\tAvr: " <<
		((double)time_stat.t_MPPC_curved_baseline_v8*coeff) / time_stat.n_MPPC_curved_baseline_v8 << std::endl;
	std::cout << "   | " << "     |" << std::endl;*/
	std::cout << "   | " << "     |_" << "Curved baseline's baseline:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_curved_baseline_baseline*coeff << "\tN: "
		<< time_stat.n_MPPC_curved_baseline_baseline << "\tAvr: " << 
		((double)time_stat.t_MPPC_curved_baseline_baseline*coeff) / time_stat.n_MPPC_curved_baseline_baseline << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Baseline's substraction:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_baseline_substraction*coeff << "\tN: "
		<< time_stat.n_MPPC_baseline_substraction << "\tAvr: " << 
		((double)time_stat.t_MPPC_baseline_substraction*coeff) / time_stat.n_MPPC_baseline_substraction << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Peaks finding:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_peaks_finding*coeff << "\tN: "
		<< time_stat.n_MPPC_peaks_finding << "\tAvr: " <<
		((double)time_stat.t_MPPC_peaks_finding*coeff) / time_stat.n_MPPC_peaks_finding << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Peaks processing:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_peaks_processing*coeff << "\tN: "
		<< time_stat.n_MPPC_peaks_processing << "\tAvr: " <<
		((double)time_stat.t_MPPC_peaks_processing*coeff) / time_stat.n_MPPC_peaks_processing << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Double integral:" << std::endl;
	std::cout << "   | " << "     | " << "T:" << (double)time_stat.t_MPPC_double_I*coeff << "\tN: "
		<< time_stat.n_MPPC_double_I << "\tAvr: " << ((double)time_stat.t_MPPC_double_I*coeff) / time_stat.n_MPPC_double_I << std::endl;
	std::cout << "   | " << "     |" << std::endl;
	std::cout << "   | " << "     |_" << "Residual time:" << std::endl;
	std::cout << "   | " << "     | " << coeff*(double)(time_stat.t_MPPC_proc
		- time_stat.t_MPPC_file_reading - time_stat.t_MPPC_filtering - time_stat.t_MPPC_threshold_and_first_baseline
		- time_stat.t_MPPC_curved_baseline - time_stat.t_MPPC_curved_baseline_baseline - time_stat.t_MPPC_baseline_substraction
		- time_stat.t_MPPC_peaks_finding - time_stat.t_MPPC_peaks_processing - time_stat.t_MPPC_double_I) << std::endl;
	std::cout << "________________________________________________" << std::endl;
}
#endif
