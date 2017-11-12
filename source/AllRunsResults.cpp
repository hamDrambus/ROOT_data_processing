#include "AllRunsResults.h"

AllRunsResults::AllRunsResults(ParameterPile::experiment_area experiment)
{
	_exp = experiment;
	S_peaks_cutoff = ParameterPile::PMT_SArea_peaks_acceptance;
	S_peaks_max_cutoff = S_peaks_cutoff - 1;
	N_peaks_cutoff = ParameterPile::PMT_N_peaks_acceptance;
	N_of_runs = 0;
	Iteration_N = 0;
}

void AllRunsResults::processAllRuns(STD_CONT<SingleRunResults> &single_results)
{
	//this->Clear() called before in the AnalysisManager. 
	int n_valid_runs = 0;
	N_of_runs = single_results.size();
	for (auto i = single_results.begin(); i != single_results.end(); ++i){
		if (i->isValid()) { //has meaning only at the secondary function call
			n_valid_runs++;
			_Ss.push_back((*i).PMT3_summed_peaks_area);
			_ns.push_back((*i).PMT3_n_peaks);
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
				|| (mppc_peaks_in_S2_area.size() != mppc_all_peaks_Ss.size()) || (mppc_peaks_in_S2_area.size()!=mppc_channels.size()))
				valid = false;
			if ((i->mppc_S2_peaks_area.size() != i->mppc_S2_start_t.size()) || (i->mppc_S2_peaks_area.size() != i->mppc_S2_finish_t.size())
				|| (i->mppc_S2_peaks_area.size() != i->mppc_peaks.size()) || (i->mppc_S2_peaks_area.size() != i->mppc_channels.size()))
				valid = false;
			if (!valid) {
				std::cout << "SingleRunResulst MPPC channels size mismatch" << std::endl;
				continue; 
			}
			int sz = i->mppc_S2_peaks_area.size();
			if (first){
				mppc_peaks_in_S2_area.resize(sz, DVECTOR());
				mppc_S2_start_time.resize(sz, DVECTOR());
				mppc_S2_finish_time.resize(sz, DVECTOR());
				mppc_all_peaks_Ss.resize(sz, DVECTOR());
				mppc_double_Is.resize(sz, DVECTOR());
#ifndef _USE_DEQUE
				for (int ch = 0; ch != sz; ++ch) {
					mppc_peaks_in_S2_area[ch].reserve(single_results.size());
					mppc_S2_start_time[ch].reserve(single_results.size());
					mppc_S2_finish_time[ch].reserve(single_results.size());
				}
#endif
				mppc_channels = i->mppc_channels;
			}
			for (int ch = 0; ch != sz; ++ch) {
				mppc_peaks_in_S2_area[ch].push_back(i->mppc_S2_peaks_area[ch]);
				mppc_S2_start_time[ch].push_back(i->mppc_S2_start_t[ch]);
				mppc_S2_finish_time[ch].push_back(i->mppc_S2_finish_t[ch]);
				mppc_double_Is[ch].push_back(i->mppc_double_I[ch]);
				for (auto peaks = i->mppc_peaks[ch].begin(); peaks != i->mppc_peaks[ch].end(); ++peaks){
					mppc_all_peaks_Ss[ch].push_back(peaks->S);
				}
			}
		}
#ifdef _HOTFIX_CLEAR_MEMORY
		STD_CONT<DVECTOR>().swap(i->mppc_baseline_xs);
		STD_CONT<DVECTOR>().swap(i->mppc_baseline_ys);
		STD_CONT<int>().swap(i->mppc_channels);
		STD_CONT<STD_CONT<peak>>().swap(i->mppc_peaks);
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
	DITERATOR i_mean = SignalOperations::find_x_iterator_by_value(areas.begin(), areas.end() - 1, G_mean);
	DITERATOR i_cutoff = SignalOperations::find_x_iterator_by_value(areas.begin(), areas.end() - 1, cutoff);
	if (G_mean < cutoff)
		if ((i_cutoff - i_mean) >= N_above_global_mean_acceptable)
			cutoff = G_mean;//or set to S_peaks_cutoff ?
	S_peaks_cutoff = std::max(cutoff, S_peaks_cutoff);
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
	_Ss.insert(_Ss.end(), with->_Ss.begin(), with->_Ss.end());
	_ns.insert(_ns.end(), with->_ns.begin(), with->_ns.end());
	N_peaks_cutoff = with->N_peaks_cutoff;
	S_peaks_cutoff = with->S_peaks_cutoff;
	N_of_runs += with->N_of_runs;
	//Iteration_N = with->Iteration_N;
	if (_xs_GEM_sum.empty()) {
		_xs_GEM_sum = with->_xs_GEM_sum;
		_ys_GEM_sum = with->_ys_GEM_sum;
	} else {
		for (auto i = _ys_GEM_sum.begin(), j = with->_ys_GEM_sum.begin(); (i != _ys_GEM_sum.end()) && (j != with->_ys_GEM_sum.end()); ++j, ++i)
			*i += *j;
	}

	bool empty = false, valid = true;
	if (mppc_peaks_in_S2_area.empty())
		empty = true;
	if ((mppc_peaks_in_S2_area.size() != mppc_S2_start_time.size()) || (mppc_peaks_in_S2_area.size() != mppc_S2_finish_time.size())
		|| (mppc_peaks_in_S2_area.size() != mppc_all_peaks_Ss.size()) || (mppc_peaks_in_S2_area.size() != mppc_channels.size())
		|| (mppc_peaks_in_S2_area.size() != mppc_double_Is.size()))
		valid = false;
	if ((with->mppc_peaks_in_S2_area.size() != with->mppc_S2_start_time.size()) || (with->mppc_peaks_in_S2_area.size() != with->mppc_S2_finish_time.size())
		|| (with->mppc_peaks_in_S2_area.size() != with->mppc_all_peaks_Ss.size()) || (with->mppc_peaks_in_S2_area.size() != with->mppc_channels.size())
		|| (mppc_peaks_in_S2_area.size() != mppc_double_Is.size()))
		valid = false;
	if (!valid) {
		std::cout << "AllRunsResults contains invalid MPPC data: channels' size mismatches" << std::endl;
		return;
	}
	if (!empty &&(mppc_peaks_in_S2_area.size() != with->mppc_peaks_in_S2_area.size())){
		std::cout << "WARNING Two AllRunsResults have MPPC channels' size mismatches: not Merging" << std::endl;
		return;
	}
	if (empty) {
		mppc_peaks_in_S2_area=with->mppc_peaks_in_S2_area;
		mppc_S2_start_time = with->mppc_S2_start_time;
		mppc_S2_finish_time = with->mppc_S2_finish_time;
		mppc_all_peaks_Ss = with->mppc_all_peaks_Ss;
		mppc_channels = with->mppc_channels;
		mppc_double_Is = with->mppc_double_Is;
	} else {
		for (int ch = 0; ch < mppc_peaks_in_S2_area.size();++ch) {
			mppc_peaks_in_S2_area[ch].insert(mppc_peaks_in_S2_area[ch].end(), with->mppc_peaks_in_S2_area[ch].begin(), with->mppc_peaks_in_S2_area[ch].end());
			mppc_S2_start_time[ch].insert(mppc_S2_start_time[ch].end(), with->mppc_S2_start_time[ch].begin(), with->mppc_S2_start_time[ch].end());
			mppc_S2_finish_time[ch].insert(mppc_S2_finish_time[ch].end(), with->mppc_S2_finish_time[ch].begin(), with->mppc_S2_finish_time[ch].end());
			mppc_all_peaks_Ss[ch].insert(mppc_all_peaks_Ss[ch].end(), with->mppc_all_peaks_Ss[ch].begin(), with->mppc_all_peaks_Ss[ch].end());
			mppc_double_Is[ch].insert(mppc_double_Is[ch].end(), with->mppc_double_Is[ch].begin(), with->mppc_double_Is[ch].end());
		}
	}
}

void AllRunsResults::Merged(void)
{
	if (!_Ss.empty() && (0==Iteration_N)) {//_xs_GEM is empty at the first run
		find_S_cutoff();

		DITERATOR S_max = std::max_element(_Ss.begin(), _Ss.end());
		double hist_S_max = (S_peaks_max_cutoff > S_peaks_cutoff) ? std::min(S_peaks_max_cutoff, *S_max) : *S_max;
		TH1D *hist_S = new TH1D("PMT_S_peaks", "PMT_S_peaks", 60, 0, hist_S_max);
		/*TH1I *hist_n = new TH1I("PMT_N_peaks", "PMT_N_peaks", 30, 0, 30);*/

		/*double *NNs = new double[_ns.size()];
		double *SSs = new double[_Ss.size()];*/
		for (int _n = 0, _s = 0; (_n < _ns.size()) && (_s < _Ss.size()); ++_s, ++_n){
			if ((S_peaks_max_cutoff > S_peaks_cutoff) && (_Ss[_s] > S_peaks_max_cutoff))
				continue; //do not draw if maximum limit is imposed && element exceeds it
			/*NNs[_n] = _ns[_n];
			SSs[_s] = _Ss[_s];*/
			/*hist_n->Fill(_ns[_n]);*/
			hist_S->Fill(_Ss[_s]);
		}
		TCanvas *c1 = new TCanvas(("S_peaks_distribution " + _exp.experiments.back()).c_str(),
			("S_peaks_distribution " + _exp.experiments.back()).c_str());
		c1->cd();
		hist_S->Draw();
		c1->Update();
		TLine *cutoff = new TLine(S_peaks_cutoff, c1->GetUymin(), S_peaks_cutoff, c1->GetUymax());
		cutoff->SetLineColor(kRed);
		cutoff->Draw();
		c1->Update();
		/*TCanvas *c2 = new TCanvas(("n_peaks_distribution " + _exp.experiments.back()).c_str(),
		("n_peaks_distribution " + _exp.experiments.back()).c_str());
		c2->cd();
		hist_n->Draw();
		c2->Update();

		TCanvas *c3 = new TCanvas(("N-S Scatter " + _exp.experiments.back()).c_str(), ("N-S Scatter " + _exp.experiments.back()).c_str());
		TGraph *gr = new TGraph(std::min(_ns.size(), _Ss.size()), SSs, NNs);
		c3->cd();
		gr->Draw("ap");
		c3->Update();*/
		/*delete[] NNs;
		delete[] SSs;*/
	}
	int n_valid_runs = _Ss.size();
#ifdef _PROCESS_GEMS
	if (!_xs_GEM_sum.empty() && (0 != n_valid_runs)&&(1==Iteration_N)){
		//x_GEM is not empty only at the second processing (when proper cutoffs were applied)
		for (auto i = _ys_GEM_sum.begin(); i != _ys_GEM_sum.end(); ++i)
			*i /= n_valid_runs;
		ParameterPile::experiment_area area_ = ParameterPile::areas_to_draw.back().to_point();
		area_.experiments = _exp.experiments;
		if (ParameterPile::draw_required(area_)){
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
			Drawing *dr = graph_manager.GetDrawing("GEM_" + area_.experiments.back(), 0, ParameterPile::DrawEngine::Gnuplot);
			dr->AddToDraw(_xs_GEM_sum, _ys_GEM_sum, "GEM\\\_" + exp, "", 0);
			DVECTOR GEM_int;
			SignalOperations::integrate(_xs_GEM_sum, _ys_GEM_sum, GEM_int);
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
				std::cout << "Experiment " << area_.experiments.back() << " processed" << std::endl;
				std::cout << "# of runs " << N_of_runs << std::endl;
				std::cout << "# of valid runs " << n_valid_runs << std::endl;
				//std::cout << "# of runs with empty PMT signal" << runs_no_PMT << std::endl;
				if (x_finish != temp_xs.end() && x_start != _xs_GEM_sum.end()){
					x_finish = SignalOperations::find_x_iterator_by_value(_xs_GEM_sum.begin(), _xs_GEM_sum.end() - 1, *x_finish);
					double Integ = *(GEM_int.begin() + (x_finish - _xs_GEM_sum.begin())) - *(GEM_int.begin() + (x_start - _xs_GEM_sum.begin()));
					std::cout << "GEM integral is " << Integ << std::endl;

					std::ofstream GEM_results_file;
					GEM_results_file.open(std::string(OUTPUT_DIR) + OUTPUT_GEMS, std::ios_base::app);
					GEM_results_file << area_.experiments.back() << "\t" << Integ << "\t" << *x_start << "\t" << *x_finish << "\t" << n_valid_runs << "\t"
						<< N_of_runs << "\t" << S_peaks_cutoff << std::endl;
					GEM_results_file.close();

					std::cout << "GEM timelimits are: [" << *x_start << ";" << *x_finish << "]" << std::endl;
				} else {
					std::cout << "GEM time boundaries are invalid" << std::endl;
				}
			} else {
				std::cout << "GEM time boundaries are invalid" << std::endl;
			}
		}
	}
#endif
	bool valid = true;
	if ((mppc_peaks_in_S2_area.size() != mppc_S2_start_time.size()) || (mppc_peaks_in_S2_area.size() != mppc_S2_finish_time.size())
		|| (mppc_peaks_in_S2_area.size() != mppc_all_peaks_Ss.size())||(mppc_peaks_in_S2_area.size()!=mppc_channels.size())
		|| (mppc_peaks_in_S2_area.size() != mppc_double_Is.size()))
		valid = false;
	if (valid && !(mppc_peaks_in_S2_area.empty())) {
		ParameterPile::experiment_area area_ = ParameterPile::areas_to_draw.back().to_point();
		area_.experiments = _exp.experiments; //TODO: account for MPPC channel
		if (ParameterPile::draw_required(area_)) {
			std::string fname = std::string(OUTPUT_DIR) + OUTPUT_MPPCS + area_.experiments.back();
			std::ofstream output;
			open_output_file(fname, output);
			output << "MPPC#\tS_1pe_avr\tS_1pe_sigma\tS_2pe_avr\tS_2pe_sigma\tS2_S_avr\tS2_S_sigma\ttime_left_avr\ttime_right_avr\tdouble_I_avr\tdouble_I_sigma" << std::endl;
			for (int ch = 0; ch < mppc_peaks_in_S2_area.size(); ++ch) {
				std::string Ss_name = area_.experiments.back()+"MPPC#" + std::to_string(mppc_channels[ch]) + "_peaks_S";
				std::string S2_S_name = area_.experiments.back() + "MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_S";
				std::string S2_start_t_name = area_.experiments.back() + "MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_start_t";
				std::string S2_finish_t_name = area_.experiments.back() + "MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_finish_t";
				std::string double_I_name = area_.experiments.back() + "MPPC#" + std::to_string(mppc_channels[ch]) + "_double_I";
				
				//done: TODO: Figure out MPPC channel. And add historam creating function
				//TODO: ParameterPile
				TH1D *hist_S = createMPPCHist(mppc_all_peaks_Ss[ch], Ss_name, 0, 4.0, 60);
				TH1D *hist_S2_S = createMPPCHist(mppc_peaks_in_S2_area[ch], S2_S_name, 0, 4.0, 60);
				TH1D *hist_S2_start_t = createMPPCHist(mppc_S2_start_time[ch], S2_S_name, ParameterPile::S2_start_time, 4.0, 60);
				TH1D *hist_S2_finish_t = createMPPCHist(mppc_S2_finish_time[ch], S2_S_name, ParameterPile::S2_start_time, 6.0, 60);
				TH1D *hist_double_I = createMPPCHist(mppc_double_Is[ch], double_I_name, 0, 6.0, 60);

				/*TF1 *g1 = new TF1("m1", "gaus", hist_S->GetMinimum(), hist_S->GetMaximum());
				TF1 *g2 = new TF1("m2", "gaus", hist_S->GetMinimum(), hist_S->GetMaximum());*/
				TF1 *_S_fit = new TF1(Ss_name.c_str(), "gaus(0)+gaus(3)", hist_S->GetXaxis()->GetXmin(), hist_S->GetXaxis()->GetXmax());
				_S_fit->SetParLimits(0, 0, hist_S->GetMaximum()*1.2);//TODO: ParameterPile
				_S_fit->SetParLimits(1, hist_S->GetXaxis()->GetXmin(), hist_S->GetXaxis()->GetXmax());
				_S_fit->SetParLimits(3, 0, hist_S->GetMaximum()*1);//TODO: ? ParameterPile
				_S_fit->SetParLimits(4, hist_S->GetBinCenter(hist_S->GetMaximumBin()), hist_S->GetXaxis()->GetXmax());
				
				_S_fit->SetParameter(1, hist_S->GetBinCenter(hist_S->GetMaximumBin()));
				_S_fit->SetParameter(2, hist_S->GetBinWidth(0));
				_S_fit->SetParameter(4, 2*hist_S->GetBinCenter(hist_S->GetMaximumBin())); //'2*' because 2 photoelectron peak
				_S_fit->SetParameter(5, hist_S->GetBinWidth(0));

				TF1 *_S2_S_fit = createMPPCFitFunc(hist_S2_S, S2_S_name);
				TF1 *_S2_start_t_fit = createMPPCFitFunc(hist_S2_start_t, S2_start_t_name);
				TF1 *_S2_finish_t_fit = createMPPCFitFunc(hist_S2_finish_t, S2_finish_t_name);
				TF1 *_double_I_fit = createMPPCFitFunc(hist_double_I, double_I_name);

				TCanvas *c1 = new TCanvas(Ss_name.c_str(), Ss_name.c_str());
				c1->cd();
				hist_S->Fit(_S_fit);
				std::cout << "GaussSum pars: " << _S_fit->GetParameter(0) << "; "<<_S_fit->GetParameter(1) << "; "<<_S_fit->GetParameter(2) << std::endl;
				hist_S->Draw();
				//_S_fit->Draw();
				c1->Update();

				TCanvas *c2 = new TCanvas(S2_S_name.c_str(), S2_S_name.c_str());
				c2->cd();
				hist_S2_S->Fit(_S2_S_fit);
				hist_S2_S->Draw();
				//_S2_S_fit->Draw();
				c2->Update();

				TCanvas *c3 = new TCanvas(S2_start_t_name.c_str(), S2_start_t_name.c_str());
				c3->cd();
				hist_S2_start_t->Fit(_S2_start_t_fit);
				hist_S2_start_t->Draw();
				//_S2_start_t_fit->Draw();
				c3->Update();

				TCanvas *c4 = new TCanvas(S2_finish_t_name.c_str(), S2_finish_t_name.c_str());
				c4->cd();
				hist_S2_finish_t->Fit(_S2_finish_t_fit);
				hist_S2_finish_t->Draw();
				//_S2_finish_t_fit->Draw();
				c4->Update();

				TCanvas *c5 = new TCanvas(double_I_name.c_str(), double_I_name.c_str());
				c5->cd();
				hist_double_I->Fit(_double_I_fit);
				hist_double_I->Draw();
				//_S2_finish_t_fit->Draw();
				c5->Update();

				//output << "MPPC#\tS_1pe_avr\tS_1pe_sigma\tS_2pe_avr\tS_2pe_sigma\tS2_S_avr\tS2_S_sigma\ttime_left_avr\ttime_right_avr\tdouble_I_avr\tdouble_I_sigma" << std::endl;
				output << mppc_channels[ch] << "\t" << _S_fit->GetParameter(1) << "\t" << _S_fit->GetParameter(2)
					<< "\t" << _S_fit->GetParameter(4) << "\t" << _S_fit->GetParameter(5) << "\t" << _S2_S_fit->GetParameter(1)
					<< "\t" << _S2_S_fit->GetParameter(2) << "\t" << _S2_start_t_fit->GetParameter(1) << "\t" << _S2_finish_t_fit->GetParameter(1)
					<<"\t"<<_double_I_fit->GetParameter(1)<<"\t"<<_double_I_fit->GetParameter(2)<< std::endl;
				std::cout << "MPPC " << area_.experiments.back() << " processed" << std::endl;
				std::cout << "# of runs " << N_of_runs << std::endl;
				std::cout << "# of valid runs " << n_valid_runs << std::endl;
			}
			output.close();
		}
	}
	Iteration_N++;
	graph_manager.Draw();
}

TF1* AllRunsResults::createMPPCFitFunc(TH1D* hist, std::string name)
{
	TF1 *out = new TF1(name.c_str(), "gaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
	out->SetParLimits(0, 0, hist->GetMaximum()*1.1);//TODO: ParameterPile
	out->SetParLimits(1, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
	out->SetParLimits(2, hist->GetBinWidth(0), hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin());
	out->SetParameter(2, 2 * hist->GetBinWidth(0));
	return out;
}

TH1D* AllRunsResults::createMPPCHist(DVECTOR &what, std::string name, double left_cutoff, double right_cutoff_from_RMS, int N_bins)
{
	if (what.empty())
		return NULL;
	DITERATOR V_max = std::max_element(what.begin(), what.end());
	DITERATOR V_min = std::min_element(what.begin(), what.end());
	//TODO: Figure out MPPC channel. And add historam creating function
	double V_mean = TMath::Mean(what.begin(), what.end());
	double V_RMS = TMath::RMS(what.begin(), what.end());
	double V_left_cutoff = std::max(left_cutoff, *V_min);
	double V_right_cutoff = std::min(V_mean + V_RMS * right_cutoff_from_RMS,*V_max);
	
	if (N_bins <= 0) {
		N_bins = 0;
		for (int _n = 0; _n < what.size(); ++_n)
			if ((what[_n]>V_left_cutoff) && (V_right_cutoff > what[_n]))
				N_bins++;
		if (N_bins <= 0)
			N_bins = 1;
		N_bins = std::sqrt(N_bins);
	}
	TH1D *hist = new TH1D(name.c_str(), name.c_str(), N_bins, V_left_cutoff, V_right_cutoff);
	
	for (int _n = 0; _n < what.size(); ++_n)
		if ((what[_n]>V_left_cutoff)&&(V_right_cutoff>what[_n]))
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
	DVECTOR().swap(_Ss);
	STD_CONT<DVECTOR>().swap(mppc_peaks_in_S2_area);
	STD_CONT<DVECTOR>().swap(mppc_S2_start_time);
	STD_CONT<DVECTOR>().swap(mppc_S2_finish_time);
	STD_CONT<DVECTOR>().swap(mppc_all_peaks_Ss);
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
	mppc_all_peaks_Ss.clear();
	mppc_double_Is.clear();
	mppc_channels.clear();
#endif
	//N_peaks_cutoff preserve
	//S_peaks_cutoff preserve
	N_of_runs = 0;
	//Iteration_N preserve;
}
