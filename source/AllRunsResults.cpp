#include "AllRunsResults.h"

AllRunsResults::AllRunsResults(ParameterPile::experiment_area experiment, ParameterPile::area_vector chs_to_average)
{
	_to_average = chs_to_average;
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

void AllRunsResults::find_S_cutoff(void)//TODO: explain the algorithm and maybe test it
//"maybe" because the actual test is S histograms pictures with cutoffs. If pictures are ok, there is no need to
//test the algorithm excessively
{
	if (_Ss.size() < ParameterPile::PMT_min_statistics) //not statistically significant. TODO: figure out the number,
		//and ParameterPile, however this is not important
		return;
#if defined(_NO_AUTO_PMT_SELECTION)||defined(_NO_PMT_SELECTION)
	S_peaks_cutoff = ParameterPile::PMT_SArea_peaks_acceptance;
	S_peaks_max_cutoff = S_peaks_cutoff - 1;
	return;
#endif
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

}

void AllRunsResults::Merge(AllRunsResults* with)
{
	N_peaks_cutoff = with->N_peaks_cutoff;
	S_peaks_cutoff = with->S_peaks_cutoff;
	if (0 == Iteration_N) {
		_Ss.insert(_Ss.end(), with->_Ss.begin(), with->_Ss.end());
		_ns.insert(_ns.end(), with->_ns.begin(), with->_ns.end());
		bool empty = false, valid = true;
		if (pmt_peaks.empty())
			empty = true;
		if (pmt_channels.size() != pmt_peaks.size())
			valid = false;
		if (with->pmt_peaks.size() != with->pmt_channels.size())
			valid = false;
		if (!valid) {
			std::cout << "AllRunsResults contains invalid PMT data: channels' size mismatches" << std::endl;
			return;
		}
		if (!empty && (pmt_channels.size() != with->pmt_channels.size())){
			std::cout << "WARNING Two AllRunsResults have PMT channels' size mismatches: not Merging" << std::endl;
			return;
		}
		if (empty) {
			pmt_channels = with->pmt_channels;
			pmt_peaks = with->pmt_peaks;
			pmt_S2_integral = with->pmt_S2_integral;
		} else {
			for (int ch = 0; ch < pmt_channels.size(); ++ch){
				pmt_peaks[ch].insert(pmt_peaks[ch].end(), with->pmt_peaks[ch].begin(), with->pmt_peaks[ch].end());
				pmt_S2_integral[ch].insert(pmt_S2_integral[ch].end(), with->pmt_S2_integral[ch].begin(), with->pmt_S2_integral[ch].end());
			}
		}
	}
	if (1 == Iteration_N) {
		bool empty = false, valid = true;
		if (avr_channels.empty())
			empty = true;
		if ((avr_channels.size() != _xs_sum.size())||(avr_channels.size() != _ys_sum.size())||(avr_channels.size() != _ys_disp.size()))
			valid = false;
		if ((with->avr_channels.size() != with->_xs_sum.size())||(with->avr_channels.size() != with->_ys_sum.size())||(with->avr_channels.size() != with->_ys_disp.size()))
			valid = false;
		if (!empty && (avr_channels.size() != with->avr_channels.size())){
			std::cout << "WARNING Two AllRunsResults have averaging channels' size mismatches: not Merging" << std::endl;
			return;
		}
		if (empty) {
			avr_channels = with->avr_channels;
			_xs_sum = with->_xs_sum;
			_ys_sum = with->_ys_sum;
			_ys_disp = with->_ys_disp;
		} else {
			for (int ch = 0; ch < avr_channels.size(); ++ch) {
				if (_xs_sum[ch].size()!=with->_xs_sum[ch].size()) {
					std::cout << "WARNING Two average signals (ch "<<avr_channels[ch]<<" in different AllRunsResults have different size: not Merging" << std::endl;
					return;
				}
				for (auto i=_ys_sum[ch].begin(), j=with->_ys_sum[ch].begin(), i_end_=_ys_sum[ch].end(), j_end_=with->_ys_sum[ch].end(); (i!=i_end_)&&(j!=j_end_); ++i,++j) {
					*i+=*j;
				}
			}
		}

		empty = false, valid = true;
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
				mppc_peaks[ch].insert(mppc_peaks[ch].end(), with->mppc_peaks[ch].begin(), with->mppc_peaks[ch].end());
				mppc_double_Is[ch].insert(mppc_double_Is[ch].end(), with->mppc_double_Is[ch].begin(), with->mppc_double_Is[ch].end());
			}
		}
	}
	if (2 == Iteration_N) {
		bool empty = false, valid = true;
		if (avr_channels.empty())
			empty = true;
		if ((avr_channels.size() != _xs_sum.size())||(avr_channels.size() != _ys_sum.size())||(avr_channels.size() != _ys_disp.size()))
			valid = false;
		if ((with->avr_channels.size() != with->_xs_sum.size())||(with->avr_channels.size() != with->_ys_sum.size())||(with->avr_channels.size() != with->_ys_disp.size()))
			valid = false;
		if (!empty && (avr_channels.size() != with->avr_channels.size())){
			std::cout << "WARNING Two AllRunsResults have averaging channels' size mismatches: not Merging" << std::endl;
			return;
		}
		if (empty) {
			avr_channels = with->avr_channels;
			_xs_sum = with->_xs_sum;
			_ys_sum = with->_ys_sum;
			_ys_disp = with->_ys_disp;
		} else {
			for (int ch = 0; ch < avr_channels.size(); ++ch) {
				if (_xs_sum[ch].size()!=with->_xs_sum[ch].size()) {
					std::cout << "WARNING Two average signals (ch "<<avr_channels[ch]<<" in different AllRunsResults have different size: not Merging" << std::endl;
					return;
				}
				for (auto i=_ys_disp[ch].begin(), j=with->_ys_disp[ch].begin(), i_end_=_ys_disp[ch].end(), j_end_=with->_ys_disp[ch].end(); (i!=i_end_)&&(j!=j_end_); ++i,++j) {
					*i+=*j;
				}
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
	N_of_runs += with->N_of_runs;
	N_of_valid_runs += with->N_of_valid_runs;
	with->Clear();
}

void AllRunsResults::Merged(void)
{
	std::cout << "Iteration " << Iteration_N << std::endl;
	std::cout << "N of runs " << N_of_runs << std::endl;
	std::cout << "N of valid runs " << N_of_valid_runs << std::endl;
	if (0==Iteration_N) {//only PMT are processed. No average signals yet
		if (!_Ss.empty())
			find_S_cutoff();
		for (std::size_t ch=0; ch<pmt_channels.size();++ch) {
			std::string PMT_output_prefix = std::string(ParameterPile::this_path) + "/" + std::string(OUTPUT_DIR) + OUTPUT_PMTS + _exp.experiments.back()
						+ "/" + "PMT_" + std::to_string(pmt_channels[ch]) + "/" + "PMT_" + std::to_string(pmt_channels[ch]) +"_";
			vector_to_file(pmt_peaks[ch], PMT_output_prefix + "peaks.dat");
			vector_to_file(pmt_S2_integral[ch], PMT_output_prefix + "S2_int.dat");
			if (0==pmt_channels[ch]) {
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
				for (auto i= pmt_peaks[ch].begin(),_end_=pmt_peaks[ch].end();i!=_end_;++i){
					double val=0;
					for (auto j= i->begin(),_end1_=i->end();j!=_end1_;++j)
						val+= j->S > 0 ? j->S : 0;
					S_tot_max = std::max(S_tot_max,val);
				}
				TH1D *hist_S_tot = new TH1D("3PMT_S_tot_peaks", "3PMT_S_tot_peaks", 60, 0, S_tot_max);

				for (auto i= pmt_peaks[ch].begin(),_end_=pmt_peaks[ch].end();i!=_end_;++i){
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
			if (1==pmt_channels[ch]) {
				double S2_max =0;
				for (auto i= pmt_peaks[ch].begin(),_end_=pmt_peaks[ch].end();i!=_end_;++i){
					double val=0;
					for (auto j= i->begin(),_end1_=i->end();j!=_end1_;++j){
						if ((j->t>=ParameterPile::S2_start_time.find(_exp.experiments.back())->second)&&(j->t<=ParameterPile::S2_finish_time.find(_exp.experiments.back())->second))
							val+= j->S > 0 ? j->S : 0;
					}
					S2_max = std::max(S2_max,val);
				}
				TH1D *hist_S2 = new TH1D("PMT#1_S_tot_peaks", "PMT#1_S2_tot_peaks", 60, 0, S2_max);
				for (auto i= pmt_peaks[ch].begin(),_end_=pmt_peaks[ch].end();i!=_end_;++i){
					double val=0;
					for (auto j= i->begin(),_end1_=i->end();j!=_end1_;++j){
						if ((j->t>=ParameterPile::S2_start_time.find(_exp.experiments.back())->second)&&(j->t<=ParameterPile::S2_finish_time.find(_exp.experiments.back())->second))
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

				double S_tot_max = 0;
				for (auto i= pmt_peaks[ch].begin(),_end_=pmt_peaks[ch].end();i!=_end_;++i){
					double val=0;
					for (auto j= i->begin(),_end1_=i->end();j!=_end1_;++j)
						val+= j->S > 0 ? j->S : 0;
					S_tot_max = std::max(S_tot_max,val);
				}
				TH1D *hist_S_tot = new TH1D("PMT#1_S_tot_peaks", "PMT#1_S_tot_peaks", 60, 0, S_tot_max);

				for (auto i= pmt_peaks[ch].begin(),_end_=pmt_peaks[ch].end();i!=_end_;++i){
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
	}
	if(1==Iteration_N) {
		for (std::size_t ch = 0, ch_end_=avr_channels.size(); ch!=ch_end_; ++ch)
			for (auto i = _ys_sum[ch].begin(), i_end_ = _ys_sum[ch].end(); i!=i_end_; ++i)
				*i/=(double)N_of_valid_runs;
		bool valid = true;
		if ((mppc_peaks_in_S2_area.size() != mppc_S2_start_time.size()) || (mppc_peaks_in_S2_area.size() != mppc_S2_finish_time.size())
				|| (mppc_peaks_in_S2_area.size()!=mppc_channels.size()) || (mppc_peaks_in_S2_area.size() != mppc_double_Is.size())
				|| (mppc_peaks_in_S2_area.size() != mppc_peaks.size()))
			valid = false;
		if (valid && !(mppc_peaks_in_S2_area.empty())) {
			for (int ch = 0; ch < mppc_peaks_in_S2_area.size(); ++ch) {
				std::string Ss_name = _exp.experiments.back()+"_MPPC#" + std::to_string(mppc_channels[ch]) + "_peaks_S";
				std::string S2_S_name = _exp.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_S";
				std::string S2_start_t_name = _exp.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_start_t";
				std::string S2_finish_t_name = _exp.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_finish_t";
				std::string double_I_name = _exp.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_double_I";
				
				//TODO: ParameterPile
				TH1D *hist_S = createMPPCHist_peaks_S(mppc_peaks[ch], Ss_name, 0, 4.0, 100);
				TH1D *hist_S2_S = createMPPCHist(mppc_peaks_in_S2_area[ch], S2_S_name, 0, 4.0, 100);
				TH1D *hist_S2_start_t = createMPPCHist(mppc_S2_start_time[ch], S2_start_t_name, ParameterPile::S2_start_time.find(_exp.experiments.back())->second, 4.0, 100);
				TH1D *hist_S2_finish_t = createMPPCHist(mppc_S2_finish_time[ch], S2_finish_t_name, ParameterPile::S2_start_time.find(_exp.experiments.back())->second, 6.0, 100);
				TH1D *hist_double_I = createMPPCHist(mppc_double_Is[ch], double_I_name, -1.1, 6.0, 100);

				TCanvas *c1 = new TCanvas(Ss_name.c_str(), Ss_name.c_str());
				c1->cd();
				if (hist_S)
					hist_S->Draw();
				c1->Update();

				TCanvas *c2 = new TCanvas(S2_S_name.c_str(), S2_S_name.c_str());
				c2->cd();
				if (hist_S2_S)
					hist_S2_S->Draw();
				c2->Update();

				TCanvas *c3 = new TCanvas(S2_start_t_name.c_str(), S2_start_t_name.c_str());
				c3->cd();
				if (hist_S2_start_t)
					hist_S2_start_t->Draw();
				c3->Update();

				TCanvas *c4 = new TCanvas(S2_finish_t_name.c_str(), S2_finish_t_name.c_str());
				c4->cd();
				if (hist_S2_finish_t)
					hist_S2_finish_t->Draw();
				c4->Update();

				TCanvas *c5 = new TCanvas(double_I_name.c_str(), double_I_name.c_str());
				c5->cd();
				if (hist_double_I)
					hist_double_I->Draw();
				c5->Update();
				std::string output_prefix =std::string(ParameterPile::this_path)+"/"+ std::string(OUTPUT_DIR) + OUTPUT_MPPCS_PICS + _exp.experiments.back()
					+ "/" + OUTPUT_MPPCS + std::to_string(mppc_channels[ch]) + "/" + OUTPUT_MPPCS + std::to_string(mppc_channels[ch]) + "_";
				vector_to_file(mppc_peaks_in_S2_area[ch], output_prefix + "S2_S.dat","MPPC_S2");
				vector_to_file(mppc_S2_start_time[ch], output_prefix + "S2_start_t.dat", "MPPC_st");
				vector_to_file(mppc_S2_finish_time[ch], output_prefix + "S2_finish_t.dat", "MPPC_fin");
				vector_to_file(mppc_double_Is[ch], output_prefix + "double_I.dat", "MPPC_II");
				vector_to_file(mppc_peaks[ch], output_prefix + "peaks.dat","MPPC_peaks");
#ifdef OUTPUT_MPPCS_PICS
				c1->SaveAs((output_prefix + "Ss.png").c_str(),"png");
				c2->SaveAs((output_prefix + "S2_S.png").c_str(),"png");
				c3->SaveAs((output_prefix + "S2_stat_t.png").c_str(),"png");
				c4->SaveAs((output_prefix + "S2_finish_t.png").c_str(),"png");
				c5->SaveAs((output_prefix + "double_I.png").c_str(),"png");
#endif
				c1->Close();
				c2->Close();
				c3->Close();
				c4->Close();
				c5->Close();
			}
		}
	}
	if(2==Iteration_N) {
		for (std::size_t ind = 0, ind_end_ = avr_channels.size(); ind!=ind_end_; ++ind) {
			int ch = avr_channels[ind];
			for (auto i = _ys_disp[ind].begin(), i_end_ = _ys_disp[ind].end(); i!=i_end_; ++i) {
				*i/=(double)N_of_valid_runs*(N_of_valid_runs-1);
				*i = sqrt(*i);
			}
			if (2==ch) {
				std::string GEM_output_prefix = std::string(ParameterPile::this_path)+"/" + std::string(OUTPUT_DIR) + OUTPUT_GEMS + "/" + OUTPUT_GEMS +"_"
						+ _exp.experiments.back()+"_AVR";
				STD_CONT<peak> no_peaks;
				no_peaks.push_back(peak());
				no_peaks.back().left = ParameterPile::S1_start_time;
				no_peaks.back().right = _xs_sum[ind].back();
				double baseline = SignalOperations::find_baseline_by_median(0, _xs_sum[ind], _ys_sum[ind], no_peaks);
				SignalOperations::substract_baseline(_ys_sum[ind], baseline);

				Drawing *dr = graph_manager.GetDrawing("GEM_" + _exp.experiments.back(), 2, ParameterPile::DrawEngine::Gnuplot);
				dr->AddToDraw(_xs_sum[ind], _ys_sum[ind], "GEM_" + _exp.experiments.back(), "", 0);
				DVECTOR GEM_int;
				SignalOperations::integrate(_xs_sum[ind], _ys_sum[ind], GEM_int,0/*baseline*/);
				dr->AddToDraw(_xs_sum[ind], GEM_int, "GEM_Int_" + _exp.experiments.back(), "axes x1y2", 0);

				//find start GEM time
				DITERATOR x_start = _xs_sum[ind].begin();
				find_GEM_start_time(_xs_sum[ind], _ys_sum[ind], x_start, ParameterPile::GEM_N_of_averaging, graph_manager);
				if (x_start != _xs_sum[ind].end()) {
					dr->AddToDraw_vertical(*x_start, -1, 1, "lc rgb \"#FF0000\"", 0);
					//find finish GEM time
					DITERATOR x_finish = _xs_sum[ind].begin();
					double temp_y_max;
					DVECTOR temp_xs = _xs_sum[ind], temp_ys_I = GEM_int;
					SignalOperations::apply_time_limits(temp_xs, temp_ys_I, *x_start, _xs_sum[ind].back());
					SignalOperations::get_max(temp_xs, temp_ys_I, x_finish, temp_y_max, ParameterPile::GEM_N_of_averaging);
					if (x_finish != temp_xs.end())
						dr->AddToDraw_vertical(*x_finish, -1, 1, "lc rgb \"#FF0000\"", 0);

					if (x_finish != temp_xs.end() && x_start != _xs_sum[ind].end()){
						DITERATOR _begin_ =_xs_sum[ind].begin();
						DITERATOR _end_ = _xs_sum[ind].end() - 1;
						x_finish = SignalOperations::find_x_iterator_by_value(_begin_,_end_, *x_finish);
						double Integ = *(GEM_int.begin() + (x_finish - _xs_sum[ind].begin())) - *(GEM_int.begin() + (x_start - _xs_sum[ind].begin()));
						std::cout << "GEM integral is " << Integ << std::endl;
						std::cout << "GEM time limits are: [" << *x_start << ";" << *x_finish << "]" << std::endl;
					} else {
						std::cout << "GEM time boundaries are invalid" << std::endl;
					}
				} else {
					std::cout << "GEM x-y size mismatch" << std::endl;
				}
				vector_to_file(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], GEM_output_prefix+".txt", std::string(OUTPUT_GEMS) +"_" + _exp.experiments.back()+"_AVR");
				//TODO: save integral as well + its error.
				continue;
			}
			if (ch>=32) {
				std::string MPPC_output_prefix = std::string(ParameterPile::this_path)+"/"+ std::string(OUTPUT_DIR) + OUTPUT_MPPCS_PICS + _exp.experiments.back()
						+ "_ch_" + std::to_string(ch) + "_AVR";
				STD_CONT<peak> no_peaks;
				no_peaks.push_back(peak());
				no_peaks.back().left = ParameterPile::S1_start_time;
				no_peaks.back().right = _xs_sum[ind].back();
				double baseline = SignalOperations::find_baseline_by_integral(0, _xs_sum[ind], _ys_sum[ind], no_peaks);
				SignalOperations::substract_baseline(_ys_sum[ind], baseline);
				vector_to_file(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], MPPC_output_prefix+".txt", std::string(OUTPUT_MPPCS_PICS) +"_" + _exp.experiments.back()+"_AVR");
				continue;
			}
			std::string PMT_output_prefix = std::string(ParameterPile::this_path) + "/" + std::string(OUTPUT_DIR) + OUTPUT_PMTS + _exp.experiments.back()
				+ "_ch_" + std::to_string(ch) + "_AVR";
			STD_CONT<peak> no_peaks;
			no_peaks.push_back(peak());
			no_peaks.back().left = ParameterPile::S1_start_time;
			no_peaks.back().right = _xs_sum[ind].back();
			double baseline = SignalOperations::find_baseline_by_integral(0, _xs_sum[ind], _ys_sum[ind], no_peaks);
			SignalOperations::substract_baseline(_ys_sum[ind], baseline);
			vector_to_file(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], PMT_output_prefix+".txt", std::string(OUTPUT_PMTS) +"_" + _exp.experiments.back()+"_AVR");
		}
		double S2_st = ParameterPile::S2_start_time.find(_exp.experiments.back())->second;
		double S2_ft = ParameterPile::S2_finish_time.find(_exp.experiments.back())->second;
		if (!_xs_PMT3_sum.empty()){
			for (auto i= _ys_PMT3_sum.begin(),_end_=_ys_PMT3_sum.end();i!=_end_;++i)
				*i/=N_of_valid_runs;
			DVECTOR xs_before_S1 = _xs_PMT3_sum, ys_before_S1=_ys_PMT3_sum;
			SignalOperations::apply_time_limits(xs_before_S1,ys_before_S1,*xs_before_S1.begin(),ParameterPile::S1_start_time);
			double baseline = SignalOperations::find_baseline_by_median(0,xs_before_S1,ys_before_S1);
			SignalOperations::substract_baseline(_ys_PMT3_sum, baseline);
			SignalOperations::integrate(_xs_PMT3_sum,_ys_PMT3_sum,ys_before_S1, 0);
			Drawing* dr = graph_manager.GetDrawing("3PMT_"+_exp.experiments.back()+"\\_AVR\\_",0,ParameterPile::DrawEngine::Gnuplot);
			dr->AddToDraw(_xs_PMT3_sum,_ys_PMT3_sum,"3PMT average signal "+_exp.experiments.back());
			dr->AddToDraw(_xs_PMT3_sum,ys_before_S1,"3PMT I of average signal "+_exp.experiments.back(),"axes x1y2");
			dr->AddToDraw_vertical(S2_st, -1, 1, "lc rgb \"#0000FF\"");
			dr->AddToDraw_vertical(S2_ft, -1, 1, "lc rgb \"#0000FF\"");
		}
		if (!_xs_PMT1_sum.empty()){
			for (auto i= _ys_PMT1_sum.begin(),_end_=_ys_PMT1_sum.end();i!=_end_;++i)
				*i/=N_of_valid_runs;
			DVECTOR xs_before_S1 = _xs_PMT1_sum, ys_before_S1=_ys_PMT1_sum;
			SignalOperations::apply_time_limits(xs_before_S1,ys_before_S1,*xs_before_S1.begin(),ParameterPile::S1_start_time);
			double baseline = SignalOperations::find_baseline_by_median(0,xs_before_S1,ys_before_S1);
			SignalOperations::substract_baseline(_ys_PMT1_sum, baseline);
			SignalOperations::integrate(_xs_PMT1_sum,_ys_PMT1_sum,ys_before_S1, 0);
			Drawing* dr = graph_manager.GetDrawing("PMT#1_"+_exp.experiments.back()+"\\\_AVR\\\_",1,ParameterPile::DrawEngine::Gnuplot);
			dr->AddToDraw(_xs_PMT1_sum,_ys_PMT1_sum,"PMT#1 average signal "+_exp.experiments.back());
			dr->AddToDraw(_xs_PMT1_sum,ys_before_S1,"PMT#1 I of average signal "+_exp.experiments.back(),"axes x1y2");
			dr->AddToDraw_vertical(S2_st, -1, 1, "lc rgb \"#0000FF\"");
			dr->AddToDraw_vertical(S2_ft, -1, 1, "lc rgb \"#0000FF\"");
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
	graph_manager.Draw();
	Clear();
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
			str.write((char*)&pp->t,sizeof(double));
		}
	}
	str.close();
}

void AllRunsResults::vector_to_file(DVECTOR &what, std::string fname, std::string title)
{
	//ROOT's TTree fails for unknown reason (not dictionary-related I guess). Crashes at Fill() (libRIO.dll) but not 100% of the time,
	//so it is somehow has undefined behavior. gDebug results in no useful info. Plus the files that were written were significantly
	//larger than mine, so I opted against TTree usage.

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

void vector_to_file(DVECTOR &xs, DVECTOR &ys, DVECTOR &ydisps, std::string fname, std::string title)
{
	std::ofstream str;
	open_output_file(fname, str, std::ios_base::trunc);
	str<<"//"<<title<<std::endl;
	str<<"//t [us]\tsignal\terror"<<std::endl;
	for (auto x = xs.begin(), y = ys.begin(), d = ydisps.begin(), x_end_ = xs.end(), y_end_=ys.end(), d_end_ = ydisps.end(); (x != x_end_)&&(y!=y_end_)&&(d!=d_end_); ++x,++y,++d)
		str<<*x<<'\t'<<str<<*y<<'\t'<<*d<<std::endl;
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
	return hist;
}

int AllRunsResults::Iteration(void) const
{	return Iteration_N;}

void AllRunsResults::Clear(void)
{
	N_of_runs = 0;
	N_of_valid_runs = 0;
	//N_peaks_cutoff preserve
	//S_peaks_cutoff preserve
	//Iteration_N preserve;
	graph_manager.Clear();
#ifdef _HOTFIX_CLEAR_MEMORY
	DVECTOR().swap(_Ss);
	DVECTOR().swap(_ns);
	if (1!=Iteration()) {
		STD_CONT<DVECTOR>().swap(_xs_sum);
		STD_CONT<DVECTOR>().swap(_ys_sum);
		STD_CONT<DVECTOR>().swap(_ys_disp);
		STD_CONT<int>().swap(avr_channels);
	}
	STD_CONT<DVECTOR>().swap(mppc_peaks_in_S2_area);
	STD_CONT<DVECTOR>().swap(mppc_S2_start_time);
	STD_CONT<DVECTOR>().swap(mppc_S2_finish_time);
	STD_CONT<STD_CONT<STD_CONT<peak>>>().swap(mppc_peaks);
	STD_CONT<DVECTOR>().swap(mppc_double_Is);
	STD_CONT<int>().swap(mppc_channels);
#else
	_Ss.clear();
	_ns.clear();
	if (1!=Iteration()) {
		_xs_sum.clear();
		_ys_sum.clear();
		_ys_disp.clear();
		avr_channels.clear();
	}
	mppc_peaks_in_S2_area.clear();
	mppc_S2_start_time.clear();
	mppc_S2_finish_time.clear();
	//mppc_all_peaks_Ss.clear();
	mppc_peaks.clear();
	mppc_double_Is.clear();
	mppc_channels.clear();
#endif
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
