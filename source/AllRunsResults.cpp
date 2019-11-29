#include "AllRunsResults.h"

AllRunsResults::AllRunsResults(ParameterPile::experiment_area experiment)
{
	_exp = experiment;
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

void AllRunsResults::find_GEM_start_time(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start, int N_trust, GraphCollection &man)
{
	DVECTOR xs_before_S1 = xs, ys_before_S1 = ys;
	SignalOperations::apply_time_limits(xs_before_S1, ys_before_S1, *(xs.begin()), ParameterPile::S1_start_time);
	DITERATOR x_befS1_max;
	double noize_amp;
	SignalOperations::get_max(xs_before_S1, ys_before_S1, x_befS1_max, noize_amp, N_trust);
	GnuplotDrawing *dr = man.GetDrawing(0);
	if (dr)
		dr->AddToDraw_vertical(*x_befS1_max, -0.4, 0.4, "lc rgb \"#000000\"");
	noize_amp *= ParameterPile::GEM_threshold_to_noise;
	if (dr)
		dr->AddToDraw_baseline(noize_amp, "threshold", "lc rgb \"#000000\"");
	DITERATOR x_S1_left = xs.begin();
	DITERATOR x_S1_right = xs.begin();
	SignalOperations::find_next_peak(xs, ys, x_S1_left, x_S1_right, noize_amp, N_trust);
	if (x_S1_left == xs.end()) {
		std::cout << "1st peak not found, threshold " << noize_amp << std::endl;
		x_start = xs.end();
		return;
	}
	if (dr)
		dr->AddToDraw_vertical(*x_S1_right, -0.4, 0.4, "lc rgb \"#00AA00\"");
	//x_S1_right += N_trust;
	SignalOperations::find_next_extremum(xs, ys, x_S1_right, N_trust);
	if (dr)
		dr->AddToDraw_vertical(*x_S1_right, -0.4, 0.4, "lc rgb \"#0000AA\"");
	x_S1_right += N_trust;
	SignalOperations::find_next_extremum(xs, ys, x_S1_right, N_trust);
	if (dr)
		dr->AddToDraw_vertical(*x_S1_right, -0.4, 0.4, "lc rgb \"#00AAAA\"");
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

void AllRunsResults::Merge(AllRunsResults* with)
{
	_valid.insert(_valid.end(), with->_valid.begin(), with->_valid.end());
	_status.insert(_status.end(), with->_status.begin(), with->_status.end());
	N_of_runs += with->N_of_runs;
	N_of_valid_runs += with->N_of_valid_runs;
	if (pictures.empty()) {
		pictures = with->pictures;
	} else {
		for (std::size_t j = 0, j_end_ = with->pictures.size(); j!= j_end_; ++j) {
			GraphCollection *pics = pictures.info(with->pictures.index(j));
			if (NULL == pics) {
				pictures.push(with->pictures.index(j), with->pictures[j]);
				break;
			} else {
				pics->Merge(with->pictures[j]);
				break;
			}
		}
	}
	if (pictures.isSameIndices(with->pictures)) {
		for (std::size_t i = 0, i_end_ = pictures.size(); i!= i_end_; ++i)
			pictures[i].Merge(with->pictures[i]);
	} else {
		if (pictures.empty()) {
			pictures = with->pictures;
		} else {
			if (!with->pictures.empty())
				std::cout << "AllRunsResults::Merge: Warning: Two AllRunsResults picture collections' channel mismatch: not merging pictures" << std::endl;
		}
	}
	if (0 == Iteration_N) {
		_Ss.insert(_Ss.end(), with->_Ss.begin(), with->_Ss.end());
		_ns.insert(_ns.end(), with->_ns.begin(), with->_ns.end());
		bool empty = false, valid = true;
		if (pmt_channels.empty())
			empty = true;
		valid = isSameChannels(pmt_channels, with->pmt_channels);
		valid = valid && isSameChannels(pmt_integrated_channels, with->pmt_integrated_channels);
		if (!empty && !valid){
			std::cout << "WARNING Two AllRunsResults have PMT channels' size mismatches: not Merging" << std::endl;
			return;
		}
		if (empty) {
			pmt_channels = with->pmt_channels;
			pmt_peaks = with->pmt_peaks;
			pmt_S2_integral = with->pmt_S2_integral;
			pmt_integrated_channels = with->pmt_integrated_channels;
		} else {
			pmt_peaks.insert(pmt_peaks.end(), with->pmt_peaks.begin(), with->pmt_peaks.end());
			pmt_S2_integral.insert(pmt_S2_integral.end(), with->pmt_S2_integral.begin(), with->pmt_S2_integral.end());
		}
	}
	if (1 == Iteration_N) {
		bool empty = false, valid = true;
		if (avr_channels.empty())
			empty = true;
		valid = isSameChannels(avr_channels, with->avr_channels);
		if (!empty && !valid){
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
		valid = isSameChannels(mppc_channels, with->mppc_channels);
		if (mppc_peaks.empty())
			empty = true;
		if (!empty && !valid){
			std::cout << "WARNING Two AllRunsResults have MPPC channels' size mismatches: not Merging" << std::endl;
			return;
		}
		if (empty) {
			mppc_channels = with->mppc_channels;
			mppc_double_Is = with->mppc_double_Is;
			mppc_peaks = with->mppc_peaks;
		} else {
			mppc_peaks.insert(mppc_peaks.end(), with->mppc_peaks.begin(), with->mppc_peaks.end());
			mppc_double_Is.insert(mppc_double_Is.end(), with->mppc_double_Is.begin(), with->mppc_double_Is.end());
		}
	}
	if (2 == Iteration_N) {
		bool empty = false, valid = true;
		if (avr_channels.empty())
			empty = true;
		valid = isSameChannels(avr_channels, with->avr_channels);
		if (!empty && !valid){
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
				if (_ys_disp[ch].size()!=with->_ys_disp[ch].size()) {
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
		time_stat.t_PMT_proc += with->time_stat.t_PMT_proc;
		time_stat.n_PMT_proc += with->time_stat.n_PMT_proc
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
	if(1==Iteration_N) {
		for (std::size_t ch = 0, ch_end_=avr_channels.size(); ch!=ch_end_; ++ch)
			for (auto i = _ys_sum[ch].begin(), i_end_ = _ys_sum[ch].end(); i!=i_end_; ++i)
				*i/=(double)N_of_valid_runs;
	}

	if (2==Iteration_N) {//Save PMT only at the end (with final _valid vector)
		for (std::size_t ch=0; ch<pmt_channels.size(); ++ch) {
			std::stringstream ch_str;
			ch_str<<pmt_channels[ch];
			int ch_int = getIndex(pmt_integrated_channels, pmt_channels[ch]);
			std::string PMT_output_prefix = std::string(ParameterPile::this_path) + std::string(OUTPUT_DIR) + OUTPUT_PMTS + _exp.experiments.back()
						+ "/" + "PMT_" + ch_str.str() + "/" + "PMT_" + ch_str.str() +"_";
			if (!ParameterPile::draw_only)
				vector_to_file(pmt_peaks, ch, PMT_output_prefix + "peaks.dat");
			if (ch_int >= 0) {
				if (!ParameterPile::draw_only)
					vector_to_file(pmt_S2_integral, ch_int, PMT_output_prefix + "S2_int.dat");
				double S2int_min = 0, S2int_max = 0;
				for (std::size_t run = 0, run_end_ = pmt_S2_integral.size(); run!=run_end_; ++run) {
					if (_valid[run]) {
						S2int_min = std::min(pmt_S2_integral[run][ch_int], S2int_min);
						S2int_max = std::max(pmt_S2_integral[run][ch_int], S2int_max);
					}
				}
				unsigned int N_bins = 80;
				S2int_max += (S2int_max-S2int_min)/N_bins;
				std::string hist_name = "PMT"+ch_str.str()+"_S2_integral " + _exp.experiments.back();
				TH1D *hist_S2int = new TH1D(hist_name.c_str(), hist_name.c_str(), N_bins, S2int_min, S2int_max);
				for (std::size_t run = 0, run_end_ = pmt_S2_integral.size(); run!=run_end_; ++run) {
					if (_valid[run]) {
						hist_S2int->Fill(pmt_S2_integral[run][ch_int]);
					}
				}
				std::string can_name = "PMT"+ch_str.str()+"_S2_integral_distribution " + _exp.experiments.back();
				TCanvas *c1 = new TCanvas(can_name.c_str(), can_name.c_str());
				c1->cd();
				hist_S2int->Draw();
				c1->Update();
				if (!ParameterPile::draw_only)
					c1->SaveAs((PMT_output_prefix+"S2_int.png").c_str(),"png");
			}
			if (0==pmt_channels[ch]) {
				DITERATOR S2_max = std::max_element(_Ss.begin(), _Ss.end());
				double hist_S2_max = *S2_max;
				TH1D *hist_S2 = new TH1D("3PMT_S2_peaks", "3PMT_S2_peaks", 60, 0, hist_S2_max);

				for (std::size_t _n = 0, _s = 0; (_n < _ns.size()) && (_s < _Ss.size()); ++_s, ++_n) {
					if (!_valid[_s])
						continue; //do not draw if maximum limit is imposed && element exceeds it
					hist_S2->Fill(_Ss[_s]);
				}

				TCanvas *c1 = new TCanvas(("3PMT_S2_peaks_distribution " + _exp.experiments.back()).c_str(),
					("3PMT_S2_peaks_distribution " + _exp.experiments.back()).c_str());
				c1->cd();
				hist_S2->Draw();
				c1->Update();
				if (!ParameterPile::draw_only)
					c1->SaveAs((PMT_output_prefix+"S2s.png").c_str(),"png");

				double S_tot_max =0;
				for (std::size_t i=0, i_end_ = pmt_peaks.size(); i!=i_end_; ++i) {
					if (!_valid[i])
						continue;
					double val=0;
					for (auto j= pmt_peaks[i][ch].begin(),_end1_= pmt_peaks[i][ch].end(); j!=_end1_; ++j)
						val+= j->S > 0 ? j->S : 0;
					S_tot_max = std::max(S_tot_max, val);
				}
				TH1D *hist_S_tot = new TH1D("3PMT_S_tot_peaks", "3PMT_S_tot_peaks", 60, 0, S_tot_max);

				for (std::size_t i=0, i_end_ = pmt_peaks.size();i!=i_end_;++i) {
					if (!_valid[i])
						continue;
					double val=0;
					for (auto j= pmt_peaks[i][ch].begin(), _end1_=pmt_peaks[i][ch].end(); j!=_end1_; ++j)
						val+= j->S > 0 ? j->S : 0;
					hist_S_tot->Fill(val);
				}
				TCanvas *c2 = new TCanvas(("3PMT_S_total_peaks_distribution " + _exp.experiments.back()).c_str(),
					("3PMT_S_total_peaks_distribution " + _exp.experiments.back()).c_str());
				c2->cd();
				hist_S_tot->Draw();
				c2->Update();
				if (!ParameterPile::draw_only)
					c2->SaveAs((PMT_output_prefix+"S_tot.png").c_str(),"png");
			}
			if (1==pmt_channels[ch]) {
				double S2_max =0;
				for (std::size_t i= 0, i_end_ = pmt_peaks.size(); i!=i_end_; ++i) {
					if (!_valid[i])
						continue;
					double val=0;
					for (auto j= pmt_peaks[i][ch].begin(),_end1_=pmt_peaks[i][ch].end(); j!=_end1_; ++j) {
						if ((j->t>=ParameterPile::S2_start_time.find(_exp.experiments.back())->second)&&(j->t<=ParameterPile::S2_finish_time.find(_exp.experiments.back())->second))
							val+= j->S > 0 ? j->S : 0;
					}
					S2_max = std::max(S2_max,val);
				}
				TH1D *hist_S2 = new TH1D("PMT#1_S_tot_peaks", "PMT#1_S2_tot_peaks", 60, 0, S2_max);
				for (std::size_t i= 0, i_end_ = pmt_peaks.size(); i!=i_end_; ++i) {
					if (!_valid[i])
						continue;
					double val=0;
					for (auto j= pmt_peaks[i][ch].begin(),_end1_=pmt_peaks[i][ch].end(); j!=_end1_; ++j) {
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
				if (!ParameterPile::draw_only)
					c1->SaveAs((PMT_output_prefix+"S2.png").c_str(),"png");

				double S_tot_max = 0;
				for (std::size_t i= 0, i_end_ = pmt_peaks.size(); i!=i_end_; ++i) {
					if (!_valid[i])
						continue;
					double val=0;
					for (auto j= pmt_peaks[i][ch].begin(),_end1_=pmt_peaks[i][ch].end(); j!=_end1_; ++j)
						val+= j->S > 0 ? j->S : 0;
					S_tot_max = std::max(S_tot_max,val);
				}
				TH1D *hist_S_tot = new TH1D("PMT#1_S_tot_peaks", "PMT#1_S_tot_peaks", 60, 0, S_tot_max);

				for (std::size_t i= 0, i_end_ = pmt_peaks.size(); i!=i_end_; ++i) {
					if (!_valid[i])
						continue;
					double val=0;
					for (auto j= pmt_peaks[i][ch].begin(),_end1_=pmt_peaks[i][ch].end(); j!=_end1_; ++j)
						val+= j->S > 0 ? j->S : 0;
					hist_S_tot->Fill(val);
				}
				TCanvas *c2 = new TCanvas(("PMT#1_S_total_peaks_distribution " + _exp.experiments.back()).c_str(),
					("PMT#1_S_total_peaks_distribution " + _exp.experiments.back()).c_str());
				c2->cd();
				hist_S_tot->Draw();
				c2->Update();
				if (!ParameterPile::draw_only)
					c2->SaveAs((PMT_output_prefix+"S_tot.png").c_str(),"png");
			}
		}
	}
	if(2==Iteration_N) { //Save MPPC only at the end (with final _valid vector)
		for (int ch = 0; ch < mppc_channels.size(); ++ch) {
			std::string Ss_name = _exp.experiments.back()+"_MPPC#" + std::to_string(mppc_channels[ch]) + "_peaks_S";
			std::string S2_S_name = _exp.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_S";
			std::string S2_start_t_name = _exp.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_start_t";
			std::string S2_finish_t_name = _exp.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_S2_finish_t";
			std::string double_I_name = _exp.experiments.back() + "_MPPC#" + std::to_string(mppc_channels[ch]) + "_double_I";

			//TODO: ParameterPile
			TH1D *hist_S = createMPPCHist_peaks_S(mppc_peaks, ch, Ss_name, 0, 4.0, 100);
			TH1D *hist_double_I = createMPPCHist(mppc_double_Is, ch, double_I_name, -1.1, 6.0, 100);

			TCanvas *c1 = new TCanvas(Ss_name.c_str(), Ss_name.c_str());
			c1->cd();
			if (hist_S)
				hist_S->Draw();
			c1->Update();

			TCanvas *c5 = new TCanvas(double_I_name.c_str(), double_I_name.c_str());
			c5->cd();
			if (hist_double_I)
				hist_double_I->Draw();
			c5->Update();
			std::string output_prefix =std::string(ParameterPile::this_path)+ std::string(OUTPUT_DIR) + OUTPUT_MPPCS_PICS + _exp.experiments.back()
				+ "/" + OUTPUT_MPPCS + std::to_string(mppc_channels[ch]) + "/" + OUTPUT_MPPCS + std::to_string(mppc_channels[ch]) + "_";
			if (!ParameterPile::draw_only) {
				vector_to_file(mppc_double_Is, ch, output_prefix + "double_I.dat", "MPPC_II");
				vector_to_file(mppc_peaks, ch, output_prefix + "peaks.dat","MPPC_peaks");
			}
	#ifdef OUTPUT_MPPCS_PICS
			if (!ParameterPile::draw_only) {
				c1->SaveAs((output_prefix + "Ss.png").c_str(),"png");
				c5->SaveAs((output_prefix + "double_I.png").c_str(),"png");
			}
	#endif
			c1->Close();
			c5->Close();
		}
	}
	if(2==Iteration_N) {
		for (std::size_t i = 0, i_end_ = pictures.size(); i!=i_end_; ++i) {
			std::pair<double, double> y_range = pictures[i].get_y_limits();
			pictures[i].SetYrange(y_range.first, y_range.second);
			pictures[i].SetXrange(ParameterPile::pics_t_zoom.first, ParameterPile::pics_t_zoom.second);
			if (ParameterPile::pics_trigger_position>=ParameterPile::pics_t_zoom.first && ParameterPile::pics_trigger_position<=ParameterPile::pics_t_zoom.second)
				for (std::size_t p = 0, p_end_ = pictures[i].size(); p!=p_end_; ++p)
					pictures[i].GetDrawing(p)->AddToDraw_vertical(ParameterPile::pics_trigger_position, y_range.first, y_range.second, "lc rgb \"#FF0000\"");
			pictures[i].Draw();
		}
		for (std::size_t ind = 0, ind_end_ = avr_channels.size(); ind!=ind_end_; ++ind) {
			int ch = avr_channels[ind];
			for (auto i = _ys_disp[ind].begin(), i_end_ = _ys_disp[ind].end(); i!=i_end_; ++i) {
				*i/=(double)N_of_valid_runs*(N_of_valid_runs-1);
				*i = sqrt(*i);
			}
			if (GEM_CH_==ch) {
				std::string GEM_output_prefix = std::string(ParameterPile::this_path)+ std::string(OUTPUT_DIR) + OUTPUT_GEMS + "/" + OUTPUT_GEMS +"_"
						+ _exp.experiments.back()+"_AVR";
				STD_CONT<peak> no_peaks;
				no_peaks.push_back(peak());
				no_peaks.back().left = ParameterPile::S1_start_time;
				no_peaks.back().right = _xs_sum[ind].back();
				double baseline = SignalOperations::find_baseline_by_median(0, _xs_sum[ind], _ys_sum[ind], no_peaks);
				SignalOperations::substract_baseline(_ys_sum[ind], baseline);

				GnuplotDrawing *dr = graph_manager.GetDrawing("GEM_" + _exp.experiments.back()+"_AVR");
				dr->SetGnuplotDirectory(OUTPUT_DIR+OUTPUT_GEMS +"/" +_exp.experiments.back() + "/");
				dr->AddToDraw(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], "GEM_" + _exp.experiments.back()+"_AVR", "");
				DVECTOR GEM_int, integral_variance;
				SignalOperations::integrate_with_variance(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], GEM_int, integral_variance, 0/*baseline*/);
				dr->AddToDraw(_xs_sum[ind], GEM_int, integral_variance, "GEM_Int_" + _exp.experiments.back(), "axes x1y2");

				//find start GEM time
				DITERATOR x_start = _xs_sum[ind].begin();
				find_GEM_start_time(_xs_sum[ind], _ys_sum[ind], x_start, ParameterPile::GEM_N_of_averaging, graph_manager);
				if (x_start != _xs_sum[ind].end()) {
					dr->AddToDraw_vertical(*x_start, -1, 1, "lc rgb \"#FF0000\"");
					//find finish GEM time
					DITERATOR x_finish = _xs_sum[ind].begin();
					double temp_y_max;
					DVECTOR temp_xs = _xs_sum[ind], temp_ys_I = GEM_int;
					SignalOperations::apply_time_limits(temp_xs, temp_ys_I, *x_start, _xs_sum[ind].back());
					SignalOperations::get_max(temp_xs, temp_ys_I, x_finish, temp_y_max, ParameterPile::GEM_N_of_averaging);
					if (x_finish != temp_xs.end())
						dr->AddToDraw_vertical(*x_finish, -1, 1, "lc rgb \"#FF0000\"");

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
				if (!ParameterPile::draw_only) {
					vector_to_file(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], GEM_output_prefix+".txt", std::string(OUTPUT_GEMS) +"_" + _exp.experiments.back()+"_AVR");
					vector_to_file(_xs_sum[ind], GEM_int, integral_variance, GEM_output_prefix+"_INT.txt", std::string(OUTPUT_GEMS) +"_" + _exp.experiments.back()+"_AVR_INT");
				}
				continue;
			}
			if (ch>=32) {
				std::string MPPC_output_prefix = std::string(ParameterPile::this_path)+ std::string(OUTPUT_DIR) + OUTPUT_MPPCS_PICS + _exp.experiments.back()
						+ "_ch_" + std::to_string(ch) + "_AVR";
				STD_CONT<peak> no_peaks;
				no_peaks.push_back(peak());
				no_peaks.back().left = ParameterPile::S1_start_time;
				no_peaks.back().right = _xs_sum[ind].back();
				double baseline = SignalOperations::find_baseline_by_integral(0, _xs_sum[ind], _ys_sum[ind], no_peaks);
				SignalOperations::substract_baseline(_ys_sum[ind], baseline);
				DVECTOR integral, integral_variance;
				SignalOperations::integrate_with_variance(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], integral, integral_variance, 0);
				if (!ParameterPile::draw_only) {
					vector_to_file(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], MPPC_output_prefix+".txt", std::string(OUTPUT_MPPCS_PICS) +"_" + _exp.experiments.back()+"_AVR");
					vector_to_file(_xs_sum[ind], integral, integral_variance, MPPC_output_prefix+"_INT.txt", std::string(OUTPUT_MPPCS_PICS) +"_" + _exp.experiments.back()+"_AVR_INT");
				}
				continue;
			}
			std::string PMT_output_prefix = std::string(ParameterPile::this_path) + std::string(OUTPUT_DIR) + OUTPUT_PMTS + _exp.experiments.back()
				+ "_ch_" + std::to_string(ch) + "_AVR";
			STD_CONT<peak> no_peaks;
			no_peaks.push_back(peak());
			no_peaks.back().left = ParameterPile::S1_start_time;
			no_peaks.back().right = _xs_sum[ind].back();
			double baseline = SignalOperations::find_baseline_by_median(0, _xs_sum[ind], _ys_sum[ind], no_peaks);
			SignalOperations::substract_baseline(_ys_sum[ind], baseline);
			DVECTOR integral, integral_variance;
			SignalOperations::integrate_with_variance(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], integral, integral_variance, 0);
			if (!ParameterPile::draw_only) {
				vector_to_file(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], PMT_output_prefix+".txt", std::string(OUTPUT_PMTS) +"_" + _exp.experiments.back()+"_AVR");
				vector_to_file(_xs_sum[ind], integral, integral_variance, PMT_output_prefix+"_INT.txt", std::string(OUTPUT_PMTS) +"_" + _exp.experiments.back()+"_AVR_INT");
			}
			double S2_st = ParameterPile::S2_start_time.find(_exp.experiments.back())->second;
			double S2_ft = ParameterPile::S2_finish_time.find(_exp.experiments.back())->second;
			GnuplotDrawing* dr = graph_manager.GetDrawing("PMT_"+_exp.experiments.back()+"_ch_"+std::to_string(ch)+"_AVR");
			dr->SetGnuplotDirectory(OUTPUT_DIR+OUTPUT_PMTS + _exp.experiments.back() + "/PMT_"  + std::to_string(ch) + "/");
			dr->AddToDraw(_xs_sum[ind], _ys_sum[ind], _ys_disp[ind], "PMT_"+_exp.experiments.back()+"_ch_"+std::to_string(ch)+"_AVR");
			dr->AddToDraw(_xs_sum[ind], integral, integral_variance, "PMT_"+_exp.experiments.back()+"_ch_"+std::to_string(ch)+"_Int_AVR","axes x1y2");
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
	graph_manager.Draw();
	Clear();
	++Iteration_N;
}

void vector_to_file(DVECTOR &xs, DVECTOR &ys, DVECTOR &ydisps, std::string fname, std::string title)
{
	std::ofstream str;
	open_output_file(fname, str, std::ios_base::trunc);
	str<<"//"<<title<<std::endl;
	if ((xs.size()!=ys.size())||(xs.size()!=ydisps.size()))
		str<<"//x-y-y_disp size mismatch"<<std::endl;
	else
		for (std::size_t i = 0, i_end_ = xs.size(); i != i_end_; ++i)
			str<<xs[i]<<"\t"<<ys[i]<<"\t"<<ydisps[i]<<std::endl;
	str.close();
}

void AllRunsResults::vector_to_file(STD_CONT<STD_CONT<STD_CONT<peak> > > &pks, int ch_ind, std::string fname, std::string title)
{
	std::ofstream str;
	open_output_file(fname, str, std::ios_base::trunc | std::ios_base::binary);
	std::size_t real_size= 0;
	for (std::size_t run = 0, run_end_ = pks.size(); run != run_end_; ++run)
		if (_valid[run])
			++real_size;
	str.write((char*)&real_size, sizeof(std::size_t));
	for (std::size_t run = 0, run_end_ = pks.size(); run != run_end_; ++run) {
		if (!_valid[run])
			continue;
		std::size_t pk_end_ = pks[run][ch_ind].size();
		str.write((char*)&pk_end_, sizeof(std::size_t));
		for (std::size_t p = 0; p != pk_end_; ++p) {
			str.write((char*)&pks[run][ch_ind][p].left, sizeof(double));
			str.write((char*)&pks[run][ch_ind][p].right, sizeof(double));
			str.write((char*)&pks[run][ch_ind][p].S, sizeof(double));
			str.write((char*)&pks[run][ch_ind][p].A,sizeof(double));
			str.write((char*)&pks[run][ch_ind][p].t,sizeof(double));
		}
	}
	str.close();
}

void AllRunsResults::vector_to_file(STD_CONT<STD_CONT<double> > &what, int ch_ind, std::string fname, std::string title)
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
	std::size_t real_size= 0;
	for (std::size_t run = 0, run_end_= what.size(); run != run_end_; ++run)
		if (_valid[run])
			++real_size;
	str.write((char*)&real_size, sizeof(std::size_t));
	for (std::size_t run = 0, run_end_= what.size(); run != run_end_; ++run)
		if (_valid[run])
			str.write((char*)&what[run][ch_ind], sizeof(double));
	str.close();
}

void AllRunsResults::vector_to_file(DVECTOR &xs, DVECTOR &ys, DVECTOR &ydisps, std::string fname, std::string title)
{
	std::ofstream str;
	open_output_file(fname, str, std::ios_base::trunc);
	str<<"//"<<title<<std::endl;
	str<<"//t [us]\tsignal\terror"<<std::endl;
	for (auto x = xs.begin(), y = ys.begin(), d = ydisps.begin(), x_end_ = xs.end(), y_end_=ys.end(), d_end_ = ydisps.end(); (x != x_end_)&&(y!=y_end_)&&(d!=d_end_); ++x,++y,++d)
		str<<*x<<'\t'<<*y<<'\t'<<*d<<std::endl;
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

TH1D* AllRunsResults::createMPPCHist_peaks_S(STD_CONT<STD_CONT<STD_CONT<peak> > > &what, int ch_ind, std::string name, double left_cutoff, double right_cutoff_from_RMS, int N_bins)
{
	if (what.empty())
		return NULL;
	double Val_max, Val_min;
	std::size_t first_run = 0;
	for (std::size_t run = 0, run_end_ = what.size(); run!=run_end_; ++run) {
		if (!_valid[run]) {
			++first_run;
			continue;
		}
		auto _end_ = what[run][ch_ind].end();
		STD_CONT<peak>::iterator V_max = std::max_element(what[run][ch_ind].begin(), what[run][ch_ind].end(), [](peak a, peak b) {
			return a.S < b.S;
		});
		STD_CONT<peak>::iterator V_min = std::min_element(what[run][ch_ind].begin(), what[run][ch_ind].end(), [left_cutoff](peak a, peak b) {
			if (a.S <= left_cutoff)
				return false;
			if (b.S <= left_cutoff)
				return true;
			return a.S < b.S;
		});
		if (first_run == run) { //initializing Val_max/min. run may contain no peaks so it must be accounted for
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
	if (first_run == what.size()) {
		Val_max = 10;
		Val_min = 0;
	}
	
	std::function<double(peak&)> picker = [](peak &pk){return pk.S; };
	double V_mean = Mean(what, ch_ind, picker);
	double V_RMS = RMS(what, ch_ind, picker);
	
	double V_left_cutoff = std::max(left_cutoff, Val_min);
	double V_right_cutoff = std::min(V_mean + V_RMS * right_cutoff_from_RMS, Val_max);

	if (N_bins <= 0) {
		N_bins = 0;
		for (std::size_t run = 0, run_end_ = what.size(); run < run_end_; ++run) {
			if (!_valid[run])
				continue;
			for (std::size_t n = 0, n_end_ = what[run][ch_ind].size(); n < n_end_; ++n)
				if ((what[run][ch_ind][n].S >= V_left_cutoff) && (V_right_cutoff > what[run][ch_ind][n].S))
					N_bins++;
		}
		if (N_bins <= 0)
			N_bins = 1;
		N_bins = std::sqrt(N_bins);
	}
	TH1D *hist = new TH1D(name.c_str(), name.c_str(), N_bins, V_left_cutoff, V_right_cutoff);
	hist->SetTitle(name.c_str());
	for (std::size_t run = 0, run_end_ = what.size(); run < run_end_; ++run) {
		if (!_valid[run])
			continue;
		for (std::size_t n = 0, n_end_ = what[run][ch_ind].size(); n < n_end_; ++n)
			if ((what[run][ch_ind][n].S >= V_left_cutoff) && (V_right_cutoff > what[run][ch_ind][n].S))
				hist->Fill(what[run][ch_ind][n].S);
	}
	return hist;
}

TH1D* AllRunsResults::createMPPCHist(STD_CONT<STD_CONT<double> > &what, int ch_ind, std::string name, double left_cutoff, double right_cutoff_from_RMS, int N_bins)
{
	if (what.empty())
		return NULL;
	STD_CONT<STD_CONT<double> >::iterator it_V_max = std::max_element(what.begin(), what.end(), [ch_ind](STD_CONT<double> &a, STD_CONT<double> &b){
		return a[ch_ind] < b[ch_ind];
	});
	STD_CONT<STD_CONT<double> >::iterator it_V_min = std::min_element(what.begin(), what.end(), [ch_ind, left_cutoff](STD_CONT<double> &a, STD_CONT<double> &b) {
		if (a[ch_ind] <= left_cutoff)
			return false;
		if (b[ch_ind] <= left_cutoff)
			return true;
		return a[ch_ind] < b[ch_ind];
	});
	double V_min = (*it_V_min)[ch_ind];
	double V_max = (*it_V_max)[ch_ind];
	double V_mean = Mean(what, ch_ind);
	double V_RMS = RMS(what, ch_ind);
	double V_left_cutoff = std::max(left_cutoff, V_min);
	double V_right_cutoff = std::min(V_mean + V_RMS * right_cutoff_from_RMS, V_max);
	if (0==V_RMS)
		V_right_cutoff = 2*V_mean;
	if (N_bins <= 0) {
		N_bins = 0;
		for (int run = 0, run_end_ = what.size(); run!=run_end_; ++run) {
			if (!_valid[run])
				continue;
			if ((what[run][ch_ind]>=V_left_cutoff) && (V_right_cutoff > what[run][ch_ind]))
				N_bins++;
		}
		if (N_bins <= 0)
			N_bins = 1;
		N_bins = std::sqrt(N_bins);
	}
	TH1D *hist = new TH1D(name.c_str(), name.c_str(), N_bins, V_left_cutoff, V_right_cutoff);
	hist->SetTitle(name.c_str());
	for (int run = 0, run_end_ = what.size(); run!=run_end_; ++run) {
		if (!_valid[run])
			continue;
		if ((what[run][ch_ind]>=V_left_cutoff) && (V_right_cutoff > what[run][ch_ind]))
			hist->Fill(what[run][ch_ind]);
	}
	return hist;
}

double AllRunsResults::Mean(STD_CONT<STD_CONT<STD_CONT<peak> > > &peaks, int ch_ind, std::function<double (peak& pk)> &value_picker)
{
	// Return the weighted mean of an array defined by the iterators.
	double sum = 0;
	double sumw = 0;
	for (std::size_t run = 0, run_end_=peaks.size(); run!=run_end_; ++run) {
		if (!_valid[run])
			continue;
		for (std::size_t pk=0, pk_end_=peaks[run][ch_ind].size(); pk!=pk_end_; ++pk) {
			sum += value_picker(peaks[run][ch_ind][pk]);
			sumw += 1;
		}
	}
	return sum / sumw;
}

double AllRunsResults::RMS(STD_CONT<STD_CONT<STD_CONT<peak> > > &peaks, int ch_ind, std::function<double (peak& pk)> &value_picker)
{
	// Return the Standard Deviation of an array defined by the iterators.
	// Note that this function returns the sigma(standard deviation) and
	// not the root mean square of the array.

	// Use the two pass algorithm, which is slower (! a factor of 2) but much more
	// precise.  Since we have a vector the 2 pass algorithm is still faster than the
	// Welford algorithm. (See also ROOT-5545)

	double n = 0;
	double tot = 0;
	double mean = Mean(peaks, ch_ind, value_picker);
	for (std::size_t run = 0, run_end_=peaks.size(); run!=run_end_; ++run) {
		if (!_valid[run])
			continue;
		for (std::size_t pk=0, pk_end_=peaks[run][ch_ind].size(); pk!=pk_end_; ++pk) {
			double x = value_picker(peaks[run][ch_ind][pk]);
			tot += (x - mean)*(x - mean);
			n += 1;
		}
	}
	double rms = (n > 1) ? TMath::Sqrt(tot / (n - 1)) : 0.0;
	return rms;
}

double AllRunsResults::Mean(STD_CONT<STD_CONT<double> > &vals, int ch_ind)
{
	double sum = 0;
	double sumw = 0;
	for (std::size_t run = 0, run_end_=vals.size(); run!=run_end_; ++run) {
		if (!_valid[run])
			continue;
		sum += vals[run][ch_ind];
		sumw += 1;
	}
	return sum / sumw;
}

double AllRunsResults::RMS(STD_CONT<STD_CONT<double> > &vals, int ch_ind)
{
	double n = 0;
	double tot = 0;
	double mean = Mean(vals, ch_ind);
	for (std::size_t run = 0, run_end_=vals.size(); run!=run_end_; ++run) {
		if (!_valid[run])
			continue;
		double x = vals[run][ch_ind];
		tot += (x - mean)*(x - mean);
		n += 1;
	}
	double rms = (n > 1) ? TMath::Sqrt(tot / (n - 1)) : 0.0;
	return rms;
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

	if (ParameterPile::Max_iteration_N == Iteration()) {
		graph_manager.Clear();
		pictures.clear();
		std::vector<bool>().swap(_valid);
		std::vector<Status>().swap(_status);

		STD_CONT<DVECTOR>().swap(_xs_sum);
		STD_CONT<DVECTOR>().swap(_ys_sum);
		STD_CONT<DVECTOR>().swap(_ys_disp);
		STD_CONT<int>().swap(avr_channels);
		STD_CONT<STD_CONT<STD_CONT<peak>>>().swap(mppc_peaks);
		STD_CONT<STD_CONT<double> >().swap(mppc_double_Is);
		STD_CONT<int>().swap(mppc_channels);

		DVECTOR().swap(_Ss);
		DVECTOR().swap(_ns);
		STD_CONT<STD_CONT<STD_CONT<peak> > >().swap(pmt_peaks);
		STD_CONT<STD_CONT<double> >().swap(pmt_S2_integral);
		STD_CONT<int>().swap(pmt_channels);
	}

#ifdef _USE_TIME_STATISTICS
	time_stat.t_RUN_proc_single_iter=0;
	time_stat.n_RUN_proc_single_iter=0;
#endif
}

void AllRunsResults::ClearMerged(void)
{
	std::vector<bool>().swap(_valid);
	std::vector<Status>().swap(_status);
}

AllRunsResults& AllRunsResults::operator=(const AllRunsResults& right)
{
	//Only these are passed from merged AllRuns to thread-specific ones.
	N_of_runs = right.N_of_runs;
	N_of_valid_runs = right.N_of_valid_runs ;
	Iteration_N = right.Iteration_N;
	_exp = right._exp;
	_ys_sum = right._ys_sum;
	_xs_sum = right._xs_sum;
	_ys_disp = right._ys_disp;
	avr_channels = right.avr_channels;
	return *this;
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
