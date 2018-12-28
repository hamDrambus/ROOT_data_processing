#include "SingleRunData.h"
#include <string>
#include "Savitzky_Golay_filter.h"

SingleRunData::SingleRunData(ParameterPile::experiment_area area)
{
	curr_area.channels = area.channels;
	curr_area.experiments = area.experiments;
	curr_area.runs.push_back(area.runs.back());//only one run is in this class, so no pairs
	curr_area.sub_runs.push_back(area.sub_runs.back()); //only one subrun is in this class

	//curr_area.channels.clear(); - not really required unless break is used in the cycle with following for(;;)
	for (int ch = curr_area.channels.get_next_index(); !(ch<0);ch = curr_area.channels.get_next_index()){
		found_base_lines.push_back(0);
		xs_channels.push_back(DVECTOR());
		ys_channels.push_back(DVECTOR());
	}
	PMT3_summed_peaks_area = 0;
	PMT3_n_peaks = 0;
}

void SingleRunData::readOneRun(AllRunsResults *results, int channel)
{
	int ind = curr_area.channels.get_order_index_by_index(channel);
	if (ind < 0)
		return;
	std::string path = DATA_PREFIX;
	path += curr_area.experiments.back() + "/";
	path += "run_" + std::to_string(curr_area.runs.back()) + "__ch_" + std::to_string(channel) + ".dat";
	file_to_vector(path, xs_channels[ind], ys_channels[ind], curr_area.sub_runs.back());
	if (!xs_channels[ind].empty()&&(AllRunsResults::Status::Empty==results->_status[results->N_of_runs]))
		results->_status[results->N_of_runs] = AllRunsResults::Status::Ok;
}

void SingleRunData::clearOneRun(int channel)
{
	int ind = curr_area.channels.get_order_index_by_index(channel);
	if (ind < 0)
		return;
	DVECTOR().swap(xs_channels[ind]);
	DVECTOR().swap(ys_channels[ind]);
}

void SingleRunData::file_to_vector(std::string fname, DVECTOR &xs, DVECTOR &ys, int index)
{
	DVECTOR().swap(xs);
	DVECTOR().swap(ys);
	std::ifstream file;
	file.open(fname, std::ios_base::binary | std::ios_base::ate);
	if (!file.is_open() || file.eof()){
		file.close();
		return;
	}
	long long Size = file.tellg();
#ifndef _USE_DEQUE
	xs.reserve((Size / ParameterPile::subruns_per_file)/2);
	ys.reserve((Size / ParameterPile::subruns_per_file)/2);
#endif
	long long offset = index * (Size / (ParameterPile::subruns_per_file));
	long long to_read = Size / (ParameterPile::subruns_per_file);
	file.seekg(offset, std::ios_base::beg);
	char bytes[2];
	unsigned long int read_N = 0;
	while (!file.eof() && (read_N<(to_read / 2))){
		file.read(bytes, 2);
		if (!file.eof()) {
			unsigned char b1_ = static_cast<unsigned char>(bytes[1]);
			unsigned char b2_ = static_cast<unsigned char>(bytes[0]);
			unsigned int val_ = ((b1_ << 8) | b2_);
			double val = val_;
			//depr: it is ok as it is //TODO: move to per channel processing
			val = DATA_VOLTAGE_AMPLITUDE*(val / DATA_VOLTAGE_CHANNELS) + DATA_VOLTAGE_OF_ZERO_CHANNEL;
			xs.push_back(read_N*DATA_TIME_CONSTANT);
			ys.push_back(val);
			++read_N;
		}
	}
	file.close();
	return;
}

bool SingleRunData::test_PMT_signal(int _N_threshold, double _S_threshold, double _S_max_threshold, AllRunsResults *results)
{
	int ind = curr_area.channels.get_order_index_by_index(0);
	if (ind < 0)
		return true; //can't test, so assuming the PMT signal is correct
	
#ifdef _NO_PMT_SELECTION
	return true;
#endif
	if (_S_max_threshold < _S_threshold){ //means there is no maximum limit
		if ((PMT3_summed_peaks_area >= _S_threshold) && (PMT3_n_peaks >= _N_threshold))
			return true;
	} else {
		if ((PMT3_summed_peaks_area >= _S_threshold) && (PMT3_summed_peaks_area <= _S_max_threshold) && (PMT3_n_peaks >= _N_threshold))
			return true;
	}
	return false;
}

void SingleRunData::push_event (AllRunsResults *all_runs_results)
{
	all_runs_results->_valid.push_back(true);
	all_runs_results->_status.push_back(AllRunsResults::Status::Empty);

	all_runs_results->mppc_peaks_in_S2_area.push_back(STD_CONT<double>()); //[run#][channel], size of mppc channels (depends on experiment area)
	all_runs_results->mppc_S2_start_time.push_back(STD_CONT<double>());	 //[run#][channel]
	all_runs_results->mppc_S2_finish_time.push_back(STD_CONT<double>());	 //[run#][channel]
	all_runs_results->mppc_double_Is.push_back(STD_CONT<double>());	 //[run#][channel]
	all_runs_results->mppc_peaks.push_back(STD_CONT<STD_CONT<peak> >()); //[run#][channel][peaks].

	all_runs_results->_ns.push_back(0);
	all_runs_results->_Ss.push_back(-1);
	all_runs_results->pmt_peaks.push_back(STD_CONT<STD_CONT<peak> >());	//[run#][channel][peaks]
	all_runs_results->pmt_S2_integral.push_back(STD_CONT<double>()); //[run#][channel]
}

void SingleRunData::processSingleRun_Iter_0(AllRunsResults *all_runs_results)
{
	/*	Process PMTs only. First, determine that the event is valid (not empty) by sum(PMT) signal (ch0).
	 *	This is done by test_PMT_signal function after ch0 signal is processed.
	 *	All PMT channels are processed in ensemble:
	 *		1) Inversion (if required)
	 *		2) Filter
	 *		3) Finding baseline&thresholds for finding peaks. Normally thresholds are set in ParameterPile
	 *		4) Finding peaks, taking integrals of signals in preset time windows. Draw if necessary.
	 *		5) test_PMT_signal - decide whether this run is valid
	 */
	push_event(all_runs_results);
	int run_index = all_runs_results->N_of_runs;
#ifdef _USE_TIME_STATISTICS
	auto PMT_start_timer = std::chrono::high_resolution_clock::now();
#endif

#ifdef _USE_TIME_STATISTICS
	auto PMT_read_file_start_timer = std::chrono::high_resolution_clock::
			now();
#endif
	bool first_pmt_processed = all_runs_results->pmt_channels.empty();
	int pmt_index = -1;
	for (int ch = curr_area.channels.get_next_index(); ch != -1; ch = curr_area.channels.get_next_index()) {
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if ((ch>=32)||(ind<0)||(2==ch))
			continue;
		readOneRun(all_runs_results, ch); //read all PMTs, ignore GEM and MPPCs
		if (xs_channels[ind].empty())
			continue;
		++pmt_index;
		if (first_pmt_processed) {
			all_runs_results->pmt_channels.push_back(ch);
		} else {
			if (all_runs_results->pmt_channels.size()>pmt_index)
				if (all_runs_results->pmt_channels[pmt_index]==ch)
					goto l_valid;
			all_runs_results->_status[run_index] = AllRunsResults::Status::PMT_mismatch;
			all_runs_results->_valid[run_index] = false;
			break;
			l_valid:;
		}
		all_runs_results->pmt_peaks[run_index].push_back(STD_CONT<peak>());
		all_runs_results->pmt_S2_integral[run_index].push_back(-1);

		if ((8<=ch)&&(12>=ch))
			SignalOperations::invert_y(xs_channels[ind],ys_channels[ind]); //invert fast PMT signals
		DVECTOR ys_raw = ys_channels[ind];

#ifdef _USE_TIME_STATISTICS
		auto PMT_read_file_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_PMT_file_reading++;
		all_runs_results->time_stat.t_PMT_file_reading +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(PMT_read_file_end_timer - PMT_read_file_start_timer).count();
#endif

#ifdef _USE_TIME_STATISTICS
		auto PMT_filter_start_timer = std::chrono::high_resolution_clock::now();
#endif
		//apply the filter to PMTs
		bool valid_pars = ParameterPile::filter_PMT_n_points.find(ch)!=ParameterPile::filter_PMT_n_points.end();
		valid_pars = valid_pars&&(ParameterPile::filter_PMT_order.find(ch)!=ParameterPile::filter_PMT_order.end());
		valid_pars = valid_pars&&(ParameterPile::filter_PMT_n_iterations.find(ch)!=ParameterPile::filter_PMT_n_iterations.end());
		if (valid_pars&&!xs_channels[ind].empty()) {//otherwise no filter is applied
			SavitzkyGolayFilter SGfilter(ParameterPile::filter_PMT_n_points.find(ch)->second, ParameterPile::filter_PMT_order.find(ch)->second,
					ParameterPile::filter_PMT_n_iterations.find(ch)->second);
			SGfilter(xs_channels[ind], ys_channels[ind]);
		}
#ifdef _USE_TIME_STATISTICS
		auto PMT_filter_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_PMT_filtering++;
		all_runs_results->time_stat.t_PMT_filtering +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(PMT_filter_end_timer - PMT_filter_start_timer).count();
#endif

#ifdef _USE_TIME_STATISTICS
		auto PMT_baseline_start_timer = std::chrono::high_resolution_clock::now();
#endif
		STD_CONT<peak> ppeaks;
		double threshold = 0;
		double threshold_edges = 0;
		calculate_PMT_threshold_and_baseline(xs_channels[ind], ys_channels[ind], threshold, threshold_edges, found_base_lines[ind], ppeaks, ch);

#ifdef _USE_TIME_STATISTICS
		auto PMT_baseline_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_PMT_baseline++;
		all_runs_results->time_stat.t_PMT_baseline +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(PMT_baseline_end_timer - PMT_baseline_start_timer).count();
		auto PMT_peaks_start_timer = std::chrono::high_resolution_clock::now();
#endif

		SignalOperations::find_peaks_fine(xs_channels[ind], ys_channels[ind], all_runs_results->pmt_peaks[run_index][pmt_index], found_base_lines[ind],
			threshold, threshold_edges, ParameterPile::PMT_N_of_averaging);

		if (0==ch){
			double S_sum = 0;
			for (auto i = all_runs_results->pmt_peaks[run_index][pmt_index].begin(), _end_ = all_runs_results->pmt_peaks[run_index][pmt_index].end(); i != _end_; ++i)
				if ((i->t >= ParameterPile::S2_start_time.find(curr_area.experiments.back())->second)&&(i->t <= ParameterPile::S2_finish_time.find(curr_area.experiments.back())->second)){
					S_sum += i->S >0 ? i->S : 0;
					++all_runs_results->_ns[run_index];
				}
			all_runs_results->_Ss[run_index] = S_sum;
			PMT3_summed_peaks_area = S_sum;
			PMT3_n_peaks = all_runs_results->_ns[run_index];
		}
		double S2_st = ParameterPile::S2_start_time.find(curr_area.experiments.back())->second;
		double S2_ft = ParameterPile::S2_finish_time.find(curr_area.experiments.back())->second;
		DVECTOR ys_int, xs_int;
		SignalOperations::integrate(xs_channels[ind],ys_channels[ind],xs_int, ys_int, S2_st, S2_ft,
				(*(xs_channels[ind].begin()+1)-*(xs_channels[ind].begin())), found_base_lines[ind]);
		DITERATOR x_max;
		SignalOperations::get_max(xs_int, ys_int, x_max, all_runs_results->pmt_S2_integral[run_index][pmt_index], 1);

		//=========================================================================
		ParameterPile::experiment_area area = curr_area.to_point();
		area.channels.erase();
		area.channels.push_back(ch);
		std::string plot_title = curr_area.experiments.back() + "_run" + std::to_string(curr_area.runs.back()) +
			"_ch" + std::to_string(ch) + "_sub" + std::to_string(curr_area.sub_runs.back());
		if (ParameterPile::draw_required(area)) {
			std::string plot_name = "";
			plot_name += std::string("PMT_ch") + std::to_string(ch);
			plot_name += "_run" + std::to_string(curr_area.runs.back()) + "_sub" + std::to_string(area.sub_runs.back());
			int ind = curr_area.channels.get_order_index_by_index(ch);

			Drawing *dr = graph_manager.GetDrawing(plot_name, ch, ParameterPile::DrawEngine::Gnuplot);
			if (NULL != dr) {
				dr->SetDirectory(OUTPUT_DIR+OUTPUT_PMTS + curr_area.experiments.back() + "/PMT_"  + std::to_string(ch) + "/");
				dr->AddToDraw(xs_channels[ind], ys_raw, "raw" + plot_title, "with lines", 0);
				dr->AddToDraw(xs_channels[ind], ys_channels[ind], "filtered" + plot_title, "with lines lw 2", 0);
				dr->AddToDraw_baseline(threshold, "threshold");
				if (threshold_edges!=found_base_lines[ind])
					dr->AddToDraw_baseline(threshold_edges, "threshold 2nd", "w l lc rgb \"#AC0ECD\"");
				dr->AddToDraw_baseline(found_base_lines[ind], "baseline", "w l lc rgb \"#0000FF\"");
				dr->AddToDraw_vertical(S2_st, -1, 1, "lc rgb \"#0000FF\"");
				dr->AddToDraw_vertical(S2_ft, -1, 1, "lc rgb \"#0000FF\"");
			}
		}
#ifdef _USE_TIME_STATISTICS
		auto PMT_peaks_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_PMT_peaks++;
		all_runs_results->time_stat.t_PMT_peaks +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(PMT_peaks_end_timer - PMT_peaks_start_timer).count();
#endif

	}

	if (all_runs_results->_valid[run_index]) {
		if (!test_PMT_signal(all_runs_results->N_peaks_cutoff, all_runs_results->S_peaks_cutoff, all_runs_results->S_peaks_max_cutoff, all_runs_results)) {
			all_runs_results->_status[run_index] = AllRunsResults::Status::NoPMTsignal;
			all_runs_results->_valid[run_index] = false;
		}
	}
//=============================================================================
#ifdef _USE_TIME_STATISTICS
	auto PMT_end_timer = std::chrono::high_resolution_clock::now();
	all_runs_results->time_stat.n_PMT_proc++;
	all_runs_results->time_stat.t_PMT_proc += std::chrono::duration_cast<std::chrono::nanoseconds>(PMT_end_timer - PMT_start_timer).count();
#endif
	runProcessedProc(all_runs_results);
	return;
}

//can not be channel independent. In difference to MPPC this one returns threshold shifted by baseline
void SingleRunData::calculate_PMT_threshold_and_baseline(DVECTOR &xs, DVECTOR &ys, double &threshold, double &threshold_2, double &baseline, STD_CONT<peak> &peaks_before_S1, int channel)
{
	DVECTOR xs_before_S1 = xs;
	DVECTOR ys_before_S1 = ys;

	double delta_x = *(xs.begin() + 1) - *xs.begin();
	SignalOperations::apply_time_limits(xs_before_S1, ys_before_S1, *xs_before_S1.begin(), ParameterPile::S1_start_time, delta_x);

	baseline = (0==channel) ? SignalOperations::find_baseline_by_integral(baseline, xs, ys)
		: SignalOperations::find_baseline_by_median(baseline, xs, ys);
	//calculate first-order baseline
	double noise_amp;
	noise_amp = SignalOperations::RMS(ys_before_S1.begin(), ys_before_S1.end());
	/*approx*/threshold = noise_amp*ParameterPile::PMT_run_acceptance_threshold_to_noize;

	if ((0 == channel)||(1==channel)){
		SignalOperations::find_peaks_fine(xs_before_S1, ys_before_S1, peaks_before_S1,
			baseline, threshold + baseline, 0/*noise_amp*0.5 in mppc*/ + baseline, ParameterPile::PMT_N_of_averaging);
		double exact_noise = noise_amp;
		if (!peaks_before_S1.empty()) {
			/*exact*/baseline = SignalOperations::find_baseline_by_integral(0, xs_before_S1, ys_before_S1,peaks_before_S1);
			SignalOperations::exclude_peaks(xs_before_S1, ys_before_S1, peaks_before_S1);
			exact_noise = SignalOperations::RMS(ys_before_S1.begin(), ys_before_S1.end());
		}
		/*exact*/threshold = exact_noise*ParameterPile::PMT_run_acceptance_threshold_to_noize;
	}
	bool valid_limits=(ParameterPile::PMT_maximum_thresh.find(channel)!=ParameterPile::PMT_maximum_thresh.end());
	valid_limits=valid_limits&&(ParameterPile::PMT_minimum_thresh.find(channel)!=ParameterPile::PMT_minimum_thresh.end());
	if (valid_limits){
		if (threshold < ParameterPile::PMT_minimum_thresh.find(channel)->second)
			threshold = ParameterPile::PMT_minimum_thresh.find(channel)->second;
		if (threshold > ParameterPile::PMT_maximum_thresh.find(channel)->second)
			threshold = ParameterPile::PMT_maximum_thresh.find(channel)->second;
	}
	threshold += baseline;
	threshold_2 = baseline;
	if (ParameterPile::PMT_thresh_edges.find(channel)!=ParameterPile::PMT_thresh_edges.end())
		threshold_2+=ParameterPile::PMT_thresh_edges.find(channel)->second;
}

void SingleRunData::push_average (int ch, bool is_first_call, AllRunsResults *all_runs_results)
{
	int ind =  curr_area.channels.get_order_index_by_index(ch);
	if (ind<0)
		return;
	if (is_first_call) {
		all_runs_results->avr_channels.push_back(ch);
		all_runs_results->_xs_sum.push_back(xs_channels[ind]);
		all_runs_results->_ys_sum.push_back(ys_channels[ind]);
		all_runs_results->_ys_disp.push_back(DVECTOR());
		all_runs_results->_ys_disp.back().resize(xs_channels[ind].size(), 0);
	} else {
		int avg_index = -1;
		for (int i=0, end_=all_runs_results->avr_channels.size(); i!=end_; ++i) {
			if (all_runs_results->avr_channels[i]==ch) {
				avg_index = i;
				break;
			}
		}
		if ((-1==avg_index)||(ys_channels[ind].size()!=all_runs_results->_ys_sum[avg_index].size())) {
			all_runs_results->_status[all_runs_results->N_of_runs] = AllRunsResults::Status::AVG_mismatch;
			all_runs_results->_valid[all_runs_results->N_of_runs] = false;
			return;
		}
		for (auto j=ys_channels[ind].begin(), i = all_runs_results->_ys_sum[avg_index].begin(),
				j_end_ = ys_channels[ind].end(), i_end_ = all_runs_results->_ys_sum[avg_index].end();(i!=i_end_)&&(j!=j_end_);++i,++j) {
			*i+=*j;
		}
	}
}

void SingleRunData::processSingleRun_Iter_1(AllRunsResults *all_runs_results)
{
	if (!all_runs_results->_valid[all_runs_results->N_of_runs]) {
		runProcessedProc(all_runs_results);
		return;
	}
#ifdef _USE_TIME_STATISTICS
	auto PMT_start_timer = std::chrono::high_resolution_clock::now();
#endif
	int run_index = all_runs_results->N_of_runs;
	bool is_first_call_avr = all_runs_results->avr_channels.empty();
	bool is_first_call_mppc = all_runs_results->mppc_channels.empty();
	for (int ch = curr_area.channels.get_next_index(); ch != -1; ch = curr_area.channels.get_next_index()) {
		int ind =  curr_area.channels.get_order_index_by_index(ch);
		if ((ind<0)||(ch>=32)) //process only GEM and PMTs
			continue;
		if (all_runs_results->_to_average.contains(ch)) {
			readOneRun(all_runs_results, ch);
			if ((8<=ch)&&(12>=ch))
				SignalOperations::invert_y(xs_channels[ind], ys_channels[ind]);
			push_average(ch, is_first_call_avr, all_runs_results);
			clearOneRun(ch);
		}
	}
#ifdef _USE_TIME_STATISTICS
	auto PMT_end_timer = std::chrono::high_resolution_clock::now();
	all_runs_results->time_stat.t_PMT_proc += std::chrono::duration_cast<std::chrono::nanoseconds>(PMT_end_timer - PMT_start_timer).count();
#endif
	if (!all_runs_results->_valid[run_index]) {
		runProcessedProc(all_runs_results);
		return;
	}
#ifdef _USE_TIME_STATISTICS
	auto MPPC_start_timer = std::chrono::high_resolution_clock::now();
#endif

	SavitzkyGolayFilter SGfilter(ParameterPile::filter_MPPC_n_points, ParameterPile::filter_MPPC_order, ParameterPile::filter_MPPC_n_iterations);
	int mppc_ind = -1;
	curr_area.channels.reset();
	for (int ch = curr_area.channels.get_next_index(); ch != -1; ch = curr_area.channels.get_next_index()) {
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if ((ind < 0)||(ch<32)) //process only mppc
			continue;
		++mppc_ind;
#ifdef _USE_TIME_STATISTICS
		auto MPPC_read_start_timer = std::chrono::high_resolution_clock::now();
#endif
		readOneRun(all_runs_results, ch);
#ifdef _USE_TIME_STATISTICS
		auto MPPC_read_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_file_reading++;
		all_runs_results->time_stat.t_MPPC_file_reading+=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_read_end_timer - MPPC_read_start_timer).count();
#endif
		if (xs_channels[ind].empty())
			continue;
		if (is_first_call_mppc) {
			all_runs_results->mppc_channels.push_back(ch);
		} else {
			if (all_runs_results->mppc_channels.size()>mppc_ind)
				if (all_runs_results->mppc_channels[mppc_ind]==ch)
					goto l_valid;
			all_runs_results->_status[run_index] = AllRunsResults::Status::MPPC_mismatch;
			all_runs_results->_valid[run_index] = false;
			clearOneRun(ch);
			break;
			l_valid:;
		}
		all_runs_results->mppc_peaks[run_index].push_back(STD_CONT<peak>());
		all_runs_results->mppc_peaks_in_S2_area[run_index].push_back(0);
		all_runs_results->mppc_S2_start_time[run_index].push_back(0);
		all_runs_results->mppc_S2_finish_time[run_index].push_back(0);
		all_runs_results->mppc_double_Is[run_index].push_back(-1);
		DVECTOR mppc_baseline_xs = xs_channels[ind], mppc_baseline_ys;

//==========================================================================
		// invert MPPCs
		SignalOperations::invert_y(xs_channels[ind], ys_channels[ind]);
		if (all_runs_results->_to_average.contains(ch)){
			push_average(ch, is_first_call_avr, all_runs_results);
		}
#ifdef _USE_TIME_STATISTICS
		//auto MPPC_filter_start_timer = std::chrono::high_resolution_clock::now();
#endif

		SGfilter(xs_channels[ind], ys_channels[ind]);
		
#ifdef _USE_TIME_STATISTICS
		/*auto MPPC_filter_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_filtering++;
		all_runs_results->time_stat.t_MPPC_filtering +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_filter_end_timer - MPPC_filter_start_timer).count();*/
#endif
		DVECTOR xs_raw = xs_channels[ind];
		DVECTOR ys_raw = ys_channels[ind];

//================================================================================
#ifdef _USE_TIME_STATISTICS
		auto MPPC_threshold_and_first_baseline_start_timer = std::chrono::high_resolution_clock::now();
#endif
		double global_threshold;
		double raw_signal_baseline = 0;
		STD_CONT<peak> peaks_before_S1;
		calculate_MPPC_threshold_and_baseline(xs_channels[ind], ys_channels[ind], global_threshold, raw_signal_baseline, peaks_before_S1);
#ifdef _USE_TIME_STATISTICS
		auto MPPC_threshold_and_first_baseline_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_threshold_and_first_baseline++;
		all_runs_results->time_stat.t_MPPC_threshold_and_first_baseline +=
std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_threshold_and_first_baseline_end_timer - MPPC_threshold_and_first_baseline_start_timer).count();
#endif
//================================================================================
		DVECTOR ys_cut = ys_channels[ind];
#ifndef _TEMP_CODE_ //use optimal x region
#ifdef _USE_TIME_STATISTICS
		auto MPPC_curved_baseline_start_timer = std::chrono::high_resolution_clock::now();
		auto MPPC_filter_start_timer = std::chrono::high_resolution_clock::now();
#endif
		double delta_x = *(mppc_baseline_xs.begin() + 1) - *(mppc_baseline_xs.begin());
		DITERATOR _end_ = mppc_baseline_xs.end() - 1;
		DITERATOR _begin_ = mppc_baseline_xs.begin();
		DITERATOR S2_start_i = SignalOperations::find_x_iterator_by_value(_begin_,
			_end_, ParameterPile::S2_start_time.find(curr_area.experiments.back())->second, delta_x);
		DITERATOR S2_finish_i = SignalOperations::find_x_iterator_by_value(S2_start_i,
			_end_, ParameterPile::S2_finish_time.find(curr_area.experiments.back())->second, delta_x);
		DITERATOR max_signal_pos;
		double max_signal_val;
		//TODO: ParameterPile
		++S2_finish_i;
		SignalOperations::get_max(mppc_baseline_xs, ys_cut, S2_start_i, S2_finish_i, max_signal_pos, max_signal_val, 1);
		double root_start_t = *max_signal_pos - ParameterPile::MPPC_ROOTs_bl_from_max_left;
		double root_finish_t = *max_signal_pos + ParameterPile::MPPC_ROOTs_bl_from_max_right;
		SignalOperations::apply_time_limits(mppc_baseline_xs, ys_cut, root_start_t, root_finish_t,delta_x);
#ifdef _USE_TIME_STATISTICS
		auto MPPC_filter_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_filtering++;
		all_runs_results->time_stat.t_MPPC_filtering +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_filter_end_timer - MPPC_filter_start_timer).count();
#endif //_USE_TIME_STATISTICS
#else //do not use optimal x region
		double root_start_t = ParameterPile::S2_start_time;
		double root_finish_t = ParameterPile::S2_finish_time;
		SignalOperations::apply_time_limits(_result.mppc_baseline_xs[mppc_ind], ys_cut, root_start_t, root_finish_t, delta_x);
#ifdef _USE_TIME_STATISTICS
		auto MPPC_curved_baseline_start_timer = std::chrono::high_resolution_clock::now();
#endif
#endif

#ifdef _USE_TIME_STATISTICS
		/*DVECTOR root_baseline_v2;
		DVECTOR root_baseline_v3;
		DVECTOR root_baseline_v4;*/
		//DVECTOR root_baseline_v5;
		//DVECTOR root_baseline_v6;
		//DVECTOR root_baseline_v7;
		//DVECTOR root_baseline_v8;
		//for (int a = 0; a < 10; a++) {
#endif
			SignalOperations::find_baseline_by_ROOT(mppc_baseline_xs, ys_cut, mppc_baseline_ys);
#ifdef _USE_TIME_STATISTICS
		//}
		auto MPPC_curved_baseline_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_curved_baseline++;
		all_runs_results->time_stat.t_MPPC_curved_baseline +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_curved_baseline_end_timer - MPPC_curved_baseline_start_timer).count();

		/*MPPC_curved_baseline_start_timer = std::chrono::high_resolution_clock::now();
		for (int a = 0; a < 10; a++) {
			SignalOperations::find_baseline_by_ROOT_v2(_result.mppc_baseline_xs[mppc_ind], ys_cut, root_baseline_v2);
		}
		MPPC_curved_baseline_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_curved_baseline_v2++;
		all_runs_results->time_stat.t_MPPC_curved_baseline_v2 +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_curved_baseline_end_timer - MPPC_curved_baseline_start_timer).count();
		MPPC_curved_baseline_start_timer = std::chrono::high_resolution_clock::now();
		for (int a = 0; a < 10; a++) {
			SignalOperations::find_baseline_by_ROOT_v3(_result.mppc_baseline_xs[mppc_ind], ys_cut, root_baseline_v3);
		}
		MPPC_curved_baseline_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_curved_baseline_v3++;
		all_runs_results->time_stat.t_MPPC_curved_baseline_v3 +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_curved_baseline_end_timer - MPPC_curved_baseline_start_timer).count();
		MPPC_curved_baseline_start_timer = std::chrono::high_resolution_clock::now();
		for (int a = 0; a < 10; a++) {
			SignalOperations::find_baseline_by_ROOT_v4(_result.mppc_baseline_xs[mppc_ind], ys_cut, root_baseline_v4);
		}
		MPPC_curved_baseline_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_curved_baseline_v4++;
		all_runs_results->time_stat.t_MPPC_curved_baseline_v4 +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_curved_baseline_end_timer - MPPC_curved_baseline_start_timer).count();*/
		//MPPC_curved_baseline_start_timer = std::chrono::high_resolution_clock::now();
		//for (int a = 0; a < 10; a++) {
		//	SignalOperations::find_baseline_by_ROOT_v5(_result.mppc_baseline_xs[mppc_ind], ys_cut, root_baseline_v5);
		//}
		//MPPC_curved_baseline_end_timer = std::chrono::high_resolution_clock::now();
		//all_runs_results->time_stat.n_MPPC_curved_baseline_v5++;
		//all_runs_results->time_stat.t_MPPC_curved_baseline_v5 +=
		//	std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_curved_baseline_end_timer - MPPC_curved_baseline_start_timer).count();
		//MPPC_curved_baseline_start_timer = std::chrono::high_resolution_clock::now();
		//for (int a = 0; a < 10; a++) {
		//	SignalOperations::find_baseline_by_ROOT_v6(_result.mppc_baseline_xs[mppc_ind], ys_cut, root_baseline_v6);
		//}
		//MPPC_curved_baseline_end_timer = std::chrono::high_resolution_clock::now();
		//all_runs_results->time_stat.n_MPPC_curved_baseline_v6++;
		//all_runs_results->time_stat.t_MPPC_curved_baseline_v6 +=
		//	std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_curved_baseline_end_timer - MPPC_curved_baseline_start_timer).count();
		//MPPC_curved_baseline_start_timer = std::chrono::high_resolution_clock::now();
		//for (int a = 0; a < 10; a++) {
		//	SignalOperations::find_baseline_by_ROOT_v7(_result.mppc_baseline_xs[mppc_ind], ys_cut, root_baseline_v7);
		//}
		//MPPC_curved_baseline_end_timer = std::chrono::high_resolution_clock::now();
		//all_runs_results->time_stat.n_MPPC_curved_baseline_v7++;
		//all_runs_results->time_stat.t_MPPC_curved_baseline_v7 +=
		//	std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_curved_baseline_end_timer - MPPC_curved_baseline_start_timer).count();
		//MPPC_curved_baseline_start_timer = std::chrono::high_resolution_clock::now();
		//for (int a = 0; a < 10; a++) {
		//	SignalOperations::find_baseline_by_ROOT_v8(_result.mppc_baseline_xs[mppc_ind], ys_cut, root_baseline_v8);
		//}
		//MPPC_curved_baseline_end_timer = std::chrono::high_resolution_clock::now();
		//all_runs_results->time_stat.n_MPPC_curved_baseline_v8++;
		//all_runs_results->time_stat.t_MPPC_curved_baseline_v8 +=
		//	std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_curved_baseline_end_timer - MPPC_curved_baseline_start_timer).count();
#endif
//================================================================================
		double ROOTs_baseline_baseline = 0;
#ifdef _USE_TIME_STATISTICS
		auto MPPC_baseline_baseline_start_timer = std::chrono::high_resolution_clock::now();
#endif
		STD_CONT<peak> exclude_middle;
		peak middle;
		middle.left = root_start_t + ParameterPile::MPPC_ROOTs_bl_left_offset;
		middle.right = root_finish_t - ParameterPile::MPPC_ROOTs_bl_right_offset;
		exclude_middle.push_back(middle);
		ROOTs_baseline_baseline = SignalOperations::find_baseline_by_integral(0, mppc_baseline_xs, mppc_baseline_ys, exclude_middle);
#ifdef _USE_TIME_STATISTICS
		auto MPPC_baseline_baseline_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_curved_baseline_baseline++;
		all_runs_results->time_stat.t_MPPC_curved_baseline_baseline +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_baseline_baseline_end_timer - MPPC_baseline_baseline_start_timer).count();
#endif
//================================================================================
#ifdef _USE_TIME_STATISTICS
		auto MPPC_substruction_start_timer = std::chrono::high_resolution_clock::now();
#endif

		//So in result I have to do 3 baseline substractions no matter what. 
		//One of them is hidden in SignalOperations::substract_baseline(DV&,DV&,DV&,DV&,double), so it doesn't cost anything
		SignalOperations::substract_baseline(ys_channels[ind], raw_signal_baseline);
		DVECTOR xs_raw_0bl = xs_channels[ind];
		DVECTOR ys_raw_0bl = ys_channels[ind];//I need this one for double integral.
		SignalOperations::substract_baseline(xs_channels[ind], ys_channels[ind], mppc_baseline_xs,
			mppc_baseline_ys, ROOTs_baseline_baseline); //handles that 2 signals have different x spans
		//TODO: optimize this function, using the fact that the ROOT's baseline has definite x points [from, to] (- [xs.begin(),xs.end())
		//and rework _result.mppc_baseline_xs[mppc_ind]
		found_base_lines[ind] = 0;
#ifdef _USE_TIME_STATISTICS
		auto MPPC_substruction_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_baseline_substraction++;
		all_runs_results->time_stat.t_MPPC_baseline_substraction +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_substruction_end_timer - MPPC_substruction_start_timer).count();
#endif
//================================================================================
#ifdef _USE_TIME_STATISTICS
		auto MPPC_peaks_finding_start_timer = std::chrono::high_resolution_clock::now();
#endif
		SignalOperations::find_peaks_fine(xs_channels[ind], ys_channels[ind], all_runs_results->mppc_peaks[run_index][mppc_ind],
				0, global_threshold, 0/*edge_threshold*/, ParameterPile::MPPC_N_trust);
		//^TODO: ParameterPile
#ifdef _USE_TIME_STATISTICS
		auto MPPC_peaks_finding_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_peaks_finding++;
		all_runs_results->time_stat.t_MPPC_peaks_finding +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_peaks_finding_end_timer - MPPC_peaks_finding_start_timer).count();
#endif
//================================================================================
#ifdef _USE_TIME_STATISTICS
		auto MPPC_peaks_proc_start_timer = std::chrono::high_resolution_clock::now();
#endif
		DVECTOR x_peaks, y_peaks;
		SignalOperations::peaks_to_yx(*(xs_channels[ind].begin()), xs_channels[ind].back(), all_runs_results->mppc_peaks[run_index][mppc_ind], x_peaks, y_peaks);
		DVECTOR x_peaks_spreaded_v2, y_peaks_spreaded_v2;
		SignalOperations::spread_peaks_v2(*(xs_channels[ind].begin()), xs_channels[ind].back(), all_runs_results->mppc_peaks[run_index][mppc_ind],
			x_peaks_spreaded_v2, y_peaks_spreaded_v2, ParameterPile::MPPC_peaks_smoothing_time);

		//TODO: actually at first I have to remove peaks and recalculate baseline by averaging over dt. Then peaks again.
		//And this is only valid for small fields. At first I'll reproduce previous reseachers' results and only the will try to
		//predict baseline behaviour at large field.
		//At low field and for at least some first histograms I'll just use ROOT's baseline.
		//Also TODO: dynamic determination of how strongly baseline is shifted at S2 region and how close MPPC peaks to each other in
		//order to classify experiments on the run. Classification is: {LowField, MediumField, LargeField};
		//Criterions:	1) how strongly min of (ROOT's) baseline is shifted
		//				2) min is in S2 region (approx), or cheaply calculated or even ParameterPile
		//				3) maybe also Summed Peaks in S2 (approximate ones, which should be cheap to calculate)
		//The main problem some experiments may turn out to be 50/50 (every run is different).
		//I must see how it actually is, but maybe I'll split 50/50 experiments in two subexperiments in AllRunsResult (and later in AllExperimentsResults)
		//"MediumField" is exactly for this purpose. If, say, 10% of runs fall into LargeField, use whole experiment as 'LowField' anyway.
		double sp_min;
		double S2_finish_t = ParameterPile::S2_finish_time.find(curr_area.experiments.back())->second,
				S2_start_t = ParameterPile::S2_start_time.find(curr_area.experiments.back())->second;
		sp_min = *std::min_element(y_peaks_spreaded_v2.begin(), y_peaks_spreaded_v2.end());
		double sp_approx_thr;
		double sp_threshold = find_spreaded_peaks_threshold(x_peaks_spreaded_v2, y_peaks_spreaded_v2, sp_approx_thr);
		STD_CONT<peak> sp_peaks;
		std::string __exp = curr_area.experiments.back();
		STD_CONT<peak>::iterator max_intersection; //search for largest cluster which is above threshold. In 99% if not 100% I
		//suspect there will be only single sp_peaks i.e. intersection, but just in case.
		if (sp_threshold < 0)
			goto l_quit_time;
		//no fine required here
		SignalOperations::find_peaks(x_peaks_spreaded_v2, y_peaks_spreaded_v2, sp_peaks, sp_min, sp_threshold, 1);
		//TODO: don't forget to test this^ place
		if (sp_peaks.empty())
			goto l_quit_time;
		max_intersection = std::max_element(sp_peaks.begin(), sp_peaks.end(), 
			[__exp](peak a, peak b) {
			if ((b.left > ParameterPile::S2_finish_time.find(__exp)->second) || (b.right < ParameterPile::S2_start_time.find(__exp)->second))
				return false;
			if ((a.left > ParameterPile::S2_finish_time.find(__exp)->second) || (a.right < ParameterPile::S2_start_time.find(__exp)->second))
				return true;
			return a.S < b.S; 
		});
		//TODO: ParameterPile
		S2_finish_t = std::min(max_intersection->right + ParameterPile::MPPC_peaks_smoothing_time, S2_finish_t);
		S2_start_t = std::max(max_intersection->left - ParameterPile::MPPC_peaks_smoothing_time, S2_start_t);
		if (S2_finish_t < S2_start_t) {
			S2_finish_t = 0;
			S2_start_t = 0;
		}
		l_quit_time:
		all_runs_results->mppc_S2_finish_time[run_index][mppc_ind] = S2_finish_t;
		all_runs_results->mppc_S2_start_time[run_index][mppc_ind] = S2_start_t;
		//TODO: decide whether this is required. ParameterPile
		for (auto pp = all_runs_results->mppc_peaks[run_index][mppc_ind].begin(); pp != all_runs_results->mppc_peaks[run_index][mppc_ind].end(); ++pp) {
			if ((pp->left >= S2_start_t) && (pp->right <= S2_finish_t)) {
				all_runs_results->mppc_peaks_in_S2_area[run_index][mppc_ind]+= pp->S;
			}
		}

#ifdef _USE_TIME_STATISTICS
		auto MPPC_peaks_proc_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_peaks_processing++;
		all_runs_results->time_stat.t_MPPC_peaks_processing +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_peaks_proc_end_timer - MPPC_peaks_proc_start_timer).count();
#endif
//================================================================================
		//Now calculate double integral for signal without subtracting ROOT's baseline
#ifdef _USE_TIME_STATISTICS
		auto MPPC_double_I_start_timer = std::chrono::high_resolution_clock::now();
#endif
		SignalOperations::apply_time_limits(xs_raw_0bl, ys_raw_0bl, S2_start_t, S2_finish_t, delta_x);
		DVECTOR ys_raw_i, ys_raw_ii;
		DITERATOR x_ii_max = xs_raw_0bl.end();
		double y_ii_max, y_ii_min;
		if (xs_raw_0bl.size()>=2) {
			SignalOperations::integrate(xs_raw_0bl, ys_raw_0bl, ys_raw_i, delta_x, 0);
			SignalOperations::integrate(xs_raw_0bl, ys_raw_i, ys_raw_ii, delta_x, 0);
			//TODO: ParameterPile
			SignalOperations::get_max(xs_raw_0bl, ys_raw_ii, x_ii_max, y_ii_max, ParameterPile::MPPC_double_I_N_trust);
			y_ii_min = ys_raw_ii.front();
			if (x_ii_max == xs_raw_0bl.begin())
				all_runs_results->mppc_double_Is[run_index][mppc_ind] = ys_raw_ii.back() - y_ii_min;
			else
				all_runs_results->mppc_double_Is[run_index][mppc_ind] = y_ii_max - y_ii_min;
		}
#ifdef _USE_TIME_STATISTICS
		auto MPPC_double_I_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_double_I++;
		all_runs_results->time_stat.t_MPPC_double_I +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_double_I_end_timer - MPPC_double_I_start_timer).count();
#endif
//================================================================================
		ParameterPile::experiment_area area = curr_area.to_point();
		area.channels.erase();
		area.channels.push_back(ch);
		if (ParameterPile::draw_required(area)) {
			std::string plot_title = curr_area.experiments.back() + "\\_run" + std::to_string(curr_area.runs.back()) +
						"\\_ch" + std::to_string(ch) + "\\_sub" + std::to_string(curr_area.sub_runs.back());
			std::string plot_name = "";
			plot_name += std::string("MPPC_ch") + std::to_string(ch);
			plot_name += "_run" + std::to_string(curr_area.runs.back()) + "_sub" + std::to_string(area.sub_runs.back());
			int ind = curr_area.channels.get_order_index_by_index(ch);

			Drawing *dr = graph_manager.GetDrawing(plot_name, ch, ParameterPile::DrawEngine::Gnuplot);
			if (NULL != dr) {
				dr->SetDirectory(OUTPUT_DIR+OUTPUT_MPPCS_PICS + curr_area.experiments.back() + "/" + OUTPUT_MPPCS + std::to_string(ch) + "/");
				dr->AddToDraw(xs_raw, ys_raw, "raw" + plot_title, "w l lc rgb \"#0000FF\"", 0);
				dr->AddToDraw(mppc_baseline_xs, mppc_baseline_ys, "ROOT baseline V1", "w l lw 2 lc rgb \"#00FF00\"", 0);
#ifdef _USE_TIME_STATISTICS
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v2, "ROOT baseline V2", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v3, "ROOT baseline V3", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v4, "ROOT baseline V4", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v5, "ROOT baseline V5", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v6, "ROOT baseline V6", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v7, "ROOT baseline V7", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v8, "ROOT baseline V8", "w l lw 2", 0);
#endif
				dr->AddToDraw(xs_channels[ind], ys_channels[ind], "without baseline " + std::to_string(curr_area.runs.back()), "with lines lc rgb \"#000000\"", 0);
				dr->AddToDraw_baseline(global_threshold, "threshold", "w l lc rgb \"#00FF00\"");
				dr->AddToDraw_baseline(ROOTs_baseline_baseline, "ROOTs baseline", "w l lc rgb \"#FF33FF\"");
				dr->AddToDraw_baseline(0/*edge_threshold*/, "threshold\\_edges");
				dr->AddToDraw_vertical(S2_start_t, -1, 1, "lc rgb \"#FF0000\"");
				dr->AddToDraw_vertical(S2_finish_t, -1, 1, "lc rgb \"#FF0000\"");
				dr->AddToDraw_vertical(middle.left, -1, 1, "lc rgb \"#0000FF\"");
				dr->AddToDraw_vertical(middle.right, -1, 1, "lc rgb \"#0000FF\"");
			}
			#ifdef _DRAW_CLUSTER_FINDING
			dr = graph_manager.GetDrawing(plot_name + "peaks", ch + 100, ParameterPile::DrawEngine::Gnuplot);
			if (NULL != dr) {
				dr->AddToDraw(x_peaks, y_peaks, "peaks "+ curr_area.experiments.back() + "\\_" + std::to_string(curr_area.runs.back())
					+ "\\_sub\\_" + std::to_string(curr_area.sub_runs.back()) + "\\_ch\\_" + std::to_string(ch), "w l", 0);
				dr->AddToDraw(x_peaks_spreaded_v2, y_peaks_spreaded_v2, "peaks spereaded v2 I = " + std::to_string(_result.mppc_S2_peaks_area.back()), "w l", 0);
				dr->AddToDraw_baseline(sp_approx_thr, "approx\\_threshold (median)", "w l");
				dr->AddToDraw_baseline(sp_threshold, "exact\\_threshold", "w l lc rgb \"#FF0000\"");
				dr->AddToDraw_vertical(S2_start_t, -1, 1, "lc rgb \"#FF0000\"");
				dr->AddToDraw_vertical(S2_finish_t, -1, 1, "lc rgb \"#FF0000\"");
			}
			#endif
			/*if (x_ii_max != xs_raw.end()) {
				dr = graph_manager.GetDrawing(plot_name + "integral", ch + 200, ParameterPile::DrawEngine::Gnuplot);
				if (NULL != dr) {
					dr->AddToDraw(xs_raw, ys_raw, "signal" + curr_area.experiments.back() + "\\_" + std::to_string(curr_area.runs.back())
						+ "\\_sub\\_" + std::to_string(curr_area.sub_runs.back()) + "\\_ch\\_" + std::to_string(ch), "w l", 0);
					dr->AddToDraw(xs_raw, ys_raw_i, "first integral", "w l", 0);
					dr->AddToDraw(xs_raw, ys_raw_ii, "second integral = "+std::to_string(_result.mppc_double_I.back()), "w l", 0);
					dr->AddToDraw_vertical(*x_ii_max, -1, 1, "lc rgb \"#FF0000\"");
				}
			}*/
		}
		clearOneRun(ch);
	}
#ifdef _USE_TIME_STATISTICS
	auto MPPC_end_timer = std::chrono::high_resolution_clock::now();
	all_runs_results->time_stat.n_MPPC_proc+=(mppc_ind+1);
	all_runs_results->time_stat.t_MPPC_proc += std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_end_timer - MPPC_start_timer).count();
#endif
	l_quit:
	runProcessedProc(all_runs_results);
	return;
}

void SingleRunData::calculate_MPPC_threshold_and_baseline(DVECTOR &xs, DVECTOR &ys, double &threshold, double &baseline, STD_CONT<peak> &peaks_before_S1)
{
	baseline = 0;
	peaks_before_S1.clear();
	DVECTOR before_S1_x = xs;
	DVECTOR before_S1_y = ys;
	double delta_x = *(xs.begin() + 1) - *xs.begin();
	SignalOperations::apply_time_limits(before_S1_x, before_S1_y, *before_S1_x.begin(), ParameterPile::S1_start_time, delta_x);
	baseline = SignalOperations::find_baseline_by_integral(baseline, before_S1_x, before_S1_y);
	//^not many peaks expected, so integral is used, not a median
	//SignalOperations::substract_baseline(before_S1_y, baseline);
	//SignalOperations::substract_baseline(ys, baseline);
	double approx_noise = SignalOperations::RMS(before_S1_y.begin(), before_S1_y.end());
	double approx_thresh_beforeS1 = approx_noise * ParameterPile::MPPC_threshold_to_noise;
	double approx_edge_thresh = approx_noise*0.5; //TODO: ParameterPile
	SignalOperations::find_peaks_fine(before_S1_x, before_S1_y, peaks_before_S1,
		baseline, approx_thresh_beforeS1 + baseline, approx_edge_thresh + baseline, ParameterPile::MPPC_N_trust);
	double exact_noise = approx_noise;
	if (!peaks_before_S1.empty()) {
		/*exact*/baseline = SignalOperations::find_baseline_by_integral(0, before_S1_x, before_S1_y, peaks_before_S1);
		SignalOperations::exclude_peaks(before_S1_x, before_S1_y, peaks_before_S1);
		exact_noise = SignalOperations::RMS(before_S1_y.begin(), before_S1_y.end());
	}
	threshold = exact_noise * ParameterPile::MPPC_threshold_to_noise;//global means for all times, not for all channels and runs
	if (threshold < ParameterPile::MPPC_minimum_peak_A)
		threshold = ParameterPile::MPPC_minimum_peak_A;
	if (threshold > ParameterPile::MPPC_maximum_peak_A)
		threshold = ParameterPile::MPPC_maximum_peak_A;
}

double SingleRunData::find_spreaded_peaks_threshold(DVECTOR &x_peaks_spreaded, DVECTOR &y_peaks_spreaded, double &apr_thresh)
{
	//variant0: appprox_thresh===(max+min)/2 -> {mean,RMS above ; mean,RMS below} -> spreaded peaks threshold -> S2 time
	//variant1: appprox_thresh===median -> {mean,RMS above ; mean,RMS below} -> spreaded peaks threshold -> S2 time
	//variant2: spereaded peaks threshold === median //I think it's bad choice
	STD_CONT<peak> peaks = STD_CONT<peak>();
	double sp_approx_thr = SignalOperations::find_baseline_by_median(0, x_peaks_spreaded, y_peaks_spreaded, peaks);
	apr_thresh = sp_approx_thr;
	double n_below = 0, n_above = 0;
	double mean_below = 0, mean_above = 0;
	for (auto i = y_peaks_spreaded.begin(); i != y_peaks_spreaded.end(); ++i) {
		if (*i < sp_approx_thr){
			n_below++;
			mean_below += *i;
		} else {
			n_above++;
			mean_above += *i;
		}
	}
	if ((n_below == 0) || (n_above == 0)) {
		return -1;
	}
	mean_above /= n_above;
	mean_below /= n_below;
	double RMS_below = 0, RMS_above = 0;
	for (auto i = y_peaks_spreaded.begin(); i != y_peaks_spreaded.end(); ++i)
		if (*i <= sp_approx_thr) //its better to use <= than <
			RMS_below += (*i - mean_below)*(*i - mean_below);
		else
			RMS_above += (*i - mean_above)*(*i - mean_above);
	RMS_below = sqrt(RMS_below) / n_below; //maybe sqrt(RMS_below/(n_below*(n_below-1))); not that is matters much
	RMS_above = sqrt(RMS_above) / n_above;

	double above_weight = std::sqrt(RMS_below) / std::pow(n_below,1.5);//no real explanation for this, just figuring out algorithm by trial and error.
	double below_weight = std::sqrt(RMS_above) / std::pow(n_above,1.5); //But in the case when statistics (n of points) is small, the dispersion is forcibly increased
	//TODO: maybe use this threshold only for S2 finish time and some other for S2 start time
	//TODO: try other weightings.
	return (mean_above*above_weight + mean_below * below_weight) / (above_weight + below_weight);
}

void SingleRunData::push_dispersion (int ch, AllRunsResults *all_runs_results)
{
	int ind =  curr_area.channels.get_order_index_by_index(ch);
	if ((ind<0)||(all_runs_results->avr_channels.empty()))
		return;
	int avg_index = -1;
	for (int i=0, end_=all_runs_results->avr_channels.size(); i!=end_; ++i) {
		if (all_runs_results->avr_channels[i]==ch) {
			avg_index = i;
			break;
		}
	}
	if (-1==avg_index) {
		all_runs_results->_status[all_runs_results->N_of_runs] = AllRunsResults::Status::AVG_mismatch;
		all_runs_results->_valid[all_runs_results->N_of_runs] = false;
		return;
	}
	for (auto j=ys_channels[ind].begin(), a = all_runs_results->_ys_sum[avg_index].begin(), d = all_runs_results->_ys_disp[avg_index].begin(),
			j_end_ = ys_channels[ind].end(), a_end_ = all_runs_results->_ys_sum[avg_index].end(), d_end_ = all_runs_results->_ys_disp[avg_index].end();
			(d!=d_end_)&&(a!=a_end_)&&(j!=j_end_);++a,++j,++d) {
		*d+=std::pow(*j-*a, 2);
	}
}

void SingleRunData::processSingleRun_Iter_2(AllRunsResults *all_runs_results)
{
	if (!all_runs_results->_valid[all_runs_results->N_of_runs])
		goto l_quit;
	for (int ch = curr_area.channels.get_next_index(); ch != -1; ch = curr_area.channels.get_next_index()) {
		int ind =  curr_area.channels.get_order_index_by_index(ch);
		if (ind<0)
			continue;
		if (all_runs_results->_to_average.contains(ch)) {
			readOneRun(all_runs_results, ch);
			if (ch>=3)
				SignalOperations::invert_y(xs_channels[ind], ys_channels[ind]);
			push_dispersion(ch, all_runs_results);
			clearOneRun(ch);
		}
	}
	l_quit:
	runProcessedProc(all_runs_results);
	return;
}

void SingleRunData::processSingleRun(AllRunsResults *all_runs_results)
{
	if (0 == all_runs_results->Iteration())
		processSingleRun_Iter_0(all_runs_results);
	if (1 == all_runs_results->Iteration())
		processSingleRun_Iter_1(all_runs_results);
	if (2 == all_runs_results->Iteration())
		processSingleRun_Iter_2(all_runs_results);
}

void SingleRunData::runProcessedProc(AllRunsResults *all_runs_results)
{
	graph_manager.Draw();
	graph_manager.Clear();
	if (all_runs_results->_valid[all_runs_results->N_of_runs])
		++(all_runs_results->N_of_valid_runs);
	++(all_runs_results->N_of_runs);
	clear_memory();
}

void SingleRunData::clear_memory(void) //clears only 'input' data, preserves processing results
{
	for (int ind =0,_end_ = xs_channels.size();ind!=_end_;++ind){
		DVECTOR().swap(xs_channels[ind]);
		DVECTOR().swap(ys_channels[ind]);
	}
}

ParameterPile::experiment_area SingleRunData::getArea(void) const {	return curr_area;}

std::size_t SingleRunData::real_size(void)
{
	std::size_t rsize = sizeof(*this);
	rsize+=xs_channels.size()*sizeof(DVECTOR);
	for (std::size_t i=0, _end_ = xs_channels.size();i!=_end_;++i)
		rsize+=xs_channels[i].capacity()*sizeof(double);
	rsize+=ys_channels.size()*sizeof(DVECTOR);
	for (std::size_t i=0, _end_ = ys_channels.size();i!=_end_;++i)
		rsize+=ys_channels[i].capacity()*sizeof(double);

	rsize+=found_base_lines.capacity()*sizeof(double);
	rsize+=curr_area.real_size();
	//rsize+=graph_manager.real_size();
	return rsize;
}


