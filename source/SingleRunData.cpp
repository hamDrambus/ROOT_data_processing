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
	}
	PMT3_summed_peaks_area = 0;
	PMT3_n_peaks = 0;
	PMT_peaks_analysed = false;
}

void SingleRunData::readOneRun(SingleRunResults &_result)
{
#ifdef _HOTFIX_CLEAR_MEMORY
	STD_CONT<DVECTOR>().swap(xs_channels);
	STD_CONT<DVECTOR>().swap(ys_channels);
#else
	xs_channels.clear();
	ys_channels.clear();
#endif
	bool empty_run = true;
	for (int ch = curr_area.channels.get_next_index(); !(ch < 0); ch = curr_area.channels.get_next_index()){
		xs_channels.push_back(DVECTOR());
		ys_channels.push_back(DVECTOR());
#ifdef _HOTFIX_DECREASE_MPPC_MEMORY_USAGE
		if ((ch >= 32 && ch <= 63)||(2==ch)){
			//empty_run = false;
			continue;
		}
#endif
		std::string path = DATA_PREFIX;
		path += curr_area.experiments.back() + "/";
		path += "run_" + std::to_string(curr_area.runs.back()) + "__ch_" + std::to_string(ch) + ".dat";
		file_to_vector(path, xs_channels.back(), ys_channels.back(), curr_area.sub_runs.back());
		if (!xs_channels.back().empty() && !ys_channels.back().empty())
			empty_run = false;
	}
	if (empty_run || xs_channels.empty() || ys_channels.empty() || xs_channels.size() != ys_channels.size()){
		_result._current_status = SingleRunResults::Status::NotLoaded;
		_result.setValid(false);
	}
	else {
		_result._current_status = SingleRunResults::Status::Ok;
		_result.setValid(true);
	}
}

#ifdef _HOTFIX_DECREASE_MPPC_MEMORY_USAGE
void SingleRunData::readOneRunMPPCs(SingleRunResults &results, int channel)
{
	int ind = curr_area.channels.get_order_index_by_index(channel);
	if (ind < 0)
		return;
	std::string path = DATA_PREFIX;
	path += curr_area.experiments.back() + "/";
	path += "run_" + std::to_string(curr_area.runs.back()) + "__ch_" + std::to_string(channel) + ".dat";
	file_to_vector(path, xs_channels[ind], ys_channels[ind], curr_area.sub_runs.back());
}

void SingleRunData::clearOneRunMPPCs(int channel)
{
	int ind = curr_area.channels.get_order_index_by_index(channel);
	if (ind < 0)
		return;
	DVECTOR().swap(xs_channels[ind]);
	DVECTOR().swap(ys_channels[ind]);
}
#endif //_HOTFIX_DECREASE_MPPC_MEMORY_USAGE

void SingleRunData::file_to_vector(std::string fname, DVECTOR &xs, DVECTOR &ys, int index)
{
#ifdef _HOTFIX_CLEAR_MEMORY
	DVECTOR().swap(xs);
	DVECTOR().swap(ys);
#else
	xs.clear();
	ys.clear();
#endif
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

void SingleRunData::find_time_limits(void)
{

}

void SingleRunData::add_draw_baselines(std::string prefix, GraphicOutputManager& graph_manager, ParameterPile::DrawEngine de)
{
	for (int ch = curr_area.channels.get_next_index(); !(ch < 0); ch = curr_area.channels.get_next_index()) {
		ParameterPile::experiment_area area = curr_area.to_point(); //TODO: add constructor of point from area
		area.channels.erase();
		area.channels.push_back(ch);
		if (ParameterPile::draw_required(area)){
			std::string plot_name = "";
			plot_name += curr_area.experiments.back() + "_";
			plot_name += "run_" + std::to_string(curr_area.runs.back()) + "_ch_" + std::to_string(ch);
			int ind = curr_area.channels.get_order_index_by_index(ch);

			Drawing *dr = graph_manager.GetDrawing(plot_name, ch, de);
			if (NULL == dr)
				return;
			dr->AddToDraw_baseline(found_base_lines[ind], prefix, "", 0);
			//dr->AddToDraw(found_base_lines[ind], prefix, 1);
		}
	}
}

void SingleRunData::add_draw_data(std::string prefix, GraphicOutputManager& graph_manager, ParameterPile::DrawEngine de, int channel)
{
	if (channel >= 0){
		int ind = curr_area.channels.get_order_index_by_index(channel);
		if (ind < 0)
			return;
		ParameterPile::experiment_area area = curr_area.to_point();
		area.channels.erase();
		area.channels.push_back(channel);
		if (ParameterPile::draw_required(area)){
			std::string plot_name = "";
			plot_name += curr_area.experiments.back() + "_";
			plot_name += "run" + std::to_string(curr_area.runs.back()) + "_ch" + std::to_string(channel)
				+ "_sub" + std::to_string(area.sub_runs.back());
			ind = curr_area.channels.get_order_index_by_index(channel);
			if (!xs_channels[ind].empty() && !ys_channels[ind].empty()){
				Drawing *dr = graph_manager.GetDrawing(plot_name, channel, de);
				if (NULL == dr)
					return;
				dr->AddToDraw(xs_channels[ind], ys_channels[ind], prefix, "with points pt 1", 0);
			}
		}
		return;
	}
	for (int ch = curr_area.channels.get_next_index(); !(ch < 0); ch = curr_area.channels.get_next_index()) {
		ParameterPile::experiment_area area = curr_area.to_point();
		area.channels.erase();
		area.channels.push_back(ch);
		if (ParameterPile::draw_required(area)){
			std::string plot_name = "";
			plot_name += curr_area.experiments.back() + "_";
			plot_name += "run_" + std::to_string(curr_area.runs.back()) + "_ch_" + std::to_string(ch) + "_sub_" + std::to_string(area.sub_runs.back());
			int ind = curr_area.channels.get_order_index_by_index(ch);
			if (!xs_channels[ind].empty() && !ys_channels[ind].empty()){
				Drawing *dr = graph_manager.GetDrawing(plot_name, ch, de);
				if (NULL == dr)
					return;
				dr->AddToDraw(xs_channels[ind], ys_channels[ind], prefix, "with points pt 1", 0);
			}
		}
	}
}

bool SingleRunData::test_PMT_signal(int _N_threshold, double _S_threshold, double _S_max_threshold, SingleRunResults &results)
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

SingleRunResults SingleRunData::processSingleRun_Iter_0(AllRunsResults *all_runs_results)
{
#ifdef _USE_TIME_STATISTICS
	auto PMT_start_timer = std::chrono::high_resolution_clock::now();
#endif
	SingleRunResults _result(this);
#ifdef _USE_TIME_STATISTICS
	auto PMT_read_file_start_timer = std::chrono::high_resolution_clock::now();
#endif
	readOneRun(_result);
	for (int ch = curr_area.channels.get_next_index(); ch != -1; ch = curr_area.channels.get_next_index()){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if ((ch==8)||(12==ch))
		SignalOperations::invert_y(xs_channels[ind],ys_channels[ind]);
	}
#ifdef _USE_TIME_STATISTICS
	auto PMT_read_file_end_timer = std::chrono::high_resolution_clock::now();
	all_runs_results->time_stat.n_PMT_file_reading++;
	all_runs_results->time_stat.t_PMT_file_reading +=
		std::chrono::duration_cast<std::chrono::nanoseconds>(PMT_read_file_end_timer - PMT_read_file_start_timer).count();
#endif

	add_draw_data("PMT_raw" + curr_area.experiments.back() + "\\_" + std::to_string(curr_area.runs.back()) + "\\_sub" +
		std::to_string(curr_area.sub_runs.back()) + "\\_ch" + std::to_string(8), graph_manager, ParameterPile::DrawEngine::Gnuplot, 8);
	add_draw_data("PMT_raw" + curr_area.experiments.back() + "\\_" + std::to_string(curr_area.runs.back()) + "\\_sub" +
			std::to_string(curr_area.sub_runs.back()) + "\\_ch" + std::to_string(12), graph_manager, ParameterPile::DrawEngine::Gnuplot, 12);

	if (!_result.isValid()){
		runProcessedProc();
		return _result;
	}
#ifdef _USE_TIME_STATISTICS
	auto PMT_filter_start_timer = std::chrono::high_resolution_clock::now();
#endif
	curr_area.channels.reset();
	for (int ch = curr_area.channels.get_next_index(); ch != -1; ch = curr_area.channels.get_next_index()){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if ((ind < 0)||(ch>=32)||(2==ch))
			continue;
		bool valid_pars = ParameterPile::filter_PMT_n_points.find(ch)!=ParameterPile::filter_PMT_n_points.end();
		valid_pars = valid_pars&&(ParameterPile::filter_PMT_order.find(ch)!=ParameterPile::filter_PMT_order.end());
		valid_pars = valid_pars&&(ParameterPile::filter_PMT_n_iterations.find(ch)!=ParameterPile::filter_PMT_n_iterations.end());
		if (!valid_pars) //no filter applied
			continue;
		SavitzkyGolayFilter SGfilter(ParameterPile::filter_PMT_n_points.find(ch)->second, ParameterPile::filter_PMT_order.find(ch)->second,
				ParameterPile::filter_PMT_n_iterations.find(ch)->second);
		SGfilter(xs_channels[ind], ys_channels[ind]);
		/*if (8==ch){
			int sz = ys_channels[ind].size();
			double *in_p = new double [sz];
			for (int i=0;i<sz;++i)
				in_p[i] = ys_channels[ind][i];
			DVECTOR out_r(sz);
			DVECTOR out_i(sz);
			DVECTOR out_m(sz);

			TVirtualFFT * fft = TVirtualFFT::FFT(1, &sz,"R2C");
			fft->SetPoints(in_p);
			fft->Transform();

			for (int i=0;i<sz;++i){
				fft->GetPointComplex(i,out_r[i],out_i[i]);
				out_m[i] = std::sqrt(out_r[i]*out_r[i]+out_i[i]*out_i[i]);
			}

			ParameterPile::experiment_area area = curr_area.to_point();
			area.channels.erase();
			area.channels.push_back(ch);
			std::string plot_title = curr_area.experiments.back() + "\\_run" + std::to_string(curr_area.runs.back()) +
				"\\_ch" + std::to_string(ch) + "\\_sub" + std::to_string(curr_area.sub_runs.back());
			if (ParameterPile::draw_required(area)) {
				std::string plot_name = "";
				plot_name += curr_area.experiments.back() + "_";
				plot_name += "run" + std::to_string(curr_area.runs.back()) + "_ch" + std::to_string(ch) + "_sub" + std::to_string(area.sub_runs.back());
				Drawing *dr = graph_manager.GetDrawing(plot_name, ch+100, ParameterPile::DrawEngine::Gnuplot);
				if (NULL != dr) {
					dr->AddToDraw(xs_channels[ind], out_r, "FFT::Re" + plot_title, "with points pt 1", 0);
					dr->AddToDraw(xs_channels[ind], out_i, "FFT::Im" + plot_title, "with points pt 2", 0);
					dr->AddToDraw(xs_channels[ind], out_m, "FFT::|.|" + plot_title, "with points pt 3", 0);
				}
			}

			delete [] in_p;
		}*/
	} //apply the filter to PMTs
#ifdef _USE_TIME_STATISTICS
	auto PMT_filter_end_timer = std::chrono::high_resolution_clock::now();
	all_runs_results->time_stat.n_PMT_filtering++;
	all_runs_results->time_stat.t_PMT_filtering +=
		std::chrono::duration_cast<std::chrono::nanoseconds>(PMT_filter_end_timer - PMT_filter_start_timer).count();
#endif
//#ifndef _TEMP_CODE
	/*add_draw_data("filtered_PMT" + curr_area.experiments.back() + "\\_" + std::to_string(curr_area.runs.back()) + "\\_sub" +
		std::to_string(curr_area.sub_runs.back()) + "\\_ch" + std::to_string(0), graph_manager, ParameterPile::DrawEngine::Gnuplot, 0);
	add_draw_data("filtered_PMT" + curr_area.experiments.back() + "\\_" + std::to_string(curr_area.runs.back()) + "\\_sub" +
		std::to_string(curr_area.sub_runs.back()) + "\\_ch" + std::to_string(1), graph_manager, ParameterPile::DrawEngine::Gnuplot, 1);
	add_draw_data("filtered_PMT" + curr_area.experiments.back() + "\\_" + std::to_string(curr_area.runs.back()) + "\\_sub" +
		std::to_string(curr_area.sub_runs.back()) + "\\_ch" + std::to_string(8), graph_manager, ParameterPile::DrawEngine::Gnuplot, 8);
	add_draw_data("filtered_PMT" + curr_area.experiments.back() + "\\_" + std::to_string(curr_area.runs.back()) + "\\_sub" +
		std::to_string(curr_area.sub_runs.back()) + "\\_ch" + std::to_string(12), graph_manager, ParameterPile::DrawEngine::Gnuplot, 12);*/
	//#endif

	curr_area.channels.reset();
	for (int ch = curr_area.channels.get_next_index(); ch != -1; ch = curr_area.channels.get_next_index()){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if ((ind < 0)||(ch>=32)||(2==ch))
			continue;
		_result.pmt_peaks.push_back(STD_CONT<peak>());
		_result.pmt_channels.push_back(ch);
		_result.PMT_S2_integral.push_back(-1);
#ifdef _USE_TIME_STATISTICS
		auto PMT_baseline_start_timer = std::chrono::high_resolution_clock::now();
#endif
		STD_CONT<peak> ppeaks;
		double threshold = 0;
		calculate_PMT_threshold_and_baseline(xs_channels[ind], ys_channels[ind], threshold, found_base_lines[ind], ppeaks, ch);

#ifdef _USE_TIME_STATISTICS
		auto PMT_baseline_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_PMT_baseline++;
		all_runs_results->time_stat.t_PMT_baseline +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(PMT_baseline_end_timer - PMT_baseline_start_timer).count();
		auto PMT_peaks_start_timer = std::chrono::high_resolution_clock::now();
#endif

		SignalOperations::find_peaks_fine(xs_channels[ind], ys_channels[ind], _result.pmt_peaks.back(), found_base_lines[ind],
			threshold,found_base_lines[ind], ParameterPile::PMT_N_of_averaging);

		if (0==ch){
			double S_sum = 0;
			for (auto i = _result.pmt_peaks.back().begin(), _end_ = _result.pmt_peaks.back().end(); i != _end_; ++i)
				if ((i->t >= ParameterPile::S2_start_time.find(curr_area.experiments.back())->second)&&(i->t <= ParameterPile::S2_finish_time.find(curr_area.experiments.back())->second)){
					S_sum += i->S >0 ? i->S : 0;
					++_result.PMT3_n_peaks;
				}
			_result.PMT3_summed_peaks_area = S_sum;
			PMT3_summed_peaks_area = S_sum;
			PMT3_n_peaks = _result.PMT3_n_peaks;
		}
		double S2_st = ParameterPile::S2_start_time.find(curr_area.experiments.back())->second;
		double S2_ft = ParameterPile::S2_finish_time.find(curr_area.experiments.back())->second;
		DVECTOR ys_int, xs_int;
		SignalOperations::integrate(xs_channels[ind],ys_channels[ind],xs_int, ys_int, S2_st, S2_ft,
				(*(xs_channels[ind].begin()+1)-*(xs_channels[ind].begin())), found_base_lines[ind]);
		DITERATOR x_max;
		SignalOperations::get_max(xs_int, ys_int, x_max, _result.PMT_S2_integral.back(), 1);

		//=========================================================================
		ParameterPile::experiment_area area = curr_area.to_point();
		area.channels.erase();
		area.channels.push_back(ch);
		std::string plot_title = curr_area.experiments.back() + "\\_run" + std::to_string(curr_area.runs.back()) +
			"\\_ch" + std::to_string(ch) + "\\_sub" + std::to_string(curr_area.sub_runs.back());
		if (ParameterPile::draw_required(area)) {
			std::string plot_name = "";
			plot_name += curr_area.experiments.back() + "_";
			plot_name += "run" + std::to_string(curr_area.runs.back()) + "_ch" + std::to_string(ch) + "_sub" + std::to_string(area.sub_runs.back());
			Drawing *dr = graph_manager.GetDrawing(plot_name, ch, ParameterPile::DrawEngine::Gnuplot);
			if (NULL != dr) {
				dr->AddToDraw(xs_channels[ind], ys_channels[ind], "filtered" + plot_title, "with points pt 1", 0);
				dr->AddToDraw_baseline(threshold, "threshold");
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
	PMT_peaks_analysed = true;

	if (!test_PMT_signal(all_runs_results->N_peaks_cutoff, all_runs_results->S_peaks_cutoff, all_runs_results->S_peaks_max_cutoff, _result)) {
		_result._current_status = SingleRunResults::Status::NoPMTsignal;
		_result.setValid(false);
	}
//=============================================================================
#ifdef _USE_TIME_STATISTICS
	auto PMT_end_timer = std::chrono::high_resolution_clock::now();
	all_runs_results->time_stat.n_PMT_proc++;
	all_runs_results->time_stat.t_PMT_proc += std::chrono::duration_cast<std::chrono::nanoseconds>(PMT_end_timer - PMT_start_timer).count();
#endif
	runProcessedProc();
	return _result;
}

//can not be channle independent. In defference to MPPC this one returns threshold shifter by baseline
void SingleRunData::calculate_PMT_threshold_and_baseline(DVECTOR &xs, DVECTOR &ys, double &threshold, double &baseline, STD_CONT<peak> &peaks_before_S1, int channel)
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
}

SingleRunResults SingleRunData::processSingleRun_Iter_1(AllRunsResults *all_runs_results)
{
#ifdef _USE_TIME_STATISTICS
	auto MPPC_start_timer = std::chrono::high_resolution_clock::now();
#endif
	SingleRunResults _result(this);
	_result.setValid(true);
	_result._current_status = SingleRunResults::Status::Ok;
	//#ifdef _TEMP_CODE
	//return _result;
	//#endif
//#ifndef _TEMP_CODE
	SavitzkyGolayFilter SGfilter(ParameterPile::filter_MPPC_n_points, ParameterPile::filter_MPPC_order, ParameterPile::filter_MPPC_n_iterations);
//#endif
	if (!test_PMT_signal(all_runs_results->N_peaks_cutoff, all_runs_results->S_peaks_cutoff, all_runs_results->S_peaks_max_cutoff, _result)) {
		_result._current_status = SingleRunResults::Status::NoPMTsignal;
		_result.setValid(false);
		runProcessedProc();
		return _result;
	}
	for (int ch = 0; ch < 2; ch++) {
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if (ind < 0)
			continue;
		readOneRunMPPCs(_result, ch);
		if (ParameterPile::PMT_use_average.find(ch)!=ParameterPile::PMT_use_average.end()){
			if (ParameterPile::PMT_use_average.find(ch)->second){
				if (0==ch){
					_result.xs_PMT3=xs_channels[ind];
					_result.ys_PMT3=ys_channels[ind];
				} else {
					_result.xs_PMT1=xs_channels[ind];
					_result.ys_PMT1=ys_channels[ind];
				}
			}
		}
		clearOneRunMPPCs(ch);
	}

	int mppc_ind = -1;
	for (int ch = 32; ch < 64; ch++){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if (ind < 0)
			continue;
		++mppc_ind;
#ifdef _USE_TIME_STATISTICS
		auto MPPC_read_start_timer = std::chrono::high_resolution_clock::now();
#endif
		readOneRunMPPCs(_result, ch);
#ifdef _USE_TIME_STATISTICS
		auto MPPC_read_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_file_reading++;
		all_runs_results->time_stat.t_MPPC_file_reading+=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_read_end_timer - MPPC_read_start_timer).count();
#endif
		_result.mppc_baseline_xs.push_back(xs_channels[ind]); //TODO: remove, make x_start, x_finish and dx or Nx.
		_result.mppc_baseline_ys.push_back(DVECTOR());
		_result.mppc_peaks.push_back(STD_CONT<peak>());
		_result.mppc_S2_peaks_area.push_back(0);
		_result.mppc_S2_start_t.push_back(0);
		_result.mppc_S2_finish_t.push_back(0);
		//_result.mppc_sum_peaks_area.push_back(0);
		_result.mppc_channels.push_back(ch);
		_result.mppc_double_I.push_back(-1);

		std::string plot_title = curr_area.experiments.back() + "\\_run" + std::to_string(curr_area.runs.back()) +
			"\\_ch" + std::to_string(ch) + "\\_sub" + std::to_string(curr_area.sub_runs.back());
//==========================================================================
		// invert MPPCs
		SignalOperations::invert_y(xs_channels[ind], ys_channels[ind]);
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
		//TODO: add methods get{noise,baseline,threshold,edge_threshold} (ind, mppc_ind);
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
		double delta_x = *(_result.mppc_baseline_xs[mppc_ind].begin() + 1) - *(_result.mppc_baseline_xs[mppc_ind].begin());
		DITERATOR _end_ = _result.mppc_baseline_xs[mppc_ind].end() - 1;
		DITERATOR _begin_ = _result.mppc_baseline_xs[mppc_ind].begin();
		DITERATOR S2_start_i = SignalOperations::find_x_iterator_by_value(_begin_,
			_end_, ParameterPile::S2_start_time.find(curr_area.experiments.back())->second, delta_x);
		DITERATOR S2_finish_i = SignalOperations::find_x_iterator_by_value(S2_start_i,
			_end_, ParameterPile::S2_finish_time.find(curr_area.experiments.back())->second, delta_x);
		DITERATOR max_signal_pos;
		double max_sgnal_val;
		//TODO: ParameterPile
		++S2_finish_i;
		SignalOperations::get_max(_result.mppc_baseline_xs[mppc_ind], ys_cut, S2_start_i, S2_finish_i, max_signal_pos, max_sgnal_val, 1);
		double root_start_t = *max_signal_pos - ParameterPile::MPPC_ROOTs_bl_from_max_left;
		double root_finish_t = *max_signal_pos + ParameterPile::MPPC_ROOTs_bl_from_max_right;
		SignalOperations::apply_time_limits(_result.mppc_baseline_xs[mppc_ind], ys_cut, root_start_t, root_finish_t,delta_x);
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
			SignalOperations::find_baseline_by_ROOT(_result.mppc_baseline_xs[mppc_ind], ys_cut, _result.mppc_baseline_ys[mppc_ind]);
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
		ROOTs_baseline_baseline = SignalOperations::find_baseline_by_integral(0, _result.mppc_baseline_xs[mppc_ind], _result.mppc_baseline_ys[mppc_ind], exclude_middle);
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
		SignalOperations::substract_baseline(xs_channels[ind], ys_channels[ind], _result.mppc_baseline_xs[mppc_ind], 
			_result.mppc_baseline_ys[mppc_ind], ROOTs_baseline_baseline); //handles that 2 signals have different x spans
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
		SignalOperations::find_peaks_fine(xs_channels[ind], ys_channels[ind], _result.mppc_peaks.back(), 0, global_threshold, 0/*edge_threshold*/, ParameterPile::MPPC_N_trust);
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
		SignalOperations::peaks_to_yx(*(xs_channels[ind].begin()), xs_channels[ind].back(), _result.mppc_peaks.back(), x_peaks, y_peaks);
		DVECTOR x_peaks_spreaded_v2, y_peaks_spreaded_v2;
		SignalOperations::spread_peaks_v2(*(xs_channels[ind].begin()), xs_channels[ind].back(), _result.mppc_peaks.back(),
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
		sp_min = *std::min_element(y_peaks_spreaded_v2.begin(), y_peaks_spreaded_v2.end());
		double sp_approx_thr;
		double sp_threshold = find_spreaded_peaks_threshold(x_peaks_spreaded_v2, y_peaks_spreaded_v2, sp_approx_thr);
		if (sp_threshold < 0){
			clearOneRunMPPCs(ch);
			continue;
		}
		STD_CONT<peak> sp_peaks;
		//no fine required here
		SignalOperations::find_peaks(x_peaks_spreaded_v2, y_peaks_spreaded_v2, sp_peaks, sp_min, sp_threshold, 1);
		//TODO: don't forget to test this^ place
		if (sp_peaks.empty()){
			clearOneRunMPPCs(ch);
			continue;
		}
		STD_CONT<peak>::iterator max_intersection; //search for largest cluster which is above treshold. In 99% if not 100% I 
		//suspect there will be only single sp_peaks i.e. intersection, but just in case.
		std::string __exp = curr_area.experiments.back();
		max_intersection = std::max_element(sp_peaks.begin(), sp_peaks.end(), 
			[__exp](peak a, peak b) {
			if ((b.left > ParameterPile::S2_finish_time.find(__exp)->second) || (b.right < ParameterPile::S2_start_time.find(__exp)->second))
				return false;
			if ((a.left > ParameterPile::S2_finish_time.find(__exp)->second) || (a.right < ParameterPile::S2_start_time.find(__exp)->second))
				return true;
			return a.S < b.S; 
		});
		double S2_finish_t_v2 = std::min(max_intersection->right, ParameterPile::S2_finish_time.find(curr_area.experiments.back())->second);
		double S2_start_t_v2 = std::max(max_intersection->left, ParameterPile::S2_start_time.find(curr_area.experiments.back())->second);
		//TODO: ParameterPile
		double S2_finish_t = std::min(max_intersection->right + ParameterPile::MPPC_peaks_smoothing_time,
				ParameterPile::S2_finish_time.find(curr_area.experiments.back())->second);
		double S2_start_t = std::max(max_intersection->left - ParameterPile::MPPC_peaks_smoothing_time,
				ParameterPile::S2_start_time.find(curr_area.experiments.back())->second);
		if (S2_finish_t < S2_start_t) {
			S2_finish_t = 0;
			S2_start_t = 0;
		}
		_result.mppc_S2_finish_t.back() = S2_finish_t;
		_result.mppc_S2_start_t.back() = S2_start_t;
		//TODO: decide whether this is required. ParameterPile
		for (auto pp = _result.mppc_peaks.back().begin(); pp != _result.mppc_peaks.back().end(); ++pp){
			if ((pp->left >= S2_start_t) && (pp->right <= S2_finish_t)){
				_result.mppc_S2_peaks_area.back() += pp->S;
			}
		}

#ifdef _USE_TIME_STATISTICS
		auto MPPC_peaks_proc_end_timer = std::chrono::high_resolution_clock::now();
		all_runs_results->time_stat.n_MPPC_peaks_processing++;
		all_runs_results->time_stat.t_MPPC_peaks_processing +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_peaks_proc_end_timer - MPPC_peaks_proc_start_timer).count();
#endif
//================================================================================
		//Now calculate double integral for signal without substracting ROOT's baseline
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
				_result.mppc_double_I.back() = ys_raw_ii.back() - y_ii_min;
			else
				_result.mppc_double_I.back() = y_ii_max - y_ii_min;
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
			std::string plot_name = "";
			plot_name += curr_area.experiments.back() + "_";
			plot_name += "run" + std::to_string(curr_area.runs.back()) + "_ch" + std::to_string(ch) + "_sub" + std::to_string(area.sub_runs.back());
			int ind = curr_area.channels.get_order_index_by_index(ch);

			Drawing *dr = graph_manager.GetDrawing(plot_name, ch, ParameterPile::DrawEngine::Gnuplot);
			if (NULL != dr) {
				dr->AddToDraw(xs_raw, ys_raw, "raw" + plot_title, "with points pt 1", 0);
				dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], _result.mppc_baseline_ys[mppc_ind], "ROOT baseline V1", "w l lw 2", 0);
#ifdef _USE_TIME_STATISTICS
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v2, "ROOT baseline V2", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v3, "ROOT baseline V3", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v4, "ROOT baseline V4", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v5, "ROOT baseline V5", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v6, "ROOT baseline V6", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v7, "ROOT baseline V7", "w l lw 2", 0);
				//dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], root_baseline_v8, "ROOT baseline V8", "w l lw 2", 0);
#endif
				dr->AddToDraw(xs_channels[ind], ys_channels[ind], "without baseline " + std::to_string(curr_area.runs.back()), "with lines lc rgb \"#0000FF\"", 0);
				dr->AddToDraw_baseline(global_threshold, "threshold", "w l lc rgb \"#00FF00\"");
				dr->AddToDraw_baseline(ROOTs_baseline_baseline, "ROOTs baseline", "w l lc rgb \"#0000FF\"");
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
				dr->AddToDraw_vertical(S2_start_t_v2, -1, 1, "lc rgb \"#EE82EE\"");
				dr->AddToDraw_vertical(S2_finish_t_v2, -1, 1, "lc rgb \"#EE82EE\"");
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
		clearOneRunMPPCs(ch);
	}
	STD_CONT<DVECTOR>().swap(_result.mppc_baseline_xs);
	STD_CONT<DVECTOR>().swap(_result.mppc_baseline_ys);
#ifdef _USE_TIME_STATISTICS
	auto MPPC_end_timer = std::chrono::high_resolution_clock::now();
	all_runs_results->time_stat.n_MPPC_proc+=(mppc_ind+1);
	all_runs_results->time_stat.t_MPPC_proc += std::chrono::duration_cast<std::chrono::nanoseconds>(MPPC_end_timer - MPPC_start_timer).count();
#endif

#ifdef _PROCESS_GEMS
	int ind = curr_area.channels.get_order_index_by_index(2);
	if (ind >= 0){
		readOneRunMPPCs(_result, 2);
#ifdef GEM_V1_
		STD_CONT<peak> no_peaks;
		no_peaks.push_back(peak());
		no_peaks.back().left = ParameterPile::S1_start_time;
		no_peaks.back().right = xs_channels[ind].back();
		found_base_lines[ind] =
			SignalOperations::find_baseline_by_integral(found_base_lines[ind], xs_channels[ind], ys_channels[ind], no_peaks);
		//memory for found_base_lines is allocated in the constructor
		SignalOperations::substract_baseline(ys_channels[ind], found_base_lines[ind]);
		found_base_lines[ind] = 0; //since it is substructed already
#endif //GEM_V1_
		_result.xs_GEM = xs_channels[ind];
		_result.ys_GEM = ys_channels[ind];
		clearOneRunMPPCs(2);
	}
	if (_result.xs_GEM.empty() || _result.ys_GEM.empty() || _result.xs_GEM.size() != _result.ys_GEM.size()){
		_result._current_status = SingleRunResults::Status::NoGEMsignal;
		_result.setValid(false);
		runProcessedProc();
		return _result;
	}
#endif

	runProcessedProc();
	return _result;
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
	//TODO: move code below to separate method as well
	//variant0: appprox_thresh===(max+min)/2 -> {mean,RMS above ; mean,RMS below} -> spereaded peaks threshold -> S2 time
	//variant1: appprox_thresh===meidan -> {mean,RMS above ; mean,RMS below} -> spereaded peaks threshold -> S2 time
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

	double above_weight = std::sqrt(RMS_below) / std::pow(n_below,1.5);//no real explanation for this, just figuring out algorithm by trieal and error.
	double below_weight = std::sqrt(RMS_above) / std::pow(n_above,1.5); //But in the case when statistics (n of points) is small, the dispersion is forcebly increased
	//TODO: maybe use this threshold only for S2 finish time and some other for S2 start time
	//TODO: try other weightings.
	return (mean_above*above_weight + mean_below * below_weight) / (above_weight + below_weight);
}

SingleRunResults SingleRunData::processSingleRun(AllRunsResults *all_runs_results)
{
	if (0 == all_runs_results->Iteration())
		return processSingleRun_Iter_0(all_runs_results);
	if (1 == all_runs_results->Iteration())
		return processSingleRun_Iter_1(all_runs_results);
	return SingleRunResults(this);
}

void SingleRunData::runProcessedProc(void)
{
	graph_manager.Draw();
	graph_manager.Clear();
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


