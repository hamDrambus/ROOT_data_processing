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
		found_base_lines.push_back(ParameterPile::baseline_approx_value[curr_area.channels.get_order_index_by_index(ch)]);
	}
	PMT3_summed_peaks_area = 0;
	PMT3_n_peaks = 0;
	PMT_peaks_analysed = false;
}

void SingleRunData::readOneRun(SingleRunResults &_result)
{
	xs_channels.clear();
	ys_channels.clear();
	bool empty_run = true;
	for (int ch = curr_area.channels.get_next_index(); !(ch < 0); ch = curr_area.channels.get_next_index()){
		xs_channels.push_back(std::vector<double>());
		ys_channels.push_back(std::vector<double>());
		std::string path = DATA_PREFIX;
		path += "event_x-ray_" + curr_area.experiments.back() + "\\";
		path += "run_" + std::to_string(curr_area.runs.back()) + "__ch_" + std::to_string(ch) + ".dat";
		file_to_vector(path, xs_channels.back(), ys_channels.back(), curr_area.sub_runs.back());
		if (!xs_channels.back().empty() && !ys_channels.empty())
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

void SingleRunData::file_to_vector(std::string fname, std::vector<double> &xs, std::vector<double> &ys, int index)
{
	xs.clear();
	ys.clear();
	std::ifstream file;
	file.open(fname, std::ios_base::binary | std::ios_base::ate);
	if (!file.is_open() || file.eof()){
		file.close();
		return;
	}
	int Size = file.tellg();
	int offset = index * Size / (ParameterPile::subruns_per_file);
	int to_read = Size / (ParameterPile::subruns_per_file);
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
			read_N++;
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

void SingleRunData::add_draw_data(std::string prefix, GraphicOutputManager& graph_manager, ParameterPile::DrawEngine de)
{
	for (int ch = curr_area.channels.get_next_index(); !(ch < 0); ch = curr_area.channels.get_next_index()) {
		ParameterPile::experiment_area area = curr_area.to_point();
		area.channels.erase();
		area.channels.push_back(ch);
		if (ParameterPile::draw_required(area)){
			std::string plot_name = "";
			plot_name += curr_area.experiments.back() + "_";
			plot_name += "run_" + std::to_string(curr_area.runs.back()) + "_ch_" + std::to_string(ch) + "_sub_" + std::to_string(area.sub_runs.back());
			int ind = curr_area.channels.get_order_index_by_index(ch);

			Drawing *dr = graph_manager.GetDrawing(plot_name, ch, de);
			if (NULL == dr)
				return;
			/*if (prefix == "in_")
				dr->AddToDraw(xs_channels[ind], ys_channels[ind], prefix, 1);
				dr->AddToDraw(xs_channels[ind], ys_channels[ind], prefix, prefix=="in_" ? 0:1);*/
			dr->AddToDraw(xs_channels[ind], ys_channels[ind], prefix, "", 0);
		}
	}
}

bool SingleRunData::test_PMT_signal(int _N_threshold, double _S_threshold, double _S_max_threshold, SingleRunResults &results)
{
	int ind = curr_area.channels.get_order_index_by_index(0);
	if (ind < 0)
		return true; //can't test, so assuming the PMT signal is correct
	
	if (!PMT_peaks_analysed){ //if threshold is constant, then there is no need to recalculate this code
		DVECTOR xs_before_S1 = xs_channels[ind], ys_before_S1 = ys_channels[ind];
		SignalOperations::apply_time_limits(xs_before_S1, ys_before_S1, *(xs_before_S1.begin()), ParameterPile::S1_start_time);
		double noise_amp;
		double PMT_baseline = found_base_lines[ind];
		DITERATOR x_max_noise;
		SignalOperations::get_max(xs_before_S1, ys_before_S1, x_max_noise, noise_amp, ParameterPile::PMT_N_of_averaging);
		noise_amp -= PMT_baseline;
		double threshold = noise_amp*ParameterPile::PMT_run_acceptance_threshold_to_noize;
		threshold += PMT_baseline;

		std::vector<peak> peaks;
		SignalOperations::find_peaks(xs_channels[ind], ys_channels[ind], peaks, PMT_baseline, threshold, ParameterPile::PMT_N_of_averaging);
		double S_sum = 0;
		for (auto i = peaks.begin(); i != peaks.end(); ++i)
			S_sum += (*i).S;
		results.PMT3_summed_peaks_area = S_sum;
		results.PMT3_n_peaks = peaks.size();
		PMT3_summed_peaks_area = S_sum;
		PMT3_n_peaks = peaks.size();
		PMT_peaks_analysed = true;
	}
	if (_S_max_threshold < _S_threshold){ //means there is no maximum limit
		if ((PMT3_summed_peaks_area >= _S_threshold) && (PMT3_n_peaks >= _N_threshold))
			return true;
	} else {
		if ((PMT3_summed_peaks_area >= _S_threshold) && (PMT3_summed_peaks_area <= _S_max_threshold) && (PMT3_n_peaks >= _N_threshold))
			return true;
	}
	return false;
}

SingleRunResults SingleRunData::processSingleRun_Iter_0(const AllRunsResults *all_runs_results)
{
	SingleRunResults _result(this);
	readOneRun(_result);

	if (!_result.isValid())
		return _result;

	SavitzkyGolayFilter SGfilter(ParameterPile::filter_PMT_n_points, ParameterPile::filter_PMT_order, ParameterPile::filter_PMT_n_iterations);
#ifdef _TEMP_CODE
	SavitzkyGolayFilter MPPC_SGfilter(ParameterPile::filter_MPPC_n_points, ParameterPile::filter_MPPC_order, ParameterPile::filter_MPPC_n_iterations);
#endif
	for (int ch = 0; ch < 2; ch++){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if (ind < 0)
			continue;
		std::vector<double> xs, ys;
		SGfilter(xs_channels[ind], ys_channels[ind], xs, ys);
		xs_channels[ind] = xs;
		ys_channels[ind] = ys;
	} //apply the filter to PMTs
#ifndef _TEMP_CODE
	add_draw_data("filtered_",graph_manager);
#endif

#ifdef _TEMP_CODE
	for (int ch = 32; ch < 64; ch++){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if (ind >= 0)
			SignalOperations::invert_y(xs_channels[ind], ys_channels[ind]);
	} // invert MPPCs
	//add_draw_data("in_",graph_manager);
	for (int ch = 32; ch < 64; ch++){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if (ind < 0)
			continue;
		std::vector<double> xs, ys;
		SGfilter(xs_channels[ind], ys_channels[ind], xs, ys);
		xs_channels[ind] = xs;
		ys_channels[ind] = ys;
		_result.mppc_baseline_xs.push_back(xs_channels[ind]);
		_result.mppc_baseline_ys.push_back(DVECTOR());
	}
	add_draw_data("filtered_", graph_manager);
#endif
	for (int ch = 32, mppc_ind = 0; ch < 64; ch++){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if (ind < 0)
			continue;
		_result.mppc_peaks.push_back(std::vector<peak>());
		_result.mppc_S2_peaks_area.push_back(0);
		_result.mppc_S2_start_t.push_back(0);
		_result.mppc_S2_finish_t.push_back(0);
		_result.mppc_sum_peaks_area.push_back(0);
		//TODO: add methods get{noise,baseline,threshold,edge_threshold} (ind, mppc_ind);
		DVECTOR before_S1_x = xs_channels[ind];
		DVECTOR before_S1_y = ys_channels[ind];
		SignalOperations::apply_time_limits(before_S1_x, before_S1_y, *before_S1_x.begin(), ParameterPile::S1_start_time);
		found_base_lines[ind] = SignalOperations::find_baseline_by_integral(found_base_lines[ind], before_S1_x, before_S1_y);
		//^not many peaks expected, so integral is used, not a median
		SignalOperations::substract_baseline(before_S1_y, found_base_lines[ind]);
		SignalOperations::substract_baseline(ys_channels[ind], found_base_lines[ind]);
		double approx_noise = TMath::RMS(before_S1_y.begin(), before_S1_y.end());
		double approx_thresh_beforeS1 = approx_noise * 5;	//TODO: ParameterPile
		double approx_edge_thresh = approx_noise*0.5;
		std::vector<peak> before_S1_peaks;
		//TODO: ParameterPile N_trust;
		SignalOperations::find_peaks_fine(before_S1_x, before_S1_y, before_S1_peaks,
			0, approx_thresh_beforeS1, approx_edge_thresh, 10);
		/*exact*/found_base_lines[ind] = SignalOperations::find_baseline_by_integral(0, before_S1_x, before_S1_y, before_S1_peaks);
		SignalOperations::exclude_peaks(before_S1_x, before_S1_y, before_S1_peaks);
		double exact_noise = TMath::RMS(before_S1_y.begin(), before_S1_y.end());
		double global_threshold = exact_noise * 5;//global means for all times, not for all channels and runs //TODO: ParameterPile
		//double global_ROOT_threshold = exact_noise * 3.5;	//TODO: ParameterPile. ROOT's baseline is significantly lower that the actual one
		double edge_threshold = exact_noise*0;	//TODO: ParameterPile

		SignalOperations::substract_baseline(ys_channels[ind], found_base_lines[ind]);
		found_base_lines[ind] = 0;

		SignalOperations::find_baseline_by_ROOT(xs_channels[ind], ys_channels[ind], _result.mppc_baseline_ys[mppc_ind]);
		before_S1_x = xs_channels[ind];
		before_S1_y = _result.mppc_baseline_ys[mppc_ind];
		SignalOperations::apply_time_limits(before_S1_x, before_S1_y, *before_S1_x.begin(), ParameterPile::S1_start_time);
		found_base_lines[ind] = SignalOperations::find_baseline_by_integral(0, before_S1_x, before_S1_y, before_S1_peaks);
		//now accounting for ROOT's baseline shift
		SignalOperations::substract_baseline(_result.mppc_baseline_ys[mppc_ind], found_base_lines[ind]);
		found_base_lines[ind] = 0;
		SignalOperations::substract_baseline(ys_channels[ind], _result.mppc_baseline_ys[mppc_ind]);//both functions have baseline at 0 for the first 32us


		SignalOperations::find_peaks_fine(xs_channels[ind], ys_channels[ind], _result.mppc_peaks.back(), 0, global_threshold, edge_threshold, 10);
		//^TODO: ParameterPile
		DVECTOR x_peaks, y_peaks;
		SignalOperations::peaks_to_yx(*(xs_channels[ind].begin()), xs_channels[ind].back(), _result.mppc_peaks.back(), x_peaks, y_peaks);
		DVECTOR x_peaks_spreaded, y_peaks_spreaded;
		DVECTOR x_peaks_spreaded_v2, y_peaks_spreaded_v2;
		SignalOperations::spread_peaks(*(xs_channels[ind].begin()), xs_channels[ind].back(), _result.mppc_peaks.back(),x_peaks_spreaded, y_peaks_spreaded);
		SignalOperations::spread_peaks_v2(*(xs_channels[ind].begin()), xs_channels[ind].back(), _result.mppc_peaks.back(), x_peaks_spreaded_v2, y_peaks_spreaded_v2,5);
		
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

		//TODO: move code below to separate method as well
		//variant0: appprox_thresh===(max+min)/2 -> {mean,RMS above ; mean,RMS below} -> spereaded peaks threshold -> S2 time
		//variant1: appprox_thresh===meidan -> {mean,RMS above ; mean,RMS below} -> spereaded peaks threshold -> S2 time
		//variant2: spereaded peaks threshold === median //I think it's bad choice
		double sp_max, sp_min;
		sp_max = *std::max_element(y_peaks_spreaded_v2.begin(), y_peaks_spreaded_v2.end());
		sp_min = *std::min_element(y_peaks_spreaded_v2.begin(), y_peaks_spreaded_v2.end());
		double sp_approx_thr = (sp_max + sp_min) / 2;
		double n_below=0, n_above=0;
		double mean_below=0, mean_above=0;
		for (auto i = y_peaks_spreaded_v2.begin(); i != y_peaks_spreaded_v2.end(); ++i) {
			if (*i < sp_approx_thr){
				n_below++;
				mean_below += *i;
			} else {
				n_above++;
				mean_above += *i;
			}
		}
		if ((n_below == 0) || (n_above == 0)) {
			continue;
		}
		mean_above /= n_above;
		mean_below /= n_below;
		double RMS_below = 0, RMS_above = 0;
		for (auto i = y_peaks_spreaded_v2.begin(); i != y_peaks_spreaded_v2.end(); ++i)
			if (*i < sp_approx_thr)
				RMS_below += (*i - mean_below)*(*i - mean_below);
			else
				RMS_above += (*i - mean_above)*(*i - mean_above);
		RMS_below = sqrt(RMS_below) / n_below; //maybe sqrt(RMS_below/(n_below*(n_below-1))); not that is matters much
		RMS_above = sqrt(RMS_above) / n_above;
		
		//TODO: maybe use this threshold only for S2 finish time and some other for S2 start time
		double sp_threshold = (mean_above*RMS_below + mean_below * RMS_above) / (RMS_above + RMS_below);
		std::vector<peak> sp_peaks;
		//no fine required here
		SignalOperations::find_peaks(x_peaks_spreaded_v2, y_peaks_spreaded_v2, sp_peaks, sp_min, sp_threshold, 1);
		//TODO: don't forget to test this^ place
		if (sp_peaks.empty())
			continue;
		std::vector<peak>::iterator max_intersection; //search for largest cluster which is above treshold. In 99% if not 100% I 
		//suspect there will be only single sp_peaks i.e. intersection, but just in case.
		max_intersection = std::max_element(sp_peaks.begin(), sp_peaks.end(), [](peak a, peak b){return a.S < b.S; });
		double S2_finish_t = std::min(max_intersection->right, ParameterPile::S2_finish_time);
		double S2_start_t = std::max(max_intersection->left, ParameterPile::S2_start_time);
		if (S2_finish_t < S2_start_t) {
			S2_finish_t = 0;
			S2_start_t = 0;
		}
		_result.mppc_S2_finish_t.back() = S2_finish_t;
		_result.mppc_S2_start_t.back() = S2_start_t;
		//TODO: decide whether this is required. ParameterPile
		double S2_finish_t_v2 = S2_finish_t + 5;
		double S2_start_t_v2 = S2_start_t - 5;
		for (auto pp = _result.mppc_peaks.back().begin(); pp != _result.mppc_peaks.back().end(); ++pp){
			if ((pp->left >= S2_start_t) && (pp->right <= S2_finish_t)){
				_result.mppc_S2_peaks_area.back() += pp->S;
			}
		}


#ifdef _TEMP_CODE
		ParameterPile::experiment_area area = curr_area.to_point();
		area.channels.erase();
		area.channels.push_back(ch);
		if (ParameterPile::draw_required(area)){
			std::string plot_name = "";
			plot_name += curr_area.experiments.back() + "_";
			plot_name += "run_" + std::to_string(curr_area.runs.back()) + "_ch_" + std::to_string(ch) + "_sub_" + std::to_string(area.sub_runs.back());
			int ind = curr_area.channels.get_order_index_by_index(ch);

			Drawing *dr = graph_manager.GetDrawing(plot_name, ch, ParameterPile::DrawEngine::Gnuplot);
			if (NULL != dr) {
				dr->AddToDraw(_result.mppc_baseline_xs[mppc_ind], _result.mppc_baseline_ys[mppc_ind], "ROOT baseline", "w l", 0);
				dr->AddToDraw(xs_channels[ind], ys_channels[ind], "without baseline "+std::to_string(curr_area.runs.back()), "", 0);
				dr->AddToDraw_baseline(global_threshold, "threshold");
				dr->AddToDraw_baseline(edge_threshold, "threshold\\_edges");
				dr->AddToDraw_vertical(S2_start_t, -1, 1, "lc rgb \"#FF0000\"");
				dr->AddToDraw_vertical(S2_finish_t, -1, 1, "lc rgb \"#FF0000\"");
			}
			dr = graph_manager.GetDrawing(plot_name + "peaks", ch + 100, ParameterPile::DrawEngine::Gnuplot);
			if (NULL != dr) {
				dr->AddToDraw(x_peaks, y_peaks, "peaks I = " + std::to_string(_result.mppc_S2_peaks_area.back()), "w l", 0);
				dr->AddToDraw(x_peaks_spreaded, y_peaks_spreaded, "peaks spereaded " + std::to_string(curr_area.runs.back()), "w l", 0);
				dr->AddToDraw(x_peaks_spreaded_v2, y_peaks_spreaded_v2, "peaks spereaded v2 " + std::to_string(curr_area.runs.back()), "w l", 0);
				dr->AddToDraw_baseline(sp_approx_thr, "approx\\_threshold", "w l");
				dr->AddToDraw_baseline(sp_threshold, "exact\\_threshold", "w l");
				dr->AddToDraw_vertical(S2_start_t, -1, 1, "lc rgb \"#FF0000\"");
				dr->AddToDraw_vertical(S2_finish_t, -1, 1, "lc rgb \"#FF0000\"");
				dr->AddToDraw_vertical(S2_start_t_v2, -1, 1, "lc rgb \"#EE82EE\"");
				dr->AddToDraw_vertical(S2_finish_t_v2, -1, 1, "lc rgb \"#EE82EE\"");
			}
		}
#endif
		mppc_ind++;
	}
	//int ind = curr_area.channels.get_order_index_by_index(0);
	//if (!(ind < 0)){ //find baseline only for the sum of PMTs for the first run
	//	std::vector<peak> no_peaks;
	//	found_base_lines[ind] = SignalOperations::find_baseline_by_median(found_base_lines[ind],
	//		xs_channels[ind], ys_channels[ind], no_peaks);
	//}
	//add_draw_baselines("base filtered", graph_manager);

	if (!test_PMT_signal(all_runs_results->N_peaks_cutoff, all_runs_results->S_peaks_cutoff, all_runs_results->S_peaks_max_cutoff, _result)) {
		_result._current_status = SingleRunResults::Status::NoPMTsignal;
		_result.setValid(false);
	}

	runProcessedProc();
	return _result;
}

SingleRunResults SingleRunData::processSingleRun_Iter_1(const AllRunsResults *all_runs_results)
{
	SingleRunResults _result(this);
	_result.setValid(true);
	_result._current_status = SingleRunResults::Status::Ok;
#ifndef _TEMP_CODE
	SavitzkyGolayFilter SGfilter(ParameterPile::filter_MPPC_n_points, ParameterPile::filter_MPPC_order, ParameterPile::filter_MPPC_n_iterations);
#endif
	if (!test_PMT_signal(all_runs_results->N_peaks_cutoff, all_runs_results->S_peaks_cutoff, all_runs_results->S_peaks_max_cutoff, _result)) {
		_result._current_status = SingleRunResults::Status::NoPMTsignal;
		_result.setValid(false);
		goto end_proc;
	}
#ifndef _TEMP_CODE
	for (int ch = 32; ch < 64; ch++){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if (ind >= 0)
			SignalOperations::invert_y(xs_channels[ind], ys_channels[ind]);
	} // invert MPPCs
	//add_draw_data("in_",graph_manager);
	for (int ch = 32; ch < 64; ch++){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if (ind < 0)
			continue;
		std::vector<double> xs, ys;
		SGfilter(xs_channels[ind], ys_channels[ind], xs, ys);
		xs_channels[ind] = xs;
		ys_channels[ind] = ys;
	}
	add_draw_data("filtered_", graph_manager);
#endif

	//get base lines
	for (auto i = xs_channels.begin(), j = ys_channels.begin(); (i != xs_channels.end()) && (j != ys_channels.end()); i++, j++){
		std::vector<peak> no_peaks;
		if ((i - xs_channels.begin()) == curr_area.channels.get_order_index_by_index(2)){//GEM
			peak before_S1;
			before_S1.left = ParameterPile::S1_start_time;
			before_S1.right = (*i).back();
			before_S1.S = 0;
			no_peaks.push_back(before_S1);
			found_base_lines[i - xs_channels.begin()] =
				SignalOperations::find_baseline_by_integral(found_base_lines[i - xs_channels.begin()], *i, *j, no_peaks);
			//memory for found_base_lines is allocated in the constructor
			for (auto g = (*j).begin(); g != (*j).end(); ++g)
				*g -= found_base_lines[i - xs_channels.begin()];
			found_base_lines[i - xs_channels.begin()] = 0; //since it is substructed already
			continue;
		}
		found_base_lines[i - xs_channels.begin()] = SignalOperations::find_baseline_by_median(found_base_lines[i - xs_channels.begin()],
			*i, *j, no_peaks);//memory is allocated in the constructor
	}
#ifndef _TEMP_CODE
	add_draw_baselines("base filtered", graph_manager);
#endif

#ifdef _PROCESS_GEMS
	int ind = curr_area.channels.get_order_index_by_index(2);
	if (ind >= 0){
		_result.xs_GEM = xs_channels[ind];
		_result.ys_GEM = ys_channels[ind];
	}
	if (_result.xs_GEM.empty() || _result.ys_GEM.empty() || _result.xs_GEM.size() != _result.ys_GEM.size()){
		_result._current_status = SingleRunResults::Status::NoGEMsignal;
		_result.setValid(false);
		goto end_proc;
	}
#endif

end_proc:
	runProcessedProc();
	return _result;
}

SingleRunResults SingleRunData::processSingleRun(const AllRunsResults *all_runs_results)
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
	//clear_memory();
}

void SingleRunData::clear_memory(void) //clears only 'input' data, preserves processing results
{
	xs_channels.clear();
	ys_channels.clear();
}

ParameterPile::experiment_area SingleRunData::getArea(void) const {	return curr_area;}
