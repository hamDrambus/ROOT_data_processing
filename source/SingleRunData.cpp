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
		ParameterPile::experiment_area area(ParameterPile::experiment_area::Point);
		area.experiments = curr_area.experiments;
		area.runs = curr_area.runs;
		area.sub_runs = curr_area.sub_runs;
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
		SignalOperations::apply_time_limits(xs_before_S1, ys_before_S1, *(xs_before_S1.begin()), ParameterPile::S1_time);
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
	for (int ch = 0; ch < 2; ch++){
		int ind = curr_area.channels.get_order_index_by_index(ch);
		if (ind < 0)
			continue;
		std::vector<double> xs, ys;
		SGfilter(xs_channels[ind], ys_channels[ind], xs, ys);
		xs_channels[ind] = xs;
		ys_channels[ind] = ys;
	} //apply the filter to PMTs
	add_draw_data("filtered_",graph_manager);

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
	SavitzkyGolayFilter SGfilter(ParameterPile::filter_MPPC_n_points, ParameterPile::filter_MPPC_order, ParameterPile::filter_MPPC_n_iterations);

	if (!test_PMT_signal(all_runs_results->N_peaks_cutoff, all_runs_results->S_peaks_cutoff, all_runs_results->S_peaks_max_cutoff, _result)) {
		_result._current_status = SingleRunResults::Status::NoPMTsignal;
		_result.setValid(false);
		goto end_proc;
	}

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

	//get base lines
	for (auto i = xs_channels.begin(), j = ys_channels.begin(); (i != xs_channels.end()) && (j != ys_channels.end()); i++, j++){
		std::vector<peak> no_peaks;
		if ((i - xs_channels.begin()) == curr_area.channels.get_order_index_by_index(2)){//GEM
			peak before_S1;
			before_S1.left = ParameterPile::S1_time;
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
	add_draw_baselines("base filtered", graph_manager);

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
	//clear_memory();
}

void SingleRunData::clear_memory(void) //clears only 'input' data, preserves processing results
{
	xs_channels.clear();
	ys_channels.clear();
}

ParameterPile::experiment_area SingleRunData::getArea(void) const {	return curr_area;}
