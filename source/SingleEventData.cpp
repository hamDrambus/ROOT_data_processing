#include "SingleEventData.h"
#include <string>
#include "Savitzky_Golay_filter.h"

SingleEventData::SingleEventData(const ParameterPile::experiment_manifest *to_process) : manifest(to_process), status(Status::Empty)
{
	for (std::size_t i = 0, i_end_ = manifest->channels.size(); i != i_end_; ++i) {
		channel_data.push(manifest->channels.index(i), ChannelData());
	}
	index.run = manifest->runs.back();
	index.subrun = manifest->sub_runs.back();
	
	trigger_offset = 0;
}

void SingleEventData::readOneRun(AllEventsResults *results, int channel)
{
	ChannelData *data = channel_data.info(channel);
	if (NULL == data)
		return;
	std::string path = manifest->in_folder;
	path += "run_" + std::to_string(index.run) + "__ch_" + std::to_string(channel) + ".dat";
	file_to_vector(path, data->xs, data->ys, index.subrun);
	if (!data->xs.empty()&&(Status::Empty==status))
		status = Status::Ok;
}

void SingleEventData::clearOneRun(int channel)
{
	ChannelData *data = channel_data.info(channel);
	if (NULL == data)
		return;
	DVECTOR().swap(data->xs);
	DVECTOR().swap(data->ys);
	DVECTOR().swap(data->curved_bl_xs);
	DVECTOR().swap(data->curved_bl_ys);
}

void SingleEventData::file_to_vector(std::string fname, DVECTOR &xs, DVECTOR &ys, int index)
{
	DVECTOR().swap(xs);
	DVECTOR().swap(ys);
	std::ifstream file;
	file.open(fname, std::ios_base::binary | std::ios_base::ate);
	if (!file.is_open() || file.eof()) {
		file.close();
		return;
	}
	long long Size = file.tellg();
#ifndef _USE_DEQUE
	xs.reserve((Size / manifest->subruns_per_file)/2);
	ys.reserve((Size / manifest->subruns_per_file)/2);
#endif
	long long offset = index * (Size / (manifest->subruns_per_file));
	long long to_read = Size / (manifest->subruns_per_file);
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
			//TODO: use ParameterPile instead of #defines
			val = manifest->data_voltage_amplitude*(val / manifest->data_voltage_channels) + manifest->data_voltage_of_zero_channel;
			xs.push_back(read_N*manifest->data_time_constant + trigger_offset);
			ys.push_back(val);
			++read_N;
		}
	}
	file.close();
	return;
}

void SingleEventData::push_event (AllEventsResults *all_runs_results)
{
	all_runs_results->events_data.push_back(*this); //TODO: really stupid architechture! Have to rework AnalysisManager's handling of AllEventsResults and SingleEventData.
}

void SingleEventData::processSingleEvent_Iter_0(AllEventsResults *all_runs_results)
{
	/*	Process PMTs and SiPMs only.
	 *		1) Inversion (if required)
	 *		2) Filter
	 *		3) Finding baseline&thresholds for finding peaks. Thresholds are set in manifest by hand
	 *		4) Restore curved baseline if necessary
	 *		5) Finding peaks, taking integrals of signals in preset time windows. Add to draw later if necessary.
	 */
	push_event(all_runs_results);
	if (!manifest->accepted_events_data.empty()) {
		trigger_offset = manifest->accepted_events_data(index.run, index.subrun);
		if (std::numeric_limits<double>::max() == trigger_offset) {
			status = Status::ExternalRejected;
		}
	}
	for (std::size_t ch_ind = 0, ch_end_ = manifest->channels.size(); ch_ind != ch_end_ && isValid(); ++ch_ind) { //corresponding channel data was initialized accordingly. Also both must be sorted by channel
		int ch = manifest->channels.index(ch_ind);
		readOneRun(all_runs_results, ch); //read all PMTs, ignore GEM and MPPCs
		if (channel_data[ch_ind].xs.empty())
			continue;
		double delta_x = *(channel_data[ch_ind].xs.begin() + 1) - *(channel_data[ch_ind].xs.begin());
		DVECTOR ys_raw, ys_filtered;
		double middle_left, middle_right;
#ifdef _TEMP_CODE
		DVECTOR pmt_baseline_ys1, pmt_baseline_ys2, pmt_baseline_ys3, pmt_baseline_ys4, pmt_baseline_ys5, pmt_baseline_ys6, pmt_baseline_ys7;
#endif
		if (manifest->channels[ch_ind].invert)
			SignalOperations::invert_y(channel_data[ch_ind].xs, channel_data[ch_ind].ys); //invert fast PMT and SiPM signals
		ys_raw = channel_data[ch_ind].ys;

		if (manifest->channels[ch_ind].find_average)
			push_average(ch, all_runs_results);

		//apply the filter
		{
			SavitzkyGolayFilter SGfilter(manifest->channels[ch_ind].filter.n_points, manifest->channels[ch_ind].filter.order,
				manifest->channels[ch_ind].filter.n_iterations);
			if (SGfilter.getNIter() > 0) {
				SGfilter(channel_data[ch_ind].xs, channel_data[ch_ind].ys);
				ys_filtered = channel_data[ch_ind].ys;
			}
		}
		STD_CONT<peak> peaks_before_S1;
		double threshold = 0, threshold_edges = 0;
		calculate_threshold_and_baseline(channel_data[ch_ind].xs, channel_data[ch_ind].ys, threshold, threshold_edges, channel_data[ch_ind].found_baseline, peaks_before_S1, ch);

		//Find integral if necessary
		if (manifest->channels[ch_ind].find_integral) {
			DVECTOR ys_int, xs_int;
			SignalOperations::integrate(channel_data[ch_ind].xs, channel_data[ch_ind].ys, xs_int, ys_int,
				manifest->channels[ch_ind].integral_range.first, manifest->channels[ch_ind].integral_range.second,
				delta_x, channel_data[ch_ind].found_baseline);
			DITERATOR x_max;
			SignalOperations::get_max(xs_int, ys_int, x_max, channel_data[ch_ind].integral, manifest->channels[ch_ind].N_extrapolation);
		}
		//Find double integral if necessary (can have another range entirely)
		if (manifest->channels[ch_ind].find_double_integral) {
			DVECTOR xs_int = channel_data[ch_ind].xs, ys_int = channel_data[ch_ind].ys, ys_i, ys_ii;
			SignalOperations::apply_time_limits(xs_int, ys_int, manifest->channels[ch_ind].double_integral_range.first, manifest->channels[ch_ind].double_integral_range.second, delta_x);
			DITERATOR x_ii_max = xs_int.end();
			double y_ii_max, y_ii_min;
			if (xs_int.size() >= 2) {
				SignalOperations::integrate(xs_int, ys_int, ys_i, delta_x, channel_data[ch_ind].found_baseline);
				SignalOperations::integrate(xs_int, ys_i, ys_ii, delta_x, 0);
				//TODO: ParameterPile
				SignalOperations::get_max(xs_int, ys_ii, x_ii_max, y_ii_max, manifest->channels[ch_ind].N_extrapolation);
				y_ii_min = ys_ii.front();
				channel_data[ch_ind].double_integral = y_ii_max - y_ii_min;
			}
		}
		//Find curved baseline (conduct baseline restoration):
		if (manifest->channels[ch_ind].baseline.do_find_curved) {
			channel_data[ch_ind].curved_bl_xs = channel_data[ch_ind].xs;
			DVECTOR ys_cut = channel_data[ch_ind].ys;
			double root_start_t = manifest->channels[ch_ind].baseline.curved_range.first; //TODO: add check for validity of parameters (automatic search for invalid?)
			double root_finish_t = manifest->channels[ch_ind].baseline.curved_range.second;
			SignalOperations::apply_time_limits(channel_data[ch_ind].curved_bl_xs, ys_cut, root_start_t, root_finish_t, delta_x);

			SignalOperations::find_baseline_by_ROOT(channel_data[ch_ind].curved_bl_xs, ys_cut, channel_data[ch_ind].curved_bl_ys,
				manifest->channels[ch_ind].baseline.curved_numberIterations, manifest->channels[ch_ind].baseline.curved_direction,
				manifest->channels[ch_ind].baseline.curved_filterOrder, manifest->channels[ch_ind].baseline.curved_smoothing,
				manifest->channels[ch_ind].baseline.curved_smoothWindow, manifest->channels[ch_ind].baseline.curved_compton,
				manifest->channels[ch_ind].baseline.curved_sparse);
			root_start_t = std::max(root_start_t, *channel_data[ch_ind].curved_bl_xs.begin());
			root_finish_t = std::min(root_finish_t, *(channel_data[ch_ind].curved_bl_xs.end()-1));
			root_start_t += manifest->channels[ch_ind].baseline.curved_trim.first;
			root_finish_t -= manifest->channels[ch_ind].baseline.curved_trim.second;
			SignalOperations::apply_time_limits(channel_data[ch_ind].curved_bl_xs, channel_data[ch_ind].curved_bl_ys, root_start_t, root_finish_t, delta_x);
		//================================================================================
			STD_CONT<peak> exclude_middle(1, peak());
			middle_left = exclude_middle.back().left = manifest->channels[ch_ind].baseline.curved_center.first;
			middle_right = exclude_middle.back().right = manifest->channels[ch_ind].baseline.curved_center.second;
			channel_data[ch_ind].curved_bl_baseline = SignalOperations::find_baseline_by_integral(0, channel_data[ch_ind].curved_bl_xs, channel_data[ch_ind].curved_bl_ys, exclude_middle);
		//================================================================================
			//This fixes baseline of signal when ROOT's baseline is subtracted before the start of strong curvature
			DITERATOR min_signal_pos;
			double min_signal_val;
			SignalOperations::get_min(channel_data[ch_ind].curved_bl_xs, channel_data[ch_ind].curved_bl_ys,
				channel_data[ch_ind].curved_bl_xs.begin(), channel_data[ch_ind].curved_bl_xs.end(), min_signal_pos, min_signal_val, 1);
			if (channel_data[ch_ind].curved_bl_xs.end() != min_signal_pos) {
				double x_min = std::min(*min_signal_pos, middle_right);
				for (std::size_t i = 0, i_end_ = std::min(channel_data[ch_ind].curved_bl_ys.size(), channel_data[ch_ind].curved_bl_xs.size());
					i != i_end_ && channel_data[ch_ind].curved_bl_xs[i] < x_min; ++i) {
					if (channel_data[ch_ind].curved_bl_ys[i] > channel_data[ch_ind].curved_bl_baseline)
						channel_data[ch_ind].curved_bl_ys[i] = channel_data[ch_ind].curved_bl_baseline;
				}
			}
		//================================================================================
			//So in result I have to do 2 baseline subtractions no matter what.
			//One of them is hidden in SignalOperations::substract_baseline(DV&,DV&,DV&,DV&,double), so it doesn't cost anything
			SignalOperations::substract_baseline(channel_data[ch_ind].xs, channel_data[ch_ind].ys, channel_data[ch_ind].curved_bl_xs,
				channel_data[ch_ind].curved_bl_ys, channel_data[ch_ind].curved_bl_baseline); //handles that 2 signals have different x spans
			//TODO: optimize this function, using the fact that the ROOT's baseline has definite x points [from, to] (- [xs.begin(),xs.end())
		//================================================================================
		}
		SignalOperations::find_peaks_fine(channel_data[ch_ind].xs, channel_data[ch_ind].ys, channel_data[ch_ind].peaks, channel_data[ch_ind].found_baseline,
			threshold, threshold_edges, manifest->channels[ch_ind].N_extrapolation);
		//=========================================================================
		std::string plot_title = manifest->name + "_run" + std::to_string(index.run) +
			"_ch" + std::to_string(ch) + "_sub" + std::to_string(index.subrun);
		if (manifest->do_draw(index.run, index.subrun) && manifest->channels[ch_ind].display.do_draw) {
			std::string plot_name = "";
			plot_name += manifest->channels[ch_ind].device + std::to_string(ch);
			plot_name += "_run" + std::to_string(index.run) + "_sub" + std::to_string(index.subrun);

			GraphCollection *pics = all_runs_results->pictures.info(ch);
			if (NULL == pics) {
				all_runs_results->pictures.push(ch, GraphCollection());
				pics = all_runs_results->pictures.info(ch);
				pics->SetGnuplotDirectory(manifest->out_gnuplot_folder);
				pics->SetPngDirectory(manifest->out_picture_folder);
			}
			GnuplotDrawing *dr = pics->GetDrawing(plot_name);
			if (NULL != dr) {
				if (channel_data[ch_ind].curved_bl_xs.empty()) { //simple baseline case
					dr->AddToDraw(channel_data[ch_ind].xs, ys_raw, "raw_" + plot_title, "with lines");
					if (!ys_filtered.empty())
						dr->AddToDraw(channel_data[ch_ind].xs, ys_filtered, "filtered_" + plot_title, "with lines lw 2");
					dr->AddToDraw_baseline(threshold, "threshold", "w l lc rgb \"#FF0000\"");
					if (threshold_edges != channel_data[ch_ind].found_baseline)
						dr->AddToDraw_baseline(threshold_edges, "threshold 2nd", "w l lc rgb \"#AC0ECD\"");
					dr->AddToDraw_baseline(channel_data[ch_ind].found_baseline, "baseline", "w l lc rgb \"#0000FF\"");
					//dr->AddToDraw_vertical(S2_st, -1, 1, "lc rgb \"#0000FF\"");
					//dr->AddToDraw_vertical(S2_ft, -1, 1, "lc rgb \"#0000FF\"");
					//if (!ys_int.empty()) {
					//	dr->AddToDraw(xs_int, ys_int, "integral", "axis x1y2 with lines lw 2 lc rgb \"#FF00FF\"");
					//}
				} else { //curved baseline case
					if (!ys_filtered.empty()) {
						dr->AddToDraw(channel_data[ch_ind].xs, ys_filtered, "filtered_" + plot_title, "w l lc rgb \"#0000FF\"");
						//dr->AddToDraw(xs_channels[ind], ys_raw, "raw_" + plot_title, "w l lc rgb \"#FF00FF\"");
					} else
						dr->AddToDraw(channel_data[ch_ind].xs, ys_raw, "raw_" + plot_title, "w l lc rgb \"#0000FF\"");
					dr->AddToDraw(channel_data[ch_ind].curved_bl_xs, channel_data[ch_ind].curved_bl_ys, "ROOT baseline V0", "w l lw 2 lc rgb \"#00FF00\"");
					dr->AddToDraw(channel_data[ch_ind].xs, channel_data[ch_ind].ys, "without baseline " + std::to_string(index.run), "with lines lc rgb \"#000000\"");
					dr->AddToDraw_baseline(threshold, "threshold", "w l lc rgb \"#FF0000\"");
					dr->AddToDraw_baseline(channel_data[ch_ind].curved_bl_baseline, "ROOTs baseline", "w l lc rgb \"#FF33FF\"");
					dr->AddToDraw_baseline(channel_data[ch_ind].found_baseline, "baseline", "w l lc rgb \"#0000FF\"");
					if (threshold_edges != channel_data[ch_ind].found_baseline)
						dr->AddToDraw_baseline(threshold_edges, "threshold\\_2nd");
					//dr->AddToDraw_vertical(S2_st, -1, 1, "lc rgb \"#FF0000\"");
					//dr->AddToDraw_vertical(S2_ft, -1, 1, "lc rgb \"#FF0000\"");
					dr->AddToDraw_vertical(middle_left, -1, 1, "lc rgb \"#0000FF\"");
					dr->AddToDraw_vertical(middle_right, -1, 1, "lc rgb \"#0000FF\"");
					//if (!ys_int.empty()) {
					//	dr->AddToDraw(xs_int, ys_int, "integral", "axis x1y2 with lines lw 2 lc rgb \"#FF00FF\"");
					//}
				}
			}
		}
		clearOneRun(ch);
	}

//=============================================================================

	eventProcessedProc(all_runs_results);
	return;
}

//can not be channel independent. In difference to MPPC this one returns threshold shifted by baseline
void SingleEventData::calculate_threshold_and_baseline(DVECTOR &xs, DVECTOR &ys, double &threshold, double &threshold_2, double &baseline, STD_CONT<peak> &peaks_before_S2, int channel)
{
	DVECTOR xs_before_S2 = xs;
	DVECTOR ys_before_S2 = ys;

	const ParameterPile::channel_manifest *man = manifest->channels.info(channel); //guaranteed to be not NULL

	double delta_x = *(xs.begin() + 1) - *xs.begin();
	SignalOperations::apply_time_limits(xs_before_S2, ys_before_S2, man->baseline.baseline_range.first, man->baseline.baseline_range.second, delta_x);

	/*approx*/baseline = man->baseline.baseline_by_average ? SignalOperations::find_baseline_by_integral(0, xs_before_S2, ys_before_S2)
		: SignalOperations::find_baseline_by_median(0, xs_before_S2, ys_before_S2);
	//calculate first-order baseline
	threshold = man->peaks.threshold;
	SignalOperations::find_peaks_fine(xs_before_S2, ys_before_S2, peaks_before_S2,
		baseline, threshold + baseline, baseline, man->N_extrapolation);
	if (!peaks_before_S2.empty()) {
		/*exact*/baseline = SignalOperations::find_baseline_by_integral(0, xs_before_S2, ys_before_S2, peaks_before_S2);
	}
	threshold += baseline;
	threshold_2 = baseline + man->peaks.threshold_cutoff;
}

void SingleEventData::push_average (int ch, AllEventsResults *all_events_results)
{
	ChannelData *ch_data = channel_data.info(ch); //guaranteed to be not NULL
	AllEventsResults::AverageData* avr_data = all_events_results->averages.info(ch);
	if (NULL == avr_data) {
		all_events_results->averages.push(ch, AllEventsResults::AverageData());
		avr_data = all_events_results->averages.info(ch);
	}
	if (avr_data->xs_sum.empty()) {
		avr_data->xs_sum = ch_data->xs;
		avr_data->ys_sum = ch_data->ys;
	} else {
		for (auto j = avr_data->ys_sum.begin(), i = ch_data->ys.begin(), j_end_ = avr_data->ys_sum.end(), i_end_ = ch_data->ys.end();
			(i != i_end_) && (j != j_end_); ++i, ++j) {
			*i += *j;
		}
	}
	++avr_data->average_event_n;
}

void SingleEventData::push_dispersion (int ch, AllEventsResults *all_events_results)
{
	ChannelData *ch_data = channel_data.info(ch); //guaranteed to be not NULL
	AllEventsResults::AverageData* avr_data = all_events_results->averages.info(ch); //guaranteed to be not NULL
	if (avr_data->ys_disp.empty())
		avr_data->ys_disp.resize(avr_data->ys_sum.size(), 0.0);
	for (auto j = ch_data->ys.begin(), a = avr_data->ys_sum.begin(), d = avr_data->ys_disp.begin(),
		j_end_ = ch_data->ys.end(), a_end_ = avr_data->ys_sum.end(), d_end_ = avr_data->ys_disp.end();
		(d != d_end_) && (a != a_end_) && (j != j_end_); ++a, ++j, ++d) {
		*d += std::pow(*j - *a, 2);
	}
	++avr_data->dispersion_event_n;
}

void SingleEventData::processSingleEvent_Iter_1(AllEventsResults *all_runs_results)
{
	if (!isValid())
		goto l_quit;
	for (std::size_t ch_ind = 0, ch_end_ = manifest->channels.size(); ch_ind != ch_end_ && isValid(); ++ch_ind) { //corresponding channel data was initialized accordingly. Also both must be sorted by channel
		int ch = manifest->channels.index(ch_ind);
		readOneRun(all_runs_results, ch);
		if (channel_data[ch_ind].xs.empty())
			continue;

		if (manifest->channels[ch_ind].invert)
			SignalOperations::invert_y(channel_data[ch_ind].xs, channel_data[ch_ind].ys); //invert fast PMT and SiPM signals

		if (manifest->channels[ch_ind].find_average)
			push_dispersion(ch, all_runs_results);
		clearOneRun(ch);
	}
	l_quit:
	eventProcessedProc(all_runs_results);
	return;
}

void SingleEventData::processSingleEvent(AllEventsResults *all_runs_results)
{
	if (0 == all_runs_results->Iteration())
		processSingleEvent_Iter_0(all_runs_results);
	if (1 == all_runs_results->Iteration())
		processSingleEvent_Iter_1(all_runs_results);
}

void SingleEventData::eventProcessedProc(AllEventsResults *all_runs_results)
{
	clear_memory();
}

void SingleEventData::clear_memory(void) //clears only 'input' data, preserves processing results
{
	for (std::size_t ch_ind = 0, ch_end_ = channel_data.size(); ch_ind != ch_end_; ++ch_ind) {
		DVECTOR().swap(channel_data[ch_ind].xs);
		DVECTOR().swap(channel_data[ch_ind].ys);
		DVECTOR().swap(channel_data[ch_ind].curved_bl_xs);
		DVECTOR().swap(channel_data[ch_ind].curved_bl_ys);
	}
}

std::size_t SingleEventData::real_size(void)
{
	std::size_t rsize = sizeof(*this);
	rsize += sizeof(indexed_info<ChannelData>);
	rsize += channel_data.size()*sizeof(ChannelData);
	for (std::size_t ch_ind = 0, ch_end_ = channel_data.size(); ch_ind != ch_end_; ++ch_ind) {
		rsize += channel_data[ch_ind].xs.capacity() * sizeof(double);
		rsize += channel_data[ch_ind].ys.capacity() * sizeof(double);
		rsize += channel_data[ch_ind].curved_bl_xs.capacity() * sizeof(double);
		rsize += channel_data[ch_ind].curved_bl_ys.capacity() * sizeof(double);
	}
	return rsize;
}


