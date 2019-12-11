#include "AllEventsResults.h"

AllEventsResults::AllEventsResults(const ParameterPile::experiment_manifest* to_process) : processing_manifest(to_process)
{
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

void AllEventsResults::Merge(AllEventsResults* with)
{
	if (ParameterPile::Max_iteration_N == Iteration_N) {//Merge per-event data only after the last iteration (was stored in individual AllEventsResults's before)
		events_data.insert(events_data.end(), with->events_data.begin(), with->events_data.end());
	}
	if (0 == Iteration_N) {
		if (averages.empty()) {
			averages = with->averages;
		} else {
			if (!averages.isSameIndices(with->averages)) {
				std::cerr << "AllEventsResults::Merge: Warning: average data channel mismatch!" << std::endl;
			}
			for (std::size_t i = 0, i_end_ = with->averages.size(); i != i_end_; ++i) {
				int channel = with->averages.index(i);
				AverageData* data = averages.info(channel);
				if (NULL == data) {
					averages.push(channel, with->averages[i]);
				} else {
					data->average_event_n += with->averages[i].average_event_n;
					if (data->xs_sum.size() != with->averages[i].xs_sum.size()) {
						std::cerr << "AllEventsResults::Merge: Warning: average data size mismatch for" << std::endl
							<< "\tchannel " << channel << " (" << data->xs_sum.size() << " vs " << with->averages[i].xs_sum.size() << ")!" << std::endl;
					}
					for (auto y1 = data->ys_sum.begin(), y1_end_ = data->ys_sum.end(), y2 = with->averages[i].ys_sum.begin(), y2_end_ = with->averages[i].ys_sum.end();
						(y1 != y1_end_) && (y2 != y2_end_); ++y1, ++y2)
						*y1 += *y2;
				}
			}
		}
	}
	
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
	
	if (1 == Iteration_N) {
		if (averages.empty()) {
			averages = with->averages;
		} else {
			if (!averages.isSameIndices(with->averages)) {
				std::cerr << "AllEventsResults::Merge: Warning: average data channel mismatch!" << std::endl;
			}
			for (std::size_t i = 0, i_end_ = with->averages.size(); i != i_end_; ++i) {
				int channel = with->averages.index(i);
				AverageData* data = averages.info(channel);
				if (NULL == data) {
					averages.push(channel, with->averages[i]);
				} else {
					data->dispersion_event_n += with->averages[i].dispersion_event_n;
					if (data->xs_sum.size() != with->averages[i].xs_sum.size()) {
						std::cerr << "AllEventsResults::Merge: Warning: average data size mismatch for" << std::endl
							<< "\tchannel " << channel << " (" << data->xs_sum.size() << " vs " << with->averages[i].xs_sum.size() << ")!" << std::endl;
					}
					for (auto y1 = data->ys_disp.begin(), y1_end_ = data->ys_disp.end(), y2 = with->averages[i].ys_disp.begin(), y2_end_ = with->averages[i].ys_disp.end();
						(y1 != y1_end_) && (y2 != y2_end_); ++y1, ++y2)
						*y1 += *y2;
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

void AllEventsResults::Merged(void)
{
	std::cout << "Iteration " << Iteration_N << std::endl;
	std::cout << "N of events " << events_data.size() << std::endl;
	std::size_t N_of_valid_runs = 0;
	for (std::size_t i = 0, i_end_ = events_data.size(); i != i_end_; ++i)
		N_of_valid_runs += events_data[i].isValid() ? 1 : 0;
	std::cout << "N of valid events " << N_of_valid_runs << std::endl;
	if (0 == Iteration_N) {
		for (std::size_t i = 0, i_end_ = averages.size(); i != i_end_; ++i) 
			for (auto y = averages[i].ys_sum.begin(), y_end_ = averages[i].ys_sum.end(); y!=y_end_; ++y)
				*y/=(double)averages[i].average_event_n;
	}
	if (1 == Iteration_N) {
		for (std::size_t i = 0, i_end_ = averages.size(); i != i_end_; ++i)
			for (auto y = averages[i].ys_disp.begin(), y_end_ = averages[i].ys_disp.end(); y != y_end_; ++y)
				*y /= (double)averages[i].dispersion_event_n;
	}

	if (ParameterPile::Max_iteration_N == Iteration_N) {
		for (std::size_t i = 0, i_end_ = pictures.size(); i!=i_end_; ++i) {
			int channel = pictures.index(i);
			std::pair<double, double> y_range = pictures[i].get_y_limits();
			const ParameterPile::channel_manifest *man = processing_manifest->channels.info(channel);
			if (man!=NULL) {
				pictures[i].SetXrange(man->display.X_limits.first, man->display.X_limits.second);
				y_range.first = (-DBL_MAX!=man->display.Y_limits.first ? man->display.Y_limits.first: y_range.first);
				y_range.second = (DBL_MAX!=man->display.Y_limits.second ? man->display.Y_limits.second: y_range.second);
				pictures[i].SetYrange(y_range.first, y_range.second);
			} else {
				pictures[i].SetYrange(y_range.first, y_range.second);
			}
			std::pair<double, double> x_range = pictures[i].get_x_limits();
			if (processing_manifest->trigger_at >= x_range.first && processing_manifest->trigger_at <= x_range.second) {
				for (std::size_t p = 0, p_end_ = pictures[i].size(); p!=p_end_; ++p)
					pictures[i].GetDrawing(p)->AddToDraw_vertical(processing_manifest->trigger_at, -DBL_MAX, DBL_MAX, "lc rgb \"#FF0000\"");
			}
			pictures[i].Draw();
		}

		for (std::size_t ind = 0, ind_end_ = averages.size(); ind!=ind_end_; ++ind) {
			int channel = averages.index(ind);
			const ParameterPile::channel_manifest *man = processing_manifest->channels.info(channel);
			if (NULL == man) {
				std::cout<<"AllEventsResults::Merged:averages: Error:"<<std::endl;
				std::cout<<"\tNo manifest found, skipping channel "<<channel<<std::endl;
				continue;
			}
			std::string output_prefix = std::string(ParameterPile::this_path) + processing_manifest->out_folder +  man->device + "_" + std::to_string(channel) + "/";
			DVECTOR xs_before_S2 = averages[ind].xs_sum;
			DVECTOR ys_before_S2 = averages[ind].ys_sum;
			double delta_x = *(xs_before_S2.begin() + 1) - *xs_before_S2.begin();
			SignalOperations::apply_time_limits(xs_before_S2, ys_before_S2, man->baseline.baseline_range.first, man->baseline.baseline_range.second, delta_x);
			double baseline = SignalOperations::find_baseline_by_integral(0, xs_before_S2, ys_before_S2);
			SignalOperations::substract_baseline(averages[ind].ys_sum, baseline);
			DVECTOR integral, integral_variance;
			SignalOperations::integrate_with_variance(averages[ind].xs_sum, averages[ind].ys_sum, averages[ind].ys_disp, integral, integral_variance, 0);

			GnuplotDrawing* dr = graph_manager.GetDrawing(man->device + processing_manifest->name+"_ch_"+std::to_string(channel)+"_AVR");
			if (processing_manifest->draw_only)
				dr->SetGnuplotDirectory(processing_manifest->out_gnuplot_folder);
			else
				dr->SetGnuplotDirectory(output_prefix);
			dr->SetPngDirectory(processing_manifest->out_picture_folder);
			dr->AddToDraw(averages[ind].xs_sum, averages[ind].ys_sum, averages[ind].ys_disp, man->device + processing_manifest->name+"_ch_"+std::to_string(channel)+"_AVR");
			dr->AddToDraw(averages[ind].xs_sum, integral, integral_variance, man->device + processing_manifest->name+"_ch_"+std::to_string(channel)+"_Int_AVR", "axes x1y2");
			dr->SetXrange(man->display.X_limits.first, man->display.X_limits.second);


		}
	}

	if (ParameterPile::Max_iteration_N == Iteration_N) {//Merged per-event data only after the last iteration (was stored in individual AllEventsResults's before)
		std::cout<<"Merged "<<processing_manifest->in_folder<<std::endl;
		if (processing_manifest->draw_only) {
			std::cout<<"(Draw only mode)"<<std::endl;
		}
		auto max = std::max_element(events_data.begin(), events_data.end(), [](const SingleEventData &a, const SingleEventData &b)->bool{
			if (!a.isValid())
				return true;
			if (!b.isValid())
				return false;
			return a.channel_data.size()<b.channel_data.size();
		});
		std::cout<<"\tReference event: run_"<<max->index.run<<"_sub_"<<max->index.subrun<<std::endl;
		indexed_info<int> all_channels;
		for (std::size_t c = 0, c_end_ = max->channel_data.size(); c!=c_end_; ++c)
			all_channels.push(max->channel_data.index(c), 0);

		for (auto e = events_data.begin(), e_end_ = events_data.end(); e!=e_end_; ++e) {
			if (!e->isValid() && e->status != SingleEventData::ExternalRejected)
				continue;
			if (!max->channel_data.isSameIndices(e->channel_data)) {
				std::cerr<<"\tWarning! Channels mismatch with reference for event run_"<<e->index.run<<"_sub_"<<e->index.subrun<<std::endl;
				for (std::size_t c = 0, c_end_ = e->channel_data.size(); c!=c_end_; ++c)
					all_channels.push(e->channel_data.index(c), 0);
			}
		}

		for (std::size_t c = 0, c_end_ = all_channels.size(); c!=c_end_; ++c) {
			int channel = all_channels.index(c);
			std::cout<<"\tWriting channel "<<channel<<"... ";
			const ParameterPile::channel_manifest *man = processing_manifest->channels.info(channel);
			if (NULL == man) {
				std::cout<<"ERROR"<<std::endl;
				std::cout<<"\tNo manifest found, skipping channel "<<channel<<std::endl;
				continue;
			}
			bool all_ok = true;
			std::string error = "";
			std::size_t peaks_n = 0;
			std::size_t integral_n = 0;
			std::size_t double_integral_n = 0;
			for (auto e = events_data.begin(), e_end_ = events_data.end(); e!=e_end_; ++e) {
				if (!e->isValid() && e->status != SingleEventData::ExternalRejected)
					continue;
				SingleEventData::ChannelData* ch_data = e->channel_data.info(channel);
				if (NULL!=ch_data) {
					if (man->peaks.do_find)
						++peaks_n;
					if (man->find_integral)
						++integral_n;
					if (man->find_double_integral)
						++double_integral_n;
				}
			}
			if (peaks_n!=0 && integral_n!=0)
				if (peaks_n!=integral_n) {
					all_ok = false;
					error += "\tpeaks and integral events size mismatch\n";
				}
			if (peaks_n!=0 && double_integral_n!=0)
				if (peaks_n!=double_integral_n) {
					all_ok = false;
					error = "\tpeaks and double integral events size mismatch\n";
				}
			if (integral_n!=0 && double_integral_n!=0)
				if (integral_n!=double_integral_n) {
					all_ok = false;
					error = "\tintegral and double integral events size mismatch\n";
				}
			if (!all_ok) {
				std::cout<<"WARNING"<<std::endl;
				std::cout<<error;
			} else {
				if (processing_manifest->draw_only)
					std::cout<<"DRAW ONLY MODE"<<std::endl;
				else if(0==peaks_n && 0==integral_n && 0==double_integral_n)
						std::cout<<"NONE"<<std::endl;
					else
						std::cout<<"OK"<<std::endl;
			}
			if (processing_manifest->draw_only) {
				continue;
			}
			std::string output_prefix = std::string(ParameterPile::this_path) + processing_manifest->out_folder +  man->device + "_" + std::to_string(channel) + "/";
			if (peaks_n!=0) {
				std::string fname = output_prefix + man->device + "_" + std::to_string(channel)+"_peaks.dat";
				std::ofstream str;
				open_output_file(fname, str, std::ios_base::trunc | std::ios_base::binary);
				str.write((char*)&peaks_n, sizeof(std::size_t));
				for (auto e = events_data.begin(), e_end_ = events_data.end(); e!=e_end_; ++e) {
					if (!e->isValid() && e->status != SingleEventData::ExternalRejected)
						continue;
					SingleEventData::ChannelData* ch_data = e->channel_data.info(channel);
					if (NULL!=ch_data) {
						if (man->peaks.do_find) {
							if (processing_manifest->write_event_indices) {
								str.write((char*)&e->index.run, sizeof(int));
								str.write((char*)&e->index.subrun, sizeof(int));
							}
							std::size_t pk_end_ = ch_data->peaks.size();
							str.write((char*)&pk_end_, sizeof(std::size_t));
							for (std::size_t p = 0; p != pk_end_; ++p) {
								str.write((char*)&ch_data->peaks[p].left, sizeof(double));
								str.write((char*)&ch_data->peaks[p].right, sizeof(double));
								str.write((char*)&ch_data->peaks[p].S, sizeof(double));
								str.write((char*)&ch_data->peaks[p].A,sizeof(double));
								str.write((char*)&ch_data->peaks[p].t,sizeof(double));
							}
						}
					}
				}
				str.close();
			}
			if (integral_n!=0) {
				std::string fname = output_prefix + man->device + "_" + std::to_string(channel)+"_I.dat";
				std::ofstream str;
				open_output_file(fname, str, std::ios_base::trunc | std::ios_base::binary);
				str.write((char*)&integral_n, sizeof(std::size_t));
				for (auto e = events_data.begin(), e_end_ = events_data.end(); e!=e_end_; ++e) {
					if (!e->isValid() && e->status != SingleEventData::ExternalRejected)
						continue;
					SingleEventData::ChannelData* ch_data = e->channel_data.info(channel);
					if (NULL!=ch_data) {
						if (man->find_integral) {
							if (processing_manifest->write_event_indices) {
								str.write((char*)&e->index.run, sizeof(int));
								str.write((char*)&e->index.subrun, sizeof(int));
							}
							str.write((char*)&ch_data->integral, sizeof(double));
						}
					}
				}
				str.close();
			}
			if (double_integral_n!=0) {
				std::string fname = output_prefix + man->device + "_" + std::to_string(channel)+"_double_I.dat";
				std::ofstream str;
				open_output_file(fname, str, std::ios_base::trunc | std::ios_base::binary);
				str.write((char*)&double_integral_n, sizeof(std::size_t));
				for (auto e = events_data.begin(), e_end_ = events_data.end(); e!=e_end_; ++e) {
					if (!e->isValid() && e->status != SingleEventData::ExternalRejected)
						continue;
					SingleEventData::ChannelData* ch_data = e->channel_data.info(channel);
					if (NULL!=ch_data) {
						if (man->find_double_integral) {
							if (processing_manifest->write_event_indices) {
								str.write((char*)&e->index.run, sizeof(int));
								str.write((char*)&e->index.subrun, sizeof(int));
							}
							str.write((char*)&ch_data->double_integral, sizeof(double));
						}
					}
				}
				str.close();
			}
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

int AllEventsResults::Iteration(void) const
{	return Iteration_N;}

void AllEventsResults::Clear(void)
{
	//events_data preserve
	//averages preserve
	//Iteration_N preserve;

	if (ParameterPile::Max_iteration_N == Iteration()) {
		graph_manager.Clear();
		pictures.clear();
		STD_CONT<SingleEventData>().swap(events_data);
		averages.clear();
	}

#ifdef _USE_TIME_STATISTICS
	time_stat.t_RUN_proc_single_iter=0;
	time_stat.n_RUN_proc_single_iter=0;
#endif
}

void AllEventsResults::ClearMerged(void)
{
}

AllEventsResults& AllEventsResults::operator=(const AllEventsResults& right)
{
	//Only these are passed from merged AllRuns to thread-specific ones.
	Iteration_N = right.Iteration_N;
	processing_manifest = right.processing_manifest;
	averages = right.averages;
	return *this;
}

#ifdef _USE_TIME_STATISTICS
void AllEventsResults::report_time_statistics()
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
