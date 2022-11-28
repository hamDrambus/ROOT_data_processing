#include "GlobalParameters.h"

namespace ParameterPile
{
	//Use this function to play with/display different configurations
	//The correct parameters used in the analysis are in Init190404(analysis_manifest& manifest)
	bool Init221020_tests(analysis_manifest& manifest)
	{
		std::cout << "Initializing TEST 221020 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/221020/";
		default_exp_manifest.out_folder = "../Data/221020/results_vt/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init221020_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Data/221020/results_vt/gnuplot/";
		default_exp_manifest.out_picture_folder = "";
		default_exp_manifest.draw_only = true; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs): ~6 mm EL gap at 1.00 atm pressure
		//default_exp_manifest.runs.push(775, 789); //f13 18.5kV
		//default_exp_manifest.runs.push(754, 773); //f12 18.0kV
		//default_exp_manifest.runs.push(732, 752); //f11 17.5kV
		//default_exp_manifest.runs.push(708, 730); //f10 17.0kV
		//default_exp_manifest.runs.push(676, 706); //f9 16.5kV
		//default_exp_manifest.runs.push(638, 674); //f8 16.0kV
		//default_exp_manifest.runs.push(597, 636); //f7 15.5kV
		//default_exp_manifest.runs.push(534, 595); //f6 15.0kV
		//default_exp_manifest.runs.push(473, 532); //f5 14.0kV
		//default_exp_manifest.runs.push(357, 471); //f4 13.0kV
		//default_exp_manifest.runs.push(206, 355); //f3 12.0kV
		//default_exp_manifest.runs.push(5, 204); //f2 11.0kV

		default_exp_manifest.runs.push(775, 775);
		default_exp_manifest.sub_runs.push(0, 19);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		chs_to_draw.push(43); //38, 44, 42, 39
		//chs_to_draw.push(50);
		//END OF MODIFY ONLY THIS BLOCK

		channel_manifest ch_manifest;
		ch_manifest.peaks.do_find = true;
		ch_manifest.save_with_indices = false;
		ch_manifest.display.draw_peaks = true;
		ch_manifest.display.X_limits = PAIR(0, 160);
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);

		ch_manifest.baseline.baseline_by_average = true;
		ch_manifest.baseline.baseline_range = PAIR(2, 50);
		ch_manifest.baseline.do_find_curved = false;
		//The optimal parameters are the same for all channels given discretization frequency (62.5 MHz for this case)
		ch_manifest.baseline.curved_numberIterations = 90;
		ch_manifest.baseline.curved_direction = TSpectrum::kBackDecreasingWindow;
		ch_manifest.baseline.curved_filterOrder = TSpectrum::kBackOrder2;
		ch_manifest.baseline.curved_smoothing = true;
		ch_manifest.baseline.curved_smoothWindow = TSpectrum::kBackSmoothing3;
		ch_manifest.baseline.curved_compton = false;
		ch_manifest.baseline.curved_sparse = 2;

		ch_manifest.find_average = false;
		ch_manifest.find_integral = false;
		ch_manifest.find_double_integral = false;


		ch_manifest.peaks.threshold_cutoff = 0;

		ch_manifest.N_extrapolation = 1;
		//----------------------------------------------
		//Set channel specifics for all experiments (kV)
		//----------------------------------------------
		//slow PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = false;
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 32;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(0, 160);
		ch_manifest.baseline.curved_center = PAIR(50, 150);
		ch_manifest.baseline.curved_trim = PAIR(3, 3);
		ch_manifest.baseline.curved_numberIterations = 90;
		ch_manifest.baseline.curved_direction = TSpectrum::kBackDecreasingWindow;
		ch_manifest.baseline.curved_filterOrder = TSpectrum::kBackOrder2;
		ch_manifest.baseline.curved_smoothWindow = TSpectrum::kBackSmoothing7;
		ch_manifest.baseline.curved_use_unfiltered = true;

		ch_manifest.peaks.threshold = 0.042;
		ch_manifest.peaks.threshold_cutoff = 0.012;
		default_exp_manifest.channels.push(0, ch_manifest);

		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(1, ch_manifest);

		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(2, ch_manifest);

		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.006;
		default_exp_manifest.channels.push(3, ch_manifest);

		ch_manifest.peaks.threshold = 0.019;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(4, ch_manifest);

		ch_manifest.filter.n_iterations = 0;

		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold_cutoff = 0.0;

		ch_manifest.peaks.threshold = 0.0075;
		default_exp_manifest.channels.push(5, ch_manifest);

		ch_manifest.peaks.threshold = 0.0053;
		default_exp_manifest.channels.push(6, ch_manifest);

		ch_manifest.peaks.threshold = 0.0050;
		default_exp_manifest.channels.push(7, ch_manifest);

		ch_manifest.peaks.threshold = 0.0040;
		default_exp_manifest.channels.push(8, ch_manifest);

		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.invert = true;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, 0.1);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 32;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(0, 160);
		ch_manifest.baseline.curved_center = PAIR(50, 150);
		ch_manifest.baseline.curved_trim = PAIR(3, 3);
		ch_manifest.baseline.curved_numberIterations = 70;
		ch_manifest.baseline.curved_direction = TSpectrum::kBackDecreasingWindow;
		ch_manifest.baseline.curved_filterOrder = TSpectrum::kBackOrder2;
		ch_manifest.baseline.curved_smoothWindow = TSpectrum::kBackSmoothing9;
		ch_manifest.baseline.curved_use_unfiltered = true;
		ch_manifest.peaks.threshold = 0.0045;
		ch_manifest.peaks.threshold_cutoff = 0.0008;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			default_exp_manifest.channels.push(ch, ch_manifest);
		}
		default_exp_manifest.channels.info(44)->baseline.curved_use_unfiltered = false;

		for (std::size_t i = 0, i_end_ = default_exp_manifest.channels.size(); i!=i_end_; ++i) {
			int ch = default_exp_manifest.channels.index(i);
			default_exp_manifest.channels[i].display.do_draw = chs_to_draw.contains(ch);
		}

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 221020 DATA!
		//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_18.5kV_0V") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_18.0kV_338V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_17.5kV_563V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_17.0kV_917V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_16.5kV_1238V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_16.0kV_1238V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_15.5kV_1238V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_15.0kV_1238V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_14.0kV_1688V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_13.0kV_2250V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_12.0kV_2250V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_11.0kV_2025V") + "/");
		manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded TEST 221020 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init221020(analysis_manifest& manifest)
	{
		std::cout << "Initializing 221020 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
			SiPM_channels.push(32, 44);
			SiPM_channels.push(48, 59);
			experiment_manifest default_exp_manifest;
			default_exp_manifest.data_time_constant = 1.6e-2;
			default_exp_manifest.data_voltage_channels = 4095;
			default_exp_manifest.data_voltage_amplitude = 2.0;
			default_exp_manifest.data_voltage_of_zero_channel = -1.0;
			default_exp_manifest.in_folder = "../Data/221020/";
			default_exp_manifest.out_folder = "../Data/221020/results_v1/";
			default_exp_manifest.write_event_indices = false;
			default_exp_manifest.accepted_events_fname = "";
			if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
				std::cout << "Init221020:: No event selection - processing everything" << std::endl;
				default_exp_manifest.accepted_events_data.clear();
			}
			default_exp_manifest.subruns_per_file = 1000;
			default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
			default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

			//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
			default_exp_manifest.out_gnuplot_folder = "";
			default_exp_manifest.out_picture_folder = "";
			default_exp_manifest.draw_only = false; //if set to true, no data is written to output
			default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
			default_exp_manifest.sub_runs.push(0, default_exp_manifest.subruns_per_file - 1);
			default_exp_manifest.trigger_at = -32;

			area_vector chs_to_draw;	 //DRAW only these channels
			//END OF MODIFY ONLY THIS BLOCK

			channel_manifest ch_manifest;
			ch_manifest.peaks.do_find = true;
			ch_manifest.save_with_indices = false;
			ch_manifest.display.draw_peaks = true;
			ch_manifest.display.X_limits = PAIR(0, 160);
			ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);

			ch_manifest.baseline.baseline_by_average = true;
			ch_manifest.baseline.baseline_range = PAIR(2, 50);
			ch_manifest.baseline.do_find_curved = false;
			//The optimal parameters are the same for all channels given discretization frequency (62.5 MHz for this case)
			ch_manifest.baseline.curved_numberIterations = 90;
			ch_manifest.baseline.curved_direction = TSpectrum::kBackDecreasingWindow;
			ch_manifest.baseline.curved_filterOrder = TSpectrum::kBackOrder2;
			ch_manifest.baseline.curved_smoothing = true;
			ch_manifest.baseline.curved_smoothWindow = TSpectrum::kBackSmoothing3;
			ch_manifest.baseline.curved_compton = false;
			ch_manifest.baseline.curved_sparse = 2;

			ch_manifest.find_average = false;
			ch_manifest.find_integral = false;
			ch_manifest.find_double_integral = false;


			ch_manifest.peaks.threshold_cutoff = 0;

			ch_manifest.N_extrapolation = 1;
			//----------------------------------------------
			//Set channel specifics for all experiments (kV)
			//----------------------------------------------
			//slow PMTs
			ch_manifest.device = "PMT";
			ch_manifest.invert = false;
			ch_manifest.filter.order = 4;
			ch_manifest.filter.n_points = 32;
			ch_manifest.filter.n_iterations = 1;
			ch_manifest.baseline.do_find_curved = true;
			ch_manifest.baseline.curved_range = PAIR(0, 160);
			ch_manifest.baseline.curved_center = PAIR(50, 150);
			ch_manifest.baseline.curved_trim = PAIR(3, 3);
			ch_manifest.baseline.curved_numberIterations = 90;
			ch_manifest.baseline.curved_direction = TSpectrum::kBackDecreasingWindow;
			ch_manifest.baseline.curved_filterOrder = TSpectrum::kBackOrder2;
			ch_manifest.baseline.curved_smoothWindow = TSpectrum::kBackSmoothing7;
			ch_manifest.baseline.curved_use_unfiltered = true;

			ch_manifest.peaks.threshold = 0.042;
			ch_manifest.peaks.threshold_cutoff = 0.012;
			default_exp_manifest.channels.push(0, ch_manifest);

			ch_manifest.peaks.threshold = 0.022;
			ch_manifest.peaks.threshold_cutoff = 0.005;
			default_exp_manifest.channels.push(1, ch_manifest);

			ch_manifest.peaks.threshold = 0.022;
			ch_manifest.peaks.threshold_cutoff = 0.005;
			default_exp_manifest.channels.push(2, ch_manifest);

			ch_manifest.peaks.threshold = 0.022;
			ch_manifest.peaks.threshold_cutoff = 0.006;
			default_exp_manifest.channels.push(3, ch_manifest);

			ch_manifest.peaks.threshold = 0.019;
			ch_manifest.peaks.threshold_cutoff = 0.005;
			default_exp_manifest.channels.push(4, ch_manifest);

			ch_manifest.filter.n_iterations = 0;

			//fast PMTs
			ch_manifest.device = "PMT";
			ch_manifest.invert = true;
			ch_manifest.baseline.do_find_curved = false;
			ch_manifest.peaks.threshold_cutoff = 0.0;

			ch_manifest.peaks.threshold = 0.0075;
			default_exp_manifest.channels.push(5, ch_manifest);

			ch_manifest.peaks.threshold = 0.0053;
			default_exp_manifest.channels.push(6, ch_manifest);

			ch_manifest.peaks.threshold = 0.0050;
			default_exp_manifest.channels.push(7, ch_manifest);

			ch_manifest.peaks.threshold = 0.0040;
			default_exp_manifest.channels.push(8, ch_manifest);

			//SiPMs
			ch_manifest.device = "SiPM";
			ch_manifest.invert = true;
			ch_manifest.display.Y_limits = PAIR(-DBL_MAX, 0.1);
			ch_manifest.filter.order = 4;
			ch_manifest.filter.n_points = 32;
			ch_manifest.filter.n_iterations = 1;
			ch_manifest.baseline.do_find_curved = true;
			ch_manifest.baseline.curved_range = PAIR(0, 160);
			ch_manifest.baseline.curved_center = PAIR(50, 150);
			ch_manifest.baseline.curved_trim = PAIR(3, 3);
			ch_manifest.baseline.curved_numberIterations = 70;
			ch_manifest.baseline.curved_direction = TSpectrum::kBackDecreasingWindow;
			ch_manifest.baseline.curved_filterOrder = TSpectrum::kBackOrder2;
			ch_manifest.baseline.curved_smoothWindow = TSpectrum::kBackSmoothing9;
			ch_manifest.baseline.curved_use_unfiltered = true;
			ch_manifest.peaks.threshold = 0.0045;
			ch_manifest.peaks.threshold_cutoff = 0.0008;
			for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
				default_exp_manifest.channels.push(ch, ch_manifest);
			}
			default_exp_manifest.channels.info(44)->baseline.curved_use_unfiltered = false;

			for (std::size_t i = 0, i_end_ = default_exp_manifest.channels.size(); i!=i_end_; ++i) {
				int ch = default_exp_manifest.channels.index(i);
				default_exp_manifest.channels[i].display.do_draw = chs_to_draw.contains(ch);
			}

			//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 221020 DATA!
			//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
			experiment_manifest new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_18.5kV_0V") + "/");
			//this is a place to tweak individual channel for specific fields
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_18.0kV_338V") + "/");
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_17.5kV_563V") + "/");
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_17.0kV_917V") + "/");
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_16.5kV_1238V") + "/");
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_16.0kV_1238V") + "/");
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_15.5kV_1238V") + "/");
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_15.0kV_1238V") + "/");
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_14.0kV_1688V") + "/");
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_13.0kV_2250V") + "/");
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_12.0kV_2250V") + "/");
			manifest.manifests.push_back(new_manifest);

			new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = "221020_S2_LAr_Pu_WLS_11.0kV_2025V") + "/");
			manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 221020 manifests" << std::endl;
		return true;
	}
};

