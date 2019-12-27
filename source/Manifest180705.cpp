#include "GlobalParameters.h"

namespace ParameterPile
{
	//Use this function to play with/display different configurations
	//The correct parameters used in the analysis are in Init190404(analysis_manifest& manifest)
	bool Init180705_tests(analysis_manifest& manifest)
	{
		std::cout << "Initializing TEST 180705 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/180705/";
		default_exp_manifest.out_folder = "../Data/180705/results_vt/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Post_processing/180705/results_v4/Cd_48V_8kV_800V/rejected_events/gnuplot/";
		default_exp_manifest.out_picture_folder = "../Post_processing/180705/results_v4/Cd_48V_8kV_800V/rejected_events/";//Do not save pictures
		default_exp_manifest.draw_only = true; //if set to true, no data is written to output
		default_exp_manifest.accepted_events_fname = "../Post_processing/180705/results_v4/Cd_48V_8kV_800V/08_rejected_events.txt";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init180705_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(281, 303); //8kV
		//default_exp_manifest.runs.push(270, 279); //9KV
		//default_exp_manifest.runs.push(259, 268); //10kV
		//default_exp_manifest.runs.push(248, 257); //11kV
		//default_exp_manifest.runs.push(305, 314); //12kV
		//default_exp_manifest.runs.push(316, 325); //13kV
		//default_exp_manifest.runs.push(327, 336); //14kV
		//default_exp_manifest.runs.push(338, 347); //16kV
		//default_exp_manifest.runs.push(349, 358); //18kV
		//default_exp_manifest.runs.push(360, 369); //20kV
		default_exp_manifest.runs.push(281, 303);
		default_exp_manifest.sub_runs.push(0, 999);// default_exp_manifest.subruns_per_file - 1);
		default_exp_manifest.trigger_at = 32;

		area_vector chs_to_draw;	 //DRAW only these channels
		chs_to_draw.push(0, 0);
		//chs_to_draw.push(12);
		//chs_to_draw.push(42);
		//END OF MODIFY ONLY THIS BLOCK

		channel_manifest ch_manifest;
		ch_manifest.peaks.do_find = true;
		ch_manifest.save_with_indices = false;
		ch_manifest.display.draw_peaks = true;
		ch_manifest.display.X_limits = PAIR(-DBL_MAX, DBL_MAX);//PAIR(25, 90);
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);

		ch_manifest.baseline.baseline_by_average = true;
		ch_manifest.baseline.baseline_range = PAIR(0, 22);
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
		ch_manifest.integral_range = PAIR(20, 40);
		ch_manifest.find_double_integral = false;
		ch_manifest.double_integral_range = PAIR(20, 40);

		ch_manifest.filter.n_iterations = 0;
		ch_manifest.filter.n_points = 2;
		ch_manifest.filter.order = 1;

		ch_manifest.peaks.threshold_cutoff = 0;

		ch_manifest.N_extrapolation = 1;
		//----------------------------------------------
		//Set channel specifics for all experiments (kV)
		//----------------------------------------------
		//slow PMTs and trigger
		ch_manifest.invert = false;
		ch_manifest.device = "PMT";
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(17, 110);
		ch_manifest.baseline.curved_center = PAIR(27, 89);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.038;
		ch_manifest.peaks.threshold_cutoff = 0.0044;
		default_exp_manifest.channels.push(0, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.065;
		ch_manifest.peaks.threshold_cutoff = 0.008;
		default_exp_manifest.channels.push(1, ch_manifest);

		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(17, 85);
		ch_manifest.baseline.curved_center = PAIR(27, 85);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.display.Y_limits = PAIR(-0.55, -0.4);
		ch_manifest.peaks.threshold_cutoff = 0.0;
		ch_manifest.display.do_draw = chs_to_draw.contains(8);
		ch_manifest.peaks.threshold = 0.00324;
		default_exp_manifest.channels.push(8, ch_manifest);

		ch_manifest.display.do_draw = chs_to_draw.contains(9);
		ch_manifest.peaks.threshold = 0.00324;
		default_exp_manifest.channels.push(9, ch_manifest);

		ch_manifest.display.do_draw = chs_to_draw.contains(10);
		ch_manifest.peaks.threshold = 0.00324;
		default_exp_manifest.channels.push(10, ch_manifest);

		ch_manifest.display.do_draw = chs_to_draw.contains(11);
		ch_manifest.peaks.threshold = 0.0043;
		default_exp_manifest.channels.push(11, ch_manifest);

		ch_manifest.display.do_draw = chs_to_draw.contains(12);
		ch_manifest.peaks.threshold = 0.0065;
		default_exp_manifest.channels.push(12, ch_manifest);
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);

		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.display.Y_limits = PAIR(-0.04, 0.0);
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(14, 65);
		ch_manifest.baseline.curved_center = PAIR(27, 50);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.peaks.threshold = 0.0140;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			if (ch == 44 || ch == 56 || ch == 57 || ch == 41)
				ch_manifest.filter.n_iterations = 0; //1
			else
				ch_manifest.filter.n_iterations = 0;
			default_exp_manifest.channels.push(ch, ch_manifest);
		}


		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_20kV_800V_12bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_18kV_800V_12bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_16kV_800V_12bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_14kV_800V_12bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_13kV_800V_12bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		default_exp_manifest.channels.info(0)->peaks.threshold = 0.052;
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_12kV_800V_6bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_11kV_800V_6bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_10kV_800V_6bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			default_exp_manifest.channels.info(ch)->baseline.do_find_curved = false;
		}
		default_exp_manifest.channels.info(0)->peaks.threshold = 0.085;

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_9kV_800V_0bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_8kV_800V_0bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded TEST 180705 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init180705(analysis_manifest& manifest)
	{
		std::cout << "Initializing 180705 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/180705/";
		default_exp_manifest.out_folder = "../Data/180705/results_v2/";
		default_exp_manifest.write_event_indices = false;
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data))
			default_exp_manifest.accepted_events_data.clear();
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Data/180705/results_v2/temp_gnuplot/";
		default_exp_manifest.out_picture_folder = "";//Do not save pictures
		default_exp_manifest.draw_only = false; //if set to true, no data is written to output
		default_exp_manifest.accepted_events_fname = "";
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		default_exp_manifest.runs.push(281, 303);
		default_exp_manifest.runs.push(270, 279);
		default_exp_manifest.runs.push(259, 268);
		default_exp_manifest.runs.push(248, 257);
		default_exp_manifest.runs.push(305, 314);
		default_exp_manifest.runs.push(316, 325);
		default_exp_manifest.runs.push(327, 336);
		default_exp_manifest.runs.push(338, 347);
		default_exp_manifest.runs.push(349, 358);
		default_exp_manifest.runs.push(360, 369);
		default_exp_manifest.sub_runs.push(0, default_exp_manifest.subruns_per_file - 1);
		default_exp_manifest.trigger_at = 32;

		area_vector chs_to_draw;	 //DRAW only these channels
		//chs_to_draw.push(38);
		//chs_to_draw.push(44);
		//chs_to_draw.push(42);
		//END OF MODIFY ONLY THIS BLOCK

		channel_manifest ch_manifest;
		ch_manifest.peaks.do_find = true;
		ch_manifest.save_with_indices = false;
		ch_manifest.display.draw_peaks = true;
		ch_manifest.display.X_limits = PAIR(0, 160);
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);

		ch_manifest.baseline.baseline_by_average = true;
		ch_manifest.baseline.baseline_range = PAIR(0, 22);
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
		ch_manifest.integral_range = PAIR(20, 40);
		ch_manifest.find_double_integral = false;
		ch_manifest.double_integral_range = PAIR(20, 40);

		ch_manifest.filter.n_iterations = 0;
		ch_manifest.filter.n_points = 2;
		ch_manifest.filter.order = 1;

		ch_manifest.peaks.threshold_cutoff = 0;

		ch_manifest.N_extrapolation = 1;
		//----------------------------------------------
		//Set channel specifics for all experiments (kV)
		//----------------------------------------------
		//slow PMTs and trigger
		ch_manifest.device = "PMT";
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(17, 110);
		ch_manifest.baseline.curved_center = PAIR(27, 89);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.038;
		ch_manifest.peaks.threshold_cutoff = 0.0044;
		default_exp_manifest.channels.push(0, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.065;
		ch_manifest.peaks.threshold_cutoff = 0.008;
		default_exp_manifest.channels.push(1, ch_manifest);

		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(17, 85);
		ch_manifest.baseline.curved_center = PAIR(27, 85);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold_cutoff = 0.0;
		ch_manifest.display.do_draw = chs_to_draw.contains(8);
		ch_manifest.peaks.threshold = 0.00324;
		default_exp_manifest.channels.push(8, ch_manifest);

		ch_manifest.display.do_draw = chs_to_draw.contains(9);
		ch_manifest.peaks.threshold = 0.00324;
		default_exp_manifest.channels.push(9, ch_manifest);

		ch_manifest.display.do_draw = chs_to_draw.contains(10);
		ch_manifest.peaks.threshold = 0.00324;
		default_exp_manifest.channels.push(10, ch_manifest);

		ch_manifest.display.do_draw = chs_to_draw.contains(11);
		ch_manifest.peaks.threshold = 0.0043;
		default_exp_manifest.channels.push(11, ch_manifest);

		ch_manifest.display.do_draw = chs_to_draw.contains(12);
		ch_manifest.peaks.threshold = 0.0065;
		default_exp_manifest.channels.push(12, ch_manifest);

		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.display.Y_limits = PAIR(-0.04, 0.0);
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(14, 65);
		ch_manifest.baseline.curved_center = PAIR(27, 50);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.peaks.threshold = 0.0140;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			if (ch == 44 || ch == 56 || ch == 57 || ch == 41)
				ch_manifest.filter.n_iterations = 1;
			else
				ch_manifest.filter.n_iterations = 0;
			default_exp_manifest.channels.push(ch, ch_manifest);
		}


		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_20kV_800V_12bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_18kV_800V_12bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_16kV_800V_12bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_14kV_800V_12bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_13kV_800V_12bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		default_exp_manifest.channels.info(0)->peaks.threshold = 0.052;
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_12kV_800V_6bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_11kV_800V_6bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_10kV_800V_6bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			default_exp_manifest.channels.info(ch)->baseline.do_find_curved = false;
		}
		default_exp_manifest.channels.info(0)->peaks.threshold = 0.085;

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_9kV_800V_0bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "180705_Cd_8kV_800V_0bB_48V") + "/");
		manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 180705 manifests" << std::endl;
		return true;
	}

};
