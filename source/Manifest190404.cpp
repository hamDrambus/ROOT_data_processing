#include "GlobalParameters.h"

namespace ParameterPile
{
	//Use this function to play with/display different configurations
	//The correct parameters used in the analysis are in Init190404(analysis_manifest& manifest)
	bool Init190404_tests(analysis_manifest& manifest)
	{
		std::cout << "Initializing TEST 190404 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/190404/";
		default_exp_manifest.out_folder = "../Data/190404/results_vt/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "../Post_processing/190404/results_v5/Cd_46V_20kV_850V/cut_events/cut6/events.txt";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init190404_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Post_processing/190404/results_v5/Cd_46V_20kV_850V/cut_events/cut6/gnuplot/";
		default_exp_manifest.out_picture_folder = "../Post_processing/190404/results_v5/Cd_46V_20kV_850V/cut_events/cut6/";//Do not save pictures
		default_exp_manifest.draw_only = true; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(199, 227);
		//default_exp_manifest.runs.push(172, 197);
		//default_exp_manifest.runs.push(149, 170);
		//default_exp_manifest.runs.push(127, 147);
		//default_exp_manifest.runs.push(96, 125);
		//default_exp_manifest.runs.push(64, 94);
		default_exp_manifest.runs.push(33, 62);
		//default_exp_manifest.runs.push(1, 31);
		//default_exp_manifest.runs.push(0, 9999);
		default_exp_manifest.sub_runs.push(0, 9999);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		//chs_to_draw.push(44);
		chs_to_draw.push(100, 101);
		//chs_to_draw.push(38);
		//chs_to_draw.push(39);
		//chs_to_draw.push(42);
		//END OF MODIFY ONLY THIS BLOCK

		channel_manifest ch_manifest;
		ch_manifest.peaks.do_find = true;
		ch_manifest.save_with_indices = false;
		ch_manifest.display.draw_peaks = true;
		ch_manifest.display.X_limits = PAIR(0, 160);
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);

		ch_manifest.baseline.baseline_by_average = true;
		ch_manifest.baseline.baseline_range = PAIR(0, 18.5);
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
		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(8);
		ch_manifest.peaks.threshold = 0.0031;
		default_exp_manifest.channels.push(8, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(9);
		ch_manifest.peaks.threshold = 0.0040;
		default_exp_manifest.channels.push(9, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(10);
		ch_manifest.peaks.threshold = 0.0031;
		default_exp_manifest.channels.push(10, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(11);
		ch_manifest.peaks.threshold = 0.0060;
		default_exp_manifest.channels.push(11, ch_manifest);

		//slow PMTs and trigger
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 110);
		ch_manifest.baseline.curved_center = PAIR(20, 89);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.display.do_draw = chs_to_draw.contains(12);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(12, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(13);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(13, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(14);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(14, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(15);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(15, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(16);
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.080;
		default_exp_manifest.channels.push(16, ch_manifest);

		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.invert = true;
		ch_manifest.display.Y_limits = PAIR(-0.04, 0.0);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 65);
		ch_manifest.baseline.curved_center = PAIR(20, 50);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.0080;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}

		//sum of slow PMTs (for display only)
		ch_manifest.invert = false;
		ch_manifest.peaks.do_find = false;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.device = "Virtual";
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.baseline.curved_range = PAIR(10, 110);
		ch_manifest.baseline.curved_center = PAIR(20, 89);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.020; //Still required for baseline restoration
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(12, 15);
		ch_manifest.display.do_draw = chs_to_draw.contains(100);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(100, ch_manifest);

		//sum of fast PMTs (for display only)
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.0080;
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(8, 11);
		ch_manifest.display.do_draw = chs_to_draw.contains(101);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(101, ch_manifest);

		//sum of SiPM (for display only)
		ch_manifest.device = "Virtual";
		ch_manifest.invert = true;
		ch_manifest.peaks.do_find = false;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.baseline.curved_range = PAIR(10, 65);
		ch_manifest.baseline.curved_center = PAIR(20, 50);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.0160; //Still required for baseline restoration
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(32, 36);
		ch_manifest.display.do_draw = chs_to_draw.contains(102);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(102, ch_manifest);

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 190404 DATA!
		//Only SiPM thresholds are different for 46 and 48 V

		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_20kV_850V_46V_th250mV") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_20kV_850V_46V_th250mV_0") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_18kV_850V_46V_th230mV") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_16kV_850V_46V_th210mV") + "/");
		manifest.manifests.push_back(new_manifest);


		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			default_exp_manifest.channels.info(ch)->baseline.do_find_curved = false;
		}

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_14kV_850V_46V_th200mV") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_12kV_850V_46V_th160mV") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_10kV_850V_46V_th150mV") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_8kV_850V_46V_th140mV") + "/");
		manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded TEST 190404 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init190404(analysis_manifest& manifest)
	{
		std::cout << "Initializing 190404 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/190404/";
		default_exp_manifest.out_folder = "../Data/190404/results_v2/";
		default_exp_manifest.write_event_indices = false;
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init190404:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Data/190404/results_v2/temp_gnuplot/";
		default_exp_manifest.out_picture_folder = "";//Do not save pictures
		default_exp_manifest.draw_only = false; //if set to true, no data is written to output
		default_exp_manifest.accepted_events_fname = "";
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		default_exp_manifest.runs.push(199, 227);
		default_exp_manifest.runs.push(172, 197);
		default_exp_manifest.runs.push(149, 170);
		default_exp_manifest.runs.push(127, 147);
		default_exp_manifest.runs.push(96, 125);
		default_exp_manifest.runs.push(64, 94);
		default_exp_manifest.runs.push(33, 62);
		default_exp_manifest.runs.push(1, 31);
		default_exp_manifest.sub_runs.push(0, default_exp_manifest.subruns_per_file - 1);
		default_exp_manifest.trigger_at = 32;

		area_vector chs_to_draw;	 //DRAW only these channels
		//END OF MODIFY ONLY THIS BLOCK

		channel_manifest ch_manifest;
		ch_manifest.peaks.do_find = true;
		ch_manifest.save_with_indices = false;
		ch_manifest.display.draw_peaks = true;
		ch_manifest.display.X_limits = PAIR(0, 160);
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);

		ch_manifest.baseline.baseline_by_average = true;
		ch_manifest.baseline.baseline_range = PAIR(0, 18.5);
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
		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(8);
		ch_manifest.peaks.threshold = 0.0031;
		default_exp_manifest.channels.push(8, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(9);
		ch_manifest.peaks.threshold = 0.0040;
		default_exp_manifest.channels.push(9, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(10);
		ch_manifest.peaks.threshold = 0.0031;
		default_exp_manifest.channels.push(10, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(11);
		ch_manifest.peaks.threshold = 0.0060;
		default_exp_manifest.channels.push(11, ch_manifest);

		//slow PMTs and trigger
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 110);
		ch_manifest.baseline.curved_center = PAIR(20, 89);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.display.do_draw = chs_to_draw.contains(12);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(12, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(13);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(13, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(14);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(14, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(15);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(15, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(16);
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.080;
		default_exp_manifest.channels.push(16, ch_manifest);

		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.invert = true;
		ch_manifest.display.Y_limits = PAIR(-0.04, 0.0);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 65);
		ch_manifest.baseline.curved_center = PAIR(20, 50);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.0080;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}

		//sum of slow PMTs (for display only)
		ch_manifest.invert = false;
		ch_manifest.peaks.do_find = false;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.device = "Virtual";
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, -DBL_MAX);
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.baseline.curved_range = PAIR(10, 110);
		ch_manifest.baseline.curved_center = PAIR(20, 89);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.020; //Still required for baseline restoration
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(12, 15);
		ch_manifest.display.do_draw = chs_to_draw.contains(100);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(100, ch_manifest);

		//sum of fast PMTs (for display only)
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.0080;
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(8, 11);
		ch_manifest.display.do_draw = chs_to_draw.contains(101);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(101, ch_manifest);


		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 190404 DATA!
		//Only SiPM thresholds are different for 46 and 48 V

		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_20kV_850V_46V_th250mV") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_20kV_850V_46V_th250mV_0") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_18kV_850V_46V_th230mV") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_16kV_850V_46V_th210mV") + "/");
		manifest.manifests.push_back(new_manifest);


		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			default_exp_manifest.channels.info(ch)->baseline.do_find_curved = false;
		}

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_14kV_850V_46V_th200mV") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_12kV_850V_46V_th160mV") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_10kV_850V_46V_th150mV") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_8kV_850V_46V_th140mV") + "/");
		manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 190404 manifests" << std::endl;
		return true;
	}
};

