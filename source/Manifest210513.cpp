#include "GlobalParameters.h"

namespace ParameterPile
{
	//Use this function to play with/display different configurations
	//The correct parameters used in the analysis are in Init190404(analysis_manifest& manifest)
	bool Init210513_tests(analysis_manifest& manifest)
	{
		std::cout << "Initializing TEST 210513 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 42);
		SiPM_channels.push(44, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/210513/";
		default_exp_manifest.out_folder = "../Data/210513/results_vt/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init210513_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Data/210513/results_vt/gnuplot/";
		default_exp_manifest.out_picture_folder = "";
		default_exp_manifest.draw_only = true; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(1, 31); //20kV, 800V
		//default_exp_manifest.runs.push(33, 62); //20kV, 850V
		//default_exp_manifest.runs.push(64, 93); //18kV
		//default_exp_manifest.runs.push(95, 124); //16kV
		//default_exp_manifest.runs.push(126, 155); //14kV
		//default_exp_manifest.runs.push(157, 196); //13kV
		//default_exp_manifest.runs.push(198, 238); //12kV
		//default_exp_manifest.runs.push(240, 289); //11kV
		//default_exp_manifest.runs.push(291, 340); //10kV
		//default_exp_manifest.runs.push(342, 391); //9kV
		//default_exp_manifest.runs.push(393, 442); //8kV

		default_exp_manifest.runs.push(342, 342);
		default_exp_manifest.sub_runs.push(0, 19);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		//chs_to_draw.push(100);
		//chs_to_draw.push(101);
		//chs_to_draw.push(102);
		//chs_to_draw.push(42); //38, 44, 42, 39
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

		ch_manifest.filter.order = 2;
		ch_manifest.filter.n_points = 17;
		ch_manifest.filter.n_iterations = 0;

		ch_manifest.peaks.threshold_cutoff = 0;

		ch_manifest.N_extrapolation = 1;
		//----------------------------------------------
		//Set channel specifics for all experiments (kV)
		//----------------------------------------------
		//slow PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 110);
		ch_manifest.baseline.curved_center = PAIR(20, 80);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.044;
		ch_manifest.peaks.threshold_cutoff = 0.010;
		default_exp_manifest.channels.push(0, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.0096;
		ch_manifest.peaks.threshold_cutoff = 0.001;
		default_exp_manifest.channels.push(1, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(2);
		ch_manifest.peaks.threshold = 0.013;
		ch_manifest.peaks.threshold_cutoff = 0.002;
		default_exp_manifest.channels.push(2, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(3);
		ch_manifest.peaks.threshold = 0.011;
		ch_manifest.peaks.threshold_cutoff = 0.001;
		default_exp_manifest.channels.push(3, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(4);
		ch_manifest.peaks.threshold = 0.012;
		ch_manifest.peaks.threshold_cutoff = 0.002;
		default_exp_manifest.channels.push(4, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold_cutoff = 0.0;
		ch_manifest.display.do_draw = chs_to_draw.contains(5);
		ch_manifest.peaks.threshold = 0.0040;
		default_exp_manifest.channels.push(5, ch_manifest); //Bad channel - periodic noise (issue with PMT#1?)
		ch_manifest.display.do_draw = chs_to_draw.contains(6);
		ch_manifest.peaks.threshold = 0.0035;
		default_exp_manifest.channels.push(6, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(7);
		ch_manifest.peaks.threshold = 0.0040;
		default_exp_manifest.channels.push(7, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(8);
		ch_manifest.peaks.threshold = 0.0028;
		default_exp_manifest.channels.push(8, ch_manifest);

		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.invert = true;
		ch_manifest.display.Y_limits = PAIR(-0.04, 0.0);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 85);
		ch_manifest.baseline.curved_center = PAIR(20, 70);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.0083;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}

		//sum of slow PMTs (for display only)
		ch_manifest.device = "Virtual";
		ch_manifest.peaks.do_find = false;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.05;
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(1, 4);
		ch_manifest.display.do_draw = chs_to_draw.contains(100);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(100, ch_manifest);

		//sum of fast PMTs (for display only)
		ch_manifest.device = "Virtual";
		ch_manifest.peaks.do_find = false;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.0080;
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(5, 8);
		ch_manifest.display.do_draw = chs_to_draw.contains(101);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(101, ch_manifest);

		//sum of SiPM (for display only)
		ch_manifest.device = "Virtual";
		ch_manifest.invert = true;
		ch_manifest.peaks.do_find = false;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.baseline.curved_range = PAIR(10, 85);
		ch_manifest.baseline.curved_center = PAIR(20, 70);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.0160; //Still required for baseline restoration
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(51); //central
		ch_manifest.summarize_channels.push(53); //central
		ch_manifest.summarize_channels.push(38, 39); //central
		ch_manifest.summarize_channels.push(41); //central
		//ch_manifest.summarize_channels.push(32, 44);
		//ch_manifest.summarize_channels.push(48, 59);
		ch_manifest.display.do_draw = chs_to_draw.contains(102);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(102, ch_manifest);

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 210128 DATA!
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_20kV_800V_46V_12dB") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		//Thresholds for 850V on PMTs
		//slow PMTs
		default_exp_manifest.channels.info(0)->peaks.threshold = 0.044;
		default_exp_manifest.channels.info(0)->peaks.threshold_cutoff = 0.010;
		default_exp_manifest.channels.info(1)->peaks.threshold = 0.011;
		default_exp_manifest.channels.info(1)->peaks.threshold_cutoff = 0.001;
		default_exp_manifest.channels.info(2)->peaks.threshold = 0.015;
		default_exp_manifest.channels.info(2)->peaks.threshold_cutoff = 0.002;
		default_exp_manifest.channels.info(3)->peaks.threshold = 0.012;
		default_exp_manifest.channels.info(3)->peaks.threshold_cutoff = 0.003;
		default_exp_manifest.channels.info(4)->peaks.threshold = 0.015;
		default_exp_manifest.channels.info(4)->peaks.threshold_cutoff = 0.003;
		//fast PMTs
		default_exp_manifest.channels.info(5)->peaks.threshold = 0.0046;
		default_exp_manifest.channels.info(5)->peaks.threshold_cutoff = 0.0;
		default_exp_manifest.channels.info(6)->peaks.threshold = 0.0035;
		default_exp_manifest.channels.info(6)->peaks.threshold_cutoff = 0.0;
		default_exp_manifest.channels.info(7)->peaks.threshold = 0.0045;
		default_exp_manifest.channels.info(7)->peaks.threshold_cutoff = 0.0;
		default_exp_manifest.channels.info(8)->peaks.threshold = 0.0030;
		default_exp_manifest.channels.info(8)->peaks.threshold_cutoff = 0.0;

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_20kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_18kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_16kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_14kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		default_exp_manifest.channels.info(1)->baseline.do_find_curved = false;
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_13kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_12kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_11kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_10kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_9kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_8kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded TEST 210513 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init210513(analysis_manifest& manifest)
	{
		std::cout << "Initializing 210513 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 42);
		SiPM_channels.push(44, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/210513/";
		default_exp_manifest.out_folder = "../Data/210513/results_v1/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init210513:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "";
		default_exp_manifest.out_picture_folder = "";
		default_exp_manifest.draw_only = false; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(1, 31); //20kV, 800V
		//default_exp_manifest.runs.push(33, 62); //20kV, 850V
		//default_exp_manifest.runs.push(64, 93); //18kV
		//default_exp_manifest.runs.push(95, 124); //16kV
		//default_exp_manifest.runs.push(126, 155); //14kV
		//default_exp_manifest.runs.push(157, 196); //13kV
		//default_exp_manifest.runs.push(198, 238); //12kV
		//default_exp_manifest.runs.push(240, 289); //11kV
		//default_exp_manifest.runs.push(291, 340); //10kV
		//default_exp_manifest.runs.push(342, 391); //9kV
		//default_exp_manifest.runs.push(393, 442); //8kV

		default_exp_manifest.runs.push(0, 9999);
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

		ch_manifest.filter.order = 2;
		ch_manifest.filter.n_points = 17;
		ch_manifest.filter.n_iterations = 0;

		ch_manifest.peaks.threshold_cutoff = 0;

		ch_manifest.N_extrapolation = 1;
		//----------------------------------------------
		//Set channel specifics for all experiments (kV)
		//----------------------------------------------
		//slow PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 110);
		ch_manifest.baseline.curved_center = PAIR(20, 80);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.044;
		ch_manifest.peaks.threshold_cutoff = 0.010;
		default_exp_manifest.channels.push(0, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.0096;
		ch_manifest.peaks.threshold_cutoff = 0.001;
		default_exp_manifest.channels.push(1, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(2);
		ch_manifest.peaks.threshold = 0.013;
		ch_manifest.peaks.threshold_cutoff = 0.002;
		default_exp_manifest.channels.push(2, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(3);
		ch_manifest.peaks.threshold = 0.011;
		ch_manifest.peaks.threshold_cutoff = 0.001;
		default_exp_manifest.channels.push(3, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(4);
		ch_manifest.peaks.threshold = 0.012;
		ch_manifest.peaks.threshold_cutoff = 0.002;
		default_exp_manifest.channels.push(4, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold_cutoff = 0.0;
		ch_manifest.display.do_draw = chs_to_draw.contains(5);
		ch_manifest.peaks.threshold = 0.0040;
		default_exp_manifest.channels.push(5, ch_manifest); //Bad channel - periodic noise (issue with PMT#1?)
		ch_manifest.display.do_draw = chs_to_draw.contains(6);
		ch_manifest.peaks.threshold = 0.0035;
		default_exp_manifest.channels.push(6, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(7);
		ch_manifest.peaks.threshold = 0.0040;
		default_exp_manifest.channels.push(7, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(8);
		ch_manifest.peaks.threshold = 0.0028;
		default_exp_manifest.channels.push(8, ch_manifest);

		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.invert = true;
		ch_manifest.display.Y_limits = PAIR(-0.04, 0.0);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 85);
		ch_manifest.baseline.curved_center = PAIR(20, 70);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.0083;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}

		//sum of slow PMTs (for display only)
		ch_manifest.device = "Virtual";
		ch_manifest.peaks.do_find = false;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.05;
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(1, 4);
		ch_manifest.display.do_draw = chs_to_draw.contains(100);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(100, ch_manifest);

		//sum of fast PMTs (for display only)
		ch_manifest.device = "Virtual";
		ch_manifest.peaks.do_find = false;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.0080;
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(5, 8);
		ch_manifest.display.do_draw = chs_to_draw.contains(101);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(101, ch_manifest);

		//sum of SiPM (for display only)
		ch_manifest.device = "Virtual";
		ch_manifest.invert = true;
		ch_manifest.peaks.do_find = false;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.baseline.curved_range = PAIR(10, 85);
		ch_manifest.baseline.curved_center = PAIR(20, 70);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.0160; //Still required for baseline restoration
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(51); //central
		ch_manifest.summarize_channels.push(53); //central
		ch_manifest.summarize_channels.push(38, 39); //central
		ch_manifest.summarize_channels.push(41); //central
		//ch_manifest.summarize_channels.push(32, 44);
		//ch_manifest.summarize_channels.push(48, 59);
		ch_manifest.display.do_draw = chs_to_draw.contains(102);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(102, ch_manifest);

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 210128 DATA!
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_20kV_800V_46V_12dB") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		//Thresholds for 850V on PMTs
		//slow PMTs
		default_exp_manifest.channels.info(0)->peaks.threshold = 0.044;
		default_exp_manifest.channels.info(0)->peaks.threshold_cutoff = 0.010;
		default_exp_manifest.channels.info(1)->peaks.threshold = 0.011;
		default_exp_manifest.channels.info(1)->peaks.threshold_cutoff = 0.001;
		default_exp_manifest.channels.info(2)->peaks.threshold = 0.015;
		default_exp_manifest.channels.info(2)->peaks.threshold_cutoff = 0.002;
		default_exp_manifest.channels.info(3)->peaks.threshold = 0.012;
		default_exp_manifest.channels.info(3)->peaks.threshold_cutoff = 0.003;
		default_exp_manifest.channels.info(4)->peaks.threshold = 0.015;
		default_exp_manifest.channels.info(4)->peaks.threshold_cutoff = 0.003;
		//fast PMTs
		default_exp_manifest.channels.info(5)->peaks.threshold = 0.0046;
		default_exp_manifest.channels.info(5)->peaks.threshold_cutoff = 0.0;
		default_exp_manifest.channels.info(6)->peaks.threshold = 0.0035;
		default_exp_manifest.channels.info(6)->peaks.threshold_cutoff = 0.0;
		default_exp_manifest.channels.info(7)->peaks.threshold = 0.0045;
		default_exp_manifest.channels.info(7)->peaks.threshold_cutoff = 0.0;
		default_exp_manifest.channels.info(8)->peaks.threshold = 0.0030;
		default_exp_manifest.channels.info(8)->peaks.threshold_cutoff = 0.0;

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_20kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_18kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_16kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_14kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		default_exp_manifest.channels.info(1)->baseline.do_find_curved = false;
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_13kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_12kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_11kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_10kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_9kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210513_Cd_8kV_850V_46V_12dB") + "/");
		manifest.manifests.push_back(new_manifest);
#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 210513 manifests" << std::endl;
		return true;
	}
};

