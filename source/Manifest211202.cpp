#include "GlobalParameters.h"

namespace ParameterPile
{
	//Use this function to play with/display different configurations
	//The correct parameters used in the analysis are in Init190404(analysis_manifest& manifest)
	bool Init211202_tests(analysis_manifest& manifest)
	{
		std::cout << "Initializing TEST 211202 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 42);
		//SiPM_channels.push(44, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/211202/";
		default_exp_manifest.out_folder = "../Data/211202/results_vt/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init211202_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Data/211202/results_vt/gnuplot/";
		default_exp_manifest.out_picture_folder = "";
		default_exp_manifest.draw_only = true; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(683, 732); //f19 7.8kV 9.6mm
		//default_exp_manifest.runs.push(632, 681); //f18 8.6kV 9.6mm
		//default_exp_manifest.runs.push(586, 630); //f17 9.5kV 9.6mm
		//default_exp_manifest.runs.push(538, 584); //f16 10.4kV 9.6mm
		//default_exp_manifest.runs.push(496, 536); //f15 11.2kV 9.6mm
		//default_exp_manifest.runs.push(457, 494); //f14 12.1kV 9.6mm
		//default_exp_manifest.runs.push(426, 455); //f13 13.8kV 9.6mm
		//default_exp_manifest.runs.push(394, 424); //f12 15.6kV 9.6mm
		//default_exp_manifest.runs.push(359, 392); //f11 17.3kV 9.6mm

		//default_exp_manifest.runs.push(328, 357); //f10 16kV 15.9mm
		//default_exp_manifest.runs.push(297, 326); //f9 15kV 15.9mm
		//default_exp_manifest.runs.push(293, 295); //f8 15kV 18.0mm
		//default_exp_manifest.runs.push(257, 291); //f7 14kV 18.0mm
		//default_exp_manifest.runs.push(216, 255); //f6 13kV 18.0mm
		//default_exp_manifest.runs.push(175, 214); //f5 12kV 18.0mm
		//default_exp_manifest.runs.push(134, 173); //f4 11kV 18.0mm
		//default_exp_manifest.runs.push(93, 132); //f3 10kV 18.0mm
		//default_exp_manifest.runs.push(52, 91); //f2 9kV 18.0mm
		//default_exp_manifest.runs.push(1, 50); //f1 8kV 18.0mm

		default_exp_manifest.runs.push(328, 328);
		default_exp_manifest.sub_runs.push(0, 19);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		chs_to_draw.push(39); //38, 44, 42, 39
		//chs_to_draw.push(50);
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
		ch_manifest.baseline.curved_range = PAIR(5, 100);
		ch_manifest.baseline.curved_center = PAIR(20, 75);
		ch_manifest.baseline.curved_trim = PAIR(2, 2);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.044;
		ch_manifest.peaks.threshold_cutoff = 0.010;
		default_exp_manifest.channels.push(0, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.011;
		ch_manifest.peaks.threshold_cutoff = 0.001;
		ch_manifest.filter.n_iterations = 1;
		default_exp_manifest.channels.push(1, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(2);
		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.002;
		default_exp_manifest.channels.push(2, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(3);
		ch_manifest.peaks.threshold = 0.013;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		ch_manifest.filter.n_iterations = 1;
		default_exp_manifest.channels.push(3, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(4);
		ch_manifest.peaks.threshold = 0.016;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		default_exp_manifest.channels.push(4, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold_cutoff = 0.0;
		ch_manifest.display.do_draw = chs_to_draw.contains(5);
		ch_manifest.peaks.threshold = 0.0046;
		default_exp_manifest.channels.push(5, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(6);
		ch_manifest.peaks.threshold = 0.0035;
		default_exp_manifest.channels.push(6, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(7);
		ch_manifest.peaks.threshold = 0.0045;
		default_exp_manifest.channels.push(7, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(8);
		ch_manifest.peaks.threshold = 0.0030;
		default_exp_manifest.channels.push(8, ch_manifest);

		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.invert = true;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, 0.1);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(5, 100);
		ch_manifest.baseline.curved_center = PAIR(20, 75);
		ch_manifest.baseline.curved_trim = PAIR(2, 2);
		ch_manifest.baseline.curved_numberIterations = 50;
		ch_manifest.peaks.threshold = 0.0090;
		ch_manifest.peaks.threshold_cutoff = 0.0018;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 211202 DATA!
		//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_16kV_850V_46V_12dB_16mm") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_15kV_850V_46V_12dB_16mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_15kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_14kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_13kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_12kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_11kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_10kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_9kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_8kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_17.3kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_15.6kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_13.8kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_12.1kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_11.2kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_10.4kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_9.5kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_8.6kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_7.8kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);
#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded TEST 211202 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init211202(analysis_manifest& manifest)
	{
		std::cout << "Initializing 211202 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 42);
		//SiPM_channels.push(44, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/211202/";
		default_exp_manifest.out_folder = "../Data/211202/results_v1/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init211202:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "";
		default_exp_manifest.out_picture_folder = "";
		default_exp_manifest.draw_only = false; //if set to true, no data is written to output

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
		ch_manifest.baseline.curved_range = PAIR(5, 100);
		ch_manifest.baseline.curved_center = PAIR(20, 75);
		ch_manifest.baseline.curved_trim = PAIR(2, 2);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.044;
		ch_manifest.peaks.threshold_cutoff = 0.010;
		default_exp_manifest.channels.push(0, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.011;
		ch_manifest.peaks.threshold_cutoff = 0.001;
		ch_manifest.filter.n_iterations = 1;
		default_exp_manifest.channels.push(1, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(2);
		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.002;
		default_exp_manifest.channels.push(2, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(3);
		ch_manifest.peaks.threshold = 0.013;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		ch_manifest.filter.n_iterations = 1;
		default_exp_manifest.channels.push(3, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(4);
		ch_manifest.peaks.threshold = 0.016;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		default_exp_manifest.channels.push(4, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold_cutoff = 0.0;
		ch_manifest.display.do_draw = chs_to_draw.contains(5);
		ch_manifest.peaks.threshold = 0.0046;
		default_exp_manifest.channels.push(5, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(6);
		ch_manifest.peaks.threshold = 0.0035;
		default_exp_manifest.channels.push(6, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(7);
		ch_manifest.peaks.threshold = 0.0045;
		default_exp_manifest.channels.push(7, ch_manifest);
		ch_manifest.display.do_draw = chs_to_draw.contains(8);
		ch_manifest.peaks.threshold = 0.0030;
		default_exp_manifest.channels.push(8, ch_manifest);

		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.invert = true;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, 0.1);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(5, 100);
		ch_manifest.baseline.curved_center = PAIR(20, 75);
		ch_manifest.baseline.curved_trim = PAIR(2, 2);
		ch_manifest.baseline.curved_numberIterations = 50;
		ch_manifest.peaks.threshold = 0.0090;
		ch_manifest.peaks.threshold_cutoff = 0.0018;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 211202 DATA!
		//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_16kV_850V_46V_12dB_16mm") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_15kV_850V_46V_12dB_16mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_15kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_14kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_13kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_12kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_11kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_10kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_9kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_8kV_850V_46V_12dB_18mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_17.3kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_15.6kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_13.8kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_12.1kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_11.2kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_10.4kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_9.5kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_8.6kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211202_Pu_7.8kV_850V_46V_12dB_10mm") + "/");
		manifest.manifests.push_back(new_manifest);
#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 211202 manifests" << std::endl;
		return true;
	}
};

