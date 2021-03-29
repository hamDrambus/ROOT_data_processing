#include "GlobalParameters.h"

namespace ParameterPile
{
	//Use this function to play with/display different configurations
	//The correct parameters used in the analysis are in Init190404(analysis_manifest& manifest)
	bool Init210218_tests(analysis_manifest& manifest)
	{
		std::cout << "Initializing TEST 210218 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/210218/";
		default_exp_manifest.out_folder = "../Data/210218/results_vt/";
		default_exp_manifest.write_event_indices = false;
		//default_exp_manifest.accepted_events_fname = "";
		std::string cut = "cut5";
		default_exp_manifest.accepted_events_fname = "../Post_processing/210218/results_v2/Pu_54V_5.5kV_850V_low_thr/sPMTs_Fr_v1/" + cut +"/events.txt";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init210218_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Post_processing/210218/results_v2/Pu_54V_5.5kV_850V_low_thr/sPMTs_Fr_v1/" + cut +"/";
		default_exp_manifest.out_picture_folder = "../Post_processing/210218/results_v2/Pu_54V_5.5kV_850V_low_thr/sPMTs_Fr_v1/" + cut +"/";
		default_exp_manifest.draw_only = true; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(1, 50); //5.5kV, 2V threshold
		//default_exp_manifest.runs.push(52, 101); //5.5kV, 250mV threshold

		default_exp_manifest.runs.push(52, 54);
		default_exp_manifest.sub_runs.push(0, default_exp_manifest.subruns_per_file - 1);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		chs_to_draw.push(100);
		chs_to_draw.push(101);
		//chs_to_draw.push(9); //38, 44, 42, 39
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
		ch_manifest.baseline.curved_range = PAIR(0, 160);
		ch_manifest.baseline.curved_center = PAIR(12, 160);
		ch_manifest.baseline.curved_trim = PAIR(2, 2);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.044;
		ch_manifest.peaks.threshold_cutoff = 0.009;
		default_exp_manifest.channels.push(0, ch_manifest); //sum of slow PMTs

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.013;
		ch_manifest.peaks.threshold_cutoff = 0.0025;
		default_exp_manifest.channels.push(1, ch_manifest); //slow PMT#1, poor channel - periodic noise - filtered (issue with PMT#1?)
		ch_manifest.filter.n_iterations = 0;

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(2);
		ch_manifest.peaks.threshold = 0.023;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		default_exp_manifest.channels.push(2, ch_manifest); //slow PMT#2

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(3);
		ch_manifest.peaks.threshold = 0.014;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		default_exp_manifest.channels.push(3, ch_manifest); //slow PMT#3

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(4);
		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(4, ch_manifest); //slow PMT#4

		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold_cutoff = 0.0;
		ch_manifest.display.do_draw = chs_to_draw.contains(5);
		ch_manifest.peaks.threshold = 0.0042;
		default_exp_manifest.channels.push(5, ch_manifest); // fast PMT#1, bad channel - periodic noise (issue with PMT#1?)
		ch_manifest.display.do_draw = chs_to_draw.contains(6);
		ch_manifest.peaks.threshold = 0.0037;
		default_exp_manifest.channels.push(6, ch_manifest); //fast PMT#2
		ch_manifest.display.do_draw = chs_to_draw.contains(7);
		ch_manifest.peaks.threshold = 0.0030;
		default_exp_manifest.channels.push(7, ch_manifest); //fast PMT#3
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(8);
		ch_manifest.peaks.threshold = 0.0030;
		default_exp_manifest.channels.push(8, ch_manifest); //fast PMT#4

		//trigger (PMTs)
		ch_manifest.device = "PMT";
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.display.do_draw = chs_to_draw.contains(9);
		ch_manifest.peaks.threshold = 0.11;
		ch_manifest.peaks.threshold_cutoff = 0.07;
		default_exp_manifest.channels.push(9, ch_manifest); //trigger (PMTs)

		//sum of slow PMTs (for display only)
		ch_manifest.device = "Virtual";
		ch_manifest.peaks.do_find = false;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.05;
		ch_manifest.peaks.threshold_cutoff = 0.0;
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

		//TRIGGER POSITIONS ARE DIFFERENT FOR TWO FOLDERS FOR 210218 DATA!
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210218_Pu_5.5kV_850V_0V_12dB_th2V") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210218_Pu_5.5kV_850V_0V_12dB_th250mV") + "/");
		manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded TEST 210218 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init210218(analysis_manifest& manifest)
	{
		std::cout << "Initializing 210218 manifests..." << std::endl;
		#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/210218/";
		default_exp_manifest.out_folder = "../Data/210218/results_v1/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init210218_tests:: No event selection - processing everything" << std::endl;
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
		//default_exp_manifest.runs.push(1, 50); //5.5kV, 2V threshold
		//default_exp_manifest.runs.push(52, 101); //5.5kV, 250mV threshold

		default_exp_manifest.runs.push(0, 9999);
		default_exp_manifest.sub_runs.push(0, default_exp_manifest.subruns_per_file - 1);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		//chs_to_draw.push(100);
		//chs_to_draw.push(101);
		//chs_to_draw.push(102);
		//chs_to_draw.push(9); //38, 44, 42, 39
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
		ch_manifest.baseline.curved_range = PAIR(0, 160);
		ch_manifest.baseline.curved_center = PAIR(12, 160);
		ch_manifest.baseline.curved_trim = PAIR(2, 2);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.044;
		ch_manifest.peaks.threshold_cutoff = 0.009;
		default_exp_manifest.channels.push(0, ch_manifest); //sum of slow PMTs

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.013;
		ch_manifest.peaks.threshold_cutoff = 0.0025;
		default_exp_manifest.channels.push(1, ch_manifest); //slow PMT#1, poor channel - periodic noise - filtered (issue with PMT#1?)
		ch_manifest.filter.n_iterations = 0;

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(2);
		ch_manifest.peaks.threshold = 0.023;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		default_exp_manifest.channels.push(2, ch_manifest); //slow PMT#2

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(3);
		ch_manifest.peaks.threshold = 0.014;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		default_exp_manifest.channels.push(3, ch_manifest); //slow PMT#3

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(4);
		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(4, ch_manifest); //slow PMT#4

		//fast PMTs
		ch_manifest.device = "PMT";
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold_cutoff = 0.0;
		ch_manifest.display.do_draw = chs_to_draw.contains(5);
		ch_manifest.peaks.threshold = 0.0042;
		default_exp_manifest.channels.push(5, ch_manifest); // fast PMT#1, bad channel - periodic noise (issue with PMT#1?)
		ch_manifest.display.do_draw = chs_to_draw.contains(6);
		ch_manifest.peaks.threshold = 0.0037;
		default_exp_manifest.channels.push(6, ch_manifest); //fast PMT#2
		ch_manifest.display.do_draw = chs_to_draw.contains(7);
		ch_manifest.peaks.threshold = 0.0030;
		default_exp_manifest.channels.push(7, ch_manifest); //fast PMT#3
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(8);
		ch_manifest.peaks.threshold = 0.0030;
		default_exp_manifest.channels.push(8, ch_manifest); //fast PMT#4

		//trigger (PMTs)
		ch_manifest.device = "PMT";
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.display.do_draw = chs_to_draw.contains(9);
		ch_manifest.peaks.threshold = 0.11;
		ch_manifest.peaks.threshold_cutoff = 0.07;
		default_exp_manifest.channels.push(9, ch_manifest); //trigger (PMTs)

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

		//TRIGGER POSITIONS ARE DIFFERENT FOR TWO FOLDERS FOR 210218 DATA!
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210218_Pu_5.5kV_850V_0V_12dB_th2V") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "210218_Pu_5.5kV_850V_0V_12dB_th250mV") + "/");
		manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 210218 manifests" << std::endl;
		return true;
	}
};

