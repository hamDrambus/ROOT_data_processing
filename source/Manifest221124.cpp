#include "GlobalParameters.h"

namespace ParameterPile
{
	//Use this function to play with/display different configurations
	//The correct parameters used in the analysis are in Init190404(analysis_manifest& manifest)
	bool Init221124_tests(analysis_manifest& manifest)
	{
		std::cout << "Initializing TEST 221124 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 43);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../hdda/Data/221124/";
		default_exp_manifest.out_folder = "../hdda/Data/221124/results_vt/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init221124_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../hdda/Data/221124/results_vt/gnuplot/";
		default_exp_manifest.out_picture_folder = "";
		default_exp_manifest.draw_only = true; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs): single phase, NBrS S2 in thin GEM1 from X-ray with 14 mm collimator and filter plate #3
		//default_exp_manifest.runs.push(357, 358); //f28 785V _2, after discharge at 831V
		//default_exp_manifest.runs.push(355, 355); //f26 831V
		//default_exp_manifest.runs.push(346, 353); //f25 785V
		//default_exp_manifest.runs.push(335, 344); //f24 739V
		//default_exp_manifest.runs.push(324, 333); //f23 693V
		//default_exp_manifest.runs.push(313, 322); //f22 646V
		//default_exp_manifest.runs.push(302, 311); //f21 600V
		//default_exp_manifest.runs.push(291, 300); //f20 554V
		//default_exp_manifest.runs.push(275, 289); //f19 508V
		//default_exp_manifest.runs.push(259, 273); //f18 462V
		//default_exp_manifest.runs.push(243, 257); //f17 416V
		//default_exp_manifest.runs.push(227, 241); //f16 369V
		//default_exp_manifest.runs.push(211, 225); //f15 323V
		//default_exp_manifest.runs.push(195, 209); //f14 277V
		//default_exp_manifest.runs.push(174, 193); //f13 231kV
		//default_exp_manifest.runs.push(153, 172); //f12 185V
		//default_exp_manifest.runs.push(132, 151); //f11 139V
		//default_exp_manifest.runs.push(111, 130); //f10 92kV
		//default_exp_manifest.runs.push(88, 109); //f9 46V
		//default_exp_manifest.runs.push(52, 86); //f8 0V

		default_exp_manifest.runs.push(346, 346);
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
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.baseline.curved_range = PAIR(0, 160);
		ch_manifest.baseline.curved_center = PAIR(50, 150);
		ch_manifest.baseline.curved_trim = PAIR(3, 3);
		ch_manifest.baseline.curved_numberIterations = 90;
		ch_manifest.baseline.curved_direction = TSpectrum::kBackDecreasingWindow;
		ch_manifest.baseline.curved_filterOrder = TSpectrum::kBackOrder2;
		ch_manifest.baseline.curved_smoothWindow = TSpectrum::kBackSmoothing7;
		ch_manifest.baseline.curved_use_unfiltered = true;

		ch_manifest.peaks.threshold = 0.038;
		ch_manifest.peaks.threshold_cutoff = 0.012;
		ch_manifest.baseline.do_find_curved = true;
		default_exp_manifest.channels.push(0, ch_manifest);

		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(1, ch_manifest);

		ch_manifest.peaks.threshold = 0.017;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(2, ch_manifest);

		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.006;
		default_exp_manifest.channels.push(3, ch_manifest);

		ch_manifest.peaks.threshold = 0.019;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(4, ch_manifest);
		ch_manifest.baseline.do_find_curved = false;

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
		ch_manifest.baseline.do_find_curved = false;
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

		for (std::size_t i = 0, i_end_ = default_exp_manifest.channels.size(); i!=i_end_; ++i) {
			int ch = default_exp_manifest.channels.index(i);
			default_exp_manifest.channels[i].display.do_draw = chs_to_draw.contains(ch);
		}

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 221124 DATA!
		std::vector<std::string> folder_list = {
			"221124_X-ray_S2_LAr_20kV_785V_850V_46V_14mm_coll_filt3_2",
			"221124_X-ray_S2_LAr_20kV_831V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_785V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_739V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_693V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_646V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_600V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_554V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_508V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_462V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_416V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_369V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_323V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_277V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_231V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_185V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_139V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_92V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_46V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_0V_850V_46V_14mm_coll_filt3",
		};
		for (auto const &f : folder_list) {
			experiment_manifest new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = f) + "/");
			//this is a place to tweak individual channel for specific fields
			//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
			manifest.manifests.push_back(new_manifest);
		}
#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded TEST 221124 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init221124(analysis_manifest& manifest)
	{
		std::cout << "Initializing 221124 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 43);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../hdda/Data/221124/";
		default_exp_manifest.out_folder = "../hdda/Data/221124/results_v1/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init221124:: No event selection - processing everything" << std::endl;
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
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.baseline.curved_range = PAIR(0, 160);
		ch_manifest.baseline.curved_center = PAIR(50, 150);
		ch_manifest.baseline.curved_trim = PAIR(3, 3);
		ch_manifest.baseline.curved_numberIterations = 90;
		ch_manifest.baseline.curved_direction = TSpectrum::kBackDecreasingWindow;
		ch_manifest.baseline.curved_filterOrder = TSpectrum::kBackOrder2;
		ch_manifest.baseline.curved_smoothWindow = TSpectrum::kBackSmoothing7;
		ch_manifest.baseline.curved_use_unfiltered = true;

		ch_manifest.peaks.threshold = 0.038;
		ch_manifest.peaks.threshold_cutoff = 0.012;
		ch_manifest.baseline.do_find_curved = true;
		default_exp_manifest.channels.push(0, ch_manifest);

		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(1, ch_manifest);

		ch_manifest.peaks.threshold = 0.017;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(2, ch_manifest);

		ch_manifest.peaks.threshold = 0.022;
		ch_manifest.peaks.threshold_cutoff = 0.006;
		default_exp_manifest.channels.push(3, ch_manifest);

		ch_manifest.peaks.threshold = 0.019;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(4, ch_manifest);
		ch_manifest.baseline.do_find_curved = false;

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
		ch_manifest.baseline.do_find_curved = false;
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

		for (std::size_t i = 0, i_end_ = default_exp_manifest.channels.size(); i!=i_end_; ++i) {
			int ch = default_exp_manifest.channels.index(i);
			default_exp_manifest.channels[i].display.do_draw = chs_to_draw.contains(ch);
		}

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 221124 DATA!
		std::vector<std::string> folder_list = {
			"221124_X-ray_S2_LAr_20kV_785V_850V_46V_14mm_coll_filt3_2",
			"221124_X-ray_S2_LAr_20kV_831V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_785V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_739V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_693V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_646V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_600V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_554V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_508V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_462V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_416V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_369V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_323V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_277V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_231V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_185V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_139V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_92V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_46V_850V_46V_14mm_coll_filt3",
			"221124_X-ray_S2_LAr_20kV_0V_850V_46V_14mm_coll_filt3",
		};
		for (auto const &f : folder_list) {
			experiment_manifest new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = f) + "/");
			//this is a place to tweak individual channel for specific fields
			//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
			manifest.manifests.push_back(new_manifest);
		}
#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 221124 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init221124_Q(analysis_manifest& manifest)
	{
		std::cout << "Initializing charge 221124 manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../hdda/Data/221124/";
		default_exp_manifest.out_folder = "../hdda/Data/221124/results_v1/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init221124_Q:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../hdda/Data/221124/results_v1/gnuplot/";
		default_exp_manifest.out_picture_folder = "";
		default_exp_manifest.draw_only = false; //if set to true, no data is written to output
		default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		default_exp_manifest.sub_runs.push(0, default_exp_manifest.subruns_per_file - 1);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		//END OF MODIFY ONLY THIS BLOCK

		channel_manifest ch_manifest;
		ch_manifest.save_with_indices = false;
		ch_manifest.display.X_limits = PAIR(0, 160);
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);

		ch_manifest.baseline.baseline_by_average = true;
		ch_manifest.baseline.baseline_range = PAIR(0, 29.0);

		ch_manifest.find_average = false;
		ch_manifest.find_integral = false;
		ch_manifest.find_double_integral = false;

		ch_manifest.filter.n_iterations = 0;
		ch_manifest.peaks.threshold_cutoff = 0;
		ch_manifest.N_extrapolation = 1;
		//----------------------------------------------
		//Set channel specifics for all experiments (kV)
		//----------------------------------------------
		//Charge signal
		ch_manifest.device = "Q";
		ch_manifest.invert = false;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.display.do_draw = chs_to_draw.contains(10);
		ch_manifest.peaks.do_find = false;
		ch_manifest.find_average = true;
		default_exp_manifest.channels.push(10, ch_manifest);


		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 221124 DATA!
		std::vector<std::string> folder_list = {
			"221124_X-ray_Q_20kV_850V_46V_14mm_coll_filt3",
			"221124_X-ray_Q_18kV_850V_46V_14mm_coll_filt3",
			"221124_X-ray_Q_16kV_850V_46V_14mm_coll_filt3",
			"221124_X-ray_Q_14kV_850V_46V_14mm_coll_filt3",
			"221124_X-ray_Q_12kV_850V_46V_14mm_coll_filt3",
			"221124_X-ray_Q_10kV_850V_46V_14mm_coll_filt3",
			"221124_X-ray_Q_8kV_800V_46V_14mm_coll_filt3"
		};
		for (auto const &f : folder_list) {
			experiment_manifest new_manifest = default_exp_manifest;
			new_manifest.append_folder((new_manifest.name = f) + "/");
			//this is a place to tweak individual channel for specific fields
			//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
			manifest.manifests.push_back(new_manifest);
		}
#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded charge 221124 manifests" << std::endl;
		return true;
	}
};

