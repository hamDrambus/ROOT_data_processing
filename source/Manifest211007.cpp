#include "GlobalParameters.h"

namespace ParameterPile
{
	//Use this function to play with/display different configurations
	//The correct parameters used in the analysis are in Init190404(analysis_manifest& manifest)
	bool Init211007_tests(analysis_manifest& manifest)
	{
		std::cout << "Initializing TEST 211007 manifests..." << std::endl;
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
		default_exp_manifest.in_folder = "../Data/211007/";
		default_exp_manifest.out_folder = "../Data/211007/results_vt/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_19.8kV_800V_2500V/events/05_events.txt";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init211007_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Post_processing/211007/results_v5/Pu_46V_19.8kV_800V_2500V/events/05_events/";
		default_exp_manifest.out_picture_folder = "../Post_processing/211007/results_v5/Pu_46V_19.8kV_800V_2500V/events/05_events/";
		default_exp_manifest.draw_only = true; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(1, 25); //f1 20.0kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(27, 51); //f2 18.8kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(53, 77); //f3 17.0kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(79, 103); //f4 15.3kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(105, 129); //f5 14.4kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(131, 160); //f6 13.6kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(162, 191); //f7 12.7kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(193, 222); //f8 11.8kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(224, 253); //f9 11.0kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(255, 279); //f10 19.9kV, 2550V on THGEM1
		//default_exp_manifest.runs.push(281, 305); //f11 19.8kV, 2500V on THGEM1
		//default_exp_manifest.runs.push(307, 331); //f12 19.9kV, 2600V on THGEM1

		default_exp_manifest.runs.push(281, 282);
		default_exp_manifest.sub_runs.push(0, 999);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		//chs_to_draw.push(100);
		//chs_to_draw.push(101);
		//chs_to_draw.push(102);
		chs_to_draw.push(38); //38, 44, 42, 39
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

		/* PMT#1 is not working after setup re-assembly (swapped signal wires)
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.011;
		ch_manifest.peaks.threshold_cutoff = 0.001;
		default_exp_manifest.channels.push(1, ch_manifest);
		*/

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
		/* PMT#1 is not working after setup re-assembly (swapped signal wires)
		ch_manifest.display.do_draw = chs_to_draw.contains(5);
		ch_manifest.peaks.threshold = 0.0046;
		default_exp_manifest.channels.push(5, ch_manifest);
		 */
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
		ch_manifest.baseline.curved_range = PAIR(10, 85);
		ch_manifest.baseline.curved_center = PAIR(20, 70);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.baseline.curved_numberIterations = 50;
		ch_manifest.peaks.threshold = 0.0083;
		ch_manifest.peaks.threshold_cutoff = 0.0018;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}

		//Charge signal
		ch_manifest.device = "Charge";
		ch_manifest.invert = false;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.baseline.baseline_range = PAIR(0, 18.5);
		ch_manifest.peaks.do_find = false;
		ch_manifest.find_average = true;
		ch_manifest.display.do_draw = false;
		default_exp_manifest.channels.push(10, ch_manifest);

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 211007 DATA!
		//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_20.0kV_800V_46V_12dB_2650V") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_18.8kV_800V_46V_12dB_2650V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_17.0kV_800V_46V_12dB_2650V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_19.8kV_800V_46V_12dB_2500V") + "/");
		manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded TEST 211007 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init211007(analysis_manifest& manifest)
	{
		std::cout << "Initializing 211007 manifests..." << std::endl;
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
		default_exp_manifest.in_folder = "../Data/211007/";
		default_exp_manifest.out_folder = "../Data/211007/results_v1/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init211007:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Data/211007/results_v1/";
		default_exp_manifest.out_picture_folder = "";
		default_exp_manifest.draw_only = false; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(1, 25); //f1 20.0kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(27, 51); //f2 18.8kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(53, 77); //f3 17.0kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(79, 103); //f4 15.3kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(105, 129); //f5 14.4kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(131, 160); //f6 13.6kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(162, 191); //f7 12.7kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(193, 222); //f8 11.8kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(224, 253); //f9 11.0kV, 2650V on THGEM1
		//default_exp_manifest.runs.push(255, 279); //f10 19.9kV, 2550V on THGEM1
		//default_exp_manifest.runs.push(281, 305); //f11 19.8kV, 2500V on THGEM1
		//default_exp_manifest.runs.push(307, 331); //f12 19.9kV, 2600V on THGEM1

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

		/* PMT#1 is not working after setup re-assembly (swapped signal wires)
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.011;
		ch_manifest.peaks.threshold_cutoff = 0.001;
		default_exp_manifest.channels.push(1, ch_manifest);
		*/

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
		/* PMT#1 is not working after setup re-assembly (swapped signal wires)
		ch_manifest.display.do_draw = chs_to_draw.contains(5);
		ch_manifest.peaks.threshold = 0.0046;
		default_exp_manifest.channels.push(5, ch_manifest);
		 */
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
		ch_manifest.baseline.curved_range = PAIR(10, 85);
		ch_manifest.baseline.curved_center = PAIR(20, 70);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.baseline.curved_numberIterations = 50;
		ch_manifest.peaks.threshold = 0.0083;
		ch_manifest.peaks.threshold_cutoff = 0.0018;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}

		//Charge signal
		ch_manifest.device = "Charge";
		ch_manifest.invert = false;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.baseline.baseline_range = PAIR(0, 18.5);
		ch_manifest.peaks.do_find = false;
		ch_manifest.find_average = true;
		ch_manifest.display.do_draw = false;
		default_exp_manifest.channels.push(10, ch_manifest);

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 211007 DATA!
		//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_20.0kV_800V_46V_12dB_2650V") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_18.8kV_800V_46V_12dB_2650V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_17.0kV_800V_46V_12dB_2650V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_15.3kV_800V_46V_12dB_2650V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_14.4kV_800V_46V_12dB_2650V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_13.6kV_800V_46V_12dB_2650V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_12.7kV_800V_46V_12dB_2650V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_11.8kV_800V_46V_12dB_2650V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_11.0kV_800V_46V_12dB_2650V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_19.9kV_800V_46V_12dB_2550V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_19.8kV_800V_46V_12dB_2500V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_19.9kV_800V_46V_12dB_2600V") + "/");
		manifest.manifests.push_back(new_manifest);
#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 211007 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init211007_charge(analysis_manifest& manifest)
	{
		std::cout << "Initializing 211007 charge only manifests..." << std::endl;
#define PAIR std::pair<double, double>
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/211007/";
		default_exp_manifest.out_folder = "../Data/211007/results_v2/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Data/211007/results_v2/";
		default_exp_manifest.out_picture_folder = "";
		default_exp_manifest.draw_only = false; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):

		default_exp_manifest.runs.push(0, 9999);
		default_exp_manifest.sub_runs.push(0, default_exp_manifest.subruns_per_file - 1);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		//END OF MODIFY ONLY THIS BLOCK

		channel_manifest ch_manifest;
		ch_manifest.save_with_indices = false;
		ch_manifest.baseline.baseline_by_average = true;
		ch_manifest.find_average = false;
		ch_manifest.find_integral = false;
		ch_manifest.find_double_integral = false;
		ch_manifest.N_extrapolation = 1;
		//----------------------------------------------
		//Set channel specifics for all experiments (kV)
		//----------------------------------------------
		//Charge signal
		ch_manifest.device = "Charge";
		ch_manifest.invert = false;
		ch_manifest.display.Y_limits = PAIR(-DBL_MAX, DBL_MAX);
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.baseline.baseline_range = PAIR(0, 18.5);
		ch_manifest.peaks.do_find = false;
		ch_manifest.find_average = true;
		ch_manifest.display.do_draw = false;
		default_exp_manifest.channels.push(10, ch_manifest);

		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_20.0kV_800V_46V_12dB_2650V") + "/");
		//this is a place to tweak individual channel for specific fields
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_20.0kV_800V_2650V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_18.8kV_800V_46V_12dB_2650V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_18.8kV_800V_2650V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_17.0kV_800V_46V_12dB_2650V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_17.0kV_800V_2650V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_15.3kV_800V_46V_12dB_2650V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_15.3kV_800V_2650V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_14.4kV_800V_46V_12dB_2650V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_14.4kV_800V_2650V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_13.6kV_800V_46V_12dB_2650V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_13.6kV_800V_2650V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_12.7kV_800V_46V_12dB_2650V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_12.7kV_800V_2650V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_11.8kV_800V_46V_12dB_2650V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_11.8kV_800V_2650V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_11.0kV_800V_46V_12dB_2650V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_11.0kV_800V_2650V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_19.9kV_800V_46V_12dB_2550V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_19.9kV_800V_2550V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_19.8kV_800V_46V_12dB_2500V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_19.8kV_800V_2500V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211007_Pu_19.9kV_800V_46V_12dB_2600V") + "/");
		new_manifest.accepted_events_fname = "../Post_processing/211007/results_v5/Pu_46V_19.9kV_800V_2600V/forms_Pu_peak/events.txt";
		if (!read_accepted_events(new_manifest.accepted_events_fname, new_manifest.accepted_events_data)) {
			std::cout << "Init211007_charge:: No event selection for "<<new_manifest.name<<"! Processing everything" << std::endl;
			new_manifest.accepted_events_data.clear();
		}
		manifest.manifests.push_back(new_manifest);
#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 211007 charge only manifests" << std::endl;
		return true;
	}
};

