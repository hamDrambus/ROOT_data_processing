#include "GlobalParameters.h"

namespace ParameterPile
{
	//Use this function to play with/display different configurations
	//The correct parameters used in the analysis are in Init190404(analysis_manifest& manifest)
	bool Init190307Xray_tests(analysis_manifest& manifest)
	{
		std::cout << "Initializing TEST 190307 (X-ray) manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/190307/";
		default_exp_manifest.out_folder = "../Data/190307/results_vt/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Data/190307/results_vt/gnuplot/";
		default_exp_manifest.out_picture_folder = "";//Do not save pictures
		default_exp_manifest.draw_only = true; //if set to true, no data is written to output
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init190307Xray_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(1, 7); //20kV
		//default_exp_manifest.runs.push(12, 13); //18kV
		//default_exp_manifest.runs.push(17, 19); //16kV
		//default_exp_manifest.runs.push(23, 25); //14kV
		//default_exp_manifest.runs.push(29, 31); //12kV
		//default_exp_manifest.runs.push(35, 37); //10kV
		//default_exp_manifest.runs.push(41, 43); //9kV
		default_exp_manifest.runs.push(1, 1);
		default_exp_manifest.sub_runs.push(0, 19);// default_exp_manifest.subruns_per_file - 1);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		chs_to_draw.push(32, 32);
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
		ch_manifest.baseline.baseline_range = PAIR(0, 25);
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
		ch_manifest.baseline.curved_range = PAIR(20, 155);
		ch_manifest.baseline.curved_center = PAIR(43, 130);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.056;
		ch_manifest.peaks.threshold_cutoff = 0.011;
		default_exp_manifest.channels.push(0, ch_manifest);  //Sum of 3 slow PMTs
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.040;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(1, ch_manifest); //Sum of 4 slow PMTs
		ch_manifest.display.do_draw = chs_to_draw.contains(2);

		ch_manifest.device = "Q";
		ch_manifest.peaks.do_find = false;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.0;
		ch_manifest.peaks.threshold_cutoff = 0.0;
		ch_manifest.find_average = true;
		default_exp_manifest.channels.push(2, ch_manifest);


		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.peaks.do_find = true;
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_numberIterations = 60;
		ch_manifest.baseline.curved_direction = TSpectrum::kBackDecreasingWindow;
		ch_manifest.baseline.curved_filterOrder = TSpectrum::kBackOrder2;
		ch_manifest.baseline.curved_smoothing = true;
		ch_manifest.baseline.curved_smoothWindow = TSpectrum::kBackSmoothing3;
		ch_manifest.baseline.curved_compton = false;
		ch_manifest.baseline.curved_sparse = 1;
		ch_manifest.baseline.curved_range = PAIR(40, 130);
		ch_manifest.baseline.curved_center = PAIR(45, 112);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.peaks.threshold = 0.0140;
		ch_manifest.peaks.threshold_cutoff = 0.004;
		ch_manifest.find_average = false;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}


		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190307_X-ray_20kV_650V_46V_coll_6mm") + "/");
		manifest.manifests.push_back(new_manifest);


#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded TEST 190307 (X-ray) manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init190307Xray(analysis_manifest& manifest)
	{
		std::cout << "Initializing 190307 (X-ray) manifests..." << std::endl;
#define PAIR std::pair<double, double>
		area_vector SiPM_channels;
		SiPM_channels.push(32, 44);
		SiPM_channels.push(48, 59);
		experiment_manifest default_exp_manifest;
		default_exp_manifest.data_time_constant = 1.6e-2;
		default_exp_manifest.data_voltage_channels = 4095;
		default_exp_manifest.data_voltage_amplitude = 2.0;
		default_exp_manifest.data_voltage_of_zero_channel = -1.0;
		default_exp_manifest.in_folder = "../Data/190307/";
		default_exp_manifest.out_folder = "../Data/190307/results_v1/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Data/190307/results_v1/gnuplot/";
		default_exp_manifest.out_picture_folder = "";//Do not save pictures
		default_exp_manifest.draw_only = false; //if set to true, no data is written to output
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init190307Xray_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(1, 7); //20kV
		//default_exp_manifest.runs.push(12, 13); //18kV
		//default_exp_manifest.runs.push(17, 19); //16kV
		//default_exp_manifest.runs.push(23, 25); //14kV
		//default_exp_manifest.runs.push(29, 31); //12kV
		//default_exp_manifest.runs.push(35, 37); //10kV
		//default_exp_manifest.runs.push(41, 43); //9kV
		default_exp_manifest.sub_runs.push(0, default_exp_manifest.subruns_per_file - 1);// default_exp_manifest.subruns_per_file - 1);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		//chs_to_draw.push(32, 32);
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
		ch_manifest.baseline.baseline_range = PAIR(0, 25);
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
		ch_manifest.baseline.curved_range = PAIR(20, 155);
		ch_manifest.baseline.curved_center = PAIR(43, 130);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.056;
		ch_manifest.peaks.threshold_cutoff = 0.011;
		default_exp_manifest.channels.push(0, ch_manifest);  //Sum of 3 slow PMTs
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.040;
		ch_manifest.peaks.threshold_cutoff = 0.005;
		default_exp_manifest.channels.push(1, ch_manifest); //Sum of 4 slow PMTs
		ch_manifest.display.do_draw = chs_to_draw.contains(2);

		ch_manifest.device = "Q";
		ch_manifest.peaks.do_find = false;
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.0;
		ch_manifest.peaks.threshold_cutoff = 0.0;
		ch_manifest.find_average = true;
		default_exp_manifest.channels.push(2, ch_manifest);


		//SiPMs
		ch_manifest.device = "SiPM";
		ch_manifest.peaks.do_find = true;
		ch_manifest.invert = true;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_numberIterations = 60;
		ch_manifest.baseline.curved_direction = TSpectrum::kBackDecreasingWindow;
		ch_manifest.baseline.curved_filterOrder = TSpectrum::kBackOrder2;
		ch_manifest.baseline.curved_smoothing = true;
		ch_manifest.baseline.curved_smoothWindow = TSpectrum::kBackSmoothing3;
		ch_manifest.baseline.curved_compton = false;
		ch_manifest.baseline.curved_sparse = 1;
		ch_manifest.baseline.curved_range = PAIR(40, 130);
		ch_manifest.baseline.curved_center = PAIR(45, 112);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.filter.order = 4;
		ch_manifest.filter.n_points = 15;
		ch_manifest.filter.n_iterations = 0;
		ch_manifest.peaks.threshold = 0.0140;
		ch_manifest.peaks.threshold_cutoff = 0.004;
		ch_manifest.find_average = false;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}


		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190307_X-ray_20kV_650V_46V_coll_6mm") + "/");
		manifest.manifests.push_back(new_manifest);
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190307_X-ray_18kV_650V_46V_coll_6mm") + "/");
		manifest.manifests.push_back(new_manifest);
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190307_X-ray_16kV_650V_46V_coll_6mm") + "/");
		manifest.manifests.push_back(new_manifest);
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190307_X-ray_14kV_650V_46V_coll_6mm") + "/");
		manifest.manifests.push_back(new_manifest);
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190307_X-ray_12kV_650V_46V_coll_6mm") + "/");
		manifest.manifests.push_back(new_manifest);
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190307_X-ray_10kV_650V_46V_coll_6mm") + "/");
		manifest.manifests.push_back(new_manifest);
		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190307_X-ray_9kV_650V_46V_coll_6mm") + "/");
		manifest.manifests.push_back(new_manifest);
#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 190307 (X-ray) manifests" << std::endl;
		return true;
	}

};
