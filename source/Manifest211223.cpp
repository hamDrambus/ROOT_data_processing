#include "GlobalParameters.h"

namespace ParameterPile
{
	//Use this function to play with/display different configurations
	//The correct parameters used in the analysis are in Init190404(analysis_manifest& manifest)
	bool Init211223_tests(analysis_manifest& manifest)
	{
		std::cout << "Initializing TEST 211223 manifests..." << std::endl;
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
		default_exp_manifest.in_folder = "../Data/211223/";
		default_exp_manifest.out_folder = "../Data/211223/results_vt/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init211223_tests:: No event selection - processing everything" << std::endl;
			default_exp_manifest.accepted_events_data.clear();
		}
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs_to_draw.push(0, 9999); //DRAW all
		default_exp_manifest.sub_runs_to_draw.push(0, default_exp_manifest.subruns_per_file - 1); //DRAW all

		//MODIFY ONLY THIS BLOCK AND DISPLAY-RELATED VALUES FOR CHANNELS
		default_exp_manifest.out_gnuplot_folder = "../Data/211223/results_vt/gnuplot/";
		default_exp_manifest.out_picture_folder = "";
		default_exp_manifest.draw_only = true; //if set to true, no data is written to output
		//default_exp_manifest.runs.push(0, 9999); //Use only when all invalid files are deleted from folders.
		//List of valid files (runs):
		//default_exp_manifest.runs.push(804, 823); //f32 11.1kV 3592V // P = 1.5 atm
		//default_exp_manifest.runs.push(788, 802); //f31 11.1kV 3233V
		//default_exp_manifest.runs.push(772, 786); //f30 11.1kV 2874V
		//default_exp_manifest.runs.push(756, 770); //f29 11.1kV 2514V
		//default_exp_manifest.runs.push(735, 754); //f28 11.1kV 2155V
		//default_exp_manifest.runs.push(714, 733); //f27 11.1kV 1796V
		//default_exp_manifest.runs.push(693, 712); //f26 11.1kV 1615V
		//default_exp_manifest.runs.push(662, 691); //f25 11.1kV 1332V
		//default_exp_manifest.runs.push(630, 660); //f24 11.1kV 1065V
		//default_exp_manifest.runs.push(598, 628); //f23 11.1kV 852V
		//default_exp_manifest.runs.push(568, 596); //f22 11.1kV 744V
		//default_exp_manifest.runs.push(537, 566); //f21 11.1kV 639V
		//default_exp_manifest.runs.push(506, 535); //f20 11.1kV 532V
		//default_exp_manifest.runs.push(475, 504); //f19 11.1kV 426V
		//default_exp_manifest.runs.push(444, 473); //f18 11.1kV 318V
		//default_exp_manifest.runs.push(412, 442); //f17 11.1kV 0V

		//default_exp_manifest.runs.push(396, 410); //f16 7.8kV 2513V // P = 1.0 atm
		//default_exp_manifest.runs.push(378, 394); //f15 7.8kV 2260V
		//default_exp_manifest.runs.push(362, 376); //f14 7.8kV 2009V
		//default_exp_manifest.runs.push(346, 360); //f13 7.8kV 1758V
		//default_exp_manifest.runs.push(325, 344); //f12 7.8kV 1507V
		//default_exp_manifest.runs.push(303, 323); //f11 7.8kV 1256V
		//default_exp_manifest.runs.push(282, 301); //f10 7.8kV 1130V
		//default_exp_manifest.runs.push(251, 280); //f9 7.8kV 931V
		//default_exp_manifest.runs.push(220, 249); //f8 7.8kV 745V
		//default_exp_manifest.runs.push(189, 218); //f7 7.8kV 596V
		//default_exp_manifest.runs.push(158, 187); //f6 7.8kV 522V
		//default_exp_manifest.runs.push(127, 156); //f5 7.8kV 448V
		//default_exp_manifest.runs.push(95, 125); //f4 7.8kV 372V
		//default_exp_manifest.runs.push(64, 93); //f3 7.8kV 298V
		//default_exp_manifest.runs.push(32, 62); //f2 7.8kV 222V
		//default_exp_manifest.runs.push(1, 30); //f1 7.8kV 0V

		default_exp_manifest.runs.push(804, 804);
		default_exp_manifest.sub_runs.push(0, 19);
		default_exp_manifest.trigger_at = -32;

		area_vector chs_to_draw;	 //DRAW only these channels
		chs_to_draw.push(38); //38, 44, 42, 39
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
		ch_manifest.baseline.curved_range = PAIR(0, 160);
		ch_manifest.baseline.curved_center = PAIR(20, 75);
		ch_manifest.baseline.curved_trim = PAIR(3, 3);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.044;
		ch_manifest.peaks.threshold_cutoff = 0.010;
		default_exp_manifest.channels.push(0, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.013;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		ch_manifest.filter.n_iterations = 1;
		default_exp_manifest.channels.push(1, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(2);
		ch_manifest.peaks.threshold = 0.013;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		ch_manifest.filter.n_iterations = 1;
		default_exp_manifest.channels.push(2, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

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
		ch_manifest.peaks.threshold_cutoff = 0.004;
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
		ch_manifest.filter.n_points = 20;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 90);
		ch_manifest.baseline.curved_center = PAIR(20, 65);
		ch_manifest.baseline.curved_trim = PAIR(3, 3);
		ch_manifest.baseline.curved_numberIterations = 50;
		ch_manifest.peaks.threshold = 0.0065;
		ch_manifest.peaks.threshold_cutoff = 0.0008;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}
		default_exp_manifest.channels.info(44)->peaks.threshold = 0.0085;

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 211223 DATA!
		//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_3592V") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_3233V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_2874V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_2514V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_2155V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_1796V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_1615V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_1332V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_1065V") + "/");
		manifest.manifests.push_back(new_manifest);


#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded TEST 211223 manifests" << std::endl;
		return true;
	}

	//Do not touch this unless intending to redo the analysis
	bool Init211223(analysis_manifest& manifest)
	{
		std::cout << "Initializing 211223 manifests..." << std::endl;
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
		default_exp_manifest.in_folder = "../Data/211223/";
		default_exp_manifest.out_folder = "../Data/211223/results_v1/";
		default_exp_manifest.write_event_indices = false;
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data)) {
			std::cout << "Init211223:: No event selection - processing everything" << std::endl;
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
		ch_manifest.baseline.curved_range = PAIR(0, 160);
		ch_manifest.baseline.curved_center = PAIR(20, 75);
		ch_manifest.baseline.curved_trim = PAIR(3, 3);
		ch_manifest.display.do_draw = chs_to_draw.contains(0);
		ch_manifest.peaks.threshold = 0.044;
		ch_manifest.peaks.threshold_cutoff = 0.010;
		default_exp_manifest.channels.push(0, ch_manifest);

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.display.do_draw = chs_to_draw.contains(1);
		ch_manifest.peaks.threshold = 0.013;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		ch_manifest.filter.n_iterations = 1;
		default_exp_manifest.channels.push(1, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.display.do_draw = chs_to_draw.contains(2);
		ch_manifest.peaks.threshold = 0.013;
		ch_manifest.peaks.threshold_cutoff = 0.003;
		ch_manifest.filter.n_iterations = 1;
		default_exp_manifest.channels.push(2, ch_manifest);
		ch_manifest.filter.n_iterations = 0;

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
		ch_manifest.peaks.threshold_cutoff = 0.004;
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
		ch_manifest.filter.n_points = 20;
		ch_manifest.filter.n_iterations = 1;
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 90);
		ch_manifest.baseline.curved_center = PAIR(20, 65);
		ch_manifest.baseline.curved_trim = PAIR(3, 3);
		ch_manifest.baseline.curved_numberIterations = 50;
		ch_manifest.peaks.threshold = 0.0065;
		ch_manifest.peaks.threshold_cutoff = 0.0008;
		for (int ch = SiPM_channels.get_next_index(); ch >= 0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			default_exp_manifest.channels.push(ch, ch_manifest);
		}
		default_exp_manifest.channels.info(44)->peaks.threshold = 0.0085;

		//ALL PARAMETERS (THRESHOLDS) ARE THE SAME FOR ALL FOLDERS FOR 211223 DATA!
		//default_exp_manifest.channels.info(3)->baseline.do_find_curved = false;
		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_3592V") + "/");
		//this is a place to tweak individual channel for specific fields
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_3233V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_2874V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_2514V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_2155V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_1796V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_1615V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_1332V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_1065V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_852V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_744V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_639V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_532V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_426V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_318V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.5atm_11.1kV_850V_46V_12dB_0V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_2513V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_2260V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_2009V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_1758V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_1507V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_1256V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_1130V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_931V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_745V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_596V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_522V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_448V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_372V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_298V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_222V") + "/");
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "211223_Pu_1.0atm_7.8kV_850V_46V_12dB_0V") + "/");
		manifest.manifests.push_back(new_manifest);
#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
		std::cout << "Loaded 211223 manifests" << std::endl;
		return true;
	}
};

