#include <stdio.h>

#include "GlobalParameters.h"
#include "AnalysisManager.h"

DITERATOR iter_add(DITERATOR& to, int what, DITERATOR& end)
{
	if (what < 0)
		return end;
	return ((int)(end - to) < what) ? end : to + what;
}

void open_output_file(std::string name, std::ofstream &str, std::ios_base::openmode _mode)
{
	ensure_file(name);
	str.open(name.c_str(), _mode);
	if (!str.is_open()){
		std::cout << "Failed to open \"" << name << "\"" << std::endl;
	}
}

void ensure_file(std::string fname)
{
	std::string folder = fname;
	while ((folder.back() != '\\') &&(folder.back()!='/') &&!folder.empty())
		folder.pop_back();
	if (!folder.empty())
		folder.pop_back();
	ensure_folder(folder);
}

void ensure_folder(std::string folder)
{
#if defined(__WIN32__)
	if (!folder.empty()) {
		DWORD ftyp = GetFileAttributesA(folder.c_str());
		if (!(ftyp & FILE_ATTRIBUTE_DIRECTORY) || ftyp == INVALID_FILE_ATTRIBUTES) {
			int code = system(("mkdir \"" + folder + "\"").c_str());
			if (code)
				std::cout << "mkdir error: " << GetLastError() << std::endl;
		}
	}
#else
	struct stat st;
	if (-1==stat(folder.c_str(), &st)) {
		int err = errno;
		switch (err) {
		case (EACCES): {
			std::cout<<"Access error"<<std::endl;
			break;
		}
		case (ENAMETOOLONG): {
			std::cout<<"Path is too long"<<std::endl;
			break;
		}
		case (ENOENT) :
		case (ENOTDIR): {
			int code = system(("mkdir -p \"" + folder + "\"").c_str());
			if (code)
				std::cout << "mkdir -p error: " << code << std::endl;
			break;
		}
		default:{
			std::cout<<"stat(\""<<folder<<"\") returned -1; errno == "<<err<<std::endl;
			break;
		}
		}
	} else {
		if (!S_ISDIR(st.st_mode)) {
			int code = system(("mkdir -p \"" + folder + "\"").c_str());
			if (code)
				std::cout << "mkdir -p error: " << code << std::endl;
		}
	}
#endif //_WIN32__
}

std::string strtoken(std::string &in, std::string break_symbs)
{
	std::string out;
	while (!in.empty()) {
		char a = in.front();
		in.erase(in.begin());
		bool break_ = false;
		for (auto h = break_symbs.begin(); h != break_symbs.end(); ++h)
			if (a == *h) {
				break_ = true;
				break;
			}
		if ((break_) && (out.empty()))
			continue;
		if (break_)
			return out;
		out.push_back(a);
	}
	return out;
}

void DrawFileData(std::string name, std::vector<double> xs, std::vector<double> ys, ParameterPile::DrawEngine de)
{
	if (xs.size() != ys.size()){
		std::cout << "DrawFileData::input data size mismatch" << std::endl;
		return;
	}
	std::string mod_name = name;
	for (int s = 0; s < mod_name.size(); s++)
		if (mod_name[s] == '\\' || mod_name[s] == '/')
			mod_name[s] = '.';
	std::ofstream file;
	open_output_file("temp_gnuplot_files/" + mod_name, file);
	std::cout << "file " << "temp_gnuplot_files/" + mod_name << ".is_open() " << file.is_open() << std::endl;
#if defined(__WIN32__)
	if (!file.is_open())
		std::cout << GetLastError() << std::endl;
#endif
	for (int h = 0; h < xs.size(); h++)
		file << xs[h] << '\t' << ys[h] << std::endl;
	file.close();
	open_output_file("temp_gnuplot_files/script.sc", file);
	file << "plot '" << ParameterPile::this_path + "temp_gnuplot_files/" + mod_name << "' u 1:2" << std::endl;
	file << "pause -1";
	file.close();
	INVOKE_GNUPLOT(ParameterPile::this_path + "temp_gnuplot_files/script.sc");
}

namespace ParameterPile
{
	analysis_manifest gManifest;
	std::string this_path;
	int threads_number; //obv. must be >=1
	std::size_t max_pics_number;
	bool gnuplot_presits;
	bool quiet_mode;

	int Max_iteration_N = 1;

	int gnuplot_pad_size = 400;
	int gnuplot_max_size = 1600;
	int gnuplot_width = 1200; //default for gnuplot is 640

	bool read_accepted_events(std::string file, accepted_events<double> &info)
	{
		std::ifstream str;
		str.open(file);
		if (!str.is_open()) {
			std::cerr << "ParameterPile::read_accepted_events: Error: Failed to open \"" << file << "\"" << std::endl;
			return false;
		}
		std::string line, word;
		int line_n = 0;
		while (!str.eof() && str.is_open()) {
			std::getline(str, line);
			++line_n;
			if (line.size() >= 2) //Ignore simple c style comment
				if ((line[0] == '/') && (line[1] == '/'))
					continue;
			try {
				word = strtoken(line, "\t ");
				int run = std::stoi(word);
				word = strtoken(line, "\t ");
				int subrun = std::stoi(word);
				word = strtoken(line, "\t ");
				double trigger = 0;
				if (!word.empty())
					trigger = std::stod(word);
				info.push(run, subrun, trigger);
			}
			catch (std::exception &e) {
				std::cerr << "ParameterPile::read_accepted_events: Error: Exception in \"" << file << "\"" << std::endl
					<< "\ton line " << line_n << std::endl;
				std::cerr << "\t" << e.what() << std::endl;
				continue;
			}
		}
		return true;
	}

	void Init_globals(void)
	{
		char path[FILENAME_MAX];
#if defined(__WIN32__)
		this_path = _getcwd(path, FILENAME_MAX);
#else
		this_path = getcwd(path, FILENAME_MAX);
#endif //__WIN32__
		if (!this_path.empty())
			if (this_path.back()!='/')
				this_path.push_back('/');

		TThread::Initialize();
		threads_number = 9; //obv. must be >=1
		max_pics_number = 25;
		gnuplot_presits = false;
		quiet_mode = true;

		Init190404(gManifest);
        Init180705(gManifest);
	}

	bool Init190404(analysis_manifest& manifest)
	{
		std::cout<<"Initializing 190404 manifests..."<<std::endl;
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
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data))
			default_exp_manifest.accepted_events_data.clear();
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
		//chs_to_draw.push(11, 12);
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
		for (int ch = SiPM_channels.get_next_index(); ch>=0; ch = SiPM_channels.get_next_index()) {
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


		for (int ch = SiPM_channels.get_next_index(); ch>=0; ch = SiPM_channels.get_next_index()) {
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
		std::cout<<"Loaded 190404 manifests"<<std::endl;
		return true;
	}


	bool Init180705(analysis_manifest& manifest)
	{
		std::cout<<"Initializing 180705 manifests..."<<std::endl;
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
		//chs_to_draw.push(38, 38);
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
		for (int ch = SiPM_channels.get_next_index(); ch>=0; ch = SiPM_channels.get_next_index()) {
			ch_manifest.display.do_draw = chs_to_draw.contains(ch);
			if (ch==44 || ch==56 || ch==57 || ch==41)
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

		for (int ch = SiPM_channels.get_next_index(); ch>=0; ch = SiPM_channels.get_next_index()) {
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
		std::cout<<"Loaded 180705 manifests"<<std::endl;
		return true;
	}
};
