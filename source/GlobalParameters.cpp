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

bool isSameChannels(const STD_CONT<int>& a, const STD_CONT<int>& b)
{
	if (a.size() != b.size())
		return false;
	else {
		for (std::size_t ind = 0, ind_end_ = a.size(); ind<ind_end_; ++ind)
			if (a[ind]!=b[ind])
				return false;
	}
	return true;
}

int getIndex(const STD_CONT<int>& channels, int ch)
{
	for (std::size_t i = 0, i_end_ = channels.size(); i!=i_end_; ++i)
		if (channels[i]==ch)
			return i;
	return -1;
}

void DrawFileData(std::string name, std::vector<double> xs, std::vector<double> ys, ParameterPile::DrawEngine de)
{
	if (xs.size() != ys.size()){
		std::cout << "DrawFileData::input data size mismatch" << std::endl;
		return;
	}
	if (de == ParameterPile::DrawEngine::ROOT){
		double *xxxs = new double[xs.size()];
		double *yyys = new double[ys.size()];
		for (int h = 0; h < xs.size(); h++){
			xxxs[h] = xs.at(h);
			yyys[h] = ys.at(h);
		}
		TGraph* gr = new TGraph(xs.size(), xxxs, yyys);
		TCanvas* can = new TCanvas(name.c_str(), name.c_str());
		can->SetWindowPosition(100, 150);
		can->Update();
		gr->Draw();
		delete[] xxxs;
		delete[] yyys;
	} else {
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
}

namespace ParameterPile
{
	analysis_manifest gManifest;
	std::string this_path;
	int threads_number = 1; //obv. must be >=1
	std::size_t max_pics_number;
	bool gnuplot_presits;

<<<<<<< HEAD
	double dt_quant = 0.1; //us

	//int filter_MPPC_n_points = 15;	//62.5 MHz
	//int filter_MPPC_order = 4;		//62.5 MHz
	int filter_MPPC_n_points = 32;	//250 MHz
	int filter_MPPC_order = 8;		//250 MHz
	int filter_MPPC_n_iterations = 0;
	std::map<int,int> filter_PMT_n_points;
	std::map<int,int> filter_PMT_order;
	std::map<int,int> filter_PMT_n_iterations;

	//these values are approximate, especially S2
	double S1_start_time = 20; //in us
	double S1_finish_time = 20.5; //in us
	std::map<std::string,double> S2_start_time; //in us
	std::map<std::string,double> S2_finish_time; //in us

	double GEM_threshold_to_noise = 1.1;
	int GEM_N_of_averaging = 30; //=== N_trust
	
	double PMT_run_acceptance_threshold_to_noise = 4;//~ 2-3
	std::map<int,double> PMT_thresh;
	std::map<int,double> PMT_thresh_edges;

	int PMT_N_of_averaging = 1; //=== N_trust. PMT signal is already smoothed by filter
	double PMT_ROOTs_bl_left_offset = 9; //for baseline's baseline
	double PMT_ROOTs_bl_right_offset = 20; //for baseline's baseline
	double PMT_ROOTs_bl_trim = 1; //remove edge artifacts
	double PMT_ROOTs_bl_from = 10;
	double PMT_ROOTs_bl_to = 110;

	area_vector ch_use_average;
	area_vector ch_integrate_S2;
	area_vector ch_use_curved_baseline;
	area_vector ch_inverse;
	area_vector ch_to_sum;

	double MPPC_peaks_smoothing_time = 5; //us
	int MPPC_N_trust = 1;
	int MPPC_double_I_N_trust = 1;
	double MPPC_ROOTs_bl_left_offset = 5; //for baseline's baseline
	double MPPC_ROOTs_bl_right_offset = 15; //for baseline's baseline
	double MPPC_ROOTs_bl_trim = 1; //remove edge artifacts
	double MPPC_ROOTs_bl_from = 10;
	double MPPC_ROOTs_bl_to = 110;
	double MPPC_threshold = 0.0120; //

	int Max_iteration_N = 2;
=======
	int Max_iteration_N = 1;
>>>>>>> 378a3efde92f0411e374bdcf346c9d3648a6059e

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
		threads_number = 1;
		max_pics_number = 50;
		gnuplot_presits = false;
	}

	bool Init190404(analysis_manifest manifest)
	{
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
		default_exp_manifest.out_folder = "../Data/190404/results/";
		default_exp_manifest.accepted_events_fname = "";
		if (!read_accepted_events(default_exp_manifest.accepted_events_fname, default_exp_manifest.accepted_events_data))
			default_exp_manifest.accepted_events_data.clear();
		default_exp_manifest.out_gnuplot_folder = "../Data/190404/results/temp_gnuplot/";
		default_exp_manifest.out_picture_folder = "";//Do not save pictures
		default_exp_manifest.subruns_per_file = 1000;
		default_exp_manifest.runs.push(0, 9999);
		default_exp_manifest.sub_runs.push(0, default_exp_manifest.subruns_per_file - 1);
		default_exp_manifest.runs_to_draw; //DRAW none
		default_exp_manifest.sub_runs_to_draw; //DRAW none

		channel_manifest ch_manifest;
		ch_manifest.peaks.do_find = true;
		ch_manifest.save_with_indices = false;
		ch_manifest.display.do_draw = false;					//DRAW
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
		ch_manifest.integral_range = PAIR(20, 40); //not tested
		ch_manifest.find_double_integral = false;
		ch_manifest.double_integral_range = PAIR(20, 40); //not tested

		ch_manifest.filter.n_iterations = 0;
		ch_manifest.filter.n_points = 2;
		ch_manifest.filter.order = 1;
		
		ch_manifest.peaks.threshold_cutoff = 0;

		ch_manifest.N_extrapolation = 1;
		//----------------------------------------------
		//Set channel specifics for all experiments (kV)
		//----------------------------------------------
		//fast PMTs
		ch_manifest.invert = true;
		ch_manifest.display.do_draw = false;					//DRAW
		ch_manifest.device = "PMT";
		ch_manifest.peaks.threshold = 0.0031;
		default_exp_manifest.channels.push(8, ch_manifest);
		ch_manifest.peaks.threshold = 0.0040;
		default_exp_manifest.channels.push(9, ch_manifest);
		ch_manifest.peaks.threshold = 0.0031;
		default_exp_manifest.channels.push(10, ch_manifest);
		ch_manifest.peaks.threshold = 0.0060;
		default_exp_manifest.channels.push(11, ch_manifest);

		//slow PMTs and trigger
		ch_manifest.invert = false;
		ch_manifest.display.do_draw = false;					//DRAW
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.baseline.curved_range = PAIR(10, 110);
		ch_manifest.baseline.curved_center = PAIR(20, 89);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(12, ch_manifest);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(13, ch_manifest);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(14, ch_manifest);
		ch_manifest.peaks.threshold = 0.020;
		default_exp_manifest.channels.push(15, ch_manifest);
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.080;
		default_exp_manifest.channels.push(16, ch_manifest);

		//SiPMs
		ch_manifest.invert = true;
		ch_manifest.display.do_draw = false;					//DRAW
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.device = "SiPM";
		ch_manifest.peaks.threshold = 0.0070;
		for (int ch = SiPM_channels.get_next_index(); ch>=0; ch = SiPM_channels.get_next_index())
			default_exp_manifest.channels.push(ch, ch_manifest);

		//sum of slow PMTs (for display only)
		ch_manifest.invert = false;
		ch_manifest.peaks.do_find = false;
		ch_manifest.display.do_draw = false;					//DRAW
		ch_manifest.baseline.do_find_curved = true;
		ch_manifest.device = "Virtual";
		ch_manifest.baseline.curved_range = PAIR(10, 110);
		ch_manifest.baseline.curved_center = PAIR(20, 89);
		ch_manifest.baseline.curved_trim = PAIR(1, 1);
		ch_manifest.peaks.threshold = 0.020; //Still required for baseline restoration
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(12, 15);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(100, ch_manifest);

		//sum of fast PMTs (for display only)
		ch_manifest.invert = true;
		ch_manifest.display.do_draw = false;					//DRAW
		ch_manifest.baseline.do_find_curved = false;
		ch_manifest.peaks.threshold = 0.0080;
		ch_manifest.summarize_channels.erase();
		ch_manifest.summarize_channels.push(8, 11);
		if (ch_manifest.display.do_draw || ch_manifest.find_average)
			default_exp_manifest.channels.push(101, ch_manifest);

		experiment_manifest new_manifest = default_exp_manifest;
		new_manifest.append_folder((new_manifest.name = "190404_Cd_20kV_850V_46V_th250mV") + "/");
		//this is a place to tweak individual channel for specific fields (e.g.:
		//new_manifest.channels.info(11)->threshold = 0.0055;
		//new_manifest.channels.info(38)->find_double_integral = true;
		//new_manifest.channels.info(38)->double_integral_range = std::pair<double, double>(25, 32);
		//) to calculate double integral for 38th SiPM and set 11 PMT threshold to 0.0055 instead of 0.0060
		manifest.manifests.push_back(new_manifest);

		new_manifest = default_exp_manifest;
		new_manifest.append_folder( (new_manifest.name = "190404_Cd_18kV_850V_46V_th240mV") + "/");
		manifest.manifests.push_back(new_manifest);

#undef PAIR
		for (std::size_t m = 0, m_end_ = manifest.manifests.size(); m != m_end_; ++m)
			manifest.manifests[m].channels.sort();
	}
};
