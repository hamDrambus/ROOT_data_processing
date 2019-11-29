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
	STD_CONT <experiment_area> areas_to_draw;
	std::string this_path;
	int subruns_per_file = 10;
	//bool override_analysis = true;
	experiment_area exp_area;
	int threads_number = 1; //obv. must be >=1
	bool draw_only = false;
	accepted_events<double> events_to_process;
	std::pair<double, double> pics_t_zoom;
	double pics_trigger_position;
	std::string save_pics_to;
	std::size_t max_pics_number;
	bool gnuplot_presits;

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

	int gnuplot_pad_size = 400;
	int gnuplot_max_size = 1600;
	int gnuplot_width = 1200; //default for gnuplot is 640

	bool draw_required(ParameterPile::experiment_area what)
	{
		for (auto i = areas_to_draw.begin(); i != areas_to_draw.end(); ++i)
			if (i->contains(what))
				return true;
		return false;
	}

	bool read_accepted_events(std::string file, ParameterPile::accepted_events<double> &info)
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
		//clear the GEM output file
#ifdef _PROCESS_GEMS
		std::ofstream file;
		open_output_file(std::string(OUTPUT_DIR) + OUTPUT_GEMS + "/" + OUTPUT_GEMS + ".txt", file);//creates folder AND trucates the file
		file << "Experiment\tIntegral[V*us]\tt_from\tt_to\tN_valid\tN\tS_cutoff" << std::endl;
		file.close();
#endif

		TThread::Initialize();
		
		subruns_per_file = 1000;
		MPPC_threshold = 0.0070; //
		threads_number = 1;
		draw_only = true;
		std::string accepted_events_fname = "../Post_processing/190404/results_v3/Cd_46V_20kV_850V/forms_Cd_peak/13_events_cuts_06+07+08+12.txt";
		save_pics_to = "../Post_processing/190404/results_v3/Cd_46V_20kV_850V/forms_Cd_peak/events/";
		pics_t_zoom = std::pair<double, double> (0, 160);
		pics_trigger_position = 32;
		max_pics_number = 50;
		gnuplot_presits = false;
		read_accepted_events(accepted_events_fname, events_to_process);
		S1_start_time = 18.5; //in us
		S1_finish_time = 19.0; //in us

		filter_PMT_n_points.insert		(std::pair<int,int>(0,50));//(channel, value)
		filter_PMT_order.insert			(std::pair<int,int>(0,8));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(0,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(1,50));
		filter_PMT_order.insert			(std::pair<int,int>(1,8));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(1,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(8,32));
		filter_PMT_order.insert			(std::pair<int,int>(8,10));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(8,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(9,32));
		filter_PMT_order.insert			(std::pair<int,int>(9,10));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(9,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(10,32));
		filter_PMT_order.insert			(std::pair<int,int>(10,10));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(10,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(11,32));
		filter_PMT_order.insert			(std::pair<int,int>(11,10));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(11,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(12,32));
		filter_PMT_order.insert			(std::pair<int,int>(12,10));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(12,0));

		PMT_thresh.insert 	(std::pair<int,double>(8,0.0031));//(channel, value)
		PMT_thresh_edges.insert		(std::pair<int,double>(8,0.0));
		PMT_thresh.insert 	(std::pair<int,double>(9,0.0040));
		PMT_thresh_edges.insert		(std::pair<int,double>(9,0.0));
		PMT_thresh.insert 	(std::pair<int,double>(10,0.0031));
		PMT_thresh_edges.insert		(std::pair<int,double>(10,0.0));
		PMT_thresh.insert 	(std::pair<int,double>(11,0.0060));
		PMT_thresh_edges.insert		(std::pair<int,double>(11,0.0));

		PMT_thresh.insert 	(std::pair<int,double>(12,0.020));
		PMT_thresh_edges.insert		(std::pair<int,double>(12,0.0));
		PMT_thresh.insert 	(std::pair<int,double>(13,0.020));
		PMT_thresh_edges.insert		(std::pair<int,double>(13,0.0));
		PMT_thresh.insert 	(std::pair<int,double>(14,0.020));
		PMT_thresh_edges.insert		(std::pair<int,double>(10,0.0));
		PMT_thresh.insert 	(std::pair<int,double>(15,0.020));
		PMT_thresh_edges.insert		(std::pair<int,double>(15,0.0));
		PMT_thresh.insert 	(std::pair<int,double>(16,0.080));
		PMT_thresh_edges.insert		(std::pair<int,double>(16,0.0));

		//ch_use_average.push_pair(GEM_CH_, GEM_CH_);

		//ch_integrate_S2.push_pair(12, 16);
		ch_use_curved_baseline.push_pair(12, 15);
		//ch_use_curved_baseline.push_pair(9, 12);
		//ch_use_curved_baseline.push_pair(32, 64);
		ch_inverse.push_pair(8, 11);
		ch_inverse.push_pair(32, 64);

		S2_start_time.insert(std::pair<std::string,double>("190404_Cd_20kV_850V_46V_th250mV", 20));

		S2_finish_time.insert(std::pair<std::string,double>("190404_Cd_20kV_850V_46V_th250mV", 70));

		areas_to_draw.push_back(experiment_area());

		areas_to_draw.back().experiments.push_back("190404_Cd_20kV_850V_46V_th250mV");

		areas_to_draw.back().runs.push_pair(0, 9999);

		areas_to_draw.back().channels.push_pair(10, 10);

		//areas_to_draw.back().channels.push_pair(32, 44); //13
		//areas_to_draw.back().channels.push_pair(48, 59); //12 =>25 channels

		areas_to_draw.back().sub_runs.push(0, subruns_per_file-1);

		exp_area.runs.push_pair(0, 34);
		exp_area.channels.push_pair(10, 10);
		//exp_area.channels.push_pair(8, 12);
		//exp_area.channels.push_pair(GEM_CH_, GEM_CH_);

		//exp_area.channels.push_pair(32, 44); //13
		//exp_area.channels.push_pair(48, 59); //12 =>25 channels

		exp_area.sub_runs.push(0, subruns_per_file-1);

		exp_area.experiments.push_back("190404_Cd_20kV_850V_46V_th250mV");
	}
};
