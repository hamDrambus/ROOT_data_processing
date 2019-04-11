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

	double dt_quant = 0.1; //us

	int filter_MPPC_n_points = 15;
	int filter_MPPC_order = 4;
	int filter_MPPC_n_iterations = 1;
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
	
	double PMT_run_acceptance_threshold_to_noize = 4;//~ 2-3
	std::map<int,double> PMT_maximum_thresh;
	std::map<int,double> PMT_minimum_thresh;
	std::map<int,double> PMT_thresh_edges;

	int PMT_N_of_averaging = 1; //=== N_trust. PMT signal is already smoothed by filter
	int PMT_N_peaks_acceptance = 0;//1 //condition is >=1, not >1
	double PMT_SArea_peaks_acceptance = 0; //V*us
	//done - for every runs S_acceptance is obtained from S distribution histogram
	//TODO: figure out the appropriate
	double PMT_min_fraction_above_cutoff = 0.2; //S cutoff is set so that at least 20% of runs are accepted (by S). See code for usage
	int PMT_min_statistics = 10; //TODO: figure out the number
	double PMT_mean_above_cutoff_acceptance = 0.05; //if cutoff is above the mean S this may be either too high cutoff
	//or there are just two large separate peaks in S histogram and it just happened that mean<cutoff
	//If there are almost no events (some small % of total) between cutoff and mean the cutoff is correct:
	//^ dN														^ dN
	//|															|			  mean
	//|wrong events     cutoff   events to be accepted			|			   |_          cutoff
	//|      _      mean  |       ___							|			 __| |           |
	//|    _| |		  |   |      |   |							|			 | | |__         |
	//|    |  |		  |   |      |   |							|			 | |    |        |
	//|    |   --	  |   |      |   |__			   VS		|_			 | |    |__      |
	//|    |    |	  |   _      |     |						| |_		_| |       |     |
	//| __ |    |	  |  | |_    |     |						|   |		|  |        |    |     _______
	//|_| |     |_	  | _|   |   |     |      S					|   |__     |  |        |    |    |       |
	//+-----------+----+------+-+------+------>					+------+----+--+--------+---------------------->
	//TODO: add illustrative pictures to the project for algorithms
	double PMT_right_cutoff_from_RMS = 4.5;//4.5
	double PMT_left_cutoff_from_RMS = 2.5;//2

	double PMT_ROOTs_bl_from_max_left = 15; //until the start of signal
	double PMT_ROOTs_bl_from_max_right = 80; //until the end of signal
	double PMT_ROOTs_bl_left_offset = 0; //for baseline's baseline
	double PMT_ROOTs_bl_right_offset = 20; //for baseline's baseline
	double PMT_ROOTs_bl_trim = 1; //remove edge artifacts

	area_vector ch_use_average;
	area_vector ch_integrate_S2;
	area_vector ch_use_curved_baseline;
	area_vector ch_inverse;

	double MPPC_peaks_smoothing_time = 5; //us
	int MPPC_N_trust = 1;
	int MPPC_double_I_N_trust = 1;
	double MPPC_ROOTs_bl_from_max_left = 30;
	double MPPC_ROOTs_bl_from_max_right = 50;
	double MPPC_ROOTs_bl_left_offset = 5; //for baseline's baseline
	double MPPC_ROOTs_bl_right_offset = 12; //for baseline's baseline
	double MPPC_ROOTs_bl_trim = 1; //remove edge artifacts
	double MPPC_threshold_to_noise = 5.5;
	double MPPC_minimum_peak_A = 0.0128; //
	double MPPC_maximum_peak_A = 0.0128; //

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
		
		PMT_SArea_peaks_acceptance = -1; //V*us
		subruns_per_file = 1000;
		MPPC_minimum_peak_A = 0.0070; //
		MPPC_maximum_peak_A = 0.0070; //
		threads_number = 4;
		S1_start_time = 18.5; //in us
		S1_finish_time = 19; //in us

		filter_PMT_n_points.insert		(std::pair<int,int>(0,50));//(channel, value)
		filter_PMT_order.insert			(std::pair<int,int>(0,8));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(0,1));
		filter_PMT_n_points.insert		(std::pair<int,int>(1,50));
		filter_PMT_order.insert			(std::pair<int,int>(1,8));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(1,1));
		filter_PMT_n_points.insert		(std::pair<int,int>(8,25));
		filter_PMT_order.insert			(std::pair<int,int>(8,10));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(8,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(9,25));
		filter_PMT_order.insert			(std::pair<int,int>(9,10));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(9,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(10,25));
		filter_PMT_order.insert			(std::pair<int,int>(10,10));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(10,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(11,25));
		filter_PMT_order.insert			(std::pair<int,int>(11,10));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(11,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(12,50));
		filter_PMT_order.insert			(std::pair<int,int>(12,8));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(12,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(13,50));
		filter_PMT_order.insert			(std::pair<int,int>(13,8));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(13,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(14,50));
		filter_PMT_order.insert			(std::pair<int,int>(14,8));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(14,0));
		filter_PMT_n_points.insert		(std::pair<int,int>(15,50));
		filter_PMT_order.insert			(std::pair<int,int>(15,8));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(15,0));

		PMT_minimum_thresh.insert 	(std::pair<int,double>(0,0.027));//(channel, value)
		PMT_maximum_thresh.insert	(std::pair<int,double>(0,0.027));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(1,0.046));
		PMT_maximum_thresh.insert	(std::pair<int,double>(1,0.046));

		PMT_minimum_thresh.insert 	(std::pair<int,double>(8,0.0031));
		PMT_maximum_thresh.insert	(std::pair<int,double>(8,0.0031));
		PMT_thresh_edges.insert		(std::pair<int,double>(8,0.0));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(9,0.0040));
		PMT_maximum_thresh.insert	(std::pair<int,double>(9,0.0040));
		PMT_thresh_edges.insert		(std::pair<int,double>(9,0.0));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(10,0.0031));
		PMT_maximum_thresh.insert	(std::pair<int,double>(10,0.0031));
		PMT_thresh_edges.insert		(std::pair<int,double>(10,0.0));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(11,0.0060));
		PMT_maximum_thresh.insert	(std::pair<int,double>(11,0.0060));
		PMT_thresh_edges.insert		(std::pair<int,double>(11,0.0));

		PMT_minimum_thresh.insert 	(std::pair<int,double>(12,0.020));
		PMT_maximum_thresh.insert	(std::pair<int,double>(12,0.020));
		PMT_thresh_edges.insert		(std::pair<int,double>(12,0.0));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(13,0.020));
		PMT_maximum_thresh.insert	(std::pair<int,double>(13,0.020));
		PMT_thresh_edges.insert		(std::pair<int,double>(13,0.0));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(14,0.020));
		PMT_maximum_thresh.insert	(std::pair<int,double>(14,0.020));
		PMT_thresh_edges.insert		(std::pair<int,double>(14,0.0));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(15,0.020));
		PMT_maximum_thresh.insert	(std::pair<int,double>(15,0.020));
		PMT_thresh_edges.insert		(std::pair<int,double>(15,0.0));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(16,0.080));
		PMT_maximum_thresh.insert	(std::pair<int,double>(16,0.080));
		PMT_thresh_edges.insert		(std::pair<int,double>(16,0.0));

		//ch_use_average.push_pair(0, 1);
		//ch_use_average.push_pair(8, 12);
		ch_use_average.push_pair(GEM_CH_, GEM_CH_);

		ch_integrate_S2.push_pair(12, 16);
		ch_use_curved_baseline.push_pair(12, 15);
		//ch_use_curved_baseline.push_pair(32, 65);
		ch_inverse.push_pair(8, 11);
		ch_inverse.push_pair(32, 65);

		S2_start_time.insert(std::pair<std::string,double>("190404_Cd_20kV_850V_46V_th250mV_0", 20));
		S2_start_time.insert(std::pair<std::string,double>("190404_Cd_20kV_850V_46V_th250mV", 20));
		S2_start_time.insert(std::pair<std::string,double>("190404_Cd_18kV_850V_46V_th230mV", 20));
		S2_start_time.insert(std::pair<std::string,double>("190404_Cd_16kV_850V_46V_th210mV", 20));

		S2_finish_time.insert(std::pair<std::string,double>("190404_Cd_20kV_850V_46V_th250mV_0",70));
		S2_finish_time.insert(std::pair<std::string,double>("190404_Cd_20kV_850V_46V_th250mV",70));
		S2_finish_time.insert(std::pair<std::string,double>("190404_Cd_18kV_850V_46V_th230mV",70));
		S2_finish_time.insert(std::pair<std::string,double>("190404_Cd_16kV_850V_46V_th210mV",70));

		areas_to_draw.push_back(experiment_area());

		//areas_to_draw.back().experiments.push_back("190404_Cd_20kV_850V_46V_th250mV_0");

		areas_to_draw.back().runs.push_pair(1, 1);

		//areas_to_draw.back().channels.push_pair(38, 38);

		//areas_to_draw.back().channels.push_pair(12, 12);
		//areas_to_draw.back().channels.push_pair(36, 36);
		//areas_to_draw.back().channels.push_pair(38, 38);
		//areas_to_draw.back().channels.push_pair(44, 44);
		areas_to_draw.back().channels.push_pair(32, 44); //13
		areas_to_draw.back().channels.push_pair(48, 59); //12 =>25 channels

		areas_to_draw.back().sub_runs.push_pair(0, 9);

		exp_area.runs.push_pair(0, 9999);
		exp_area.channels.push_pair(8, 16);
		//exp_area.channels.push_pair(GEM_CH_, GEM_CH_);

		exp_area.channels.push_pair(32, 44); //13
		exp_area.channels.push_pair(48, 59); //12 =>25 channels

		exp_area.sub_runs.push_pair(0, 999); //subruns_per_file-1);

		exp_area.experiments.push_back("190404_Cd_20kV_850V_46V_th250mV_0");
		exp_area.experiments.push_back("190404_Cd_20kV_850V_46V_th250mV");
		exp_area.experiments.push_back("190404_Cd_18kV_850V_46V_th230mV");
		exp_area.experiments.push_back("190404_Cd_16kV_850V_46V_th210mV");
	}
};
