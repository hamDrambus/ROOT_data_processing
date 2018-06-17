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
	std::string folder = name;
	while ((folder.back() != '\\') &&(folder.back()!='/') &&!folder.empty())
		folder.pop_back();
	if (!folder.empty())
		folder.pop_back();
#if defined(__WIN32__)
	if (!folder.empty()){
		DWORD ftyp = GetFileAttributesA(folder.c_str());
		if (!(ftyp & FILE_ATTRIBUTE_DIRECTORY) || ftyp == INVALID_FILE_ATTRIBUTES){
				int code = system(("mkdir \"" + folder + "\"").c_str());
				if (code)
					std::cout << "mkdir error: " << GetLastError() << std::endl;
			}
	}
#else
	struct stat st;
	stat(folder.c_str(),&st);
	if(!S_ISDIR(st.st_mode)){
		int code = system(("mkdir -p \"" + folder + "\"").c_str());
			if (code)
				std::cout << "mkdir error: " << code << std::endl;
	}
#endif //_WIN32__
	str.open(name.c_str(), std::ios_base::trunc);
	if (!str.is_open()){
		std::cout << "Failed to open \"" << name << "\"" << std::endl;
	}
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
		file << "plot '" << ParameterPile::this_path + "/temp_gnuplot_files/" + mod_name << "' u 1:2" << std::endl;
		file << "pause -1";
		file.close();
		INVOKE_GNUPLOT(ParameterPile::this_path + "/temp_gnuplot_files/script.sc");
	}
}

namespace ParameterPile
{
	STD_CONT <experiment_area> areas_to_draw;
	std::string this_path;
	int subruns_per_file = 10;
	//bool override_analysis = true;
	experiment_area exp_area;
	int threads_number = 3; //obv. must be >=1

	double dt_quant = 0.1; //us

	int filter_MPPC_n_points = 15;
	int filter_MPPC_order = 4;
	int filter_MPPC_n_iterations = 0;
	std::map<int,int> filter_PMT_n_points;
	std::map<int,int> filter_PMT_order;
	std::map<int,int> filter_PMT_n_iterations;

	//these values are approximate, especially S2
	double S1_start_time = 27; //in us
	double S1_finish_time = 35; //in us
	std::map<std::string,double> S2_start_time; //in us
	std::map<std::string,double> S2_finish_time; //in us

	double GEM_threshold_to_noise = 1.1;
	int GEM_N_of_averaging = 30; //=== N_trust
	
	double PMT_run_acceptance_threshold_to_noize = 4;//~ 2-3
	std::map<int,double> PMT_maximum_thresh;
	std::map<int,double> PMT_minimum_thresh;
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
	//TODO: add illustrative pictures to the project for algotithms
	double PMT_right_cutoff_from_RMS = 4.5;//4.5
	double PMT_left_cutoff_from_RMS = 2.5;//2
	std::map<int,bool> PMT_use_average;

	double MPPC_peaks_smoothing_time = 5; //us
	int MPPC_N_trust = 1;
	int MPPC_double_I_N_trust = 1;
	double MPPC_ROOTs_bl_from_max_left = 20;
	double MPPC_ROOTs_bl_from_max_right = 40;
	double MPPC_ROOTs_bl_left_offset = 5; //for baseline's baseline
	double MPPC_ROOTs_bl_right_offset = 12; //for baseline's baseline
	double MPPC_threshold_to_noise = 5.5;
	double MPPC_minimum_peak_A = 0.0128; //
	double MPPC_maximum_peak_A = 0.0128; //

	int Max_iteration_N = 1;

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

		//clear the GEM output file
#ifdef _PROCESS_GEMS
		std::ofstream file;
		open_output_file(std::string(OUTPUT_DIR) + OUTPUT_GEMS + "/" + OUTPUT_GEMS + ".txt", file);//creates folder AND trucates the file
		file << "Experiment\tIntegral[V*us]\tt_from\tt_to\tN_valid\tN\tS_cutoff" << std::endl;
		file.close();
#endif

		TThread::Initialize();
		
		PMT_SArea_peaks_acceptance = 0.2; //V*us
		subruns_per_file = 10;
		MPPC_minimum_peak_A = 0.0128; //
		MPPC_maximum_peak_A = 0.0128; //

		filter_PMT_n_points.insert		(std::pair<int,int>(0,50));
		filter_PMT_order.insert			(std::pair<int,int>(0,8));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(0,1));
		filter_PMT_n_points.insert		(std::pair<int,int>(1,50));
		filter_PMT_order.insert			(std::pair<int,int>(1,8));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(1,1));
		filter_PMT_n_points.insert		(std::pair<int,int>(8,50));
		filter_PMT_order.insert			(std::pair<int,int>(8,4));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(8,1));
		filter_PMT_n_points.insert		(std::pair<int,int>(12,50));
		filter_PMT_order.insert			(std::pair<int,int>(12,4));
		filter_PMT_n_iterations.insert	(std::pair<int,int>(12,1));

		PMT_minimum_thresh.insert 	(std::pair<int,double>(0,0.012));
		PMT_maximum_thresh.insert	(std::pair<int,double>(0,0.012));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(1,0.024));
		PMT_maximum_thresh.insert	(std::pair<int,double>(1,0.024));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(8,0.00046));
		PMT_maximum_thresh.insert	(std::pair<int,double>(8,0.00046));
		PMT_minimum_thresh.insert 	(std::pair<int,double>(12,0.00038));
		PMT_maximum_thresh.insert	(std::pair<int,double>(12,0.00038));

		PMT_use_average.insert(std::pair<int,bool>(0,true));
		PMT_use_average.insert(std::pair<int,bool>(1,true));


		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_4_thmV",80));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_5_thmV",78));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_6_thmV",70));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_7_thmV",63));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_8_thmV",55));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_9_thmV",47));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_10_thmV_recalib",44));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_10_thmV",44));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_12_thmV",44));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_14_thmV",44));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_16_thmV",44));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_18_thmV",44));
		S2_start_time.insert(std::pair<std::string,double>("event_x-ray_20_thmV",44));

		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_4_thmV",140));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_5_thmV",120));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_6_thmV",110));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_7_thmV",110));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_8_thmV",105));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_9_thmV",105));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_10_thmV_recalib",105));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_10_thmV",90));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_12_thmV",90));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_14_thmV",90));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_16_thmV",90));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_18_thmV",90));
		S2_finish_time.insert(std::pair<std::string,double>("event_x-ray_20_thmV",90));

		areas_to_draw.push_back(experiment_area());
		
		areas_to_draw.back().experiments.push_back("event_x-ray_4_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_5_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_6_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_7_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_8_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_9_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_10_thmV_recalib");
		areas_to_draw.back().experiments.push_back("event_x-ray_10_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_12_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_14_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_16_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_18_thmV");
		areas_to_draw.back().experiments.push_back("event_x-ray_20_thmV");

		areas_to_draw.back().runs.push_pair(0, 0);
		areas_to_draw.back().channels.push_pair(0, 1);
		areas_to_draw.back().channels.push_pair(2, 2);
		areas_to_draw.back().channels.push_pair(8, 8);
		areas_to_draw.back().channels.push_pair(12, 12);

		areas_to_draw.back().channels.push_pair(37, 38);
		areas_to_draw.back().channels.push_pair(44, 44);

		//areas_to_draw.back().channels.push_pair(32, 44); //13
		//areas_to_draw.back().channels.push_pair(48, 56); //9
		//areas_to_draw.back().channels.push_pair(57, 59); //3 =>25

		areas_to_draw.back().sub_runs.push_pair(3, 5);

		exp_area.runs.push_pair(0, 9999);
		exp_area.channels.push_pair(0, 1);
		exp_area.channels.push_pair(2, 2);
		//exp_area.channels.push_pair(8, 8);
		//exp_area.channels.push_pair(12, 12);

		//exp_area.channels.push_pair(56, 56);

		//exp_area.channels.push_pair(32, 44); //13
		//exp_area.channels.push_pair(48, 55); //8
		//exp_area.channels.push_pair(57, 59); //3 =>24 channels

		exp_area.sub_runs.push_pair(0, 9); //subruns_per_file-1);

		//exp_area.experiments.push_back("event_x-ray_4_thmV");
		//exp_area.experiments.push_back("event_x-ray_5_thmV");
		//exp_area.experiments.push_back("event_x-ray_6_thmV");
		//exp_area.experiments.push_back("event_x-ray_7_thmV");
		//exp_area.experiments.push_back("event_x-ray_8_thmV");
		//exp_area.experiments.push_back("event_x-ray_9_thmV");
		//exp_area.experiments.push_back("event_x-ray_10_thmV_recalib");
		//exp_area.experiments.push_back("event_x-ray_10_thmV");
		//exp_area.experiments.push_back("event_x-ray_12_thmV");
		//exp_area.experiments.push_back("event_x-ray_14_thmV");
		exp_area.experiments.push_back("event_x-ray_16_thmV");
		exp_area.experiments.push_back("event_x-ray_18_thmV");
		exp_area.experiments.push_back("event_x-ray_20_thmV");
	}
};
