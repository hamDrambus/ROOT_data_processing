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
	int subruns_per_file = 1000;
	//bool override_analysis = true;
	experiment_area exp_area;
	int threads_number = 1; //obv. must be >=1

	int filter_MPPC_n_points = 15;
	int filter_MPPC_order = 4;
	int filter_MPPC_n_iterations = 1;
	int filter_PMT_n_points = 50;//50
	int filter_PMT_order = 8;//8
	int filter_PMT_n_iterations = 1;

	//TODO: depr
	//std::pair<double, double> baseline_search_limits(-0.1, 0.1);// = std::pair<double, double>(-0.05, 0.05);
	//int baseline_search_max_iterations=4;

	//these values are approximate, especially S2
	double S1_start_time = 31; //in ms
	double S1_finish_time = 42; //in ms
	std::map<std::string,double> S2_start_time; //in ms
	std::map<std::string,double> S2_finish_time; //in ms

	double GEM_threshold_to_noise = 1.1;
	int GEM_N_of_averaging = 30; //=== N_trust
	
	double PMT_run_acceptance_threshold_to_noize = 4;//~ 2-3
	double PMT_minimum_thresh = 0.050;
	double PMT_maximum_thresh = 0.050;
	double PMT1_minimum_thresh = 0.051;
	double PMT1_maximum_thresh = 0.051;
	int PMT_N_of_averaging = 1; //=== N_trust. PMT signal is already smoothed by filter
	int PMT_N_peaks_acceptance = 0;//1 //condition is >=1, not >1
	double PMT_SArea_peaks_acceptance = 0.00; //V*us
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
	int gnuplot_width = 1400; //default for gnuplot is 640

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
		
		PMT_use_average.insert(std::pair<int,bool>(0,true));
		PMT_use_average.insert(std::pair<int,bool>(1,true));

		S2_start_time.insert(std::pair<std::string,double>("7kV_SiPM_46V_xray_240Hz_PMT_750V",50));
		S2_start_time.insert(std::pair<std::string,double>("8kV_SiPM_46V_xray_240Hz_PMT_750V",42));
		S2_start_time.insert(std::pair<std::string,double>("9kV_SiPM_46V_xray_240Hz_",50));
		S2_start_time.insert(std::pair<std::string,double>("10kV_SiPM_46V_xray_240Hz",42));
		S2_start_time.insert(std::pair<std::string,double>("11kV_SiPM_46V_xray_240Hz",42));
		S2_start_time.insert(std::pair<std::string,double>("12kV_SiPM_46V_xray_240Hz",42));
		S2_start_time.insert(std::pair<std::string,double>("13kV_SiPM_46V_xray_240Hz",42));
		S2_start_time.insert(std::pair<std::string,double>("14kV_SiPM_46V_xray_240Hz",42));
		S2_start_time.insert(std::pair<std::string,double>("15kV_SiPM_46V_xray_240Hz_PMT_700V_6dB",42));

		S2_finish_time.insert(std::pair<std::string,double>("7kV_SiPM_46V_xray_240Hz_PMT_750V",90));
		S2_finish_time.insert(std::pair<std::string,double>("8kV_SiPM_46V_xray_240Hz_PMT_750V",88));
		S2_finish_time.insert(std::pair<std::string,double>("9kV_SiPM_46V_xray_240Hz_",88));
		S2_finish_time.insert(std::pair<std::string,double>("10kV_SiPM_46V_xray_240Hz",10));
		S2_finish_time.insert(std::pair<std::string,double>("11kV_SiPM_46V_xray_240Hz",88));
		S2_finish_time.insert(std::pair<std::string,double>("12kV_SiPM_46V_xray_240Hz",88));
		S2_finish_time.insert(std::pair<std::string,double>("13kV_SiPM_46V_xray_240Hz",88));
		S2_finish_time.insert(std::pair<std::string,double>("14kV_SiPM_46V_xray_240Hz",88));
		S2_finish_time.insert(std::pair<std::string,double>("15kV_SiPM_46V_xray_240Hz_PMT_700V_6dB",88));

		areas_to_draw.push_back(experiment_area());
		
		areas_to_draw.back().experiments.push_back("7kV_SiPM_46V_xray_240Hz_PMT_750V");
		/*areas_to_draw.back().experiments.push_back("8kV_SiPM_46V_xray_240Hz_PMT_750V");
		areas_to_draw.back().experiments.push_back("9kV_SiPM_46V_xray_240Hz");
		areas_to_draw.back().experiments.push_back("10kV_SiPM_46V_xray_240Hz");
		areas_to_draw.back().experiments.push_back("11kV_SiPM_46V_xray_240Hz");
		areas_to_draw.back().experiments.push_back("12kV_SiPM_46V_xray_240Hz");
		areas_to_draw.back().experiments.push_back("13kV_SiPM_46V_xray_240Hz");
		areas_to_draw.back().experiments.push_back("14kV_SiPM_46V_xray_240Hz");
		areas_to_draw.back().experiments.push_back("15kV_SiPM_46V_xray_240Hz_PMT_700V_6dB");*/

		areas_to_draw.back().runs.push_pair(0, 9999);
		areas_to_draw.back().channels.push_pair(0, 1);
		areas_to_draw.back().channels.push_pair(2, 2);

		areas_to_draw.back().channels.push_pair(32, 64);

		/*areas_to_draw.back().channels.push_pair(32, 44); //13
		areas_to_draw.back().channels.push_pair(48, 55); //8
		areas_to_draw.back().channels.push_pair(57, 59);*/ //3 =>24
		areas_to_draw.back().sub_runs.push_pair(1, 1);

		exp_area.runs.push_pair(67, 67);
		exp_area.channels.push_pair(0, 1);
		exp_area.channels.push_pair(2, 2);

		exp_area.channels.push_pair(32, 44); //13
		exp_area.channels.push_pair(48, 56); //9
		exp_area.channels.push_pair(57, 59); //3 =>25 channels
		exp_area.sub_runs.push_pair(0, 9); //subruns_per_file-1);

		exp_area.experiments.push_back("7kV_SiPM_46V_xray_240Hz_PMT_750V");
		//exp_area.experiments.push_back("8kV_SiPM_46V_xray_240Hz_PMT_750V");
		//exp_area.experiments.push_back("9kV_SiPM_46V_xray_240Hz");
		//exp_area.experiments.push_back("10kV_SiPM_46V_xray_240Hz");
		//exp_area.experiments.push_back("11kV_SiPM_46V_xray_240Hz");
		//exp_area.experiments.push_back("12kV_SiPM_46V_xray_240Hz");
		//exp_area.experiments.push_back("13kV_SiPM_46V_xray_240Hz");
		//exp_area.experiments.push_back("14kV_SiPM_46V_xray_240Hz");
		//exp_area.experiments.push_back("15kV_SiPM_46V_xray_240Hz_PMT_700V_6dB");
	}
};
