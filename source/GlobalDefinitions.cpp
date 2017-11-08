#include <stdio.h>
#include <direct.h>
#include "GlobalDefinitions.h"
#include "AnalysisManager.h"

#include <windows.h>

void open_output_file(std::string name, std::ofstream &str)
{
	std::string folder = name;
	while ((folder.back() != '\\') &&(folder.back()!='/') &&!folder.empty())
		folder.pop_back();
	if (!folder.empty())
		folder.pop_back();
	if (!folder.empty()){
		DWORD ftyp = GetFileAttributesA(folder.c_str());
		if (!(ftyp & FILE_ATTRIBUTE_DIRECTORY) || ftyp == INVALID_FILE_ATTRIBUTES){
				int code = system(("mkdir \"" + folder + "\"").c_str());
				if (code)
					std::cout << "mkdir error: " << GetLastError() << std::endl;
			}
	}
	str.open(name, std::ios_base::trunc);
	if (!str.is_open())
		std::cout << "Failed to open \"" << name << "\"" << std::endl;
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
		open_output_file("temp_gnuplot_files\\" + mod_name, file);
		std::cout << "file " << "temp_gnuplot_files\\" + mod_name << ".is_open() " << file.is_open() << std::endl;
		if (!file.is_open())
			std::cout << GetLastError() << std::endl;
		for (int h = 0; h < xs.size(); h++)
			file << xs[h] << '\t' << ys[h] << std::endl;
		file.close();
		open_output_file("temp_gnuplot_files\\script.sc", file);
		file << "plot '" << ParameterPile::this_path + "\\temp_gnuplot_files\\" + mod_name << "' u 1:2" << std::endl;
		file << "pause -1";
		file.close();
		system(("start \"\" \"%GNUPLOT%\\gnuplot.exe\" -c \"" + ParameterPile::this_path + "\\temp_gnuplot_files\\script.sc\"").c_str());
		//std::cout << "Gnuplot is not supported yet" << std::endl;
	}
}

namespace ParameterPile
{
	std::vector <experiment_area> areas_to_draw;
	std::string this_path;
	int subruns_per_file = 10;
	bool override_analysis = true;
	experiment_area exp_area;
	int threads_number = 1; //obv. must be >=1

	int filter_MPPC_n_points = 15;
	int filter_MPPC_order = 4;
	int filter_MPPC_n_iterations = 1;
	int filter_PMT_n_points = 50;
	int filter_PMT_order = 8;
	int filter_PMT_n_iterations = 1;

	std::pair<double, double> baseline_search_limits(-0.1, 0.1);// = std::pair<double, double>(-0.05, 0.05);
	int baseline_search_max_iterations=4;
	std::vector<double> baseline_approx_value;

	//these values are approximate, especially S2
	double S1_start_time = 32; //in ms
	double S1_finish_time = 35; //in ms
	double S2_start_time = 45; //in ms
	double S2_finish_time = 145; //in ms

	double GEM_threshold_to_noise = 1.1;
	int GEM_N_of_averaging = 30; //=== N_trust
	
	double PMT_run_acceptance_threshold_to_noize = 2;//~ 2-3
	int PMT_N_of_averaging = 1; //=== N_trust. PMT signal is already smoothed by filter
	int PMT_N_peaks_acceptance = 1;//1 //condition is >1, not >=1
	double PMT_SArea_peaks_acceptance = 0.0; //V*ms 
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
	double PMT_right_cutoff_from_RMS = 4;
	double PMT_left_cutoff_from_RMS = 2;

	int Max_iteration_N = 1;

	int gnuplot_pad_size = 400;
	int gnuplot_max_size = 1600;
	int gnuplot_width = 900; //default for gnuplot is 640

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
		this_path = _getcwd(path, FILENAME_MAX);

		//clear the GEM output file
		std::ofstream file;
		open_output_file(std::string(OUTPUT_DIR) + OUTPUT_GEMS, file);//creates folder AND trucates the file
		file << "Experiment\tIntegral[V*us]\tt_from\tt_to\tN_valid\tN\tS_cutoff" << std::endl;
		file.close();

		TThread::Initialize();
		
		areas_to_draw.push_back(experiment_area());
		areas_to_draw.back().experiments.push_back("4_thmV");
		areas_to_draw.back().experiments.push_back("5_thmV");
		areas_to_draw.back().experiments.push_back("6_thmV");
		areas_to_draw.back().experiments.push_back("7_thmV");
		areas_to_draw.back().experiments.push_back("8_thmV");
		areas_to_draw.back().experiments.push_back("9_thmV");
		areas_to_draw.back().experiments.push_back("10_thmV");
		areas_to_draw.back().experiments.push_back("10_thmV_recalib");
		areas_to_draw.back().experiments.push_back("12_thmV");
		areas_to_draw.back().experiments.push_back("14_thmV");
		areas_to_draw.back().experiments.push_back("16_thmV");
		areas_to_draw.back().experiments.push_back("18_thmV");
		areas_to_draw.back().experiments.push_back("20_thmV");
		areas_to_draw.back().runs.push_back(3396);
		areas_to_draw.back().runs.push_back(3400);
		//areas_to_draw.back().channels.push_back(0);
		//areas_to_draw.back().channels.push_back(0);
		//areas_to_draw.back().channels.push_back(2);
		//areas_to_draw.back().channels.push_back(2);
		areas_to_draw.back().channels.push_back(40);
		areas_to_draw.back().channels.push_back(41);
		areas_to_draw.back().sub_runs.push_back(1);
		areas_to_draw.back().sub_runs.push_back(1);

		exp_area.channels.push_back(0);
		exp_area.channels.push_back(0);
		exp_area.channels.push_back(2);
		exp_area.channels.push_back(2);
		exp_area.channels.push_back(41);
		exp_area.channels.push_back(41);
		
		exp_area.runs.push_back(2000);
		exp_area.runs.push_back(5000);
		
		exp_area.sub_runs.push_back(1);
		exp_area.sub_runs.push_back(1);//subruns_per_file-1);

		exp_area.experiments.push_back("4_thmV");
		/*exp_area.experiments.push_back("5_thmV");
		exp_area.experiments.push_back("6_thmV");
		exp_area.experiments.push_back("7_thmV");
		exp_area.experiments.push_back("8_thmV");
		exp_area.experiments.push_back("9_thmV");
		exp_area.experiments.push_back("10_thmV");
		exp_area.experiments.push_back("10_thmV_recalib");
		exp_area.experiments.push_back("12_thmV");
		exp_area.experiments.push_back("14_thmV");
		exp_area.experiments.push_back("16_thmV");
		exp_area.experiments.push_back("18_thmV");*/
		//exp_area.experiments.push_back("20_thmV");

		//TODO: ? get rid of baseline_approx? 
		for (int ch = exp_area.channels.get_next_index(); !(ch < 0); ch = exp_area.channels.get_next_index()){
			if ((ch < 2) && (ch >= 0)){ //PMT
				baseline_approx_value.push_back(-0.4);
				continue;
			}
			if (ch == 2){ //GEM
				baseline_approx_value.push_back(-0.45);
				continue;
			}
			if ((ch < 64) && (ch >= 32)){ //MPPC
				baseline_approx_value.push_back(0);
				continue;
			}
			baseline_approx_value.push_back(0);
		}
	}
};
