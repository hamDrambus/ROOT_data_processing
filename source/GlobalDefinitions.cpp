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

//int get_next_index(std::vector<int> areas, int curr_index)
//{
//	bool even= true;
//	int L = -1, R = -1;
//	for (std::vector<int>::iterator sub = areas.begin(); sub != areas.end(); (sub==areas.end()?sub: sub++), even = !even){
//		if (even)
//			L = *sub;
//		else {
//			R = *sub;
//			if (!(curr_index< L) && (curr_index< R))
//				return curr_index + 1;
//			if (curr_index == R) {
//				if (++sub != areas.end()){ //go to the next area, but check it is valid
//					int out =*sub;
//					if ((sub++) != areas.end())
//						return out;
//				}
//			}
//		}
//	}
//	return -1;
//}

int get_order_index_by_index(int ind, std::vector<int> &areas)
{
	bool even = true;
	int l = -1, r = -1;
	int out = -1, N = 0;
	for (auto i = areas.begin(); i != areas.end(); i++, even = !even){
		if (even)
			l = *i;
		else {
			r = *i;
			if (!(ind<l) && !(ind>r)){
				out = N + (ind - l);
				break;
			}
			N += r - l + 1;
		}
	}
	return out;
}

namespace ParameterPile
{
	area_vector::area_vector(void) : _is_valid(false), _last_returned_index(-1)
	{
		_last_returned_index_left = _vec.end();
	}
	int area_vector::get_order_index_by_index(int ind)
	{
		if (!_is_valid)
			return -1;
		bool even = true;
		int l = -1, r = -1, N = 0;
		for (auto i = _vec.begin(); i != _vec.end(); i++, even = !even){
			if (even)
				l = *i;
			else {
				r = *i;
				if (!(ind<l) && !(ind>r))
					return  N + (ind - l);
				N += r - l + 1;
			}
		}
		return -1;
	}
	int area_vector::get_index_by_order_index(int ind)
	{
		if (!_is_valid)
			return -1;
		bool even = true;
		int l = -1, r = -1, N = 0;
		for (auto i = _vec.begin(); i != _vec.end(); i++, even = !even){
			if (even)
				l = *i;
			else {
				r = *i;
				N += r - l + 1;
				if (N>=ind)
					return r - N + ind + 1;
			}
		}
		return -1;
	}

	int area_vector::get_next_index(void) //for running through all indices
	{
		if (!_is_valid)
			return -1;
		if ((_last_returned_index < 0) || _last_returned_index_left == _vec.end()){
			_last_returned_index_left = _vec.begin();
			return _last_returned_index = *_last_returned_index_left;
		}
		std::vector<int>::iterator right = _last_returned_index_left+1;
		_last_returned_index++;
		if (_last_returned_index>*right){
			_last_returned_index_left = right++;
			if (_last_returned_index_left == _vec.end())
				return _last_returned_index = -1;
			return _last_returned_index = *_last_returned_index_left;
		}
		return _last_returned_index;
	}

	int area_vector::get_next_index(int after)
	{
		if (!_is_valid)
			return -1;
		if (-1 == after)
			return *_vec.begin();
		bool even = true;
		int L = -1, R = -1;
		for (auto sub = _vec.begin(); sub != _vec.end(); (sub == _vec.end() ? sub : sub++), even = !even) {
			if (even)
				L = *sub;
			else {
				R = *sub;
				if (!(after< L) && (after < R))
					return after + 1;
				if (after == R) {
					if (++sub != _vec.end()){ //go to the next area, but check it is valid
						int out = *sub;
						if ((sub++) != _vec.end())
							return out;
					}
				}
			}
		}
		return -1;
	}

	void area_vector::push_pair(int left, int right)
	{
		if (_vec.empty())
			_is_valid = true;
		if (left > right) {
			_vec.push_back(right);
			_vec.push_back(left);
		} else {
			_vec.push_back(right);
			_vec.push_back(left);
		}
		clear();
	}

	void area_vector::push(int val)
	{
		if (!_vec.empty())
			_is_valid = !_is_valid;
		_vec.push_back(val);
		clear();
	}

	bool area_vector::contains(int index)
	{	return !(get_order_index_by_index(index)<0);}

	bool area_vector::empty(void)
	{	return _vec.empty();}
	std::vector<area_vector> area_vector::split_area(int N)
	{	
		std::vector <area_vector> out_;
		for (int h = 0; h < N; h++){
			out_.push_back(area_vector());
		}
		int N_runs = 0;
		bool even = true;
		int l, r;
		for (auto h = _vec.begin(); h != _vec.end(); ++h, even = !even){
			if (even){
				l = *h;
			} else {
				r = *h;
				N_runs += r - l + 1;
			}
		}
		int N_accum = 0;
		int curr_index = 0;
		even = true;
		double delta = N_runs / (double)N;
		for (auto h = _vec.begin(); h != _vec.end(); ++h, even = !even) {
			if (even) {
				l = *h;
			} else {
				r = *h;
				N_accum += r - l + 1;
				int new_l = l, new_r = r;
				while (N_accum >(int)((curr_index + 1)*delta)) {
					if ((N_accum - (int)((curr_index + 1)*delta)) > (int)(delta))
						new_r = new_l + (int)(delta)-1;
					else
						new_r = N_accum - (int)((curr_index + 1)*delta) + new_l - 1;
					out_[curr_index].push_pair(new_l, new_r);
					new_l = new_r + 1;
					curr_index++;
				}
				new_r = r;
				if (N_accum <= (int)((curr_index + 1)*delta) && (new_l <= new_r))
					out_[curr_index].push_pair(new_l, new_r);
			}
		}
		return out_;
	}
	area_vector area_vector::intersect(area_vector& with)
	{
		area_vector out;
		bool even_out = true;
		int left_out, right_out;
		bool finished = false;
		for (auto g = out_area.runs.begin(); g != out_area.runs.end(); ++g, even_out = !even_out){
			if (even_out){
				left_out = *g;
			} else {
				right_out = *g;
				bool even = true;
				int left, right;
				for (auto h = area.runs.begin(); h != area.runs.end(); ++h, even = !even){
					if (even){
						left = *h;
					} else {
						right = *h;
						if ((right_out > right) && (left_out <= right)){
							*g = right;
						}
						if ((left_out < left) && (right_out >= left)){
							*(g - 1) = left;
						}
						left_out = *(g - 1);
						right_out = *(g);
					}
				}
			}
		}
		return out_area;
	}
	void area_vector::clear() //clears _last_returned_index etc.
	{
		_last_returned_index = -1;
		_last_returned_index_left = _vec.end();
	}
	void area_vector::erase() //clears vector
	{
		_vec.clear();
		clear();
	}
	//void area_vector::refine (void); //[2,3][3,4] to [2,4]|OR| [4,5] [1,7] to 


	std::vector <experiment_area> areas_to_draw;//TODO
	std::string this_path;
	int subruns_per_file = 10;
	bool override_analysis = true;
	experiment_area exp_area;
	int threads_number = 6; //obv. must be >=1

	int filter_MPPC_n_points = 15;
	int filter_MPPC_order = 4;
	int filter_MPPC_n_iterations = 1;
	int filter_PMT_n_points = 50;
	int filter_PMT_order = 8;
	int filter_PMT_n_iterations = 1;

	std::pair<double, double> baseline_search_limits(-0.1, 0.1);// = std::pair<double, double>(-0.05, 0.05);
	int baseline_search_max_iterations=4;
	std::vector<double> baseline_approx_value;

	double S1_time = 32; //in ms

	double GEM_threshold_to_noise = 1.1;
	int GEM_N_of_averaging = 30; //=== N_trust
	double PMT_run_acceptance_threshold_to_noize = 2;//~ 2-3
	int PMT_N_of_averaging = 1; //=== N_trust. PMT signal is already smoothed by filter
	int PMT_N_peaks_acceptance = 1;//1 //condition is >1, not >=1
	double PMT_SArea_peaks_acceptance = 0.0; //V*ms //done - for every runs S_acceptance is obtained from S distribution histogram
	//TODO: figure out the appropriate

	int gnuplot_pad_size = 400;
	int gnuplot_max_size = 1600;
	int gnuplot_width = 900; //default for gnuplot is 640

	bool draw_required(ParameterPile::experiment_area what) //'what' must be a single run
	{
		for (auto i = areas_to_draw.begin(); i != areas_to_draw.end(); i++) {
			for (auto exp = (*i).experiments.begin(); exp != (*i).experiments.end(); exp++)
				if ((*exp) == what.experiments.back()){
					bool even = true;
					int l = -1, r = -1;
					for (auto runs = (*i).runs.begin(); runs != (*i).runs.end();runs++,even=!even) {
						if (even)
							l = *runs;
						else {
							r = *runs;
							if (!(what.runs.back() < l) && !(what.runs.back() > r)){
								bool even_ch = true;
								int l_ch = -1, r_ch = -1;
								for (auto channels = (*i).channels.begin(); channels != (*i).channels.end(); channels++,even_ch=!even_ch) {
									if (even_ch)
										l_ch = *channels;
									else {
										r_ch = *channels;
										if (!(what.channels.back() < l_ch) && !(what.channels.back() > r_ch)){
											bool even_sub = true;
											int l_sub = -1, r_sub = -1;
											for (auto sub = (*i).sub_runs.begin(); sub != (*i).sub_runs.end(); sub++, even_sub = !even_sub) {
												if (even_sub)
													l_sub = *sub;
												else {
													r_sub = *sub;
													if (!(what.sub_runs.back() < l_sub) && !(what.sub_runs.back() > r_sub)){
														return true;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
		}
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
		areas_to_draw.back().runs.push_back(0);
		areas_to_draw.back().runs.push_back(0);
		areas_to_draw.back().channels.push_back(0);
		areas_to_draw.back().channels.push_back(0);
		areas_to_draw.back().channels.push_back(2);
		areas_to_draw.back().channels.push_back(2);
		areas_to_draw.back().channels.push_back(40);
		areas_to_draw.back().channels.push_back(41);
		areas_to_draw.back().sub_runs.push_back(1);
		areas_to_draw.back().sub_runs.push_back(1);

		exp_area.channels.push_back(0);
		exp_area.channels.push_back(0);
		exp_area.channels.push_back(2);
		exp_area.channels.push_back(2);
		/*exp_area.channels.push_back(40);
		exp_area.channels.push_back(41);*/
		
		exp_area.runs.push_back(2000);
		exp_area.runs.push_back(5000);
		
		exp_area.sub_runs.push_back(0);
		exp_area.sub_runs.push_back(9);//subruns_per_file-1);

		exp_area.experiments.push_back("4_thmV");
		exp_area.experiments.push_back("5_thmV");
		exp_area.experiments.push_back("6_thmV");
		exp_area.experiments.push_back("7_thmV");
		exp_area.experiments.push_back("8_thmV");
		exp_area.experiments.push_back("9_thmV");
		exp_area.experiments.push_back("10_thmV");
		exp_area.experiments.push_back("10_thmV_recalib");
		exp_area.experiments.push_back("12_thmV");
		exp_area.experiments.push_back("14_thmV");
		exp_area.experiments.push_back("16_thmV");
		exp_area.experiments.push_back("18_thmV");
		exp_area.experiments.push_back("20_thmV");

		int ch = exp_area.channels[0];
		for (ch = exp_area.channels[0]; !(ch < 0); ch = get_next_index(exp_area.channels, ch))
			baseline_approx_value.push_back(0);
		for (int ch = 0; ch < 2; ch++){
			int ind = get_order_index_by_index(ch, exp_area.channels);
			if (!(ind < 0)){
				baseline_approx_value[ind] = -0.4;//PMT
			}
		}
		ch = get_order_index_by_index(2, exp_area.channels);
		if (!(ch < 0)){
			baseline_approx_value[ch] = -0.45;//GEM
		}
		for (ch = 32; ch < 64; ch++){
			int ind = get_order_index_by_index(ch, exp_area.channels);
			if (!(ind < 0)){
				baseline_approx_value[ind] = 0;//MPPC
			}
		}

	}
};
