#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <sstream>

#include <TThread.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TVector.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/Point2D.h>
#include <windows.h>

#undef max
#undef min

#define DATA_PREFIX std::string("../../../Data/")
#define DATA_NAME_FORMAT "^run_\d+__ch_\d+\.dat$"
#define DATA_EXPERIMENT_FORMAT "^event_x-ray_\d{1,2}_.*$"

#define RUNS_CACHE_PATH "output/processed_runs.txt"
#define ALL_RUNS_PATH "output/processed_experiments.txt"
#define ALL_EXPERIMENTS_PATH "output/processed_summary.txt"

#define DATA_TIME_CONSTANT 1.6e-2
//^in microseconds
#define DATA_VOLTAGE_CHANNELS 4095
#define DATA_VOLTAGE_AMPLITUDE 2.0
//in volts
#define DATA_VOLTAGE_OF_ZERO_CHANNEL (-1.0)
//in volts
#define OUTPUT_DIR "results\\"
#define OUTPUT_GEMS "GEMs.txt"
#define _TEMP_CODE
#define DVECTOR std::vector<double>
#define DITERATOR std::vector<double>::iterator

void open_output_file(std::string name, std::ofstream &str);
int get_next_index(std::vector<int> areas, int curr_index);
int get_order_index_by_index(int ind, std::vector<int> &areas);

struct peak
{
	double left;
	double right;
	double S;
};

namespace ParameterPile
{
	enum DrawEngine { Gnuplot, ROOT };

	class area_vector
	{
	protected:
		bool _is_valid;
		std::vector<int> _vec;
		int _last_returned_index;
		std::vector<int>::iterator _last_returned_index_left;
	public:
		area_vector(void);
		int get_order_index_by_index(int ind);
		int get_index_by_order_index(int ind);
		int get_next_index(void); //for running through all indices
		int get_next_index(int after);
		void push_pair(int left, int right);
		void push(int val);
		bool contains(int index);
		bool empty(void);
		std::vector<area_vector> split_area(int N);
		area_vector intersect (area_vector& with);
		void clear(); //clears _last_returned_index etc.
		void erase(); //clears vector
		//void refine (void); //[2,3][3,4] to [2,4]|OR| [4,5] [1,7] to 
	};

	class experiment_area //TODO - make analysis via this class. //->NextFile?
	{
	public:
		enum Type {Area, Point};
	protected:
		Type _type;
	public:
		experiment_area(Type type = Type::Area);

		std::vector<std::string> experiments;
		std::vector<int> runs; //contains pairs [from, to]
		std::vector<int> channels; //contains pairs [from, to]
		std::vector<int> sub_runs; //contains pairs [from, to]
		
	};

	bool draw_required(ParameterPile::experiment_area what);
	extern std::vector <experiment_area> areas_to_draw;//TODO
	extern std::string this_path;
	extern int subruns_per_file;
	extern bool override_analysis;
	extern experiment_area exp_area;
	extern int threads_number;

	extern int filter_MPPC_n_points;
	extern int filter_MPPC_order;
	extern int filter_MPPC_n_iterations;
	extern int filter_PMT_n_points;
	extern int filter_PMT_order;
	extern int filter_PMT_n_iterations;

	extern std::pair<double, double> baseline_search_limits;//depricated
	extern int baseline_search_max_iterations; //depr
	extern std::vector<double> baseline_approx_value;

	extern double S1_time; //in ms

	extern double GEM_threshold_to_noise;
	extern int GEM_N_of_averaging; //=== N_trust

	extern double PMT_run_acceptance_threshold_to_noize;
	extern int PMT_N_of_averaging; //=== N_trust
	extern int PMT_N_peaks_acceptance;
	extern double PMT_SArea_peaks_acceptance; //V*ms

	extern int gnuplot_pad_size;
	extern int gnuplot_max_size;
	extern int gnuplot_width;

	void Init_globals(void);
};
void DrawFileData(std::string name, std::vector<double> xs, std::vector<double> ys, ParameterPile::DrawEngine de = ParameterPile::DrawEngine::ROOT);

#endif