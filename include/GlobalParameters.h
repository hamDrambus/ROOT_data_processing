#ifndef _GLOBAL_PARAMETERS_H
#define _GLOBAL_PARAMETERS_H

#include "GlobalDefinitions.h"
#include "ExperimentArea.h"
#ifdef _USE_TIME_STATISTICS
#include <chrono>
#endif

//unmodified algorithm parameters:
//find_background_v_raw(f_ys, ys.size(), 60,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false);
//optimized (with sparse) optimal algorithm:
//find_background_v_0(f_ys, ys.size(), 80, TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false, 2);

//for 60MHz discretization frequency (9999 points per 160 us)
#define ROOT_BL_CALL_V0 find_background_v_0(f_ys, ys.size(), 90,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,2);
#define ROOT_BL_CALL_V1 find_background_v_0(f_ys, ys.size(), 60,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,2);
#define ROOT_BL_CALL_V2 find_background_v_0(f_ys, ys.size(), 70,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,2);
#define ROOT_BL_CALL_V3 find_background_v_raw(f_ys, ys.size(), 60,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V4 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V5 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V6 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V7 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);

//for 250MHz discretization frequency (40000 points per 160 us)
//#define ROOT_BL_CALL_V0 find_background_v_0(f_ys, ys.size(), 360,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,7);
//#define ROOT_BL_CALL_V1 find_background_v_0(f_ys, ys.size(), 360,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing5, false,6);
//#define ROOT_BL_CALL_V2 find_background_v_0(f_ys, ys.size(), 500,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,6);
//#define ROOT_BL_CALL_V3 find_background_v_0(f_ys, ys.size(), 360,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,7);
//#define ROOT_BL_CALL_V4 find_background_v_0(f_ys, ys.size(), 360,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing5, false,4);
//#define ROOT_BL_CALL_V5 find_background_v_0(f_ys, ys.size(), 250,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing5, false,4);
//#define ROOT_BL_CALL_V6 find_background_v_0(f_ys, ys.size(), 500,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing5, false,4);
//#define ROOT_BL_CALL_V7 find_background_v_0(f_ys, ys.size(), 360,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing5, false,6);


namespace ParameterPile
{
	enum DrawEngine { Gnuplot, ROOT };

	bool draw_required(ParameterPile::experiment_area what);

	extern STD_CONT <experiment_area> areas_to_draw;
	extern std::string this_path;
	extern int subruns_per_file;
	//extern bool override_analysis;
	extern experiment_area exp_area;
	extern int threads_number;

	extern double dt_quant;

	extern int filter_MPPC_n_points;
	extern int filter_MPPC_order;
	extern int filter_MPPC_n_iterations;
	extern std::map<int,int> filter_PMT_n_points;
	extern std::map<int,int> filter_PMT_order;
	extern std::map<int,int> filter_PMT_n_iterations;

	//all times below are approximate
	extern double S1_start_time; //in ms
	extern double S1_finish_time; //in ms
	extern std::map<std::string,double> S2_start_time; //in ms
	extern std::map<std::string,double> S2_finish_time; //in ms

	extern double GEM_threshold_to_noise;
	extern int GEM_N_of_averaging; //=== N_trust

	extern double PMT_run_acceptance_threshold_to_noize;
	extern std::map<int,double> PMT_maximum_thresh;
	extern std::map<int,double> PMT_minimum_thresh;
	extern std::map<int,double> PMT_thresh_edges;
	extern int PMT_N_of_averaging; //=== N_trust
	extern int PMT_N_peaks_acceptance;
	extern double PMT_SArea_peaks_acceptance; //V*ms
	extern double PMT_min_fraction_above_cutoff;
	extern int PMT_min_statistics;
	extern double PMT_mean_above_cutoff_acceptance;
	extern double PMT_right_cutoff_from_RMS;
	extern double PMT_left_cutoff_from_RMS;
	extern double PMT_ROOTs_bl_from_max_left;
	extern double PMT_ROOTs_bl_from_max_right;
	extern double PMT_ROOTs_bl_left_offset; //for baseline's baseline
	extern double PMT_ROOTs_bl_right_offset; //for baseline's baseline
	extern double PMT_ROOTs_bl_trim;

	extern area_vector ch_use_average;
	extern area_vector ch_integrate_S2;
	extern area_vector ch_use_curved_baseline;
	extern area_vector ch_inverse;

	extern double MPPC_peaks_smoothing_time;
	extern int MPPC_N_trust;
	extern int MPPC_double_I_N_trust;
	extern double MPPC_ROOTs_bl_from_max_left;
	extern double MPPC_ROOTs_bl_from_max_right;
	extern double MPPC_ROOTs_bl_left_offset; //for baseline's baseline
	extern double MPPC_ROOTs_bl_right_offset; //for baseline's baseline
	extern double MPPC_ROOTs_bl_trim;
	extern double MPPC_threshold_to_noise;
	extern double MPPC_maximum_peak_A;
	extern double MPPC_minimum_peak_A;

	extern int Max_iteration_N;

	extern int gnuplot_pad_size;
	extern int gnuplot_max_size;
	extern int gnuplot_width;

	void Init_globals(void);
};
void DrawFileData(std::string name, DVECTOR xs, DVECTOR ys, ParameterPile::DrawEngine de = ParameterPile::DrawEngine::ROOT);

#endif
