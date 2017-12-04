#ifndef _GLOBAL_PARAMETERS_H
#define _GLOBAL_PARAMETERS_H

#include "GlobalDefinitions.h"
#include "ExperimentArea.h"
#ifdef _USE_TIME_STATISTICS
#include <chrono>
#endif

#define ROOT_BL_CALL_V0 spec->Background(f_ys, ys.size(), 60,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V2 spec->Background(f_ys, ys.size(), 25,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V3 spec->Background(f_ys, ys.size(), 60,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V4 spec->Background(f_ys, ys.size(), 50,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V5 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V6 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V7 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);
#define ROOT_BL_CALL_V8 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);

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

	extern int filter_MPPC_n_points;
	extern int filter_MPPC_order;
	extern int filter_MPPC_n_iterations;
	extern int filter_PMT_n_points;
	extern int filter_PMT_order;
	extern int filter_PMT_n_iterations;

	//all times below are approximate
	extern double S1_start_time; //in ms
	extern double S1_finish_time; //in ms
	extern double S2_start_time; //in ms
	extern double S2_finish_time; //in ms

	extern double GEM_threshold_to_noise;
	extern int GEM_N_of_averaging; //=== N_trust

	extern double PMT_run_acceptance_threshold_to_noize;
	extern int PMT_N_of_averaging; //=== N_trust
	extern int PMT_N_peaks_acceptance;
	extern double PMT_SArea_peaks_acceptance; //V*ms
	extern double PMT_min_fraction_above_cutoff;
	extern int PMT_min_statistics;
	extern double PMT_mean_above_cutoff_acceptance;
	extern double PMT_right_cutoff_from_RMS;
	extern double PMT_left_cutoff_from_RMS;

	extern double MPPC_peaks_smoothing_time;
	extern int MPPC_N_trust;
	extern double MPPC_ROOTs_bl_from_max_left;
	extern double MPPC_ROOTs_bl_from_max_right;
	extern double MPPC_ROOTs_bl_left_offset; //for baseline's baseline
	extern double MPPC_ROOTs_bl_right_offset; //for baseline's baseline

	extern int Max_iteration_N;

	extern int gnuplot_pad_size;
	extern int gnuplot_max_size;
	extern int gnuplot_width;

	void Init_globals(void);
};
void DrawFileData(std::string name, DVECTOR xs, DVECTOR ys, ParameterPile::DrawEngine de = ParameterPile::DrawEngine::ROOT);

#endif