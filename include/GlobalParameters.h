#ifndef _GLOBAL_PARAMETERS_H
#define _GLOBAL_PARAMETERS_H

#include "GlobalDefinitions.h"
#include "ExperimentArea.h"

//unmodified algorithm parameters:
//find_background_v_raw(f_ys, ys.size(), 60,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false);
//optimized (with sparse) optimal algorithm:
//find_background_v_0(f_ys, ys.size(), 80, TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false, 2);

//for 62.5 MHz discretization frequency (9999 points per 160 us)
//#define ROOT_BL_CALL_V0 find_background_v_0(f_ys, ys.size(), 90,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,2);
//#define ROOT_BL_CALL_V1 find_background_v_0(f_ys, ys.size(), 60,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,2);
//#define ROOT_BL_CALL_V2 find_background_v_0(f_ys, ys.size(), 70,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,2);
//#define ROOT_BL_CALL_V3 find_background_v_raw(f_ys, ys.size(), 60,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false);
//#define ROOT_BL_CALL_V4 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);
//#define ROOT_BL_CALL_V5 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);
//#define ROOT_BL_CALL_V6 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);
//#define ROOT_BL_CALL_V7 spec->Background(f_ys, ys.size(), 20,	TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, false, TSpectrum::kBackSmoothing3, false);

//for 250MHz discretization frequency (40000 points per 160 us)
//#define ROOT_BL_CALL_V0 find_background_v_0(f_ys, ys.size(), 360,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,7);
//#define ROOT_BL_CALL_V1 find_background_v_0(f_ys, ys.size(), 360,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing5, false,6);
//#define ROOT_BL_CALL_V2 find_background_v_0(f_ys, ys.size(), 500,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,6);
//#define ROOT_BL_CALL_V3 find_background_v_0(f_ys, ys.size(), 360,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false,7);
//#define ROOT_BL_CALL_V4 find_background_v_0(f_ys, ys.size(), 360,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing5, false,4);
//#define ROOT_BL_CALL_V5 find_background_v_0(f_ys, ys.size(), 250,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing5, false,4);
//#define ROOT_BL_CALL_V6 find_background_v_0(f_ys, ys.size(), 500,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing5, false,4);
//#define ROOT_BL_CALL_V7 find_background_v_0(f_ys, ys.size(), 360,	TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing5, false,6);


namespace ParameterPile //TODO: rename to Settings or gSettings. Maybe move to class instead of namespace
{
	enum DrawEngine { Gnuplot, ROOT };

	extern analysis_manifest gManifest;
	//Following functions extend manifest, not erase already present data
	bool Init190404(analysis_manifest& manifest);
	bool Init190404_tests(analysis_manifest& manifest);
	bool Init180705(analysis_manifest& manifest);
	bool Init180705_tests(analysis_manifest& manifest);
	bool Init180830Xray(analysis_manifest& manifest);
	bool Init180830Xray_tests(analysis_manifest& manifest);
	bool Init191107(analysis_manifest& manifest);
	bool Init191107_tests(analysis_manifest& manifest);
	bool Init191128(analysis_manifest& manifest);
	bool Init191128_tests(analysis_manifest& manifest);
	bool Init200116(analysis_manifest& manifest);
	bool Init200116_tests(analysis_manifest& manifest);
	bool Init200213(analysis_manifest& manifest);
	bool Init200213_tests(analysis_manifest& manifest);
	bool Init170622Cd(analysis_manifest& manifest);
	bool Init170622Cd_tests(analysis_manifest& manifest);
	bool Init200910Pu(analysis_manifest& manifest);
	bool Init200910Pu_tests(analysis_manifest& manifest);
	bool Init190307Xray(analysis_manifest& manifest);
	bool Init190307Xray_tests(analysis_manifest& manifest);
	bool Init201015(analysis_manifest& manifest);
	bool Init201015_tests(analysis_manifest& manifest);
	bool Init201015Xray(analysis_manifest& manifest);
	bool Init201015Xray_tests(analysis_manifest& manifest);
	bool Init201217(analysis_manifest& manifest);
	bool Init201217_tests(analysis_manifest& manifest);
	bool Init210121(analysis_manifest& manifest);
	bool Init210121_tests(analysis_manifest& manifest);
	bool Init210128(analysis_manifest& manifest);
	bool Init210128_tests(analysis_manifest& manifest);
	bool Init210218(analysis_manifest& manifest);
	bool Init210218_tests(analysis_manifest& manifest);
	bool Init210302(analysis_manifest& manifest);
	bool Init210302_tests(analysis_manifest& manifest);
	bool Init210311(analysis_manifest& manifest);
	bool Init210311_tests(analysis_manifest& manifest);
	bool Init210316(analysis_manifest& manifest);
	bool Init210316_tests(analysis_manifest& manifest);
	bool Init210429(analysis_manifest& manifest);
	bool Init210429_tests(analysis_manifest& manifest);
	bool Init210513(analysis_manifest& manifest);
	bool Init210513_tests(analysis_manifest& manifest);

	bool read_accepted_events(std::string file, accepted_events<double> &info); //does not erase already present data in info

	extern std::string this_path;
	extern int threads_number;
	extern std::size_t max_pics_number;
	extern bool gnuplot_presits;
	extern bool quiet_mode;

	extern int Max_iteration_N;

	extern int gnuplot_pad_size;
	extern int gnuplot_max_size;
	extern int gnuplot_width;

	void Init_globals(void);
};

#endif
