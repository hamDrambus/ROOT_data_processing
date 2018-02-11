#include "GlobalParameters.h"
#include "MTAnalysisManager.h"

//done TODO: 0 add constructor of experiment_area point
//done TODO: 1 clean up cutoff algorithms (comment them?)
//done (not tested) TODO: 2 add peak amplitude and refine peak finding operations
//TODO: 3 global: analyse MPPCs :
	// done 1) fix spread peaks function
	// done 2) improve find peak function: thresh1 for finding peak as currenly + 
	//		thresh2 for finding actual peak edges
	// 3) Tweak filter's + ROOT background function parameters
	// 4) It's time to analyse MPPC for low fields:
	//		1) find baseline by median/integral for first 32 ms
	//		2) find noise for the first 32 ms
	//		3) Find peaks for first 32 ms
	//		4) improve baseline omitting peaks
	//		5) using noise*Factor find all peaks. use ROOT's baseline for that? Then get baseline without peaks by avr over dt. update peaks.
	//			should for for small baseline deviation
	//		6) determine time windows for S1 and S2 signals using peaks (possible peak spread)
	//			For stability use limitation for that windows
	//		7) Tot. histograms of peaks outside S2 and S1 and inside S2
	//		8) Average summed area of peaks inside S2 per each MPPC
	// 5) Average field with not const baseline : 
	//		Use time step: find baseline (median/integral), exclude peaks, update baseline
	//		Same as for low fields: {S1,S2,noise} peaks ->histograms
	//		+ shaped baseline as function of peaks! determine non linearities
//TODO: 4 add description of programm structure as well as algorithms

int main(int argc, char *argv[])
{
	ParameterPile::Init_globals();
	/*std::ofstream str;
	str.open("test_bin.dat", std::ios_base::trunc | std::ios_base::binary);
	std::size_t sz = sizeof(std::size_t);
	str.write((char*)&sz, sizeof(std::size_t));
	str.close();*/
	int n_par = 0;
	char **f = NULL;
	TApplication app("test_app",&n_par,f);
	/*TCanvas* c1 = new TCanvas("test", "test_title", 800, 500);
	TF1 *func = new TF1("test_func", "sin(x)+3*x", 0, 10);
	func->Draw();*/
	
#ifdef _USE_TIME_STATISTICS
	auto timer_start = std::chrono::high_resolution_clock::now();
#endif
	if (1 >= ParameterPile::threads_number) {
		AnalysisManager man(ParameterPile::exp_area);
		man.processAllExperiments();
	} else {
		MTAnalysisManager man(ParameterPile::exp_area);
		man.processAllExperiments();
	}
	std::cout<<std::endl<<"========================================="<<std::endl << "Finished" << std::endl;
#ifdef _USE_TIME_STATISTICS
	auto timer_finish = std::chrono::high_resolution_clock::now();
	std::cout << "Time elapsed [s]: " << std::chrono::duration_cast<std::chrono::seconds>(timer_finish - timer_start).count();
#endif

	app.Run();
	//app->Delete();
	return 0;
}
