#include "GlobalParameters.h"
#include "MTAnalysisManager.h"

//TODO: add description of programm structure as well as algorithms

//09.04.2019 UPD TODO: so much rework is required
//	DONE 0) Simplify (delete) fancy unnecessary algorithms. Thresholds are selected by hand anyways, and post processing is done separately
//	!!!1) ParameterPile--->gSettings
//	2) make area_vector to be thread-save
//	DONE 3) all analysis settings are in separate locations: Init190404()
//	4) support writing data in ROOT trees, not only in my custom files (too cumbersome hierarchy)
//	DONE 5) drawing manager grouping by channel, not in order of processing
//	6) DrawSignals mode for choosing thresholds and filters and such
//	*7) Qt (Python?) interactive signal displaying and on-the run selection of processing parameters (like oscilloscope)
//	DONE 8) MPPC-PMT separation is unnecessary.
//	DONE 9) Implement pre-selection (empty PMT signal) via external file with accepted indices.

int main(int argc, char *argv[])
{
	ParameterPile::Init_globals();
	/*std::ofstream str;
	str.open("test_bin.dat", std::ios_base::trunc | std::ios_base::binary);
	std::size_t sz = sizeof(std::size_t);
	str.write((char*)&sz, sizeof(std::size_t));
	str.close();*/
	std::cout << "This: " << ParameterPile::this_path << std::endl;
	
#ifdef _USE_TIME_STATISTICS
	auto timer_start = std::chrono::high_resolution_clock::now();
#endif
	if (1 >= ParameterPile::threads_number) {
		AnalysisManager man(ParameterPile::gManifest);
		man.processAllExperiments();
	} else {
		MTAnalysisManager man(ParameterPile::gManifest);
		man.processAllExperiments();
	}
	std::cout<<std::endl<<"========================================="<<std::endl << "Finished" << std::endl;
#ifdef _USE_TIME_STATISTICS
	auto timer_finish = std::chrono::high_resolution_clock::now();
	std::cout << "Time elapsed [s]: " << std::chrono::duration_cast<std::chrono::seconds>(timer_finish - timer_start).count();
#endif

	return 0;
}
