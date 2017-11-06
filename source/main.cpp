#include "GlobalDefinitions.h"
#include "MTAnalysisManager.h"

//done TODO: 0 add constructor of experiment_area point
//done TODO: 1 clean up cutoff algorithms (comment them?)
//done (not tested) TODO: 2 add peak amplitude and refine peak finding operations
//TODO: 3 global: analyse MPPCs
//TODO: 4 add description of programm structure as well as algorithms

int main(int argc, char *argv[])
{
	ParameterPile::Init_globals();
	int n_par = 0;
	char **f = NULL;
	TApplication* app = new TApplication("test_app",&n_par,f);
	TCanvas* c1 = new TCanvas("test", "test_title", 800, 500);
	TF1 *func = new TF1("test_func", "sin(x)+5*x", 0, 10);
	func->Draw();
	
	if (1 >= ParameterPile::threads_number) {
		AnalysisManager man(ParameterPile::exp_area);
		man.processAllExperiments();
	}
	else {
		MTAnalysisManager man(ParameterPile::exp_area);
		man.processAllExperiments();
	}
	std::cout << "Finished" << std::endl;

	app->Run();
	delete app;
	return 0;
}
