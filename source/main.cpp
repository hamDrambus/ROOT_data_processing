#include "GlobalDefinitions.h"
#include "MTAnalysisManager.h"

int main(int argc, char *argv[])
{
	ParameterPile::Init_globals();
	int n_par = 0;
	char **f = NULL;
	TApplication* app = new TApplication("test_app",&n_par,f);
	TCanvas* c1 = new TCanvas("test", "test_title", 800, 500);
	TF1 *func = new TF1("test_func", "sin(x)+5*x", 0, 10);
	func->Draw();
	MTAnalysisManager man(ParameterPile::exp_area);

	man.processAllExperiments();
	std::cout << "Finished" << std::endl;

	app->Run();
	delete app;
	return 0;
}
