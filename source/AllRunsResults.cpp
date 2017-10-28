#include "AllRunsResults.h"

AllRunsResults::AllRunsResults(ParameterPile::experiment_area experiment)
{
	_exp = experiment;
	S_peaks_cutoff = ParameterPile::PMT_SArea_peaks_acceptance;
	N_peaks_cutoff = ParameterPile::PMT_N_peaks_acceptance;
}

void AllRunsResults::processAllRuns(std::vector<SingleRunResults> &single_results)
{
	_Ss.clear();
	_ns.clear();
	int n_valid_runs = 0;
	for (auto i = single_results.begin(); i != single_results.end(); ++i){
		if (i->isValid()){ //has meaning only at the secondary function call
			n_valid_runs++;
			_Ss.push_back((*i).PMT3_summed_peaks_area);
			_ns.push_back((*i).PMT3_n_peaks);
			if (_xs_GEM_sum.empty()){
				_xs_GEM_sum = (*i).xs_GEM;
				_ys_GEM_sum = (*i).ys_GEM;
			} else {
				for (auto ys1 = _ys_GEM_sum.begin(), ys2 = (*i).ys_GEM.begin();
					(ys1 != _ys_GEM_sum.end() && ys2 != (*i).ys_GEM.end()); ++ys1, ++ys2)
					*ys1 += *ys2;
			}
		}
	}
	if (!_Ss.empty() && _xs_GEM_sum.empty()) {//_xs_GEM is empty at the first run
		find_S_cutoff();
		find_S_cutoff_v2();

		DITERATOR S_max = std::max_element(_Ss.begin(), _Ss.end());
		TH1D *hist_S = new TH1D("PMT_S_peaks", "PMT_S_peaks", 60, 0, *S_max);
		/*TH1I *hist_n = new TH1I("PMT_N_peaks", "PMT_N_peaks", 30, 0, 30);*/

		/*double *NNs = new double[_ns.size()];
		double *SSs = new double[_Ss.size()];*/
		for (int _n = 0, _s = 0; (_n < _ns.size()) && (_s < _Ss.size()); ++_s, ++_n){
			/*NNs[_n] = _ns[_n];
			SSs[_s] = _Ss[_s];*/
			/*hist_n->Fill(_ns[_n]);*/
			hist_S->Fill(_Ss[_s]);
		}
		TCanvas *c1 = new TCanvas(("S_peaks_distribution " + _exp.experiments.back()).c_str(),
			("S_peaks_distribution " + _exp.experiments.back()).c_str());
		c1->cd();
		hist_S->Draw();
		c1->Update();
		TLine *cutoff = new TLine(S_peaks_cutoff, c1->GetUymin(), S_peaks_cutoff, c1->GetUymax());
		cutoff->SetLineColor(kRed);
		cutoff->Draw();
		c1->Update();
		/*TCanvas *c2 = new TCanvas(("n_peaks_distribution " + _exp.experiments.back()).c_str(),
			("n_peaks_distribution " + _exp.experiments.back()).c_str());
		c2->cd();
		hist_n->Draw();
		c2->Update();

		TCanvas *c3 = new TCanvas(("N-S Scatter " + _exp.experiments.back()).c_str(), ("N-S Scatter " + _exp.experiments.back()).c_str());
		TGraph *gr = new TGraph(std::min(_ns.size(), _Ss.size()), SSs, NNs);
		c3->cd();
		gr->Draw("ap");
		c3->Update();*/
		/*delete[] NNs;
		delete[] SSs;*/
	}
	if (!_xs_GEM_sum.empty() && (0!=n_valid_runs)){ //x_GEM is not empty only at the second processing (when proper cutoffs were applied)
		for (auto i = _ys_GEM_sum.begin(); i != _ys_GEM_sum.end(); ++i)
			*i /= n_valid_runs;
		ParameterPile::experiment_area area_ = _exp;
		area_.runs.back()=ParameterPile::areas_to_draw.back().runs.back();
		area_.sub_runs.back()=ParameterPile::areas_to_draw.back().sub_runs.back();
		area_.channels.back()=ParameterPile::areas_to_draw.back().channels.back();
		if (ParameterPile::draw_required(area_)){
			//xs_GEM.insert(xs_GEM.begin(), 0);
			//ys_GEM.insert(ys_GEM.begin(), 0);
			std::string exp = _exp.experiments.back();
			bool invalid = false;
			for (auto i = exp.begin(); i != exp.end(); invalid ? i = exp.begin() : ++i){
				invalid = false;
				if (*i == '_'){
					if (i != exp.begin())
						if (*(i - 1) == '\\')
							continue;
					exp.insert(i, '\\');
					invalid = true;
				}
			}
			Drawing *dr = graph_manager.GetDrawing("GEM_" + area_.experiments.back(), 0, ParameterPile::DrawEngine::Gnuplot);
			dr->AddToDraw(_xs_GEM_sum, _ys_GEM_sum, "GEM\\\_" + exp, "", 0);
			std::vector<double> GEM_int;
			SignalOperations::integrate(_xs_GEM_sum, _ys_GEM_sum, GEM_int);
			dr->AddToDraw(_xs_GEM_sum, GEM_int, "GEM\\\_I\\\_" + exp, "axes x1y2", 0);

			//find start GEM time 
			std::vector<double>::iterator x_start = _xs_GEM_sum.begin();
			find_GEM_start_time(_xs_GEM_sum, _ys_GEM_sum, x_start, ParameterPile::GEM_N_of_averaging, graph_manager);
			if (x_start != _xs_GEM_sum.end())
				dr->AddToDraw_vertical(*x_start, -1, 1, "lc rgb \"#FF0000\"", 0);
			//find finish GEM time
			std::vector<double>::iterator x_finish = _xs_GEM_sum.begin();
			double temp_y_max;
			DVECTOR temp_xs = _xs_GEM_sum, temp_ys_I = GEM_int;
			SignalOperations::apply_time_limits(temp_xs, temp_ys_I, *x_start, _xs_GEM_sum.back());
			SignalOperations::get_max(temp_xs, temp_ys_I, x_finish, temp_y_max, ParameterPile::GEM_N_of_averaging);
			if (x_finish != temp_xs.end())
				dr->AddToDraw_vertical(*x_finish, -1, 1, "lc rgb \"#FF0000\"", 0);
			std::cout << "Experiment " << area_.experiments.back() << " processed" << std::endl;
			std::cout << "# of runs " << single_results.size() << std::endl;
			//std::cout << "# of runs with empty PMT signal" << runs_no_PMT << std::endl;
			if (x_finish != temp_xs.end() && x_start != _xs_GEM_sum.end()){
				x_finish = SignalOperations::find_x_iterator_by_value(_xs_GEM_sum.begin(), _xs_GEM_sum.end() - 1, *x_finish);
				double Integ = *(GEM_int.begin() + (x_finish - _xs_GEM_sum.begin())) - *(GEM_int.begin() + (x_start - _xs_GEM_sum.begin()));
				std::cout << "GEM integral is " << Integ << std::endl;
#ifdef _TEMP_CODE
				std::ofstream GEM_results_file;
				GEM_results_file.open(std::string(OUTPUT_DIR) + OUTPUT_GEMS, std::ios_base::app);
				GEM_results_file << area_.experiments.back() << "\t" << Integ << "\t" << *x_start << "\t" << *x_finish << std::endl;
				GEM_results_file.close();
#endif
				std::cout << "GEM timelimits are: [" << *x_start << ";" << *x_finish << "]" << std::endl;
			} else {
				std::cout << "GEM time boundaries are invalid" << std::endl;
			}
		}
		//TODO: Add averaging over the runs etc.
		//TODO: Probably will have to use average value for per run processing again.
		//te best way I guess is to reverse nextRun (add prevRun, use reverse iterators) and 
		//call processOneRun (averaging results);
		//the thing is the same may be requred for experiment processing
		//The other way is do not clean SingleRunData, and simply rerun processOneRun for every member of one_run_data with external pars
	}
	graph_manager.Draw();
}

void AllRunsResults::find_GEM_start_time(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start, int N_trust, GraphicOutputManager &man)
{
	std::vector<double> xs_before_S1 = xs, ys_before_S1 = ys;
	SignalOperations::apply_time_limits(xs_before_S1, ys_before_S1, *(xs.begin()), ParameterPile::S1_time);
	DITERATOR x_befS1_max;
	double noize_amp;
	SignalOperations::get_max(xs_before_S1, ys_before_S1, x_befS1_max, noize_amp, N_trust);
	Drawing *dr = man.GetDrawing(0);
	if (dr)
		dr->AddToDraw_vertical(*x_befS1_max, -0.4, 0.4, "lc rgb \"#000000\"", 0);
	noize_amp *= ParameterPile::GEM_threshold_to_noise;
	if (dr)
		dr->AddToDraw_baseline(noize_amp, "threshold", "lc rgb \"#000000\"", 0);
	DITERATOR x_S1_left = xs.begin();
	DITERATOR x_S1_right = xs.begin();
	SignalOperations::find_next_peak(xs, ys, x_S1_left, x_S1_right, noize_amp, N_trust);
	if (x_S1_left == xs.end()){
		std::cout << "1st peak not found, threshold " << noize_amp << std::endl;
		x_start = xs.end();
		return;
	}
	if (dr)
		dr->AddToDraw_vertical(*x_S1_right, -0.4, 0.4, "lc rgb \"#00AA00\"", 0);
	//x_S1_right += N_trust;
	SignalOperations::find_next_extremum(xs, ys, x_S1_right, N_trust);
	if (dr)
		dr->AddToDraw_vertical(*x_S1_right, -0.4, 0.4, "lc rgb \"#0000AA\"", 0);
	x_S1_right += N_trust;
	SignalOperations::find_next_extremum(xs, ys, x_S1_right, N_trust);
	if (dr)
		dr->AddToDraw_vertical(*x_S1_right, -0.4, 0.4, "lc rgb \"#00AAAA\"", 0);
	x_S1_right += N_trust;
	SignalOperations::find_next_extremum(xs, ys, x_S1_right, N_trust);
	if (x_S1_right == xs.end()){
		x_start = xs.end();
		return;
	}
	if (*(ys.begin() + (x_S1_right - xs.begin())) > 0){
		x_start = x_S1_right;

		return;
	}
	x_S1_left = x_S1_right;
	SignalOperations::find_next_peak(xs, ys, x_S1_left, x_S1_right, 0, N_trust); //effectively finds intersection with 0
	x_start = x_S1_left;
}

void AllRunsResults::find_S_cutoff(void)
{
	if (_Ss.size() < 10) //not statistically significant. TODO: figure out the number, and ParameterPile, however this is not important
		return;
	DVECTOR areas = _Ss;
	std::sort(areas.begin(), areas.end());
	DVECTOR areas_full = areas;
	double max_ds = 0;
	DITERATOR i_max_ds = areas.begin();
	bool found_max_delta_S = false;
	while (!found_max_delta_S){
		for (auto i = areas.begin(); i != (areas.end() - 1); ++i) {
			if ((*(i + 1) - *i) > max_ds){
				max_ds = (*(i + 1) - *i);
				i_max_ds = i + 1;
			}
		}
		if ((areas.end() - i_max_ds) < 10) {//it is not the max of interest
			areas.erase(i_max_ds, areas.end());
			max_ds = 0;
			i_max_ds = areas.begin();
		} else {
			found_max_delta_S = true;
		}
		if (areas.size() < 10){
			i_max_ds = areas.begin();
			break;
		}
	}
	if (found_max_delta_S) {
		double mean = TMath::Mean(areas_full.begin() + (i_max_ds - areas.begin()), areas_full.end());
		double rms = TMath::RMS(areas_full.begin() + (i_max_ds - areas.begin()), areas_full.end());
		double cutoff = mean - rms * 2 ; //TODO: Parameter Pile
		S_peaks_cutoff = std::max(cutoff, S_peaks_cutoff);
	} else {
		double mean = TMath::Mean(areas_full.begin(), areas_full.end());
		double rms = TMath::RMS(areas_full.begin(), areas_full.end());
		double cutoff = mean - rms * 2; //TODO: Parameter Pile
		S_peaks_cutoff = std::max(cutoff, S_peaks_cutoff);
	}
}

void AllRunsResults::find_S_cutoff_v2(void)
{
	if (_Ss.size() < 10) //not statistically significant. TODO: figure out the number, and ParameterPile, however this is not important
		return;
	DVECTOR areas = _Ss;
	std::sort(areas.begin(), areas.end());
	DVECTOR ys;
	ys.resize(areas.size(), 1);
	DVECTOR x_spreaded, y_spreaded;
	SignalOperations::spread_peaks(areas, ys, x_spreaded, y_spreaded);
	Drawing *dr = graph_manager.GetDrawing(_exp.experiments.back() + " S_spreaded", 1, ParameterPile::DrawEngine::Gnuplot);
	dr->AddToDraw(x_spreaded, y_spreaded, _exp.experiments.back() + " S_spreaded", "w lines", 0);
	//dr->DrawData();
}
