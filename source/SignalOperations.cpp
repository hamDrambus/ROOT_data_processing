#include "SignalOperations.h"
#include "Math/Functor.h"
#include "Math/BrentMinimizer1D.h"
#include "Polynom2Order.h"
#include "GraphicOutputManager.h"

namespace SignalOperations {

	void invert_y(std::vector<double> &x_in_out, std::vector<double> &y_in_out)
	{
		for (auto i = y_in_out.begin(); i != y_in_out.end(); i++)
			*i = -*i;
	}

	double find_baseline_by_median(double approx, std::vector<double>&xs, std::vector<double>&ys, std::vector<peak> &peaks)
	{
		DVECTOR selected_y;
		for (auto i = xs.begin(), j = ys.begin(); (i != xs.end()) && (j != ys.end()); ++i, ++j){
			bool do_account = true;
			for (auto pp = peaks.begin(); pp != peaks.end(); pp++)
				if (!(*i<(*pp).left) && !(*i>(*pp).right)) {
					do_account = false;
					break;
				}
			if (do_account){
				selected_y.push_back(*j-approx);
			}
		}
		std::sort(selected_y.begin(), selected_y.end());
		int ind = selected_y.size();
		if (0 == ind)
			return approx;
		if (0 == ind % 2)
			return approx + selected_y[ind / 2];
		else
			return approx + 0.5*(selected_y[ind / 2] + selected_y[1 + (ind / 2)]);
	}

	double find_baseline_by_integral(double approx, std::vector<double>&xs, std::vector<double>&ys, std::vector<peak> &peaks)
	{
		std::vector<std::vector<double>> xs_cut, ys_cut;
		xs_cut.push_back(std::vector<double>());
		ys_cut.push_back(std::vector<double>());
		bool one_vector = true;
		for (auto i = xs.begin(), j = ys.begin(); (i != xs.end()) && (j != ys.end()); ++i, ++j){
			bool do_account = true;
			for (auto pp = peaks.begin(); pp != peaks.end(); pp++)
				if (!(*i<(*pp).left) && !(*i>(*pp).right)) {
					do_account = false;
					break;
				}
			if (do_account&&!one_vector){
				if (xs_cut.back().size() < 2){
					xs_cut.pop_back();
					ys_cut.pop_back();
				}
				xs_cut.push_back(std::vector<double>());
				ys_cut.push_back(std::vector<double>());
				one_vector = true;
			}
			if (do_account&&one_vector){
				xs_cut.back().push_back(*i);
				ys_cut.back().push_back((*j) - approx);
			}
			if (!do_account)
				one_vector = false;
		}
		if (xs_cut.back().size() < 2){
			xs_cut.pop_back();
			ys_cut.pop_back();
		}
		double Sum_dx = 0, Sum_int = 0;
		for (auto i = xs_cut.begin(), j = ys_cut.begin(); (i != xs_cut.end()), (j != ys_cut.end()); ++i, ++j){
			Sum_dx += ((*i).back()) - *(*i).begin();
			std::vector<double> tmp;
			integrate(*i, *j, tmp);
			Sum_int += tmp.back();
		}
		if (0 == Sum_dx)
			return approx;
		return (Sum_int / Sum_dx) + approx;
	}

	void find_peaks(DVECTOR &xs, DVECTOR &ys, std::vector<peak> &peaks, double base_line, double threshold, int N_trust)
	{
		peaks.clear();
		DITERATOR x_peak_l = xs.begin(), x_peak_r = xs.begin();
		while (x_peak_l != xs.end()){
			SignalOperations::find_next_peak(xs, ys, x_peak_l, x_peak_r, threshold, N_trust);
			if (x_peak_l != xs.end()){
				peak pk;
				pk.left = *x_peak_l;
				pk.right = *x_peak_r;
				SignalOperations::integrate(xs, ys, pk.S, x_peak_l, x_peak_r, base_line);
				peaks.push_back(pk);
				x_peak_l = x_peak_r;
				int delta_i = 1 + N_trust / 2;
				x_peak_l = ((xs.end() - x_peak_l) < delta_i) ? xs.end() : (x_peak_l + delta_i);
			}
		}
	}

	void integrate(std::vector<double>&xs, std::vector<double>&ys, std::vector<double> &y_out, double baseline)
	{
		y_out.clear();
		y_out.reserve(xs.size());
		for (auto ix = xs.begin(), iy = ys.begin(); (ix != xs.end()) && (iy != ys.end()); ix++, iy++){
			double prev = y_out.empty() ? 0 : y_out.back();
			double dx = (ix == (xs.end() - 1)) ? (*ix - *(ix - 1)) :
				(ix == xs.begin()) ? (*(ix + 1) - *ix) : (*(ix + 1) - *(ix - 1)) / 2;
			y_out.push_back(dx*(*iy - baseline) + prev);
		}
	}

	void integrate(std::vector<double>&xs, std::vector<double>&ys, std::vector<double> &x_out, std::vector<double> &y_out, double left, double right, double baseline)
	{
		for (auto ix = xs.begin(), iy = ys.begin(); (ix != xs.end()) && (iy != ys.end()); ix++, iy++){
			if (!((*ix) < left) && !((*ix)>right)){
				double prev = y_out.empty() ? 0 : y_out.back();
				double dx = (ix == (xs.end() - 1)) ? (*ix - *(ix - 1)) : (*(ix + 1) - *ix);
				y_out.push_back(dx*(*iy - baseline) + prev);
				x_out.push_back(*ix);
			}
		}
	}

	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out, DITERATOR left, DITERATOR right, double baseline)
	{
		x_out.clear();
		y_out.clear();
		if (right == xs.end()||left==xs.end()||(left>right))
			return;
		x_out.reserve(right - left + 1);
		y_out.reserve(right - left + 1);
		DITERATOR x_end = ++right;
		for (auto ix = left, iy = ys.begin() + (left-xs.begin()); (ix != x_end) && (iy != ys.end()); ++ix, ++iy){
			double prev = y_out.empty() ? 0 : y_out.back();
			double dx = (ix == (xs.end() - 1)) ? (*ix - *(ix - 1)) : (*(ix + 1) - *ix);
			y_out.push_back(dx*(*iy - baseline) + prev);
			x_out.push_back(*ix);
		}
	}

	void integrate(DVECTOR &xs, DVECTOR &ys, double &y_out, DITERATOR left, DITERATOR right, double baseline)
	{
		if (right == xs.end() || left == xs.end() || (left>right))
			return;
		DITERATOR x_end = ++right;
		y_out = 0;
		for (auto ix = left, iy = ys.begin() + (left - xs.begin()); (ix != x_end) && (iy != ys.end()); ++ix, ++iy){
			double dx = (ix == (xs.end() - 1)) ? (*ix - *(ix - 1)) : (*(ix + 1) - *ix);
			y_out+=dx*(*iy - baseline);
		}
	}

	void apply_time_limits(std::vector<double>&xs, std::vector<double>&ys, double x_left, double x_right)
	{
		std::vector<double>::iterator iterator_left = xs.begin(), iterator_right = xs.end() - 1;
		std::vector<double>::iterator iterator_left_y = ys.begin(), iterator_right_y = ys.end() - 1;
		for (auto ix = xs.begin(), iy = ys.begin(); (ix != xs.end()) && (iy != ys.end()); ix++, iy++){
			if (!((*ix) > x_left)){
				iterator_left = ix;
				iterator_left_y = iy;
			}
			if (!((*ix) > x_right)){
				iterator_right = ix;
				iterator_right_y = iy;
			}
		}
		xs.erase(iterator_right + 1, xs.end());
		xs.erase(xs.begin(), iterator_left);
		ys.erase(iterator_right_y + 1, ys.end());
		ys.erase(ys.begin(), iterator_left_y);
	}

	std::vector<double>::iterator find_x_iterator_by_value(std::vector<double>::iterator &x_left, std::vector<double>::iterator &x_right, double x)
	{
		for (auto h = x_left; h != x_right; h++) //find in which point of the vector the maximum realizes
			if (!(*h > x) && !(*(h + 1) < x)){
				if ((x - *h) > (*(h + 1) - x))
					return h + 1;
				else
					return h;
				break;
			}
		return x_left; // mustn't occur
	}

	void get_max(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_max, double &y_max, int N_trust)
	{
		if (xs.size() != ys.size()){
			x_max = xs.end();
			return;
		}
		bool use_fit = true;
		if (N_trust < 3) {//2nd order polynom
			N_trust = 1;
			use_fit = false;
		}
		int delta = N_trust / 2 + 1;
		y_max = *(ys.begin());
		x_max = xs.begin();

		if (use_fit){
			Polynom2Order fitter;
			for (auto i = xs.begin(), j = ys.begin(); (i != xs.end()) && (j != ys.end());
				((delta<(xs.end() - i)) ? i = i + delta : i = xs.end()), ((delta<(ys.end() - j)) ? j = j + delta : j = ys.end())){
				int shift = (int)(xs.size() - (i - xs.begin()) - N_trust) < 0 ? (xs.size() - (i - xs.begin()) - N_trust) : 0; //accounts for the end
				
				TVectorD coefs;
				DITERATOR x_left = i + shift;
				fitter(xs, ys, x_left-xs.begin(), N_trust, coefs, *x_left);

				double y_max_exact,x_max_exact;
				DITERATOR x_max_here;
				fitter.FindMaximum(x_max_here, x_max_exact, y_max_exact);
				if (y_max < y_max_exact){
					y_max = y_max_exact;
					x_max = x_max_here;
				}
			}
		} else //do not use the fit
			for (auto i = xs.begin(), j = ys.begin(); (i != xs.end()) && (j != ys.end()); ++i, ++j){
				if (*j > y_max){
					x_max = i;
					y_max = *j;
				}
			}
	}
	//after getting the peak, the next search must be started from (x_finish+1) or +(N_trust/2)!!!
	void find_next_peak(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start,
		DITERATOR &x_finish, double threshold, int N_trust)//TODO: does not handle double peak which middle
		//is slightly above the threshold (slightly means that 2nd order fit would intersect the threshold)
	{
		x_finish = xs.end();
		if ((xs.size() != ys.size())||((xs.end()-x_start)<N_trust)){
			x_start = xs.end();
			x_finish = xs.begin();
			return;
		}
		std::vector<double>::iterator minimal_iterator = x_start;
		bool use_fit = true;
		if (N_trust < 3) {//2nd order polynom
			N_trust = 1;
			use_fit = false;
		}
		int delta = N_trust / 2 + 1;
		std::vector<double>::iterator approx_x_left = minimal_iterator;
		std::vector<double>::iterator approx_x_right = xs.end();
		bool found_peak = false;

		if (use_fit) {
			Polynom2Order fitter;
			for (auto i = minimal_iterator, j = ys.begin() + (minimal_iterator - xs.begin()); (i != xs.end()) && (j != ys.end());
				((delta<(xs.end() - i)) ? i = i + delta : i = xs.end()), ((delta<(ys.end() - j)) ? j = j + delta : j = ys.end())){
				int shift = (int)(xs.size() - (i - xs.begin()) - N_trust) < 0 ? (xs.size() - (i - xs.begin()) - N_trust) : 0; //accounts for the end
				TVectorD coefs;
				std::vector<double>::iterator x_left = i + shift;
				fitter(xs, ys, (x_left - xs.begin()), N_trust, coefs, *x_left);

				DITERATOR x_inter1, x_inter2;
				double x_inter_exact1, x_inter_exact2;
				fitter.FindIntersection(x_inter1, x_inter2, x_inter_exact1, x_inter_exact2, threshold);
				if (x_inter2 != xs.end()){
					if (found_peak){
						if (fitter.Derivative(x_inter_exact2) <= 0.0){
							x_finish = x_inter2;
							return;
						}
					} else {
						if (fitter.Derivative(x_inter_exact2) >= 0.0){
							x_start = x_inter2;
							found_peak = true;
						}
					}
				}
				if (x_inter1 != xs.end()){
					if (found_peak){
						if (fitter.Derivative(x_inter_exact1) <= 0.0){
							x_finish = x_inter1;
							return;
						}
					} else {
						if (fitter.Derivative(x_inter_exact1) >= 0.0){
							x_start = x_inter1;
							found_peak = true;
						}
					}
				}
			}
		} else { //do not use the fit
			for (auto i = minimal_iterator, j = ys.begin() + (minimal_iterator - xs.begin()); (i != xs.end()) && (j != ys.end()); ++i, ++j){
				if (found_peak){
					if ((*j) < threshold){
						x_finish= i - 1; //must be valid
						return;
					}
				} else {
					if ((*j) >= threshold){
						x_start = i;
						found_peak = true;
					}
				}
			}
		}
		x_start = xs.end();
	}

	void find_next_extremum(std::vector<double>&xs, std::vector<double>&ys, std::vector<double>::iterator &x_start, int N_trust)
	{
		bool use_fit = true;
		if (N_trust < 3) {//2nd order polynom
			N_trust = 1;
			use_fit = false;
		}
		int delta = N_trust / 2 + 1;
		if (use_fit){
			Polynom2Order fitter;
			for (auto i = x_start, j = ys.begin()+(x_start - xs.begin()); (i != xs.end()) && (j != ys.end());
				((delta<(xs.end() - i)) ? i = i + delta : i = xs.end()), ((delta<(ys.end() - j)) ? j = j + delta : j = ys.end())){
				//TODO: is valid?
				int shift = (int)(xs.size() - (i - xs.begin()) - N_trust) < 0 ? (xs.size() - (i - xs.begin()) - N_trust) : 0; //accounts for the end
				TVectorD coefs;
				std::vector<double>::iterator x_left = i + shift;
				fitter(xs, ys, (i - xs.begin()) + shift, N_trust, coefs, *x_left);
				DITERATOR x_extr;
				double x_extr_exact, y_extr_exact;
				fitter.FindExtremum(x_extr, x_extr_exact, y_extr_exact);
				if (x_extr != xs.end()){
					if (x_extr <= x_start)
						continue; //accidently found previous extremum
					x_start = x_extr;
					return;
				}
			}
		} else{ //do not use the fit
			for (auto i = x_start, j = ys.begin() +(x_start-xs.begin()) ; (i != xs.end()) && (j != ys.end()); ++i, ++j){
				if ((i == xs.begin()) || ((i + 1) == xs.end()))
					continue;
				if (!((*i - *(i - 1))*(*(i + 1) - *i) < 0)){
					x_start = i;
					return;
				}
			}
		}
		x_start = xs.end(); //not found
	}

	//void find_GEM_start_time(std::vector<double>&xs, std::vector<double>&ys, std::vector<double>::iterator &x_start, int N_trust)
	//{
	//	DVECTOR xs_before_S1 = xs, ys_before_S1 = ys;
	//	apply_time_limits(xs_before_S1, ys_before_S1, *(xs.begin()), ParameterPile::S1_time);
	//	DITERATOR x_befS1_max;
	//	double noize_amp;
	//	get_max(xs_before_S1, ys_before_S1, x_befS1_max, noize_amp, N_trust);
	//	noize_amp *= 1.1;
	//	DITERATOR x_S1_left = xs.begin();
	//	DITERATOR x_S1_right = xs.begin();
	//	find_next_peak(xs, ys, x_S1_left, x_S1_right, noize_amp, N_trust);
	//	find_next_extremum(xs, ys, (++x_S1_right), N_trust);
	//	find_next_extremum(xs, ys, (++x_S1_right), N_trust);
	//	find_next_extremum(xs, ys, (++x_S1_right), N_trust);
	//	if (x_S1_right == xs.end()){
	//		x_start = xs.end();
	//		return;
	//	}
	//	if (*(ys.begin() + (x_S1_right - xs.begin())) > 0){
	//		x_start = x_S1_right;
	//		return;
	//	}
	//	x_S1_left = x_S1_right;
	//	find_next_peak(xs, ys, x_S1_left, x_S1_right, 0, N_trust); //effectively finds intersection with 0
	//	x_start = x_S1_left;
	//}

	void spread_peaks(double x_left, double x_right, std::vector<peak> &peaks)
	{

	}

	void spread_peaks(DVECTOR &xs_in, DVECTOR &ys_in,DVECTOR &xs_out, DVECTOR& ys_out)
	{
		xs_out.clear();
		ys_out.clear();
		xs_out.reserve(xs_in.size());
		ys_out.reserve(ys_in.size());
		bool first= true; //eliminating the same xs
		double x_prev;
		double y_val;
		for (auto i = xs_in.begin(), j = ys_in.begin(); (i != xs_in.end()) && (j != ys_in.end()); ++i, ++j) {
			if (first){
				x_prev = *i;
				y_val = *j;
				first= false;
			} else {
				if (x_prev == *i)
					y_val += *j;
				else {
					xs_out.push_back(x_prev);
					ys_out.push_back(y_val);
					x_prev = *i;
					y_val = *j;
				}
			}
		}
		xs_out.push_back(x_prev);
		ys_out.push_back(y_val);
		xs_in = xs_out;
		ys_in = ys_out; //eliminated the same xs

		if ((xs_in.size() < 2)||ys_in.size()!=xs_in.size())
			return;
		xs_out.reserve(2*xs_in.size());
		ys_out.reserve(2*ys_in.size());
		for (auto i = xs_in.begin(), j = ys_in.begin(); (i != (xs_in.end() - 1)) && (j != (ys_in.end() - 1)); ++i, ++j){
			xs_out.push_back(*i);
			xs_out.push_back(*(i+1)-1e-6); //so I can plot it with gnuplot
			ys_out.push_back((*(j + 1) + (*j)) / (*(i + 1) - *i)); //xs are now different.
			ys_out.push_back((*(j + 1) + (*j)) / (*(i + 1) - *i));
		}
	}

};
