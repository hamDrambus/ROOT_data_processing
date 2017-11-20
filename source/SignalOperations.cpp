#include "SignalOperations.h"
#include "Math/Functor.h"
#include "Math/BrentMinimizer1D.h"
#include "TSpectrum.h"
#include "Polynom2Order.h"
#include "GraphicOutputManager.h"

namespace SignalOperations {

	void invert_y(DVECTOR &x_in_out, DVECTOR &y_in_out)
	{
		for (auto i = y_in_out.begin(); i != y_in_out.end(); i++)
			*i = -*i;
	}

	double find_baseline_by_median(double approx, DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks)
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

	double find_baseline_by_integral(double approx, DVECTOR &xs, DVECTOR &ys)
	{
		if ((xs.size() <= 1) || (xs.size() != ys.size()))
			return approx;
		double val = 0;
		double dx = xs.back()-*(xs.begin());
		if (0 == dx)
			return approx;
		integrate(xs, ys, val,xs.begin(),--xs.end());
		return val/dx;
	}

	double find_baseline_by_integral(double approx, DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks)
	{
		STD_CONT<DVECTOR> xs_cut, ys_cut;
		xs_cut.push_back(DVECTOR());
		ys_cut.push_back(DVECTOR());
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
				xs_cut.push_back(DVECTOR());
				ys_cut.push_back(DVECTOR());
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
			DVECTOR tmp;
			integrate(*i, *j, tmp);
			Sum_int += tmp.back();
		}
		if (0 == Sum_dx)
			return approx;
		return (Sum_int / Sum_dx) + approx;
	}

	void find_baseline_by_ROOT(DVECTOR &xs, DVECTOR &ys, DVECTOR &ys_out)
	{
#ifdef _HOTFIX_CLEAR_MEMORY
		DVECTOR().swap(ys_out);
#else
		ys_out.clear();
#endif
		TSpectrum *spec = new TSpectrum();
		float *f_ys = new float[ys.size()];
		for (int h = 0; h != ys.size(); ++h)
			f_ys[h] = ys[h];
		//TODO: ? ParameterPile and as input parameters?
		spec->Background(f_ys, ys.size(), 50, TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false);
#ifndef _USE_DEQUE
		ys_out.reserve(ys.size());
#endif
		for (int h = 0; h != ys.size(); ++h)
			ys_out.push_back(f_ys[h]);
		delete [] f_ys;
		spec->Delete();
	}

	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &y_out, double baseline)
	{
#ifdef _HOTFIX_CLEAR_MEMORY
		DVECTOR().swap(y_out);
#else
		y_out.clear();
#endif
#ifndef _USE_DEQUE
		y_out.reserve(xs.size());
#endif
		if (xs.size() == 1){
			y_out.push_back(0);
			return;
		}
		for (auto ix = xs.begin(), iy = ys.begin(); (ix != xs.end()) && (iy != ys.end()); ix++, iy++){
			double prev = y_out.empty() ? 0 : y_out.back();
			double dx = (ix == (xs.end() - 1)) ? (*ix - *(ix - 1)) :
				(ix == xs.begin()) ? (*(ix + 1) - *ix) : (*(ix + 1) - *(ix - 1)) / 2;
			y_out.push_back(dx*(*iy - baseline) + prev);
		}
	}

	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out, double left, double right, double baseline)
	{
		if (xs.size() == 1){
			if ((xs.back() < left) || (xs.back() > right))
				return;
			y_out.push_back(0);
			x_out.push_back(xs.back());
			return;
		}
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
#ifdef _HOTFIX_CLEAR_MEMORY
		DVECTOR().swap(x_out);
		DVECTOR().swap(y_out);
#else
		x_out.clear();
		y_out.clear();
#endif
		if (right == xs.end()||left==xs.end()||(left>right))
			return;
#ifndef _USE_DEQUE
		x_out.reserve(right - left + 1);
		y_out.reserve(right - left + 1);
#endif
		if (xs.size() == 1) {
			if (left != xs.end()){
				x_out.push_back(*left);
				y_out.push_back(0);
			}
			return;
		}
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
		if (xs.size() <= 1)
			return;
		for (auto ix = left, iy = ys.begin() + (left - xs.begin()); (ix != x_end) && (iy != ys.end()); ++ix, ++iy){
			double dx = (ix == (xs.end() - 1)) ? (*ix - *(ix - 1)) : (*(ix + 1) - *ix);
			y_out+=dx*(*iy - baseline);
		}
	}

	void apply_time_limits(DVECTOR &xs, DVECTOR &ys, double x_left, double x_right)
	{
		DITERATOR iterator_left = xs.begin(), iterator_right = xs.end() - 1;
		DITERATOR iterator_left_y = ys.begin(), iterator_right_y = ys.end() - 1;
		for (auto ix = xs.begin(), iy = ys.begin(); (ix != xs.end()) && (iy != ys.end()); ix++, iy++){
			if (!((*ix) > x_left)) {
				iterator_left = ix;
				iterator_left_y = iy;
			}
			if (!((*ix) > x_right)) {
				iterator_right = ix;
				iterator_right_y = iy;
			}
		}
		xs.erase(iterator_right + 1, xs.end());
		xs.erase(xs.begin(), iterator_left);
		ys.erase(iterator_right_y + 1, ys.end());
		ys.erase(ys.begin(), iterator_left_y);
	}

	DITERATOR find_x_iterator_by_value(DITERATOR &x_left, DITERATOR &x_right, double x)
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
	void get_max(DVECTOR &xs, DVECTOR &ys, DITERATOR x_start, DITERATOR x_finish, DITERATOR &x_max, double &y_max, int N_trust)
	{
		N_trust = std::min((int)(x_finish - x_start), N_trust); //other funtions return invalid results in this case
		if (xs.size() != ys.size()){
			x_max = xs.end();
			return;
		}
		bool use_fit = true;
		if (N_trust < 3) {//2nd order polynom
			N_trust = 1;
			use_fit = false;
		}
		int delta = N_trust / 3;
		y_max = *(ys.begin()+(x_start - xs.begin()));
		x_max = x_start;

		if (use_fit){
			Polynom2Order fitter;
			for (auto i = x_start, j = ys.begin() + (x_start-xs.begin()); (i != xs.end()) && (j != ys.end())&&(i<x_finish);
				((delta<(xs.end() - i)) ? i = i + delta : i = xs.end()), ((delta<(ys.end() - j)) ? j = j + delta : j = ys.end())){
				int shift = (int)(xs.size() - (i - xs.begin()) - N_trust) < 0 ? (xs.size() - (i - xs.begin()) - N_trust) : 0; //accounts for the end

				TVectorD coefs;
				DITERATOR x_left = i + shift;
				fitter(xs, ys, x_left - xs.begin(), N_trust, coefs, *x_left);

				double y_max_exact, x_max_exact;
				DITERATOR x_max_here;
				fitter.FindMaximum(x_max_here, x_max_exact, y_max_exact);
				if (y_max < y_max_exact){
					y_max = y_max_exact;
					x_max = x_max_here;
				}
			}
		} else //do not use the fit
			for (auto i = x_start, j = ys.begin() + (x_start - xs.begin()); (i != xs.end()) && (j != ys.end()) && (i<x_finish); ++i, ++j){
				if (*j > y_max){
					x_max = i;
					y_max = *j;
				}
			}
	}
	void get_max(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_max, double &y_max, int N_trust)
	{	get_max(xs, ys, xs.begin(), xs.end(), x_max, y_max, N_trust);}
	//after getting the peak, the next search must be started from (x_finish+1) or +(N_trust/2)!!!
	void find_next_peak(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start,
		DITERATOR &x_finish, double threshold, int N_trust)//done - now every step by x uses fitting 
		//TODO: does not handle double peak which middle is slightly above
		//the threshold (slightly means that 2nd order fit would intersect the threshold)
	{
		x_finish = xs.end();
		if ((xs.size() != ys.size()) || ((xs.end() - x_start) < N_trust)){
			x_start = xs.end();
			x_finish = xs.begin();
			return;
		}
		DITERATOR minimal_iterator = x_start;
		x_start = xs.end();
		bool use_fit = true;
		if (N_trust < 3) {//2nd order polynom
			N_trust = 1;
			use_fit = false;
		}
		int delta = N_trust / 3;
		DITERATOR approx_x_left = minimal_iterator;
		DITERATOR approx_x_right = xs.end();
		bool found_peak = false;

		if (use_fit) {
			Polynom2Order fitter;
			for (auto i = minimal_iterator, j = ys.begin() + (minimal_iterator - xs.begin()); (i != xs.end()) && (j != ys.end());
				((delta < (xs.end() - i)) ? i = i + delta : i = xs.end()), ((delta < (ys.end() - j)) ? j = j + delta : j = ys.end())) {
				int shift = (int)(xs.size() - (i - xs.begin()) - N_trust) < 0 ? (xs.size() - (i - xs.begin()) - N_trust) : 0; //accounts for the end
				TVectorD coefs;
				DITERATOR x_left = i + shift;
				fitter(xs, ys, (x_left - xs.begin()), N_trust, coefs, *x_left);

				//#ifdef _TEMP_CODE
				//				if ((*i <= 32.49) && (32.49 <= *(i + N_trust - 1))) {
				//					GraphicOutputManager man;
				//					Drawing *dr = man.GetDrawing("Peak find test "+std::to_string(*j), 0, ParameterPile::DrawEngine::Gnuplot);
				//					DVECTOR tmp_x, tmp_y;
				//					for (auto ti = i, tj = j; (ti < (i + N_trust)) && (tj < (j + N_trust)); ++ti, ++tj) {
				//						tmp_x.push_back(*ti);
				//						tmp_y.push_back(*tj);
				//					}
				//					dr->AddToDraw(tmp_x, tmp_y, "peak " + std::to_string(*j));
				//					TVectorD coefs;
				//					fitter.getCoefs(coefs);
				//					double a = coefs[2];
				//					double b = coefs[1]-2**x_left*coefs[2];
				//					double c = coefs[0] - coefs[1] * *x_left + coefs[2] * *x_left * *x_left;
				//					std::stringstream aa, bb, cc;
				//					aa << std::setprecision(12) << a;
				//					bb << std::setprecision(12) << b;
				//					cc << std::setprecision(12) << c;
				//					dr->AddToDraw("a = " + aa.str() + "\nb = " + bb.str() + "\nc = " + cc.str() + "\nf(x) = a*x*x + b*x + c", "f(x)", "fit", "w l", 0);
				//					man.Draw();
				//				}
				//#endif
				DITERATOR x_inter1, x_inter2;
				double x_inter_exact1, x_inter_exact2;
				fitter.FindIntersection(x_inter1, x_inter2, x_inter_exact1, x_inter_exact2, threshold);
				if (x_inter2 != xs.end()){
					if (found_peak){
						if (fitter.Derivative(x_inter_exact2) < 0.0){
							x_finish = x_inter2;
							return;
						}
					} else {
						if (fitter.Derivative(x_inter_exact2) > 0.0){
							x_start = x_inter2;
							found_peak = true;
						}
					}
				}
				if (x_inter1 != xs.end()){
					if (found_peak){
						if (fitter.Derivative(x_inter_exact1) < 0.0){
							x_finish = x_inter1;
							return;
						}
					} else {
						if (fitter.Derivative(x_inter_exact1) > 0.0){
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
						x_finish = i - 1; //must be valid
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
		if (found_peak&&x_start != xs.end()){
			x_finish = xs.end() - 1;
			return;
		}
		x_start = xs.end();
	}

	void find_peaks(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks, double base_line, double threshold, int N_trust)
	{
#ifdef _HOTFIX_CLEAR_MEMORY
		STD_CONT<peak>().swap(peaks);
#else
		peaks.clear();
#endif
		DITERATOR x_peak_l = xs.begin(), x_peak_r = xs.begin();
		while (x_peak_l != xs.end()){
			SignalOperations::find_next_peak(xs, ys, x_peak_l, x_peak_r, threshold, N_trust);
			if (x_peak_l != xs.end()){
				peak pk;
				pk.left = *x_peak_l;
				pk.right = *x_peak_r;
				SignalOperations::integrate(xs, ys, pk.S, x_peak_l, x_peak_r, base_line);
				DITERATOR pk_max;
				SignalOperations::get_max(xs, ys, x_peak_l, x_peak_r + 1, pk_max, pk.A, N_trust);
				if ((pk_max != xs.end()) && (pk.S > 0) && (pk.right>=pk.left))
					peaks.push_back(pk);
				x_peak_l = x_peak_r;
				int delta_i = std::max(N_trust / 3, 1);
				x_peak_l = ((xs.end() - x_peak_l) < delta_i) ? xs.end() : (x_peak_l + delta_i);
			}
		}
	}

	//seaches peak from x_start, in difference to find_next peak this one first finds peak by threshold_finder and then finds its edges (wider 
	//than intersection of threshold and signal) using thresh_edges.
	void find_next_peak_fine(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start, DITERATOR &x_finish, double &Amp,
		double thresh_finder, double thresh_edges, int N_trust)
	{
		if (thresh_edges >= thresh_finder)
			thresh_edges = 0;
		x_finish = xs.end();
		DITERATOR minimal_iterator = x_start;
		DITERATOR pk_max;
		if ((xs.size() != ys.size()) || ((xs.end() - x_start)<N_trust))
			goto bad_return;
		find_next_peak(xs, ys, x_start, x_finish, thresh_finder, N_trust);
		if (x_start == xs.end())
			goto bad_return;
		SignalOperations::get_max(xs, ys, x_start, x_finish + 1, pk_max, Amp, N_trust);
		if (pk_max == xs.end())
			goto bad_return;

		bool use_fit = true;
		if (N_trust < 3) {//2nd order polynom
			N_trust = 1;
			use_fit = false;
		}
		int delta = N_trust / 3;
		if (use_fit){ //now extend edges
			Polynom2Order fitter;
			for (auto i = x_finish, j = ys.begin() + (x_finish - xs.begin()); (i != xs.end()) && (j != ys.end());
				((delta<(xs.end() - i)) ? i = i + delta : i = xs.end()), ((delta<(ys.end() - j)) ? j = j + delta : j = ys.end())){
				int shift = (int)(xs.size() - (i - xs.begin()) - N_trust) < 0 ? (xs.size() - (i - xs.begin()) - N_trust) : 0; //accounts for the end
				TVectorD coefs;
				DITERATOR x_left = i + shift;
				fitter(xs, ys, (i - xs.begin()) + shift, N_trust, coefs, *x_left);
				DITERATOR x_extr;
				DITERATOR x_intersect1, x_intersect2;
				DITERATOR x_inter = xs.end(); //of interest
				double x_extr_exact, y_extr_exact;
				double x_inter_exact1, x_inter_exact2;
				fitter.FindExtremum(x_extr, x_extr_exact, y_extr_exact);
				fitter.FindIntersection(x_intersect1, x_intersect2, x_inter_exact1, x_inter_exact2, thresh_edges);
				if (x_extr <= x_finish)
					x_extr = xs.end(); //accidently found extremum as peak maximum (in case x_start==x_finish)
				if (x_intersect2 != xs.end())
					if (fitter.Derivative(x_inter_exact2) < 0)
						x_inter = x_intersect2;
				if (x_intersect1 != xs.end()&&(x_inter==xs.end()))
					if (fitter.Derivative(x_inter_exact1) < 0)
						x_inter = x_intersect1;
				if (x_extr != xs.end()&&x_inter!=xs.end()){
					if (y_extr_exact < thresh_edges)
						x_finish = x_inter;
					else
						x_finish = x_extr;
					break;
				}
				if (x_extr != xs.end()){
					x_finish = x_extr;
					break;
				}
				if (x_inter != xs.end()){
					x_finish = x_inter;
					break;
				}
			}
			D_REV_ITERATOR x_left_peak = D_REV_ITERATOR(x_start);
			D_REV_ITERATOR x_rend = D_REV_ITERATOR(minimal_iterator);
			for (auto i = x_left_peak, j = ys.rbegin() + (x_left_peak - xs.rbegin()); (i != x_rend) && (j != ys.rend());
				((delta<(xs.rend() - i)) ? i = i + delta : i = x_rend), ((delta<(ys.rend() - j)) ? j = j + delta : j = ys.rend())){
				int shift = (int)(xs.size() - (i - xs.rbegin()) - N_trust) < 0 ? (xs.size() - (i - xs.rbegin()) - N_trust) : 0; //accounts for the rend
				TVectorD coefs;
				DITERATOR x_left = (i + shift).base();
				fitter(xs, ys, x_left - xs.begin(), N_trust, coefs, *x_left);
				DITERATOR x_extr;
				DITERATOR x_intersect1, x_intersect2;
				DITERATOR x_inter = xs.end(); //of interest
				double x_extr_exact, y_extr_exact;
				double x_inter_exact1, x_inter_exact2;
				fitter.FindExtremum(x_extr, x_extr_exact, y_extr_exact);
				fitter.FindIntersection(x_intersect1, x_intersect2, x_inter_exact1, x_inter_exact2, thresh_edges);
				if (x_extr >= x_start)
					x_extr = xs.end(); //accidently found extremum as peak maximum (in case x_start==x_finish)
				if (x_intersect2 != xs.end())
					if (fitter.Derivative(x_inter_exact2) > 0)
						x_inter = x_intersect2;
				if (x_intersect1 != xs.end() && (x_inter == xs.end()))
					if (fitter.Derivative(x_inter_exact1) > 0)
						x_inter = x_intersect1;
				if (x_extr != xs.end() && x_inter != xs.end()){
					if (y_extr_exact < thresh_edges)
						x_start = x_inter;
					else
						x_start = x_extr;
					break;
				}
				if (x_extr != xs.end()){
					x_start = x_extr;
					break;
				}
				if (x_inter != xs.end()){
					x_start = x_inter;
					break;
				}
			}
			if (x_start < minimal_iterator)
				x_start = minimal_iterator;
		} else { //do not use the fit
			for (auto i = x_finish, j = ys.begin() + (x_finish - xs.begin()); (i != xs.end()) && (j != ys.end()); ++i, ++j){
				if ((i == xs.begin()) || ((i + 1) == xs.end()))
					continue;
				if ((*j - thresh_edges)*(*(j - 1) - thresh_edges) <= 0){ //intersection
					x_finish = i; //TODO: more exact?
					break;
				}
				if (!((*j - *(j - 1))*(*(j + 1) - *j) <= 0)){ //extremum
					if (*j > thresh_edges)
						x_finish = i;
					break;
				}
			}
			D_REV_ITERATOR x_left_peak = D_REV_ITERATOR(x_start);
			D_REV_ITERATOR x_rend = D_REV_ITERATOR(minimal_iterator);
			for (auto i = x_left_peak, j = ys.rbegin() + (x_left_peak - xs.rbegin()); (i != x_rend) && (j != ys.rend()); ++i, ++j){
				if ((i == xs.rbegin()) || ((i + 1) == xs.rend()))
					continue;
				if ((*j - thresh_edges)*(*(j - 1) - thresh_edges) <= 0){
					x_finish = i.base(); //TODO: more exact?
					break;
				}
				if (!((*j - *(j - 1))*(*(j + 1) - *j) <= 0)){
					if (*j > thresh_edges) {
						x_start = i.base();
						break;
					}
				}
			}
			if (x_start < minimal_iterator)
				x_start = minimal_iterator;
		}
		return;
	bad_return:
		Amp = -1;
		x_start = xs.end();
		x_finish = xs.begin();
		return;
	}
	void find_peaks_fine(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks, double base_line, double threshold, double threshold_edges, int N_trust)
	{
		if (threshold_edges >= threshold)
			threshold_edges = 0;
#ifdef _HOTFIX_CLEAR_MEMORY
		STD_CONT<peak>().swap(peaks);
#else
		peaks.clear();
#endif
		DITERATOR x_peak_l = xs.begin(), x_peak_r = xs.begin();
		while (x_peak_l != xs.end()){
			double Amp;
			SignalOperations::find_next_peak_fine(xs, ys, x_peak_l, x_peak_r,Amp, threshold, threshold_edges, N_trust);
			if (x_peak_l != xs.end()){
				peak pk;
				pk.left = *x_peak_l;
				pk.right = *x_peak_r;
				pk.A = Amp;
				SignalOperations::integrate(xs, ys, pk.S, x_peak_l, x_peak_r, base_line);
				if ((pk.S>0) && (pk.A>0)&&(pk.right>=pk.left))
					peaks.push_back(pk);
				x_peak_l = x_peak_r;
				int delta_i = std::max(N_trust / 3, 1);
				x_peak_l = ((xs.end() - x_peak_l) < delta_i) ? xs.end() : (x_peak_l + delta_i);
			}
		}
	}

	void find_next_extremum(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start, int N_trust)
	{
		bool use_fit = true;
		if (N_trust < 3) {//2nd order polynom
			N_trust = 1;
			use_fit = false;
		}
		int delta = N_trust / 3;
		if (use_fit) {
			Polynom2Order fitter;
			for (auto i = x_start, j = ys.begin()+(x_start - xs.begin()); (i != xs.end()) && (j != ys.end());
				((delta<(xs.end() - i)) ? i = i + delta : i = xs.end()), ((delta<(ys.end() - j)) ? j = j + delta : j = ys.end())){
				int shift = (int)(xs.size() - (i - xs.begin()) - N_trust) < 0 ? (xs.size() - (i - xs.begin()) - N_trust) : 0; //accounts for the end
				TVectorD coefs;
				DITERATOR x_left = i + shift;
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
		} else { //do not use the fit
			for (auto i = x_start, j = ys.begin() +(x_start-xs.begin()) ; (i != xs.end()) && (j != ys.end()); ++i, ++j){
				if ((i == xs.begin()) || ((i + 1) == xs.end()))
					continue;
				if (!((*j - *(j - 1))*(*(j + 1) - *j) <= 0)){
					x_start = i;
					return;
				}
			}
		}
		x_start = xs.end(); //not found
	}

	void spread_peaks(double x_left, double x_right, STD_CONT<peak> &peaks, DVECTOR &xs_out, DVECTOR& ys_out)
	{
		//doesn't check whether peaks are valid (e.g. peak.A<0)
#ifdef _HOTFIX_CLEAR_MEMORY
		DVECTOR().swap(xs_out);
		DVECTOR().swap(ys_out);
#else
		xs_out.clear();
		ys_out.clear();
#endif
		for (auto pp = peaks.begin(); pp != peaks.end(); ++pp){
			bool is_first = (pp == peaks.begin());
			double x_l = is_first ? x_left - 1e-7 : ((pp - 1)->left + (pp - 1)->right) / 2;
			double x_r = (pp->left + pp->right)/2;
			double y_v = is_first ? 0.5*pp->S : 0.5*((pp - 1)->S + pp->S);
			y_v /= x_r - x_l;
			xs_out.push_back(x_l + 1e-7);
			xs_out.push_back(x_r - 1e-7);
			ys_out.push_back(y_v);
			ys_out.push_back(y_v);
		}
		double x_l = peaks.empty() ? x_left : 0.5*(peaks.back().left + peaks.back().right);
		double x_r = x_right;
		double y_v = peaks.empty() ? 0 : 0.5*(peaks.back().S);
		y_v /= x_r - x_l;
		xs_out.push_back(x_l + 1e-7);
		xs_out.push_back(x_r);
		ys_out.push_back(y_v);
		ys_out.push_back(y_v);
	}
	//in comparisson to ^ has more smooth result
	void spread_peaks_v2(double x_left, double x_right, STD_CONT<peak> &peaks, DVECTOR &xs_out, DVECTOR& ys_out, double min_dx)
	{
		//doesn't check whether peaks are valid (e.g. peak.A<0)
#ifdef _HOTFIX_CLEAR_MEMORY
		DVECTOR().swap(xs_out);
		DVECTOR().swap(ys_out);
#else
		xs_out.clear();
		ys_out.clear();
#endif
		double current_left_x;
		bool uniting_several_peaks = false;
		double y_v=0;
		for (auto pp = peaks.begin(); pp != peaks.end(); ++pp){
			bool is_first = (pp == peaks.begin());
			double x_l = uniting_several_peaks ? current_left_x :
				(is_first ? x_left - 1e-7 : ((pp - 1)->left + (pp - 1)->right) / 2);
			double x_r = (pp->left + pp->right) / 2;
			y_v += is_first ? 0.5*pp->S : 0.5*((pp - 1)->S + pp->S);
			if ((x_r - x_l) < min_dx){
				uniting_several_peaks = true;
				current_left_x = x_l;
				continue;
			} else {
				uniting_several_peaks = false;
			}
			y_v /= x_r - x_l;
			xs_out.push_back(x_l + 1e-7);
			xs_out.push_back(x_r - 1e-7);
			ys_out.push_back(y_v);
			ys_out.push_back(y_v);
			y_v = 0;
		}
		double x_l = uniting_several_peaks ? current_left_x : 
			peaks.empty() ? x_left : 0.5*(peaks.back().left + peaks.back().right);
		double x_r = x_right;
		double dx = x_r - x_l;
		if (xs_out.size() >= 2)
			dx = std::max(x_r - *(xs_out.end() - 2), std::max(min_dx,dx));
		else
			dx = std::max(min_dx,dx);
		y_v += peaks.empty() ? 0 : 0.5*(peaks.back().S);
		y_v /= dx;
		xs_out.push_back(x_l + 1e-7);
		xs_out.push_back(x_r);
		ys_out.push_back(y_v);
		ys_out.push_back(y_v);
	}

	void peaks_to_yx(double x_left, double x_right, STD_CONT<peak> &peaks, DVECTOR &xs_out, DVECTOR& ys_out)
	{
#ifdef _HOTFIX_CLEAR_MEMORY
		DVECTOR().swap(xs_out);
		DVECTOR().swap(ys_out);
#else
		xs_out.clear();
		ys_out.clear();
#endif
		xs_out.push_back(x_left);
		ys_out.push_back(0);
		for (auto pp = peaks.begin(); pp != peaks.end(); ++pp) {
			double y_v = pp->S / (pp->right - pp->left);
			xs_out.push_back(pp->left - 1e-7);
			ys_out.push_back(0);
			xs_out.push_back(pp->left + 1e-7);
			ys_out.push_back(y_v);
			xs_out.push_back(pp->right - 1e-7);
			ys_out.push_back(y_v);
			xs_out.push_back(pp->right + 1e-7);
			ys_out.push_back(0);
		}
		xs_out.push_back(x_right);
		ys_out.push_back(0);
	}

	void spread_peaks(DVECTOR &xs_in, DVECTOR &ys_in, DVECTOR &xs_out, DVECTOR& ys_out)
	{
#ifdef _HOTFIX_CLEAR_MEMORY
		DVECTOR().swap(xs_out);
		DVECTOR().swap(ys_out);
#else
		xs_out.clear();
		ys_out.clear();
#endif
#ifndef _USE_DEQUE
		xs_out.reserve(xs_in.size());
		ys_out.reserve(ys_in.size());
#endif
		bool first = true; //eliminating the same xs
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
#ifndef _USE_DEQUE
		xs_out.reserve(2*xs_in.size());
		ys_out.reserve(2*ys_in.size());
#endif
		for (auto i = xs_in.begin(), j = ys_in.begin(); (i != (xs_in.end() - 1)) && (j != (ys_in.end() - 1)); ++i, ++j){
			xs_out.push_back(*i);
			xs_out.push_back(*(i+1)-1e-7); //so I can plot it with gnuplot
			ys_out.push_back((*(j + 1) + (*j)) / (*(i + 1) - *i)); //xs are now different.
			ys_out.push_back((*(j + 1) + (*j)) / (*(i + 1) - *i));
		}
	}

	void exclude_peaks(DVECTOR &xs_in, DVECTOR &ys_in, STD_CONT<peak> &peaks)
	{
		DVECTOR x_out, y_out;
#ifndef _USE_DEQUE
		x_out.reserve(xs_in.size());
		y_out.reserve(ys_in.size());
#endif
		for (auto i = xs_in.begin(), j = ys_in.begin(); (i != xs_in.end()) && (j != ys_in.end()); ++i, ++j){
			bool do_account = true;
			for (auto pp = peaks.begin(); pp != peaks.end(); pp++)
				if (!(*i<(*pp).left) && !(*i>(*pp).right)) {
					do_account = false;
					break;
				}
			if (do_account) {
				x_out.push_back(*i);
				y_out.push_back(*j);
			}
		}
		xs_in = x_out;
		ys_in = y_out;
	}

	void substract_baseline(DVECTOR &ys_in, double base_line)
	{
		for (auto i = ys_in.begin(); i != ys_in.end(); ++i)
			*i -= base_line;
	}
	void substract_baseline(DVECTOR &ys_in, DVECTOR &base_ys)
	{
		if (ys_in.size() != base_ys.size())
			return;
		for (auto i = ys_in.begin(), j = base_ys.begin(); (i != ys_in.end())&&(j!=base_ys.end()); ++i, ++j)
			*i -= *j;
	}

};
