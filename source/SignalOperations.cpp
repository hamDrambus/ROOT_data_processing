#include "SignalOperations.h"
#ifndef _AVOID_CERN_ROOT
#include "Math/Functor.h"
#include "Math/BrentMinimizer1D.h"
#endif
#include "Polynom2Order.h"
#include "GraphCollection.h"

namespace SignalOperations {

	//Finds indices in xs bounding value val, i.e. val>xs[first] && val<xs[second].
	//If val==xs[first] returns (first, first) pair
	//If val is outside the array, returns index to the nearest element to val.
	boost::optional<std::pair<std::size_t, std::size_t> > getX_indices(const DVECTOR& xs, double val)
	{
		boost::optional<std::pair<std::size_t, std::size_t>> out;
		std::size_t sz = xs.size();
		if (0 == sz)
			return out;
		if (val <= xs.front()) {
			out = std::pair<std::size_t, std::size_t>(0, 0);
			return out;
		}
		if (val >= xs.back()) {
			out = std::pair<std::size_t, std::size_t>(sz - 1, sz - 1);
			return out;
		}
		//find first x which is not less that X_point. That is index bounding X_point: xs[first] <= X_point < xs[first + 1]
		//See std::lower_bound and std::upper_bound
		std::size_t count = sz;
		std::size_t first = 0;
		//std::lower_bound(xys.begin(), xys.end(), [](const std::pair<double, double> &a, const std::pair<double, double> &b)->bool{
		//	return a.first<b.first;
		//});
		while (count > 0) {
			std::size_t step = count / 2;
			std::size_t ind = first + step;
			if (!(val < xs[ind])) {
				first = ++ind;
				count -= step + 1;
			} else
				count = step;
		}
		//first is such, that x>=xs[first-1] and x<xs[first]
		//first always != 0 here
		--first;
		if (val == xs[first]) {
			out = std::pair<std::size_t, std::size_t>(first, first);
			return out;
		}
		out = std::pair<std::size_t, std::size_t>(first, first + 1);
		return out;
	}

	//Wrapper for std::sort. Sorts peaks by left fronts.
	//Peaks must not intersect, i.e. pk[i].right <= pk[i+1].left and pk[i].left<=pk[i].right for sorted pk array.
	void sort_peaks(STD_CONT<peak> &peaks)
	{
		std::sort(peaks.begin(), peaks.end(), [](const peak& a, const peak& b)->bool {
			return a.left < b.left;
		});
	}

	void invert_y(DVECTOR &x_in_out, DVECTOR &y_in_out)
	{
		for (auto i = y_in_out.begin(), _end_ = y_in_out.end(); i != _end_; ++i)
			*i = -*i;
	}

	//Assumes that peaks are sorted (improved performance)
	//Calculates baseline ignoring points lying inside peaks
	double find_baseline_by_median(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks)
	{
		if (xs.size() != ys.size())
			return DBL_MAX;
		DVECTOR selected_y;
#ifndef _USE_DEQUE
		selected_y.reserve(ys.size());
#endif
		std::size_t pk = 0;
		std::size_t pk_end_ = peaks.size();
		for (std::size_t i = 0, i_end_ = xs.size(); i != i_end_; ++i) {
			if ((pk >= pk_end_) || (xs[i] < peaks[pk].left)) {
				selected_y.push_back(ys[i]);
				continue;
			}
			if (xs[i] >= peaks[pk].left && xs[i] <= peaks[pk].right) {
				continue;
			}
			if (xs[i] > peaks[pk].right) { //Move to the next peak and check xs[i] against it
				++pk;
				--i;
			}
		}
		std::sort(selected_y.begin(), selected_y.end());
		int ind = selected_y.size();
		if (0 == ind)
			return DBL_MAX;
		if (0 == ind % 2)
			return selected_y[ind / 2];
		else
			return 0.5*(selected_y[ind / 2] + selected_y[1 + (ind / 2)]);
	}

	double find_baseline_by_median(DVECTOR &xs, DVECTOR &ys)
	{
		if (xs.size() != ys.size())
			return DBL_MAX;
		DVECTOR selected_y=ys;
		std::sort(selected_y.begin(), selected_y.end());
		int ind = selected_y.size();
		if (0 == ind)
			return DBL_MAX;
		if (0 == ind % 2)
			return selected_y[ind / 2];
		else
			return 0.5*(selected_y[ind / 2] + selected_y[1 + (ind / 2)]);
	}

	double find_baseline_by_integral(double approx_baseline, DVECTOR &xs, DVECTOR &ys)
	{
		if ((xs.size() <= 1) || (xs.size() != ys.size()))
			return DBL_MAX;
		double val = 0;
		double dx = xs.back()-*(xs.begin());
		if (0 == dx)
			return DBL_MAX;
		integrate(xs, ys, val, xs.begin(), --xs.end(), *(++xs.begin()) - *(xs.begin()), approx_baseline);
		return approx_baseline + val/dx;
	}

	//Assumes that peaks are sorted (improved performance)
	//Calculates baseline ignoring points lying inside peaks
	//Approximate baseline is required for enhanced precision
	double find_baseline_by_integral(double approx_baseline, DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks)
	{
		if ((xs.size() != ys.size()) || (xs.size() <= 1))
			return DBL_MAX;
		long double Sum_dx = 0, Sum_int = 0;
		bool start_point = false; //for continuous range (first accepted point)
		std::size_t pk = 0;
		std::size_t pk_end_ = peaks.size();
		for (std::size_t i = 0, i_end_ = xs.size(); i != i_end_; ++i) {
			long double dx;
			if (pk >= pk_end_ || xs[i] < peaks[pk].left) {
				if (start_point)
					dx = (i == i_end_ - 1 ? 0 : 0.5*(xs[i + 1] - xs[i]));
				else
					dx = (i == i_end_ - 1 ? 0 : 0.5*(xs[i + 1] - xs[i])) + (i == 0 ? 0 : 0.5*(xs[i] - xs[i-1]));
				start_point = false;
				Sum_dx += dx;
				Sum_int += (ys[i] - approx_baseline) * dx;
				continue;
			}
			if (xs[i] >= peaks[pk].left && xs[i] <= peaks[pk].right) {
				start_point = false;
				continue;
			}
			if (xs[i] > peaks[pk].right) { //Move to the next peak and check xs[i] against it
				start_point = true;
				++pk;
				--i;
			}
		}
		if (0 == Sum_dx)
			return DBL_MAX;
		return approx_baseline + (Sum_int / Sum_dx);
	}

	void find_baseline_by_ROOT(DVECTOR &xs, DVECTOR &ys, DVECTOR &ys_out,
		int numberIterations, int direction, int filterOrder, bool smoothing, int smoothWindow, bool compton, int sparse)
	{
		ys_out = ys;
		find_background_v_0(ys_out, numberIterations, direction, filterOrder, smoothing, smoothWindow, compton, sparse);
	}

	//My modification of CERN ROOT's TSpectrum->Background() function
	//First, it works directly on std container (vector or deque)
	//Second, sparse parameter was added to speed up the algorithm (skip every (sparse-1) points
	//in the main loop). sparse = 1 is equivalent to TSpectrum->Background()
	std::string find_background_v_0(DVECTOR& spectrum,
		int numberIterations,
		int direction, int filterOrder,
		bool smoothing, int smoothWindow,
		bool compton, int sparse)
	{

		int i, j, w, bw, b1, b2, priz;
		double a, b, c, d, e, yb1, yb2, ai, av, men, b4, c4, d4, e4, b6, c6, d6, e6, f6, g6, b8, c8, d8, e8, f8, g8, h8, i8;
		int ssize = spectrum.size();
		if (ssize <= 0)
			return "Wrong Parameters";
		if (numberIterations < 1)
			return "Width of Clipping Window Must Be Positive";
		if (ssize < 2 * numberIterations + 1)
			return "Too Large Clipping Window";
		if (smoothing == true && smoothWindow != TSpectrum::kBackSmoothing3 && smoothWindow != TSpectrum::kBackSmoothing5 && smoothWindow != TSpectrum::kBackSmoothing7 && smoothWindow != TSpectrum::kBackSmoothing9 && smoothWindow != TSpectrum::kBackSmoothing11 && smoothWindow != TSpectrum::kBackSmoothing13 && smoothWindow != TSpectrum::kBackSmoothing15)
			return "Incorrect width of smoothing window";
		if (sparse < 1)
			return "Sparse Must Be Positive";
		if (sparse >= numberIterations)
			return "Sparse must be Less than the width of clipping window";
		double *working_space = new double[2 * ssize];
		for (i = 0; i < ssize; i++){
			working_space[i] = spectrum[i];
			working_space[i + ssize] = spectrum[i];
		}
		bw = (smoothWindow - 1) / 2;
		if (direction == TSpectrum::kBackIncreasingWindow)
			i = 1;
		else if (direction == TSpectrum::kBackDecreasingWindow)
			i = numberIterations;
		if (filterOrder == TSpectrum::kBackOrder2) {
			do {
				for (j = i; j < ssize - i; ++j) {
					if (smoothing == false){
						a = working_space[ssize + j];
						b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
						if (b < a)
							a = b;
						working_space[j] = a;
					}

					else if (smoothing == true){
						a = working_space[ssize + j];
						av = 0;
						men = 0;
						for (w = j - bw; w <= j + bw; w++){
							if (w >= 0 && w < ssize){
								av += working_space[ssize + w];
								men += 1;
							}
						}
						av = av / men;
						b = 0;
						men = 0;
						for (w = j - i - bw; w <= j - i + bw; w++){
							if (w >= 0 && w < ssize){
								b += working_space[ssize + w];
								men += 1;
							}
						}
						b = b / men;
						c = 0;
						men = 0;
						for (w = j + i - bw; w <= j + i + bw; w++){
							if (w >= 0 && w < ssize){
								c += working_space[ssize + w];
								men += 1;
							}
						}
						c = c / men;
						b = (b + c) / 2;
						if (b < a)
							av = b;
						working_space[j] = av;
					}
				}
				for (j = i; j < ssize - i; j++)
					working_space[ssize + j] = working_space[j];
				if (direction == TSpectrum::kBackIncreasingWindow)
					i += sparse;
				else if (direction == TSpectrum::kBackDecreasingWindow)
					i -= sparse;
			} while ((direction == TSpectrum::kBackIncreasingWindow && i <= numberIterations) || (direction == TSpectrum::kBackDecreasingWindow && i >= 1));
		}

		else if (filterOrder == TSpectrum::kBackOrder4) {
			do{
				for (j = i; j < ssize - i; j++) {
					if (smoothing == false){
						a = working_space[ssize + j];
						b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
						c = 0;
						ai = i / 2;
						c -= working_space[ssize + j - (int)(2 * ai)] / 6;
						c += 4 * working_space[ssize + j - (int)ai] / 6;
						c += 4 * working_space[ssize + j + (int)ai] / 6;
						c -= working_space[ssize + j + (int)(2 * ai)] / 6;
						if (b < c)
							b = c;
						if (b < a)
							a = b;
						working_space[j] = a;
					}

					else if (smoothing == true){
						a = working_space[ssize + j];
						av = 0;
						men = 0;
						for (w = j - bw; w <= j + bw; w++){
							if (w >= 0 && w < ssize){
								av += working_space[ssize + w];
								men += 1;
							}
						}
						av = av / men;
						b = 0;
						men = 0;
						for (w = j - i - bw; w <= j - i + bw; w++){
							if (w >= 0 && w < ssize){
								b += working_space[ssize + w];
								men += 1;
							}
						}
						b = b / men;
						c = 0;
						men = 0;
						for (w = j + i - bw; w <= j + i + bw; w++){
							if (w >= 0 && w < ssize){
								c += working_space[ssize + w];
								men += 1;
							}
						}
						c = c / men;
						b = (b + c) / 2;
						ai = i / 2;
						b4 = 0, men = 0;
						for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								b4 += working_space[ssize + w];
								men += 1;
							}
						}
						b4 = b4 / men;
						c4 = 0, men = 0;
						for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								c4 += working_space[ssize + w];
								men += 1;
							}
						}
						c4 = c4 / men;
						d4 = 0, men = 0;
						for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								d4 += working_space[ssize + w];
								men += 1;
							}
						}
						d4 = d4 / men;
						e4 = 0, men = 0;
						for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								e4 += working_space[ssize + w];
								men += 1;
							}
						}
						e4 = e4 / men;
						b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
						if (b < b4)
							b = b4;
						if (b < a)
							av = b;
						working_space[j] = av;
					}
				}
				for (j = i; j < ssize - i; j++)
					working_space[ssize + j] = working_space[j];
				if (direction == TSpectrum::kBackIncreasingWindow)
					i += sparse;
				else if (direction == TSpectrum::kBackDecreasingWindow)
					i -= sparse;
			} while ((direction == TSpectrum::kBackIncreasingWindow && i <= numberIterations) || (direction == TSpectrum::kBackDecreasingWindow && i >= 1));
		}

		else if (filterOrder == TSpectrum::kBackOrder6) {
			do {
				for (j = i; j < ssize - i; j++) {
					if (smoothing == false){
						a = working_space[ssize + j];
						b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
						c = 0;
						ai = i / 2;
						c -= working_space[ssize + j - (int)(2 * ai)] / 6;
						c += 4 * working_space[ssize + j - (int)ai] / 6;
						c += 4 * working_space[ssize + j + (int)ai] / 6;
						c -= working_space[ssize + j + (int)(2 * ai)] / 6;
						d = 0;
						ai = i / 3;
						d += working_space[ssize + j - (int)(3 * ai)] / 20;
						d -= 6 * working_space[ssize + j - (int)(2 * ai)] / 20;
						d += 15 * working_space[ssize + j - (int)ai] / 20;
						d += 15 * working_space[ssize + j + (int)ai] / 20;
						d -= 6 * working_space[ssize + j + (int)(2 * ai)] / 20;
						d += working_space[ssize + j + (int)(3 * ai)] / 20;
						if (b < d)
							b = d;
						if (b < c)
							b = c;
						if (b < a)
							a = b;
						working_space[j] = a;
					}

					else if (smoothing == true){
						a = working_space[ssize + j];
						av = 0;
						men = 0;
						for (w = j - bw; w <= j + bw; w++){
							if (w >= 0 && w < ssize){
								av += working_space[ssize + w];
								men += 1;
							}
						}
						av = av / men;
						b = 0;
						men = 0;
						for (w = j - i - bw; w <= j - i + bw; w++){
							if (w >= 0 && w < ssize){
								b += working_space[ssize + w];
								men += 1;
							}
						}
						b = b / men;
						c = 0;
						men = 0;
						for (w = j + i - bw; w <= j + i + bw; w++){
							if (w >= 0 && w < ssize){
								c += working_space[ssize + w];
								men += 1;
							}
						}
						c = c / men;
						b = (b + c) / 2;
						ai = i / 2;
						b4 = 0, men = 0;
						for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								b4 += working_space[ssize + w];
								men += 1;
							}
						}
						b4 = b4 / men;
						c4 = 0, men = 0;
						for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								c4 += working_space[ssize + w];
								men += 1;
							}
						}
						c4 = c4 / men;
						d4 = 0, men = 0;
						for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								d4 += working_space[ssize + w];
								men += 1;
							}
						}
						d4 = d4 / men;
						e4 = 0, men = 0;
						for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								e4 += working_space[ssize + w];
								men += 1;
							}
						}
						e4 = e4 / men;
						b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
						ai = i / 3;
						b6 = 0, men = 0;
						for (w = j - (int)(3 * ai) - bw; w <= j - (int)(3 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								b6 += working_space[ssize + w];
								men += 1;
							}
						}
						b6 = b6 / men;
						c6 = 0, men = 0;
						for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								c6 += working_space[ssize + w];
								men += 1;
							}
						}
						c6 = c6 / men;
						d6 = 0, men = 0;
						for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								d6 += working_space[ssize + w];
								men += 1;
							}
						}
						d6 = d6 / men;
						e6 = 0, men = 0;
						for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								e6 += working_space[ssize + w];
								men += 1;
							}
						}
						e6 = e6 / men;
						f6 = 0, men = 0;
						for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								f6 += working_space[ssize + w];
								men += 1;
							}
						}
						f6 = f6 / men;
						g6 = 0, men = 0;
						for (w = j + (int)(3 * ai) - bw; w <= j + (int)(3 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								g6 += working_space[ssize + w];
								men += 1;
							}
						}
						g6 = g6 / men;
						b6 = (b6 - 6 * c6 + 15 * d6 + 15 * e6 - 6 * f6 + g6) / 20;
						if (b < b6)
							b = b6;
						if (b < b4)
							b = b4;
						if (b < a)
							av = b;
						working_space[j] = av;
					}
				}
				for (j = i; j < ssize - i; j++)
					working_space[ssize + j] = working_space[j];
				if (direction == TSpectrum::kBackIncreasingWindow)
					i += sparse;
				else if (direction == TSpectrum::kBackDecreasingWindow)
					i -= sparse;
			} while ((direction == TSpectrum::kBackIncreasingWindow && i <= numberIterations) || (direction == TSpectrum::kBackDecreasingWindow && i >= 1));
		}

		else if (filterOrder == TSpectrum::kBackOrder8) {
			do{
				for (j = i; j < ssize - i; j++) {
					if (smoothing == false){
						a = working_space[ssize + j];
						b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
						c = 0;
						ai = i / 2;
						c -= working_space[ssize + j - (int)(2 * ai)] / 6;
						c += 4 * working_space[ssize + j - (int)ai] / 6;
						c += 4 * working_space[ssize + j + (int)ai] / 6;
						c -= working_space[ssize + j + (int)(2 * ai)] / 6;
						d = 0;
						ai = i / 3;
						d += working_space[ssize + j - (int)(3 * ai)] / 20;
						d -= 6 * working_space[ssize + j - (int)(2 * ai)] / 20;
						d += 15 * working_space[ssize + j - (int)ai] / 20;
						d += 15 * working_space[ssize + j + (int)ai] / 20;
						d -= 6 * working_space[ssize + j + (int)(2 * ai)] / 20;
						d += working_space[ssize + j + (int)(3 * ai)] / 20;
						e = 0;
						ai = i / 4;
						e -= working_space[ssize + j - (int)(4 * ai)] / 70;
						e += 8 * working_space[ssize + j - (int)(3 * ai)] / 70;
						e -= 28 * working_space[ssize + j - (int)(2 * ai)] / 70;
						e += 56 * working_space[ssize + j - (int)ai] / 70;
						e += 56 * working_space[ssize + j + (int)ai] / 70;
						e -= 28 * working_space[ssize + j + (int)(2 * ai)] / 70;
						e += 8 * working_space[ssize + j + (int)(3 * ai)] / 70;
						e -= working_space[ssize + j + (int)(4 * ai)] / 70;
						if (b < e)
							b = e;
						if (b < d)
							b = d;
						if (b < c)
							b = c;
						if (b < a)
							a = b;
						working_space[j] = a;
					}

					else if (smoothing == true){
						a = working_space[ssize + j];
						av = 0;
						men = 0;
						for (w = j - bw; w <= j + bw; w++){
							if (w >= 0 && w < ssize){
								av += working_space[ssize + w];
								men += 1;
							}
						}
						av = av / men;
						b = 0;
						men = 0;
						for (w = j - i - bw; w <= j - i + bw; w++){
							if (w >= 0 && w < ssize){
								b += working_space[ssize + w];
								men += 1;
							}
						}
						b = b / men;
						c = 0;
						men = 0;
						for (w = j + i - bw; w <= j + i + bw; w++){
							if (w >= 0 && w < ssize){
								c += working_space[ssize + w];
								men += 1;
							}
						}
						c = c / men;
						b = (b + c) / 2;
						ai = i / 2;
						b4 = 0, men = 0;
						for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								b4 += working_space[ssize + w];
								men += 1;
							}
						}
						b4 = b4 / men;
						c4 = 0, men = 0;
						for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								c4 += working_space[ssize + w];
								men += 1;
							}
						}
						c4 = c4 / men;
						d4 = 0, men = 0;
						for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								d4 += working_space[ssize + w];
								men += 1;
							}
						}
						d4 = d4 / men;
						e4 = 0, men = 0;
						for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								e4 += working_space[ssize + w];
								men += 1;
							}
						}
						e4 = e4 / men;
						b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
						ai = i / 3;
						b6 = 0, men = 0;
						for (w = j - (int)(3 * ai) - bw; w <= j - (int)(3 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								b6 += working_space[ssize + w];
								men += 1;
							}
						}
						b6 = b6 / men;
						c6 = 0, men = 0;
						for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								c6 += working_space[ssize + w];
								men += 1;
							}
						}
						c6 = c6 / men;
						d6 = 0, men = 0;
						for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								d6 += working_space[ssize + w];
								men += 1;
							}
						}
						d6 = d6 / men;
						e6 = 0, men = 0;
						for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								e6 += working_space[ssize + w];
								men += 1;
							}
						}
						e6 = e6 / men;
						f6 = 0, men = 0;
						for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								f6 += working_space[ssize + w];
								men += 1;
							}
						}
						f6 = f6 / men;
						g6 = 0, men = 0;
						for (w = j + (int)(3 * ai) - bw; w <= j + (int)(3 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								g6 += working_space[ssize + w];
								men += 1;
							}
						}
						g6 = g6 / men;
						b6 = (b6 - 6 * c6 + 15 * d6 + 15 * e6 - 6 * f6 + g6) / 20;
						ai = i / 4;
						b8 = 0, men = 0;
						for (w = j - (int)(4 * ai) - bw; w <= j - (int)(4 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								b8 += working_space[ssize + w];
								men += 1;
							}
						}
						b8 = b8 / men;
						c8 = 0, men = 0;
						for (w = j - (int)(3 * ai) - bw; w <= j - (int)(3 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								c8 += working_space[ssize + w];
								men += 1;
							}
						}
						c8 = c8 / men;
						d8 = 0, men = 0;
						for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								d8 += working_space[ssize + w];
								men += 1;
							}
						}
						d8 = d8 / men;
						e8 = 0, men = 0;
						for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								e8 += working_space[ssize + w];
								men += 1;
							}
						}
						e8 = e8 / men;
						f8 = 0, men = 0;
						for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
							if (w >= 0 && w < ssize){
								f8 += working_space[ssize + w];
								men += 1;
							}
						}
						f8 = f8 / men;
						g8 = 0, men = 0;
						for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								g8 += working_space[ssize + w];
								men += 1;
							}
						}
						g8 = g8 / men;
						h8 = 0, men = 0;
						for (w = j + (int)(3 * ai) - bw; w <= j + (int)(3 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								h8 += working_space[ssize + w];
								men += 1;
							}
						}
						h8 = h8 / men;
						i8 = 0, men = 0;
						for (w = j + (int)(4 * ai) - bw; w <= j + (int)(4 * ai) + bw; w++){
							if (w >= 0 && w < ssize){
								i8 += working_space[ssize + w];
								men += 1;
							}
						}
						i8 = i8 / men;
						b8 = (-b8 + 8 * c8 - 28 * d8 + 56 * e8 - 56 * f8 - 28 * g8 + 8 * h8 - i8) / 70;
						if (b < b8)
							b = b8;
						if (b < b6)
							b = b6;
						if (b < b4)
							b = b4;
						if (b < a)
							av = b;
						working_space[j] = av;
					}
				}
				for (j = i; j < ssize - i; j++)
					working_space[ssize + j] = working_space[j];
				if (direction == TSpectrum::kBackIncreasingWindow)
					i += sparse;
				else if (direction == TSpectrum::kBackDecreasingWindow)
					i -= sparse;
			} while ((direction == TSpectrum::kBackIncreasingWindow && i <= numberIterations) || (direction == TSpectrum::kBackDecreasingWindow && i >= 1));
		}

		if (compton == true) {
			for (i = 0, b2 = 0; i < ssize; i++){
				b1 = b2;
				a = working_space[i], b = spectrum[i];
				j = i;
				if (std::fabs(a - b) >= 1) {
					b1 = i - 1;
					if (b1 < 0)
						b1 = 0;
					yb1 = working_space[b1];
					for (b2 = b1 + 1, c = 0, priz = 0; priz == 0 && b2 < ssize; b2++){
						a = working_space[b2], b = spectrum[b2];
						c = c + b - yb1;
						if (std::fabs(a - b) < 1) {
							priz = 1;
							yb2 = b;
						}
					}
					if (b2 == ssize)
						b2 -= 1;
					yb2 = working_space[b2];
					if (yb1 <= yb2){
						for (j = b1, c = 0; j <= b2; j++){
							b = spectrum[j];
							c = c + b - yb1;
						}
						if (c > 1){
							c = (yb2 - yb1) / c;
							for (j = b1, d = 0; j <= b2 && j < ssize; j++){
								b = spectrum[j];
								d = d + b - yb1;
								a = c * d + yb1;
								working_space[ssize + j] = a;
							}
						}
					}

					else{
						for (j = b2, c = 0; j >= b1; j--){
							b = spectrum[j];
							c = c + b - yb2;
						}
						if (c > 1){
							c = (yb1 - yb2) / c;
							for (j = b2, d = 0; j >= b1 && j >= 0; j--){
								b = spectrum[j];
								d = d + b - yb2;
								a = c * d + yb2;
								working_space[ssize + j] = a;
							}
						}
					}
					i = b2;
				}
			}
		}

		for (j = 0; j < ssize; ++j)
			spectrum[j] = working_space[ssize + j];
		delete[]working_space;
		return "";
	}

	//Does not assume that xs are equidistant
	//Returns empty y_out in case of error
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &y_out, double baseline)
	{
		if ((xs.size() != ys.size()) || (xs.size() <= 1)) {
			DVECTOR().swap(y_out);
			return;
		}
		y_out.resize(ys.size());
		double prev = 0, dx;
		for (std::size_t i = 0, i_end_ = xs.size(); (i != i_end_); ++i) {
			dx = ((i == (i_end_ - 1)) ? 0 : 0.5*(xs[i + 1] - xs[i])) +
				((i == 0) ? 0 : 0.5*(xs[i] - xs[i - 1]));
			prev = dx*(ys[i] - baseline) + prev;
			y_out[i] = prev;
		}
	}

	//Does not assume that xs are equidistant
	//Calculates integral as well as its standard deviation
	//Returns empty y_out and y_out_err in case of error
	//Output standard deviation may be wrong due to limited precision
	void integrate_with_variance(DVECTOR &xs, DVECTOR &ys, DVECTOR &ys_err,  DVECTOR &y_out, DVECTOR &y_out_err, double baseline)
	{
		if ((xs.size() != ys.size()) || (xs.size() <= 1)||(xs.size()!=ys_err.size())) {
			DVECTOR().swap(y_out);
			DVECTOR().swap(y_out_err);
			return;
		}
		y_out.resize(ys.size());
		y_out_err.resize(ys.size());
		long double variance = 0;
		double prev = 0, dx;
		for (std::size_t i = 0, i_end_ = xs.size(); i != i_end_; ++i) {
			dx = ((i == (i_end_ - 1)) ? 0 : 0.5*(xs[i+1] - xs[i])) +
				((i == 0) ? 0 : 0.5*(xs[i] - xs[i-1]));
			prev = dx*(ys[i] - baseline) + prev;
			variance += (long double) dx*dx*ys_err[i]*ys_err[i];
			y_out[i]=prev;
			y_out_err[i] = sqrt(variance);
		}
	}

	//Assumes that xs are equidistant
	//Returns empty y_out in case of error
	//Output standard deviation may be wrong due to limited precision
	void integrate(DVECTOR &ys, DVECTOR &y_out, double dx_hint, double baseline)
	{
		if (ys.size() == 0) {
			DVECTOR().swap(y_out);
			return;
		}
		y_out.resize(ys.size());
		double prev = 0;
		double dx;
		for (std::size_t i = 0, i_end_ = ys.size(); i != i_end_; ++i) {
			dx = ((i == (i_end_ - 1)) ? 0 : 0.5*dx_hint) +
				(i == 0 ? 0 : 0.5*dx_hint);
			prev = dx*(ys[i] - baseline) + prev;
			y_out[i] = prev;
		}
	}

	//Assumes that xs are equidistant
	//Carefully handles integration dx near left and right points and other special cases.
	//E.g. if (left +eps1)==(right-eps2)==xs[j] integral will be ys[j]*(eps2-eps1) assuming that each eps is far from neighbors xs
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out, double left, double right, double dx_hint, double baseline)
	{
		y_out.clear();
		x_out.clear();
		if ((xs.size() != ys.size()) || (xs.size() <= 1))
			return;
		double Left = std::min(left, right);
		double Right = std::max(left, right);
		std::size_t from, to;
		if (xs.back() < Left)
			return;
		if (xs.front() > Right)
			return;
		if (xs.front() >= Left)
			from = 0;
		else {
			from = (std::size_t)((Left - xs.front()) / dx_hint);
			if (from >= xs.size()) { //due to rounding errors
				boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Left);
				from = inds->second;
			} else
				if (xs[from] != Left) {
					if (xs[from] < Left && xs[from + 1] >= Left)
						from = from + 1;
					else {
						boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Left);
						from = inds->second;
					}
				}
		}
		if (xs.back() <= Right)
			to = xs.size();
		else {
			to = (std::size_t)((Right - xs.front()) / dx_hint);
			if (to >= xs.size()) { //due to rounding errors
				boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Right);
				to = inds->first + 1;
			} else
				if (xs[to] <= Right && xs[to + 1] > Right)
					++to;
				else {
					boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Right);
					to = inds->first + 1;
				}
		}
		long double prev = 0, dx;
		x_out.resize(to - from);
		y_out.resize(to - from);
		for (std::size_t i = from, i_end_ = xs.size(); i != to; ++i) {
			dx = (i == 0 ? 0 : (i == from ? xs[i] - Left : 0.5*dx_hint)) +
				(i == (i_end_ - 1) ? 0 : (i == (to - 1) ? Right - xs[i] : 0.5*dx_hint));
			prev = dx * (ys[i] - baseline) + prev;
			x_out[i - from] = (xs[i]);
			y_out[i - from] = prev;
		}
	}

	//Does not assume that xs are equidistant
	//Carefully handles integration dx near left and right points and other special cases.
	//E.g. if (left +eps1)==(right-eps2)==xs[j] integral will be ys[j]*(eps2-eps1) assuming that each eps is far from neighbors xs
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out, double left, double right, double baseline)
	{
		y_out.clear();
		x_out.clear();
		if ((xs.size() != ys.size()) || (xs.size() <= 1))
			return;
		double Left = std::min(left, right);
		double Right = std::max(left, right);
		std::size_t from, to;
		if (xs.back() < Left)
			return;
		if (xs.front() > Right)
			return;
		if (xs.front() >= Left)
			from = 0;
		else {
			boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Left);
			from = inds->second;
		}
		if (xs.back() <= Right)
			to = xs.size();
		else {
			boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Right);
			to = inds->first + 1;
		}
		long double prev = 0, dx;
		x_out.resize(to - from);
		y_out.resize(to - from);
		for (std::size_t i = from, i_end_ = xs.size(); i != to; ++i) {
			dx = (i == 0 ? 0 : (i == from ? xs[i] - Left : 0.5*(xs[i] - xs[i-1]))) +
				(i == (i_end_ - 1) ? 0 : (i == (to - 1) ? Right - xs[i] : 0.5*(xs[i+1] - xs[i])));
			prev = dx*(ys[i] - baseline) + prev;
			x_out[i - from] = (xs[i]);
			y_out[i - from] = prev;
		}
	}

	//Assumes that xs are equidistant
	//Integrates in the range [left, right), i.e.
	//Use xs.begin(), xs.end() to integrate whole signal
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out, DITERATOR left, DITERATOR right, double dx_hint, double baseline)
	{
		if (xs.size() != ys.size()) {
			x_out.clear();
			y_out.clear();
			return;
		}
		std::size_t from = std::min(left - xs.begin(), right - xs.begin());
		std::size_t to = std::max(left - xs.begin(), right - xs.begin());
		x_out.resize(to - from);
		y_out.resize(to - from);
		long double prev = 0, dx;
		for (std::size_t i = from, i_end_ = xs.size(); i != to; ++i) {
			dx = (i == 0 ? 0 : (i == from ? 0 : 0.5*dx_hint)) +
				(i == (i_end_ - 1) ? 0 : (i == (to - 1) ? 0 : 0.5*dx_hint));
			prev += dx * (ys[i] - baseline);
			x_out[i - from] = (xs[i]);
			y_out[i - from] = prev;
		}
	}

	//Integrates in the range [left, right), i.e.
	//Use xs.begin(), xs.end() to integrate whole signal
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out, DITERATOR left, DITERATOR right, double baseline)
	{
		if (xs.size()!=ys.size() || xs.size()<=1) {
			x_out.clear();
			y_out.clear();
			return;
		}
		std::size_t from = std::min(left - xs.begin(), right - xs.begin());
		std::size_t to = std::max(left - xs.begin(), right - xs.begin());
		x_out.resize(to - from);
		y_out.resize(to - from);
		long double prev = 0, dx;
		for (std::size_t i = from, i_end_ = xs.size(); i != to; ++i) {
			dx = (i == 0 ? 0 : (i == from ? 0 : 0.5*(xs[i] - xs[i - 1]))) +
				(i == (i_end_ - 1) ? 0 : (i == (to - 1) ? 0 : 0.5*(xs[i + 1] - xs[i])));
			prev += dx * (ys[i] - baseline);
			x_out[i - from] = (xs[i]);
			y_out[i - from] = prev;
		}
	}

	//Assumes that xs are equidistant
	//Integrates in the range [left, right), i.e.
	//Use xs.begin(), xs.end() to integrate whole signal
	void integrate(DVECTOR &xs, DVECTOR &ys, double &y_out, DITERATOR left, DITERATOR right, double dx_hint, double baseline)
	{
		y_out = 0;
		if (xs.size() != ys.size()) {
			return;
		}
		std::size_t from = std::min(left - xs.begin(), right - xs.begin());
		std::size_t to = std::max(left - xs.begin(), right - xs.begin());
		double dx;
		for (std::size_t i = from, i_end_ = xs.size(); i != to; ++i) {
			dx = (i == 0 ? 0 : (i == from ? 0 : 0.5*dx_hint)) +
				(i == (i_end_ - 1) ? 0 : (i == (to - 1) ? 0 : 0.5*dx_hint));
			y_out += dx * (ys[i] - baseline);
		}
	}

	//Does not assume that xs are equidistant
	//Integrates in the range [left, right), i.e.
	//Use xs.begin(), xs.end() to integrate whole signal
	void integrate(DVECTOR &xs, DVECTOR &ys, double &y_out, DITERATOR left, DITERATOR right, double baseline)
	{
		y_out = 0;
		if (xs.size() != ys.size()) {
			return;
		}
		std::size_t from = std::min(left - xs.begin(), right - xs.begin());
		std::size_t to = std::max(left - xs.begin(), right - xs.begin());
		double dx;
		for (std::size_t i = from, i_end_ = xs.size(); i != to; ++i) {
			dx = (i == 0 ? 0 : (i == from ? 0 : 0.5*(xs[i] - xs[i-1]))) +
				(i == (i_end_ - 1) ? 0 : (i == (to - 1) ? 0 : 0.5*(xs[i + 1] - xs[i])));
			y_out += dx * (ys[i] - baseline);
		}
	}

	//Assumes that xs are equidistant
	//Crops signal to [x_left, x_right]
	void apply_time_limits(DVECTOR &xs, DVECTOR &ys, double x_left, double x_right, double dx_hint)
	{
		if (xs.size() != ys.size())
			return;
		double Left = std::min(x_left, x_right);
		double Right = std::max(x_left, x_right);
		std::size_t from, to;
		if (xs.back() < Left) {
			xs.clear();
			ys.clear();
			return;
		}
		if (xs.front() > Right) {
			xs.clear();
			ys.clear();
			return;
		}
		if (xs.front() >= Left)
			from = 0;
		else {
			from = (std::size_t)((Left - xs.front()) / dx_hint);
			if (from >= xs.size()) { //due to rounding errors
				boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Left);
				from = inds->second;
			} else
				if (xs[from] != Left) {
					if (xs[from] < Left && xs[from + 1] >= Left)
						from = from + 1;
					else { //dx hind didn't work
						boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Left);
						from = inds->second;
					}
				}
		}
		if (xs.back() <= Right)
			to = xs.size();
		else {
			to = (std::size_t)((Right - xs.front()) / dx_hint);
			if (to >= xs.size()) { //due to rounding errors
				boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Right);
				to = inds->first + 1;
			} else
				if (xs[to] <= Right && xs[to + 1] > Right)
					++to;
				else { //dx hind didn't work
					boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Right);
					to = inds->first + 1;
				}
		}

		xs.erase(xs.begin() + to, xs.end());
		xs.erase(xs.begin(), xs.begin() + from);
		ys.erase(ys.begin() + to, ys.end());
		ys.erase(ys.begin(), ys.begin() + from);
	}

	//Crops signal to [x_left, x_right]
	void apply_time_limits(DVECTOR &xs, DVECTOR &ys, double x_left, double x_right)
	{
		if (xs.size() != ys.size())
			return;
		double Left = std::min(x_left, x_right);
		double Right = std::max(x_left, x_right);
		std::size_t from, to;
		if (xs.back() < Left) {
			xs.clear();
			ys.clear();
			return;
		}
		if (xs.front() > Right) {
			xs.clear();
			ys.clear();
			return;
		}
		if (xs.front() >= Left)
			from = 0;
		else {
			boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Left);
			from = inds->second;
		}
		if (xs.back() <= Right)
			to = xs.size();
		else {
			boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(xs, Right);
			to = inds->first + 1;
		}

		xs.erase(xs.begin() + to, xs.end());
		xs.erase(xs.begin(), xs.begin() + from);
		ys.erase(ys.begin() + to, ys.end());
		ys.erase(ys.begin(), ys.begin() + from);
	}

	//Assumes that peaks are sorted (improved performance)
	//Removes signal corresponding to peaks
	void exclude_peaks(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks)
	{
		if (xs.size() != ys.size())
			return;
		DVECTOR x_out, y_out;
#ifndef _USE_DEQUE
		x_out.reserve(xs.size());
		y_out.reserve(ys.size());
#endif
		bool start_point = false; //for continuous range (first accepted point)
		std::size_t pk = 0;
		std::size_t pk_end_ = peaks.size();
		for (std::size_t i = 0, i_end_ = xs.size(); i != i_end_; ++i) {
			if ((pk >= pk_end_) || (xs[i] < peaks[pk].left)) {
				x_out.push_back(ys[i]);
				y_out.push_back(xs[i]);
				continue;
			}
			if (xs[i] >= peaks[pk].left && xs[i] <= peaks[pk].right) {
				continue;
			}
			if (xs[i] > peaks[pk].right) { //Move to the next peak and check xs[i] against it
				++pk;
				--i;
			}
		}
		xs = x_out;
		ys = y_out;
	}

	//Splits signal into found peaks. Even peaks (counting from 0) are stored in output.first and the odd ones - in output.second.
	//Both output.first and output.second store x-y pair arrays corresponding to separate peaks.
	//That is output.first.first[0] and output.first.second[0] contain xs and ys for the first (even) peak.
	//        output.second.first[0] and output.second.second[0] contain xs and ys for the second (odd) peak.
	//        output.first.first[1] and output.first.second[1] contain xs and ys for the third (even) peak.
	//... and so forth.
	void select_peaks(const DVECTOR &xs, const DVECTOR &ys, STD_CONT<peak> &peaks, double dx_hint, double baseline,
		std::pair<std::pair<STD_CONT<DVECTOR>, STD_CONT<DVECTOR>>, std::pair<STD_CONT<DVECTOR>, STD_CONT<DVECTOR>>>  &output)
	{
		STD_CONT<DVECTOR>().swap(output.first.first);//xs for even peaks
		STD_CONT<DVECTOR>().swap(output.first.second);//ys for even peaks
		STD_CONT<DVECTOR>().swap(output.second.first);//xs for odd peaks
		STD_CONT<DVECTOR>().swap(output.second.second);//ys for odd peaks
		if (xs.size() != ys.size())
			return;
		sort_peaks(peaks);
		bool even = true;
		bool in_peak = false;
		std::size_t pk = 0;
		std::size_t pk_end_ = peaks.size();
		for (std::size_t i = 0, i_end_ = xs.size(); i != i_end_; ++i) {
			if (pk >= pk_end_)
				break;
			bool was_in_peak = in_peak;
			if (xs[i] < peaks[pk].left)
				continue;
			if (xs[i] >= peaks[pk].left && xs[i] <= peaks[pk].right) {
				in_peak = true;
				if (!was_in_peak) {
					if (even) {
						output.first.first.push_back(DVECTOR());
						output.first.second.push_back(DVECTOR());
						output.first.first.back().push_back(xs[i] - dx_hint / 20);
						output.first.second.back().push_back(baseline);
					} else {
						output.second.first.push_back(DVECTOR());
						output.second.second.push_back(DVECTOR());
						output.second.first.back().push_back(xs[i] - dx_hint / 20);
						output.second.second.back().push_back(baseline);
					}
				}
				was_in_peak = in_peak;
				if (even) {
					output.first.first.back().push_back(xs[i]);
					output.first.second.back().push_back(ys[i]);
				} else {
					output.second.first.back().push_back(xs[i]);
					output.second.second.back().push_back(ys[i]);
				}
			}
			if (xs[i] >= peaks[pk].right) {
				in_peak = false;
				if (xs[i] == peaks[pk].right) {
					if (even) {
						output.first.first.back().push_back(xs[i]);
						output.first.second.back().push_back(ys[i]);
					} else {
						output.second.first.back().push_back(xs[i]);
						output.second.second.back().push_back(ys[i]);
					}
				}
				if (was_in_peak) {
					if (even) {
						output.first.first.back().push_back(output.first.first.back().back() + dx_hint / 20);
						output.first.second.back().push_back(baseline);
					} else {
						output.second.first.back().push_back(output.second.first.back().back() + dx_hint / 20);
						output.second.second.back().push_back(baseline);
					}
				}
				++pk;
				even = !even;
				--i;//!!! for the case when peaks[pk+1].left and >= xs[i] - need to include x[i] into peak as well
			}
		}
		if (in_peak) {
			if (even) {
				output.first.first.back().push_back(output.first.first.back().back() + dx_hint / 20);
				output.first.second.back().push_back(baseline);
			} else {
				output.second.first.back().push_back(output.second.first.back().back() + dx_hint / 20);
				output.second.second.back().push_back(baseline);
			}
		}
	}

	//Searches y maximum in [x_start, x_finish) range
	std::pair<double, DITERATOR> get_max(DVECTOR &xs, DVECTOR &ys, DITERATOR x_start, DITERATOR x_finish)
	{
		std::pair<double, DITERATOR> y_max(-DBL_MAX, xs.end());
		if (xs.size() != ys.size())
			return y_max;
		DITERATOR zero = xs.begin();
		DITERATOR _begin_ = std::min(x_start, x_finish);
		DITERATOR _end_ = std::max(x_start, x_finish);
		for (auto i = _begin_; i != _end_; ++i) {
			if (*(ys.begin() + (i - zero)) > y_max.first) {
				y_max.first = *(ys.begin() + (i - zero));
				y_max.second = i;
			}
		}
		return y_max;
	}

	std::pair<double, DITERATOR> get_max(DVECTOR &xs, DVECTOR &ys)
	{
		return get_max(xs, ys, xs.begin(), xs.end());
	}

	//Searches y minimum in [x_start, x_finish) range
	std::pair<double, DITERATOR> get_min(DVECTOR &xs, DVECTOR &ys, DITERATOR x_start, DITERATOR x_finish)
	{
		std::pair<double, DITERATOR> y_min(DBL_MAX, xs.end());
		if (xs.size() != ys.size())
			return y_min;
		DITERATOR zero = xs.begin();
		DITERATOR _begin_ = std::min(x_start, x_finish);
		DITERATOR _end_ = std::max(x_start, x_finish);
		for (auto i = _begin_; i != _end_; ++i) {
			if (*(ys.begin() + (i - zero)) < y_min.first) {
				y_min.first = *(ys.begin() + (i - zero));
				y_min.second = i;
			}
		}
		return y_min;
	}

	std::pair<double, DITERATOR> get_min(DVECTOR &xs, DVECTOR &ys)
	{
		return get_min(xs, ys, xs.begin(), xs.end());
	}

	//Old version of searching peaks by threshold (before 2019.12.06) rewritten in find_peaks_fine_v2 style (more readable and faster) 
	//Everything that is above threshold is counted as peak.
	//Peak fronts (ends, edges) are determined by intersection with threshold_edges (typically == baseline)
	//Assumes xs are equidistant
	//Returns sorted peaks
	void find_peaks_fine(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks, double baseline, double threshold, double threshold_edges)
	{
		if (threshold_edges >= threshold)
			threshold_edges = baseline;
		peaks.clear();
		if (xs.size() <= 1)
			return;
		double delta_x = *(++xs.begin()) - *(xs.begin());
		std::size_t x_peak_l = 0, x_peak_r = 0;
		const std::size_t _end_ = xs.size();
		while (x_peak_l != _end_) {//when no peak found x_peak_l is set to xs.size()
			x_peak_l = _end_;
			std::size_t search_from = x_peak_r; //next peak may include this point (i.e. one point can belong to two peaks)
			//find next peak (signal > threshold)
			std::size_t x_start = _end_; //x_start!=_end_ == peak found
			for (std::size_t i = search_from; i != _end_; ++i) {
				if (ys[i] >= threshold) {
					x_start = i;
					break;
				}
			}
			if (x_start == _end_)
				continue; //no more peaks found
			//find left peak front (signal < threshold_edges, but not x>=search_from)
			for (std::size_t i = x_start; i >= search_from && i < _end_; --i) {
				if (i == 0)
					continue;
				if ((ys[i] - threshold_edges)*(ys[i - 1] - threshold_edges) <= 0) {
					x_peak_l = i;
					break;
				}
			}
			if (x_peak_l == _end_)
				x_peak_l = search_from;
			//find right peak front (either intersection with threshold_edges or minimum between this peak and the next one (above threshold).
			//Take integrals in the same loop to improve performance
			peak pk;
			double A_weight = 0;
			pk.A = -DBL_MAX;
			pk.S = 0;
			pk.t = 0;
			std::size_t edge_intersection = _end_;
			for (std::size_t i = x_start; i != _end_; ++i) {
				if ((i == 0) || (i == x_start)) {
					A_weight += ys[i] - baseline;
					pk.t += (ys[i] - baseline)*xs[i];
					pk.A = std::max(pk.A, ys[i] - baseline);
					pk.S += (ys[i] - baseline) * delta_x *((i==0 ? 0.0 : 0.5) + (i == (_end_ -1) ? 0.0 : 0.5));
					continue;
				}
				if ((ys[i] - threshold_edges) <= 0 && (ys[i - 1] - threshold_edges) >= 0) { //intersection
					edge_intersection = i - 1;
					break;
				}
				A_weight += ys[i] - baseline;
				pk.t += (ys[i] - baseline)*xs[i];
				pk.A = std::max(pk.A, ys[i] - baseline);
				pk.S += (ys[i] - baseline) * delta_x *((i == 0 ? 0.0 : 0.5) + (i == (_end_ - 1) ? 0.0 : 0.5));
			}
			
			if (edge_intersection == _end_) //The case when signal ended before peak's right front was found
				x_peak_r = _end_ - 1;
			
			pk.left = xs[x_peak_l];
			pk.right = xs[x_peak_r];
			pk.t /= A_weight;

			if ((pk.S > 0) && (pk.A > 0) && (pk.right >= pk.left) && (pk.t >= 0))
				peaks.push_back(pk);
			++x_peak_r; //otherwise loop is stuck at one-point peaks
		}
	}

	//New version of searching peaks by threshold (after 2019.12.06)
	//Everything that is above threshold is counted as peak.
	//Peak fronts (ends, edges) are determined by intersection with threshold_edges (typically == baseline)
	//If dip between two peaks is larger than threshold_edges but less than threshold
	//this new algorithm separates them. Position of the minimum between them is taken as separating point.
	//Assumes xs are equidistant
	//Returns sorted peaks
	void find_peaks_fine_v2(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks, double baseline, double threshold, double threshold_edges)
	{
		if (threshold_edges >= threshold)
			threshold_edges = baseline;
		peaks.clear();
		if (xs.size() <= 1)
			return;
		double delta_x = *(++xs.begin()) - *(xs.begin());
		std::size_t x_peak_l = 0, x_peak_r = 0;
		const std::size_t _end_ = xs.size();
		while (x_peak_l != _end_) {//when no peak found x_peak_l is set to xs.size()
			x_peak_l = _end_;
			std::size_t search_from = x_peak_r; //next peak may include this point (i.e. one point can belong to two peaks)
			//find next peak (signal > threshold)
			std::size_t x_start = _end_; //x_start!=_end_ == peak found
			for (std::size_t i = search_from; i != _end_; ++i) {
				if (ys[i] >= threshold) {
					x_start = i;
					break;
				}
			}
			if (x_start == _end_)
				continue; //no more peaks found
			//find left peak front (signal < threshold_edges, but not x>=search_from)
			for (std::size_t i = x_start; i >= search_from && i<_end_; --i) {
				if (i == 0)
					continue;
				if ((ys[i] - threshold_edges)*(ys[i - 1] - threshold_edges) <= 0) {
					x_peak_l = i;
					break;
				}
			}
			if (x_peak_l == _end_)
				x_peak_l = search_from;
			//find right peak front (either intersection with threshold_edges or minimum between this peak and the next one (above threshold).
			//Can't take integrals and find maximum in the same loop because the integral range is unknown until completion of determining the right front
			bool shared_r = false;
			bool shared_l = x_peak_r == x_peak_l; //new peak's left slope is shared with previous peak's right
			std::size_t min = _end_, thresh_intersection_1 = _end_, thresh_intersection_2 = _end_;
			std::size_t edge_intersection = _end_;
			for (std::size_t i = x_start; i != _end_; ++i) {
				if ((i == 0) || (i == x_start)) {
					continue;
				}
				if ((ys[i] - threshold_edges) <= 0 && (ys[i - 1] - threshold_edges) >= 0) { //intersection
					edge_intersection = i - 1;
					break;
				}
				if (_end_==thresh_intersection_1) {
					if ((ys[i] - threshold) <= 0 && (ys[i - 1] - threshold) > 0)
						thresh_intersection_1 = i - 1;
				} else {
					if ((ys[i] - threshold) > 0 && (ys[i - 1] - threshold) <= 0) {
						thresh_intersection_2 = i - 1;
						break;
					}
				}
				if (_end_!=thresh_intersection_1) {
					if (min == _end_)
						min = i;
					else
						min = (ys[i]<=ys[min] ? i : min);
				}
			}
			
			if ((thresh_intersection_2 == _end_) && (edge_intersection == _end_)) { //The case when signal ended before peak's right front was found
				x_peak_r = _end_ - 1;
			} else {
				if (thresh_intersection_2!=_end_) {
					x_peak_r = min;
					shared_r = true;
				}
				else
					x_peak_r = edge_intersection;
			}
			peak pk;
			double A_weight = 0;
			pk.A = -DBL_MAX;
			pk.S = 0;
			pk.t = 0;
			pk.left = xs[x_peak_l];
			pk.right = xs[x_peak_r];
			for (std::size_t i = x_peak_l; i <= x_peak_r; ++i) {
				A_weight += ys[i]-baseline;
				pk.t += (ys[i]-baseline)*xs[i];
				pk.A = std::max(pk.A, ys[i] - baseline);
				if (i == x_peak_l) {
					if (x_peak_r == x_peak_l) {
						pk.S += (ys[i] - baseline) * delta_x *((shared_l ? 0.5 : 1.0) + (shared_r ? 0.5 : 1.0));
						continue;
					}
					pk.S += (ys[i] - baseline)*delta_x*(shared_l ? 0.5 : 1.0);
					continue;
				}
				if (i == x_peak_r) {
					pk.S += (ys[i] - baseline) * delta_x *(shared_r ? 0.5 : 1.0);
					continue;
				}
				pk.S += (ys[i] - baseline)*delta_x;
			}

			pk.t /= A_weight;

			if ((pk.S>0) && (pk.A>0)&&(pk.right>=pk.left)&&(pk.t>=0))
				peaks.push_back(pk);
			if (!shared_r)
				++x_peak_r; //otherwise loop is stuck at one-point peaks
		}
	}

	void subtract_baseline(DVECTOR &ys_in, double base_line)
	{
		for (auto i = ys_in.begin(), _end_ = ys_in.end(); i != _end_; ++i)
			*i -= base_line;
	}

	//Both vectors must have the same sizes
	void subtract_baseline(DVECTOR &ys_in, DVECTOR &base_ys)
	{
		if (ys_in.size() != base_ys.size())
			return;
		for (auto i = ys_in.begin(), j = base_ys.begin(), _end_ = ys_in.end(); i != _end_; ++i, ++j)
			*i -= *j;
	}
	
	//Required for subtracting ROOT's baseline which is calculated only for some range
	//Assumes that xs_in and base_xs have the same points except for maybe some points
	//See code for details
	void subtract_baseline(DVECTOR& xs_in, DVECTOR &ys_in, DVECTOR &base_xs, DVECTOR &base_ys, double baseline_baseline)
	{
		if ((xs_in.size() != ys_in.size()) || (base_xs.size() != base_ys.size()))
			return;
		std::size_t x_ind = 0;
		std::size_t base_ind = 0;
		std::size_t x_end_ = xs_in.size();
		std::size_t base_end_ = base_xs.size();
		while (x_ind != x_end_ && base_ind != base_end_) {
			if (xs_in[x_ind] == base_xs[base_ind]) {
				ys_in[x_ind] -= (base_ys[base_ind] - baseline_baseline);
				++x_ind;
				++base_ind;
				continue;
			}
			if (xs_in[x_ind] < base_xs[base_ind]) {
				++x_ind;
				continue;
			}
			if (xs_in[x_ind] > base_xs[base_ind]) {
				++base_ind;
				continue;
			}
		}
	}

};
