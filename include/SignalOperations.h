#ifndef SIGNAL_OPERATIONS_H
#define SIGNAL_OPERATIONS_H

#include <boost/optional.hpp>
#include "GlobalParameters.h"

namespace SignalOperations
{
	//Finds indices in xs bounding value val, i.e. val>xs[first] && val<xs[second]. If val==xs[first] returns (first, first) pair.
	//If val is outside the array, returns index to the nearest element to val.
	boost::optional<std::pair<std::size_t, std::size_t> > getX_indices(const DVECTOR& xs, double val);

	//Wrapper for std::sort. Sorts peaks by left fronts.
	//Peaks must not intersect, i.e. pk[i].right <= pk[i+1].left and pk[i].left<=pk[i].right for sorted pk array.
	void sort_peaks(STD_CONT<peak> &peaks);

	void invert_y(DVECTOR &x_in_out, DVECTOR &y_in_out);

	//Assumes that peaks are sorted (improved performance)
	//Calculates baseline ignoring points lying inside peaks
	//Warning: quite expensive function
	double find_baseline_by_median(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks);

	//Warning: quite expensive function
	double find_baseline_by_median(DVECTOR &xs, DVECTOR &ys);

	//Assumes that peaks are sorted (improved performance)
	//Calculates baseline ignoring points lying inside peaks
	//Approximate baseline is required for enhanced precision
	double find_baseline_by_integral(double approx_baseline, DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks);

	//Approximate baseline is required for enhanced precision
	double find_baseline_by_integral(double approx_baseline, DVECTOR &xs, DVECTOR &ys);

	//Wrapper for find_background_v_0
	void find_baseline_by_ROOT(DVECTOR &xs, DVECTOR &ys, DVECTOR &ys_out,
		int numberIterations, int direction, int filterOrder, bool smoothing, int smoothWindow, bool compton, int sparse);

	//My modification of CERN ROOT's TSpectrum->Background() function
	//First, it works directly on std container (vector or deque)
	//Second, sparse parameter was added to speed up the algorithm (skip every (sparse-1) points
	//in the main loop). sparse = 1 is equivalent to TSpectrum->Background()
	std::string find_background_v_0(DVECTOR& spectrum,
		int numberIterations,
		int direction, int filterOrder,
		bool smoothing, int smoothWindow,
		bool compton, int sparse = 1);

	//Does not assume that xs are equidistant
	//Returns empty y_out in case of error
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &y_out, double baseline);

	//Does not assume that xs are equidistant
	//Calculates integral as well as its standard deviation
	//Returns empty y_out and y_out_err in case of error
	//Output standard deviation may be wrong due to limited precision
	void integrate_with_variance(DVECTOR &xs, DVECTOR &ys, DVECTOR &ys_err,  DVECTOR &y_out, DVECTOR &y_out_err, double baseline);

	//Assumes that xs are equidistant
	//Returns empty y_out in case of error
	//Output standard deviation may be wrong due to limited precision
	void integrate(DVECTOR &ys, DVECTOR &y_out, double dx_hint, double baseline);

	//Does not assume that xs are equidistant
	//Carefully handles integration dx near left and right points and other special cases.
	//E.g. if (left +eps1)==(right-eps2)==xs[j] integral will be ys[j]*(eps2-eps1) assuming that each eps is far from neighbors xs
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out, double left, double right, double baseline);

	//Assumes that xs are equidistant
	//Carefully handles integration dx near left and right points and other special cases.
	//E.g. if (left +eps1)==(right-eps2)==xs[j] integral will be ys[j]*(eps2-eps1) assuming that each eps is far from neighbors xs
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out, double left, double right, double dx_hint, double baseline);

	//Integrates in the range [left, right), i.e.
	//Use xs.begin(), xs.end() to integrate whole signal
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out, DITERATOR left, DITERATOR right, double baseline);

	//Assumes that xs are equidistant
	//Integrates in the range [left, right), i.e.
	//Use xs.begin(), xs.end() to integrate whole signal
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out, DITERATOR left, DITERATOR right, double dx_hint, double baseline);

	//Does not assume that xs are equidistant
	//Integrates in the range [left, right), i.e.
	//Use xs.begin(), xs.end() to integrate whole signal
	void integrate(DVECTOR &xs, DVECTOR &ys, double &y_out, DITERATOR left, DITERATOR right, double baseline);

	//Assumes that xs are equidistant
	//Integrates in the range [left, right), i.e.
	//Use xs.begin(), xs.end() to integrate whole signal
	void integrate(DVECTOR &xs, DVECTOR &ys, double &y_out, DITERATOR left, DITERATOR right, double dx_hint, double baseline);

	//Crops signal to [x_left, x_right]
	void apply_time_limits(DVECTOR &xs, DVECTOR &ys, double x_left, double x_right);

	//Assumes that xs are equidistant
	//Crops signal to [x_left, x_right]
	void apply_time_limits(DVECTOR &xs, DVECTOR &ys, double x_left, double x_right, double dx_hint);

	//Assumes that peaks are sorted (improved performance)
	//Removes signal corresponding to peaks
	void exclude_peaks(DVECTOR &xs_in, DVECTOR &ys_in, STD_CONT<peak> &peaks);

	//Splits signal into found peaks. Even peaks (counting from 0) are stored in output.first and the odd ones - in output.second.
	//Both output.first and output.second store x-y pair arrays corresponding to separate peaks.
	//That is output.first.first[0] and output.first.second[0] contain xs and ys for the first (even) peak.
	//        output.second.first[0] and output.second.second[0] contain xs and ys for the second (odd) peak.
	//        output.first.first[1] and output.first.second[1] contain xs and ys for the third (even) peak.
	//... and so forth.
	void select_peaks(const DVECTOR &xs, const DVECTOR &ys, STD_CONT<peak> &peaks, double dx_hint, double baseline,
			std::pair<std::pair<STD_CONT<DVECTOR>, STD_CONT<DVECTOR>>, std::pair<STD_CONT<DVECTOR>, STD_CONT<DVECTOR>>>  &output);

	//Searches y maximum in [x_start, x_finish) range
	//Returns value of maximum and its position (for xs)
	std::pair<double, DITERATOR> get_max(DVECTOR &xs, DVECTOR &ys, DITERATOR x_start, DITERATOR x_finish);
	//Returns value of minimum and its position (for xs)
	std::pair<double, DITERATOR> get_max(DVECTOR &xs, DVECTOR &ys);

	//Searches y minimum in [x_start, x_finish) range
	//Returns value of minimum and its position (for xs)
	std::pair<double, DITERATOR> get_min(DVECTOR &xs, DVECTOR &ys, DITERATOR x_start, DITERATOR x_finish);
	//Returns value of minimum and its position (for xs)
	std::pair<double, DITERATOR> get_min(DVECTOR &xs, DVECTOR &ys);

	//Old version of searching peaks by threshold (before 2019.12.06) rewritten in find_peaks_fine_v2 style (more readable) 
	//Everything that is above threshold is counted as peak.
	//Peak fronts (ends, edges) are determined by intersection with threshold_edges (typically == baseline)
	//Assumes xs are equidistant
	//Returns sorted peaks
	void find_peaks_fine(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks, double base_line, double threshold, double threshold_edges);

	//New version of searching peaks by threshold (after 2019.12.06)
	//Everything that is above threshold is counted as peak.
	//Peak fronts (ends, edges) are determined by intersection with threshold_edges (typically == baseline)
	//If dip between two peaks is larger than threshold_edges but less than threshold
	//this new algorithm separates them. Position of the minimum between them is taken as separating point.
	//Assumes xs are equidistant
	//Returns sorted peaks
	void find_peaks_fine_v2(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks, double baseline, double threshold, double threshold_edges);

	void subtract_baseline(DVECTOR &ys_in, double base_line);
	//Both vectors must have the same sizes
	void subtract_baseline(DVECTOR &ys_in, DVECTOR &base_ys);

	//Required for subtracting ROOT's baseline which is calculated only for some range
	//Assumes that xs_in and base_xs have the same points except for maybe some points
	//See code for details
	void subtract_baseline(DVECTOR& xs_in, DVECTOR &ys_in, DVECTOR &base_xs, DVECTOR &base_ys, double baseline_baseline);

};

#endif
