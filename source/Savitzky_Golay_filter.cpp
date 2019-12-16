#include "Savitzky_Golay_filter.h"

namespace uBLAS = boost::numeric::ublas;

double SavitzkyGolayFilter::pown(double val, unsigned int n) const
{
	if (0==val)
		return val;
	double result = 1;
	while (true) {
		if (n & 1)
			result *= val;
		n >>= 1;
		if (!n)
			break;
		val *= val;
	}
	return result;
}

SavitzkyGolayFilter::SavitzkyGolayFilter(std::size_t n_points, std::size_t order, std::size_t n_iterations) :
		_n_points(n_points), _order(order), _n_iterations(n_iterations)
{}

std::vector<double> SavitzkyGolayFilter::calculate_coefs(void) const
{
	std::vector<double> _coefs;
	if (!isValid()) {
		_coefs.clear();
		return _coefs;
	}
	if (1 == _order && 2 == _n_points ) {
		_coefs.resize(_n_points);
		_coefs[0] = 1;
		_coefs[1] = 0;
		return _coefs;
	}
	std::size_t center_point = (_n_points-1)/2;
	uBLAS::matrix<double> mat(_n_points, _order + 1);
	for (std::size_t col = 0, col_end_ = mat.size2(); col < col_end_; ++col)
		for (std::size_t row = 0, row_end_ = mat.size1(); row < row_end_; ++row)
			mat(row, col) = pown((double)((int)row - (int)center_point), col);
	uBLAS::matrix<double> matT_mat = uBLAS::prod(uBLAS::trans(mat), mat);
	int res = uBLAS::lu_factorize(matT_mat);
	if (res != 0) {
		std::cerr<<"SavitzkyGolayFilter::calculate_coefs: Error: failed lu_factorize"<<std::endl;
		_coefs.clear();
		return _coefs;
	}
	_coefs.resize(_n_points);
	for (std::size_t ny = 0, ny_end_ = _n_points; ny < ny_end_; ++ny) {
		uBLAS::vector<double> Y(_n_points);
		for (std::size_t row = 0, row_end_ = Y.size(); row < row_end_; ++row)
			Y[row] = (row == ny ? 1.0 : 0);
		//Solve the equation mat^T*mat*X = mat^T*Y for X via LU decomposition (mat is generally not diagonal)
		Y = uBLAS::prod(uBLAS::trans(mat), Y);
		uBLAS::inplace_solve(matT_mat, Y, uBLAS::unit_lower_tag());
		uBLAS::inplace_solve(matT_mat, Y, uBLAS::upper_tag());
		_coefs[ny] = Y[0];
	}
	return _coefs;
}

void SavitzkyGolayFilter::setNPoints(std::size_t n)
{
	_n_points = n;
}
void SavitzkyGolayFilter::setOrder(std::size_t n)
{
	_order = n;
}
void SavitzkyGolayFilter::setNIter(std::size_t n)
{
	_n_iterations = n;
}
void SavitzkyGolayFilter::setPars(std::size_t n_points, std::size_t order, std::size_t n_iterations)
{
	_order = order;
	_n_points = n_points;
	_n_iterations = n_iterations;
}

bool SavitzkyGolayFilter::isValid(void) const
{
	return _n_points>_order;
}

std::size_t SavitzkyGolayFilter::getNPoints(void) const
{	return _n_points;}
std::size_t SavitzkyGolayFilter::getOrder(void) const
{	return _order;}
std::size_t SavitzkyGolayFilter::getNIter(void) const
{	return _n_iterations;}
void SavitzkyGolayFilter::getPars(std::size_t &n_points, std::size_t &order, std::size_t &n_iterations) const
{
	n_points = _n_points;
	order = _order;
	n_iterations = _n_iterations;
}

void SavitzkyGolayFilter::operator ()(DVECTOR &xs_in_out, DVECTOR &ys_in_out) const
{
	if ((xs_in_out.size() != ys_in_out.size())||(xs_in_out.size() < _n_points)||(0==_n_iterations)||!isValid())
		return;
	DVECTOR ys_out (ys_in_out.size());
	PolynomialFit fit(_order);
#ifndef _AVOID_CERN_ROOT
	TVectorD A;
#else
	std::vector<double> A;
#endif
	for (int iter = 0; iter < _n_iterations; ++iter) {
		int start_index = 0;
		for (int h = 0, h_end_ = xs_in_out.size(); h <h_end_ ; ++h) {
			start_index = h - _n_points / 2;
			start_index = start_index < 0 ? 0 : start_index;
			if (start_index > xs_in_out.size() - _n_points)
				start_index = xs_in_out.size() - _n_points;
#ifndef _AVOID_CERN_ROOT
			fit(xs_in_out, ys_in_out, start_index, _n_points, A, xs_in_out[h]);
#else
			boost::optional<double> x0(xs_in_out[h]);
			A = fit(xs_in_out, ys_in_out, start_index, _n_points, x0);
#endif
			ys_out[h] = A[0]; //I moved X coordinates to the point of interest (xs_in_out[h]) in the matrix construction
			/*xs_out[h] = 0;
			for (int row = 0; row < A.GetNrows(); row++)
			xs_out[h] += A[row] * pow(xs_in[h] - xs_in[h], row);*/
		}
	}
	ys_in_out = ys_out;
}

//Assumes xs are equidistant
void SavitzkyGolayFilter::operator ()(DVECTOR &ys_in_out) const
{
	if ((ys_in_out.size() < _n_points)||(0==_n_iterations)||!isValid())
		return;
	std::vector<double> _coefs = calculate_coefs();
	if (_n_points!=_coefs.size()) {
		std::cerr<<"SavitzkyGolayFilter():Error: _coefs has wrong size"<<std::endl;
		return;
	}
	std::size_t size = ys_in_out.size();
	DVECTOR ys_out (size);
	PolynomialFit fit(_order); //for edge points;
	std::size_t half_range = (_n_points-1)/2;
#ifndef _AVOID_CERN_ROOT
	TVectorD A;
#else
	std::vector<double> A;
#endif
	std::vector<double> temp_xs(_n_points); //for calculating coefs for edge case using PolynomialFit
	for (std::size_t x=0, x_end_ = _n_points; x!=x_end_; ++x)
		temp_xs[x] = x;
	for (std::size_t iter = 0; iter < _n_iterations; ++iter) {
		std::vector<double> left_ys(_n_points);
		for (std::size_t x=0, x_end_ = _n_points; x!=x_end_; ++x)
			left_ys[x] = ys_in_out[x];
		std::vector<double> right_ys(_n_points);
		for (std::size_t x=0, x_end_ = _n_points; x!=x_end_; ++x)
			right_ys[x] = ys_in_out[size + x - _n_points];
		for (std::size_t h = 0; h <size ; ++h) {
			if (h < half_range) { //can't use _coefs for edge cases, so fall back to fit
#ifndef _AVOID_CERN_ROOT
				fit(temp_xs, left_ys, 0, _n_points, A, temp_xs[h]);
#else
				boost::optional<double> x0(temp_xs[h]);
				A = fit(temp_xs, left_ys, 0, _n_points, x0);
#endif
				ys_out[h] = A[0];
				continue;
			}
			if ((h + _n_points - half_range - 1) >= size) { //can't use _coefs for edge cases, so fall back to fit
#ifndef _AVOID_CERN_ROOT
				fit(temp_xs, right_ys, 0, _n_points, A, temp_xs[h + _n_points - size]);
#else
				boost::optional<double> x0(temp_xs[h + _n_points - size]);
				A = fit(temp_xs, right_ys, 0, _n_points, x0);
#endif
				ys_out[h] = A[0];
				continue;
			}
			ys_out[h] = 0;
			for (std::size_t c = 0; c!=_n_points; ++c)
				ys_out[h] += ys_in_out[h - half_range + c] * _coefs[c];
		}
		ys_in_out = ys_out;
	}
}
