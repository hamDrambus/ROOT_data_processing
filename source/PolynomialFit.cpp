#include "PolynomialFit.h"

double PolynomialFit::pown(double val, unsigned int n) const
{
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

#ifndef _AVOID_CERN_ROOT

PolynomialFit::PolynomialFit(int order)
{
	setOrder(order);
}

void PolynomialFit::setOrder(int n)
{
	if (n < 0)
		_order = 2;
	else
		_order = n;
}

int PolynomialFit::getOrder(void) const
{
	return _order;
}

void PolynomialFit::getCoefs(TVectorD &pars) const
{
	pars.ResizeTo(_last_coefs);
	pars = _last_coefs;
}

void PolynomialFit::operator ()(std::vector<double> &xs_in, std::vector<double> &ys_in,
	TVectorD &pars_out, double in_x0){
	return (*this)(xs_in, ys_in, 0, xs_in.size(), pars_out, in_x0);
}

void PolynomialFit::operator ()(std::vector<double> &xs_in, std::vector<double> &ys_in,
	int offset, int N_points, TVectorD &pars_out, double in_x0) //only for a part of a vector
{
	if (xs_in.size() != ys_in.size())
		return;
	if ((xs_in.size()-offset) < N_points)
		return;
	if (N_points < (_order + 1))
		return;

	TMatrixD mat(N_points, _order + 1);
    for (int col = 0, col_end_ = mat.GetNcols(); col != col_end_; ++col)
        for (int row = 0, row_end_ = mat.GetNrows(); row != row_end_; ++row)
            mat[row][col] = pown(xs_in[offset+ row] - in_x0, col);
	TVectorD Y(N_points);
	for (int row = 0; row < Y.GetNrows(); row++)
		Y[row] = ys_in[offset + row];
	TMatrixD mT(mat);
	mT.T();
	TMatrixD In(mT*mat);//because normal assignment like mat = mT*mat does not work! First resizing must be done.
	In.SetTol(1e-40);
	_last_coefs.ResizeTo(_order + 1);
	_last_coefs = (In.Invert())*mT*Y;
	pars_out.ResizeTo(_last_coefs);
	pars_out = _last_coefs;
}

#else //_AVOID_CERN_ROOT

namespace uBLAS = boost::numeric::ublas;

PolynomialFit::PolynomialFit(std::size_t order)
{
	setOrder(order);
}

PolynomialFit::~PolynomialFit()
{}

std::vector<double> PolynomialFit::operator ()(const std::vector<std::pair<double, double>> &vals_in, boost::optional<double> &in_x0) const {
	return (*this)(vals_in, 0, vals_in.size(), in_x0);
}

std::vector<double> PolynomialFit::operator ()(const std::vector<std::pair<double, double>> &vals,
	int offset, int N_points, boost::optional<double> &in_x0) const//only for a part of a vector
{
	std::vector<double> out;
	if ((vals.size() - offset) < N_points) {
		std::cout << "PolynomialFit::operator(): Error: N points is out of range:" << std::endl;
		std::cout << "\tx.size()=" << vals.size() << " offset=" << offset << " N_points=" << N_points << std::endl;
		return out;
	}
	if (offset < 0) {
		std::cout << "PolynomialFit::operator(): Error: offset is out of range:" << std::endl;
		std::cout << "\tx.size()=" << vals.size() << " offset=" << offset << " N_points=" << N_points << std::endl;
		return out;
	}
	if (N_points < (_order + 1)) {
		std::cout << "PolynomialFit::operator(): Error: no enough N points for fit:" << std::endl;
		std::cout << "\torder=" << _order << " N_points=" << N_points << std::endl;
		return out;
	}
	in_x0 = (in_x0 ? *in_x0 : vals[offset].first); //It is bad to set x0 to some fixed value (e.g. 0) because
	//interpolating too far from it will result in unstable results due to limited precision.
	//Ideally x0 should be set to the point at which we interpolate the data.
	if (1 == _order) {
		out.resize(2);
		out[1] = (vals[offset + 1].second - vals[offset].second) / (vals[offset + 1].first - vals[offset].first);
		out[0] = vals[offset].second + (*in_x0 - vals[offset].first)*out[1];
		//^value at in_x0 point
	} else {
		uBLAS::matrix<double> mat(N_points, _order + 1);
		for (int col = 0, col_end_ = mat.size2(); col < col_end_; ++col)
			for (int row = 0, row_end_ = mat.size1(); row < row_end_; ++row)
				mat(row, col) = pow(vals[offset + row].first - *in_x0, col);
		uBLAS::vector<double> Y(N_points);
		for (int row = 0, row_end_ = Y.size(); row < row_end_; ++row)
			Y[row] = vals[offset + row].second;
		//Solve the equation mat^T*mat*X = mat^T*Y for X via LU decomposition (mat is generally not diagonal)
		Y = uBLAS::prod(uBLAS::trans(mat), Y);
		mat = uBLAS::prod(uBLAS::trans(mat), mat);
		int res = uBLAS::lu_factorize(mat);
		if (res != 0)
			return out;
		uBLAS::inplace_solve(mat, Y, uBLAS::unit_lower_tag());
		uBLAS::inplace_solve(mat, Y, uBLAS::upper_tag());
		out.resize(Y.size());
		std::copy(Y.begin(), Y.end(), out.begin());
	}
	if (out.size() != (_order + 1)) {
		out.resize(0);
		return out;
	}
	return out;
}

//in_x0 - relative to what point carry out fit. Automatic value is set if boost::none is passed.
std::vector<double> PolynomialFit::operator ()(const std::vector<double> &xs, const std::vector<double> &ys, boost::optional<double> &in_x0) const
{
	return (*this)(xs, ys, 0, std::min(xs.size(), ys.size()), in_x0);
}
//Fit only part of a vector. offset+N_points-1 must be in the range of the vector
std::vector<double> PolynomialFit::operator ()(const std::vector<double> &xs, const std::vector<double> &ys, int offset, int N_points, boost::optional<double> &in_x0) const
{
	std::vector<double> out;
	if (xs.size() != ys.size()) {
		std::cout << "PolynomialFit::operator(): Error: x-y size mismatch:" << std::endl;
		std::cout << "\txs.size()=" << xs.size() << "\tys.size()=" << ys.size() << std::endl;
		return out;
	}
	if ((xs.size() - offset) < N_points) {
		std::cout << "PolynomialFit::operator(): Error: N points is out of range:" << std::endl;
		std::cout << "\tx.size()=" << xs.size() << " offset=" << offset << " N_points=" << N_points << std::endl;
		return out;
	}
	if (offset < 0) {
		std::cout << "PolynomialFit::operator(): Error: offset is out of range:" << std::endl;
		std::cout << "\tx.size()=" << xs.size() << " offset=" << offset << " N_points=" << N_points << std::endl;
		return out;
	}
	if (N_points < (_order + 1)) {
		std::cout << "PolynomialFit::operator(): Error: no enough N points for fit:" << std::endl;
		std::cout << "\torder=" << _order << " N_points=" << N_points << std::endl;
		return out;
	}
	in_x0 = (in_x0 ? *in_x0 : xs[offset]); //It is bad to set x0 to some fixed value (e.g. 0) because
	//interpolating too far from it will result in unstable results due to limited precision.
	//Ideally x0 should be set to the point at which we interpolate the data.
	if (1 == _order) {
		out.resize(2);
		out[1] = (ys[offset + 1] - ys[offset]) / (xs[offset + 1] - xs[offset]);
		out[0] = ys[offset] + (*in_x0 - xs[offset])*out[1];
		//^value at in_x0 point
	} else {
		uBLAS::matrix<double> mat(N_points, _order + 1);
		for (std::size_t col = 0, col_end_ = mat.size2(); col < col_end_; ++col)
			for (std::size_t row = 0, row_end_ = mat.size1(); row < row_end_; ++row)
				mat(row, col) = pown(xs[offset + row] - *in_x0, col);
		uBLAS::vector<double> Y(N_points);
		for (std::size_t row = 0, row_end_ = Y.size(); row < row_end_; ++row)
			Y[row] = ys[offset + row];
		//Solve the equation mat^T*mat*X = mat^T*Y for X via LU decomposition (mat is generally not diagonal)
		Y = uBLAS::prod(uBLAS::trans(mat), Y);
		mat = uBLAS::prod(uBLAS::trans(mat), mat);
		int res = uBLAS::lu_factorize(mat);
		if (res != 0)
			return out;
		uBLAS::inplace_solve(mat, Y, uBLAS::unit_lower_tag());
		uBLAS::inplace_solve(mat, Y, uBLAS::upper_tag());
		out.resize(Y.size());
		std::copy(Y.begin(), Y.end(), out.begin());
	}
	if (out.size() != (_order + 1)) {
		out.resize(0);
		return out;
	}
	return out;
}

#endif //_AVOID_CERN_ROOT