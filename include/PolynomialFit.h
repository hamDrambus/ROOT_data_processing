#ifndef POLYNOMIAL_FIT_H
#define POLYNOMIAL_FIT_H

#include "GlobalParameters.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cfloat>
#ifndef _AVOID_CERN_ROOT
#include "TMatrixD.h"
#include "TVectorD.h"

//parameters are [0]+[1]*x+[2]*x^2+...
class PolynomialFit {
protected:
	Int_t _order;
	TVectorD _last_coefs;
protected:
    double pown(double val, unsigned int n) const;
public:
	PolynomialFit(Int_t order);
	virtual void setOrder(Int_t n); //TODO: actualluy is is a bad practice to call virtual method from the constructor
	//but is is ok here, since derivative class only limits setOrder() possible values to {2}

	Int_t getOrder(void) const;
	void getCoefs(TVectorD &pars) const;

	virtual void operator ()(std::vector<double> &xs_in, std::vector<double> &ys_in,
		TVectorD &pars_out, double in_x0=0); //in_x0 - in what poInt_t set zero x (In the SG filter it is convinient to set x_in
	//to the poInt_t in which value is calculated
	virtual void operator ()(std::vector<double> &xs_in, std::vector<double> &ys_in,
		int offset, int N_points, TVectorD &pars_out, double in_x0=0); //only for a part of a vector
};

#else //_AVOID_CERN_ROOT

#include <boost/optional.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>

class PolynomialFit {
protected:
	std::size_t _order;
	double pown(double val, unsigned int n) const;
public:
	PolynomialFit(std::size_t order);
	~PolynomialFit();
	void setOrder(std::size_t n) {
		_order = n;
	}
	std::size_t getOrder(void) const {
		return _order;
	}
	//in_x0 - relative to what point carry out fit. Automatic value is set if boost::none is passed.
	std::vector<double> operator ()(const std::vector<std::pair<double, double>> &vals_in, boost::optional<double> &in_x0) const;
	//Fit only part of a vector. offset+N_points-1 must be in the range of the vector
	std::vector<double> operator ()(const std::vector<std::pair<double, double>> &vals_in, int offset, int N_points, boost::optional<double> &in_x0) const;
	//in_x0 - relative to what point carry out fit. Automatic value is set if boost::none is passed.
	std::vector<double> operator ()(const std::vector<double> &xs, const std::vector<double> &ys, boost::optional<double> &in_x0) const;
	//Fit only part of a vector. offset+N_points-1 must be in the range of the vector
	std::vector<double> operator ()(const std::vector<double> &xs, const std::vector<double> &ys, int offset, int N_points, boost::optional<double> &in_x0) const;
};

#endif //_AVOID_CERN_ROOT

#endif
