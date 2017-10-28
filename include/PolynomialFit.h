#ifndef POLYNOMIAL_FIT_H
#define POLYNOMIAL_FIT_H

#include "GlobalDefinitions.h"
#include "TMatrixD.h"
#include "TVectorD.h"

class PolynomialFit {
protected:
	int _order;
	TVectorD _last_coefs;
public:
	PolynomialFit(int order);
	virtual void setOrder(int n); //TODO: actualluy is is a bad practice to call virtual method from the constructor
	//but is is ok here, since derivative class only limits setOrder() possible values to {2}

	int getOrder(void) const;
	void getCoefs(TVectorD &pars) const;

	virtual void operator ()(const std::vector<double> &xs_in, const std::vector<double> &ys_in,
		TVectorD &pars_out, double in_x0=0); //in_x0 - in what point set zero x (In the SG filter it is convinient to set x_in
	//to the point in which value is calculated
	virtual void operator ()(const std::vector<double> &xs_in, const std::vector<double> &ys_in,
		int offset, int N_points, TVectorD &pars_out, double in_x0=0); //only for a part of a vector
};

#endif