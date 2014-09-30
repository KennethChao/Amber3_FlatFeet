// **Generated from Mathematica**

#include <cmath>

	/** \addtogroup mdefs.h
	 * @{
	 */
/**
 * \brief Modified Mathematica Definitions file
 * \author Wolfram Research Inc., Copyright 1986 through 1999
 * \author Eric Cousineau - Modified because macros suck.
 */
inline double Power(double x, double y) { return pow(x, y); }
inline double Sqrt(double x) { return sqrt(x); }

inline double Abs(double x) { return fabs(x); }

inline double Exp(double x) { return exp(x); }
inline double Log(double x) { return log(x); }

inline double Sin(double x) { return sin(x); }
inline double Cos(double x) { return cos(x); }
inline double Tan(double x) { return tan(x); }

inline double ArcSin(double x) { return asin(x); }
inline double ArcCos(double x) { return acos(x); }
inline double ArcTan(double x) { return atan(x); }

inline double Sinh(double x) { return sinh(x); }
inline double Cosh(double x) { return cosh(x); }
inline double Tanh(double x) { return tanh(x); }

const double E	= 2.71828182845904523536029;
const double Pi = 3.14159265358979323846264;
const double Degree = 0.01745329251994329576924;
    /** @} */

#include "desired.hpp"

using namespace std;
using namespace Eigen;

namespace proxi_opt {

double y_TimeExtCwf(const double &t, const Eigen::VectorXd &aRow)
{
	double value = (aRow(6) + aRow(4)*Cos(t*aRow(5)) + (aRow(0)*Cos(t*aRow(1)) + aRow(2)*Sin(t*aRow(1)))/Power(2.718281828459045,1.*t*aRow(3)) + (2.*aRow(3)*aRow(4)*aRow(5)*Sin(t*aRow(5)))/(Power(aRow(1),2) + Power(aRow(3),2) - 1.*Power(aRow(5),2)));
	return value;
}

double dy_TimeExtCwf(const double &t, const Eigen::VectorXd &aRow)
{
	double value = ((2.*aRow(3)*aRow(4)*Power(aRow(5),2)*Cos(t*aRow(5)))/(Power(aRow(1),2) + Power(aRow(3),2) - 1.*Power(aRow(5),2)) - (1.*aRow(1)*(-1.*aRow(2)*Cos(t*aRow(1)) + aRow(0)*Sin(t*aRow(1))))/Power(2.718281828459045,1.*t*aRow(3)) - (1.*aRow(3)*(aRow(0)*Cos(t*aRow(1)) + aRow(2)*Sin(t*aRow(1))))/Power(2.718281828459045,1.*t*aRow(3)) - 1.*aRow(4)*aRow(5)*Sin(t*aRow(5)));
	return value;
}

double ddy_TimeExtCwf(const double &t, const Eigen::VectorXd &aRow)
{
	double value = ((-2.*aRow(1)*aRow(2)*aRow(3)*Cos(t*aRow(1)))/Power(2.718281828459045,1.*t*aRow(3)) - 1.*aRow(4)*Power(aRow(5),2)*Cos(t*aRow(5)) - (1.*Power(aRow(1),2)*aRow(2)*Sin(t*aRow(1)))/Power(2.718281828459045,1.*t*aRow(3)) + (aRow(2)*Power(aRow(3),2)*Sin(t*aRow(1)))/Power(2.718281828459045,1.*t*aRow(3)) + (aRow(0)*(-1.*Power(aRow(1),2)*Cos(t*aRow(1)) + Power(aRow(3),2)*Cos(t*aRow(1)) + 2.*aRow(1)*aRow(3)*Sin(t*aRow(1))))/Power(2.718281828459045,1.*t*aRow(3)) - (2.*aRow(3)*aRow(4)*Power(aRow(5),3)*Sin(t*aRow(5)))/(Power(aRow(1),2) + Power(aRow(3),2) - 1.*Power(aRow(5),2)));
	return value;
}

} // namespace proxi_opt