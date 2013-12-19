/*----------------------------------------------------------------------*\
|     ____                    ______ __  __  ____  _  ________            |
|    / __ \                  /  ___ |  \/  |/ __ \| |/ /  ____|           |
|   | |  | |_ __   ___ _ __ |  (___ | \  / | |  | | ' /| |__              |
|   | |  | | '_ \ / _ \ '_ \ \___  \| |\/| | |  | |  < |  __|             |
|   | |__| | |_) |  __/ | | |____)  | |  | | |__| | . \| |____            |
|    \____/| .__/ \___|_| |_|______/|_|  |_|\____/|_|\_\______|           |
|          | |                                                            |
|          |_|                                                            |
|                                                                         |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|	Date: 07 Mar 2013                                                     |
|-------------------------------------------------------------------------|
|	License                                                               |
|                                                                         |
|   This file is part of OpenSMOKE.                                       |
|                                                                         |
|   OpenSMOKE is free software: you can redistribute it and/or modify     |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE is distributed in the hope that it will be useful,          |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE. If not, see <http://www.gnu.org/licenses/>.     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "ode_system_object_virtual_class.h"

template<typename Vector, typename Matrix>
class ODESystemObject_Test_Diffusion : public OdeSMOKE::ODESystemObjectVirtualClass<Vector, Matrix>
{
public:

	ODESystemObject_Test_Diffusion()
	{
		a_ = 5.;
		b_ = 5.;
		m_ = 0.6;
		e_ = 2.;
		f_ = -1.4;

		n_ = 100;
		L_ = 1.;
		dx_ = L_/double(n_-1);

		number_of_equations_ = 2*n_;

		p_.resize(n_);
		h_.resize(n_);
	}

	virtual void GetFunctions(const double t, const Vector& y, Vector& dy_over_dx)
	{
		const double alfa = 1.e-6;

		unsigned int count = 0;
		for(unsigned int i=0;i<n_;i++)
			p_(i) = y(count++);
		for(unsigned int i=0;i<n_;i++)
			h_(i) = y(count++);

		// Equation for p
		dy_over_dx(0) = (h_[2]-h_[1])*alfa;
		for(unsigned int i=1;i<n_-1;i++)
		{
			const double rx = e_+f_*(dx_*i);
			const double d2p_over_dx2 = ( p_(i+1) - 2.*p_(i) + p_(i-1)) / (dx_*dx_);
			dy_over_dx(i) = rx*p_(i)*(1.-p_(i)) - a_*p_(i)/(1.+b_*p_(i))*h_(i)+d_*d2p_over_dx2;
		}
		dy_over_dx(n_-1) = (h_[n_-1]-h_[n_-2])*alfa;

		// Equation for h
		dy_over_dx(n_) = (p_[2]-p_[1])*alfa;
		for(unsigned int i=1;i<n_-1;i++)
		{
			const double dhp_over_dx2 = ( h_(i+1) - 2.*h_(i) + h_(i-1)) / (dx_*dx_);
			dy_over_dx(n_+i) = a_*p_(i)/(1.+b_*p_(i))*h_(i)-m_*h_(i)+d_*dhp_over_dx2;
		}
		dy_over_dx(2*n_-1) = (p_[n_-1]-p_[n_-2])*alfa;
	}

	virtual void GetNumericalJacobian(const double x, const Vector& y, Matrix& J)
	{
	}

	virtual void GetLinearSystemSolution(const Matrix& A, const Vector& b, Vector& z)
	{
	}

	virtual void PrintStep(const double x, const Vector& y, const Vector& dy_over_dx)
	{
		unsigned int i = n_/10*8;
		std::cout << x << " " << e_+f_*(dx_*i) << " " << p_(i) << " " << h_(i) << " " << p_(i)+h_(i) << std::endl;
	}

	double final_abscissa() const 
	{
		return 100.;
	}

	const Vector initial_condition() const
	{
		Vector y0(number_of_equations_);

		unsigned int count = 0;
		for(unsigned int i=0;i<n_;i++)
			y0(count++) = 0.4;
		for(unsigned int i=0;i<n_;i++)
			y0(count++) = 0.6;

		return y0;
	}

	const Vector analytical_solution(const double x) const
	{
		Vector y(number_of_equations_);
		return y;
	}

private:

	unsigned int number_of_equations_;

	double a_;
	double b_;
	double m_;
	double d_;
	double e_;
	double f_;

	unsigned int n_;
	double L_;
	double dx_;

	Vector p_;
	Vector h_;
};