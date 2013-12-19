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

//!  Lorenz attractor
/*!
	 The Lorenz system is a system of ordinary differential equations (the Lorenz equations) 
	 first studied by Edward Lorenz. It is notable for having chaotic solutions for certain 
	 parameter values and initial conditions. 
	 In particular, the Lorenz attractor is a set of chaotic solutions of the Lorenz system 
	 which, when plotted, resemble a butterfly or figure eight.
*/

template<typename Vector>
class ODESystemObject_Test_Lorenz
{

protected: 

	unsigned int number_of_equations_;

public:

	ODESystemObject_Test_Lorenz()
	{
		number_of_equations_ = 3000;

		sigma = 10.0;
		R = 28.0;
		b = 8.0 / 3.0;
	}

	void GetFunctions(const double x, const Vector& y, Vector& dy_over_dx)
	{
		for(unsigned int i=0;i<number_of_equations_/3;i++)
		{
			const unsigned int index = 3*i;
			dy_over_dx[0+index] = sigma * ( y[1+index] - y[0+index] );
			dy_over_dx[1+index] = R * y[0+index] - y[1+index] - y[0+index] * y[2+index];
			dy_over_dx[2+index] = -b * y[2+index] + y[0+index] * y[1+index];
		}
	}

	void PrintStep(const double x, const Vector& y, const Vector& dy_over_dx)
	{
	//	std::cout << "OdeSMOKE: " << x << '\t' << y[0] << '\t' << y[1] << '\t' << y[2] << std::endl;
	}

	double final_abscissa() const 
	{
		return 25.;
	}

	const Vector initial_condition() const
	{
		Vector y0;
		OdeSMOKE::resize<Vector>(number_of_equations_,y0);
		for(unsigned int i=0;i<number_of_equations_/3;i++)
		{
			const unsigned int index = 3*i;
			y0[0+index] = 10.;
			y0[1+index] = 1.;
			y0[2+index] = 1.;
		}
		return y0;
	}

private:

	double sigma;
	double R;
	double b;
};