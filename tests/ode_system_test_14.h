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

//!  Non stiff ODE system
/*!
		This is a non stiff ODE system used for testing purposes.
		From:         Buchanan J.L., Turner P.R., Numerical Methods and Analysis,
		              McGraw-Hill Inc., New York (1992)
		Suggested by: Buzzi-Ferraris G., Metodi numerici e software in C++, Addison-Wesley Longman Italia (1998)
		              Problem #14 (page 718)

*/

template<typename Vector>
class ODESystemObject_Test_14
{
protected:

	unsigned int number_of_equations_;

public:

	ODESystemObject_Test_14()
	{
		number_of_equations_ = 4;
	}

	virtual void GetFunctions(const double x, const Vector& y, Vector& dy_over_dx)
	{
		dy_over_dx[0] = y[1];
		dy_over_dx[1] = -2.*y[0]+.5*y[2];
		dy_over_dx[2] = y[3];
		dy_over_dx[3] = 2.*y[0]-2.*y[2]+10.*cos(2.*x);
	}

	virtual void PrintStep(const double x, const Vector& y, const Vector& dy_over_dx)
	{
	}

	double final_abscissa() const 
	{
		return 20.;
	}

	const Vector initial_condition() const
	{
		Vector y0;
		OdeSMOKE::resize<Vector>(number_of_equations_,y0);
		y0[0] = 0.;
		y0[1] = 0.;
		y0[2] = 0.;
		y0[3] = 0.;
		return y0;
	}

	const Vector analytical_solution(const double x) const
	{
		Vector y;
		OdeSMOKE::resize<Vector>(number_of_equations_,y);
		y[0] = (5.*cos(x))/6. - (5.*cos(sqrt(3.)*x))/2. + (10.*cos(x)*cos(x))/3. - 5./3.;
		y[1] = (5.*sqrt(3.)*sin(sqrt(3.)*x))/2. - (5.*sin(x))/6. - (10.*sin(2.*x))/3.;
		y[2] = 5.*cos(sqrt(3.)*x) + (5.*cos(x))/3. - (40.*cos(x)*cos(x))/3. + 20./3.;
		y[3] = (40.*sin(2.*x))/3. - (5.*sin(x))/3. - 5*sqrt(3.)*sin(sqrt(3.)*x);
		return y;
	}
};