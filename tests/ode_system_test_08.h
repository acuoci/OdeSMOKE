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
		From:         Atkinson K.E., An introduction to numerical analysis,
		              John Wiley and Sons, New York (1989)
		Suggested by: Buzzi-Ferraris G., Metodi numerici e software in C++, Addison-Wesley Longman Italia (1998)
		              Problem #8 (page 717)

*/

template<typename Vector>
class ODESystemObject_Test_08
{
protected:

	unsigned int number_of_equations_;

public:

	ODESystemObject_Test_08()
	{
		number_of_equations_ = 1;
	}

	virtual void GetFunctions(const double x, const Vector& y, Vector& dy_over_dx)
	{
		dy_over_dx[0] = y[0]/4.*(1.-y[0]/20.);
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
		y0[0] = 1.;
		return y0;
	}

	const Vector analytical_solution(const double x) const
	{
		Vector y;
		OdeSMOKE::resize<Vector>(number_of_equations_,y);
		y[0] = 20./(1.+19.*exp(-x/4.));
		return y;
	}
};