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
|--------------------------------------------------------------------------
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
\*------------------------------------------------------------------/-----*/

namespace OdeSMOKE
{
	template<typename Vector, typename Solver>
	void ODESolverVirtualClass<Vector,Solver>::SetInitialConditions(const Vector& y0)
	{
		if (y0.size() != this->number_of_equations_)
			FatalError("The size of the initial condition vector is wrong");

		this->y0_ = y0;
	}

	template<typename Vector, typename Solver>
	void ODESolverVirtualClass<Vector,Solver>::SetAbsoluteTolerances(const Vector& abs_tolerances)
	{
		if (abs_tolerances.size() != this->number_of_equations_)
			FatalError("The size of the absolute tolerance vector is wrong");

		this->abs_tolerances_ = abs_tolerances;
	}

	template<typename Vector, typename Solver>
	void ODESolverVirtualClass<Vector,Solver>::SetRelativeTolerances(const Vector& rel_tolerances)
	{
		if (rel_tolerances.size() != this->number_of_equations_)
			FatalError("The size of the relative tolerance vector is wrong");

		this->rel_tolerances_ = rel_tolerances;
	}

	template<typename Vector, typename Solver>
	void ODESolverVirtualClass<Vector,Solver>::SetAbsoluteTolerances(const double abs_tolerance)
	{	
		if (abs_tolerance < 1.e-32)
			FatalError("The user-defined absolute tolerance is too small");

		this->abs_tolerances_.setConstant(abs_tolerance);
	}

	template<typename Vector, typename Solver>
	void ODESolverVirtualClass<Vector,Solver>::SetRelativeTolerances(const double rel_tolerance)
	{
		if (rel_tolerance < 1.e-32)
			FatalError("The user-defined relative tolerance is too small");

		this->rel_tolerances_.setConstant(rel_tolerance);
	}

	template<typename Vector, typename Solver>
	void ODESolverVirtualClass<Vector,Solver>::SetMaximumNumberOfSteps(const unsigned int max_number_steps)
	{
		if (max_number_steps <= 0)
			FatalError("The user-defined maximum number of steps must be at least equal to 1");

		this->max_number_steps_ = max_number_steps;
	}
		
	template<typename Vector, typename Solver>
	void ODESolverVirtualClass<Vector,Solver>::SetMaximumStep(const double max_step)
	{
		if (max_step <= 0)
			FatalError("The user-defined maximum step must be larger than 0");

		this->max_step_ = max_step;
		this->user_defined_max_step_ = true;
	}

	template<typename Vector, typename Solver>
	void ODESolverVirtualClass<Vector,Solver>::SetMinimumStep(const double min_step)
	{
		if (min_step <= 0)
			FatalError("The user-defined minimum step must be larger than 0");

		this->min_step_ = min_step;
		this->user_defined_min_step_ = true;
	}

	template<typename Vector, typename Solver>
	void ODESolverVirtualClass<Vector,Solver>::SetFirstStep(const double first_step)
	{
		if (first_step <= 0)
			FatalError("The user-defined first step must be larger than 0");

		this->first_step_ = first_step;
		this->user_defined_first_step_ = true;
	}		

	template<typename Vector, typename Solver>
	void ODESolverVirtualClass<Vector,Solver>::FatalError(const std::string message)
	{
		std::cout << "OdeSMOKE Solver Fatal Error Message" << std::endl;
		std::cout << " * " << message << std::endl;
		std::cout << " * Press enter to exit... " << std::endl;
		getchar();
		exit(-1);
	}
}