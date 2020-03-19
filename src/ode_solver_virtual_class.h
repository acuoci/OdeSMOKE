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

#ifndef OdeSMOKE_ODESOLVERVIRTUALCLASS_H
#define	OdeSMOKE_ODESOLVERVIRTUALCLASS_H

#include "utilities.h"

namespace OdeSMOKE
{
	//!  Virtual class for ODE solvers. 
	/*!
		 This class provides a common interface for implementing different ODE solvers.
	*/

	template<typename Vector, typename Solver>
	class ODESolverVirtualClass : public Solver
	{
	public:

		//! Set initial conditions
		/*!
		    Set the initial value of the dependent variable
		*/
		void SetInitialConditions(const Vector& y0);
		
		//! Set the absolute tolerances
		/*!
		    Set the absolute tolerances (default value: 1e-12)
		*/
		void SetAbsoluteTolerances(const Vector& abs_tolerances);

		//! Set the absolute tolerances
		/*!
		    Set the absolute tolerances (default value: 1e-12)
		*/
		void SetAbsoluteTolerances(const double abs_tolerance);
		
		//! Set the relative tolerances
		/*!
		    Set the relative tolerances (default value: 1e-7)
		*/
		void SetRelativeTolerances(const Vector& rel_tolerances);

		//! Set the relative tolerances
		/*!
		    Set the relative tolerances (default value: 1e-7)
		*/
		void SetRelativeTolerances(const double rel_tolerance);
		
		//! Set the maximum number of steps
		/*!
		    Set the maximum number of steps (default value: 500000)
		*/
		void SetMaximumNumberOfSteps(const unsigned int max_number_steps);
		
		//! Set the maximum allowed step
		/*!
		    Set the maximum allowed step (default value: 1e16)
		*/
		void SetMaximumStep(const double max_step);
		
		//! Set the minimum allowed step
		/*!
		    Set the minimum allowed step (default value: 1e-16)
		*/
		void SetMinimumStep(const double min_step);

		//! Set the first step
		/*!
		    Set the first step (default value: 0)
		*/
		void SetFirstStep(const double initial_step);
		
		//! Returns the current solution
		/*!
		    Returns the current solution
		*/
		const Vector& y() const { return this->y_; }

	protected:

		//! Set the default condtions
		/*!
		    Set the default conditions about the internal variables
		*/
		//void SetDefaultConditions();

		//! Fatal error message
		/*!
		    This function is called in case of a fatal error
			\param message the error message
		*/
		void FatalError(const std::string message);

	protected:


	};
}

#include "ode_solver_virtual_class.hpp"

#endif	//OdeSMOKE_DVODEPP_H