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

#ifndef ODESMOKE_RUNGEKUTTAFAMILY_H
#define	ODESMOKE_RUNGEKUTTAFAMILY_H

#include "ode_solver_virtual_class.h"

namespace OdeSMOKE
{
	//!  Virtual class for Runge-Kutta solvers. 
	/*!
		 This class provides a common interface for implementing different Runge-Kutta solvers.
	*/

	template<typename Vector, typename Method>
	class OdeRungeKuttaFamily : public Method
	{
	public:

		//! Default constructor
		/*!
		  This is the default constructor
		*/
		OdeRungeKuttaFamily();

		//! Set the safety coefficient for updating the new step size
		/*!
		  This function sets the safety coefficient for updating the new step size (default 0.70)
		  \param safety_coefficient the safety coefficient (must be larger than 0 and smaller than 1)
		*/
		void SetSafetyCoefficient(const double safety_coefficient);

		//! Set the maximum enlargement ratio for updating the new step size
		/*!
		  This function sets the maximum enlargement ratio for updating the new step size (default 4)
		  \param smax_enlargement_ratio the maximum enlargement ratio (must be strictly larger than 1)
		*/
		void SetMaximumEnlargementRatio(const double max_enlargement_ratio);

		//! Set the stabilization factor for updating the new step size
		/*!
		  This function sets the stabilization factor for updating the new step size (default 0)
		  Small values are for more aggressive policy. Larger values (<0.1) for more robust policy.
		  \param stabilization_factor the stabilization factor (must be between 0 and 0.1)
		*/
		void SetStabilizationFactor(const double stabilization_factor);

		//! Solution of the ODE system
		/*!
		  This function solves the ODE system according to the selected Runge-Kutta algorithm
		  \param x0 the starting abscissa
		  \param xF the final abscissa
		*/
		void Solve(const double x0, const double xF);

		//! Solution of the ODE system
		/*!
		  This function solves the ODE system according to the selected Runge-Kutta algorithm
		  \param x0 the starting abscissa
		  \param dx the fixed step of the requested intervals
		  \param nsteps the number of requested intervals
		*/
		void Solve(const double x0, const double dx, const unsigned int nsteps);

		//! Number of steps
		/*!
		    Returns the number of steps which were performed during the integration
		*/
		unsigned int number_steps() const;

		//! Number of accepted steps
		/*!
		    Returns the number of accepted steps which were performed during the integration
			It must be equal to the number of steps
		*/
		unsigned int number_accepted_steps() const;

		//! Number of rejected steps
		/*!
		    Returns the number of rejected steps during the integration
		*/
		unsigned int number_rejected_steps() const;

		//! Number of function evaluations
		/*!
		    Returns the number total number of function evaluations during the integration
		*/
		unsigned int number_function_evaluations() const;

	};

}

#include "runge_kutta_family.hpp"

#endif	// OpenSMOKE_RUNGEKUTTAFAMILY_H