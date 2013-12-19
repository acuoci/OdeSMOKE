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

#ifndef OpenSMOKE_RUNGEKUTTA4TH_H
#define	OpenSMOKE_RUNGEKUTTA4TH_H

#include "runge_kutta_family.h"

namespace OdeSMOKE
{
	//!  Runge-Kutta 4th Order solver. 
	/*!
		 This is the classic 4th order Runge-Kutta solver. The new step is estimated using the Richardson's extrapolation.


		 \f$
			\mathbf{y_{n+1}}=\mathbf{y_n}+\frac{\mathbf{k_1}+2\mathbf{k_2}+2\mathbf{k_3}+\mathbf{k_4}}{6} \\
			
			\left\{\begin{matrix}
				\mathbf{k_1}=h\mathbf{f}(\mathbf{y_n},t_n) \\ 
				\mathbf{k_2}=h\mathbf{f}(\mathbf{y_n}+0.5\mathbf{k_1},t_n+0.5t_n) \\ 
				\mathbf{k_3}=h\mathbf{f}(\mathbf{y_n}+0.5\mathbf{k_2},t_n+0.5t_n) \\ 
				\mathbf{k_4}=h\mathbf{f}(\mathbf{y_n}+\mathbf{k_3},t_n+h)
			\end{matrix}\right.
		\f$                                   
	*/

	template<typename ODESystemObject, typename Vector>
	class OdeRungeKutta4th : public ODESystemObject
	{
	public:

		//! Default constructor
		/*!
		  This is the default constructor
		*/
		OdeRungeKutta4th();

		//! Print useful info on the screen
		/*!
		    Print useful info on the screen
		*/
		void Info() const;

	private:

		void SetDefaultConditions();

		//! Application of the Runge-Kutta algorithm to the requested interval
		/*!
		  This function applies the Runge-Kutta algorithm to the requested integration interval
		  \param tF the final abscissa
		*/
		void AdvanceAdaptiveStepSize(const double tF);

		//! Application of the Runge-Kutta algorithm to the requested number of intervals
		/*!
		  This function applies the Runge-Kutta algorithm to the requested integration number of intervals
		  \param nsteps the number of fixed steps
		*/
		void AdvanceFixedStepSize(const unsigned int nsteps);

		//! Selection of the proper integration step
		/*!
		  This function selects the proper integration step based on the evaluation of the local error
		*/
		void AdvanceOverSingleStep();

		//! Application of the Runge-Kutta algorithm to a single step
		/*!
		  This function implements a single step of the Runge-Kutta algorithm
		  \param x_in the current abscissa
		  \param y_in the current solution (independent variable)
		  \param dy_over_dx_in the current derivatives
		  \param y_out the solution at the end of the current step
		*/
		void Step(const double x_in, const Vector& y_in, const Vector& dy_over_dx_in, Vector& y_out);

	private:

		static const double a21;
		static const double a32;
		
		static const double c2;
		static const double c3;
		static const double c4;

		static const double b1;
		static const double b2;
		static const double b3;
		static const double b4;

		Vector vb;
		Vector k2_;
		Vector k3_;
		Vector k4_;

		Vector dy_over_dx_;						/*!< current derivatives */
		Vector dy_over_dx_old_;
		Vector y_2h_;
		Vector y_h_;

		double h_;										/*!< current step */
		unsigned int number_function_evaluations_;		/*!< current number of function evaluations */
		unsigned int number_steps_;						/*!< maximum number of steps */
		unsigned int number_accepted_steps_;			/*!< number of accepted step estimations */
		unsigned int number_rejected_steps_;			/*!< number of rejected step estimations */

		double safety_coefficient_;						/*!< safety coefficient for the estimation of the new step (default 0.70) */
		double max_enlargement_ratio_;					/*!< maximum enlargement factor from one step to the next one (default 4) */
		double stabilization_factor_;					/*!< stabilization factor for stepsize control (default 0) */

		unsigned int max_number_steps_;					/*!< maximum number of steps which can be performed */
		double max_step_;								/*!< maximum step allowed during the integration */
		double min_step_;								/*!< minimum step allowed during the integration */
		double first_step_;								/*!< first step  */
		Vector abs_tolerances_;							/*!< absolute tolerances (default 1e-12) */
		Vector rel_tolerances_;							/*!< relative tolerances (default 1e-7) */
		Vector y0_;										/*!< initial conditions */
		double x0_;										/*!< initial abscissa */
		double xF_;										/*!< final abscissa */
		double x_;										/*!< current abscissa */
		Vector y_;										/*!< current solution */

		bool analytical_jacobian_;						/*!< analytical calculation of Jcobian matrix */

		bool user_defined_max_step_;					/*!< boolean variable indicating if the maximum step is defined by the user */
		bool user_defined_min_step_;					/*!< boolean variable indicating if the minimum step is defined by the user */
		bool user_defined_first_step_;					/*!< boolean variable indicating if the first step is defined by the user */

	};
}

#include "runge_kutta_4th.hpp"

#endif	// OpenSMOKE_RUNGEKUTTA4TH_H