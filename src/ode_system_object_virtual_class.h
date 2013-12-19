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

#ifndef ODESMOKE_ODESYSTEMOBJECTCLASS_H
#define	ODESMOKE_ODESYSTEMOBJECTCLASS_H

#include <iostream>

namespace OdeSMOKE
{
	//!  Virtual class for ODE systems. 
	/*!
		 This class provides a common interface for implementing a ODE system.
	*/

	template<typename Vector, typename Matrix>
	class ODESystemObjectVirtualClass
	{
	public:

		//! Default constructor
		/*!
		    This default constructor should be over-written by the derived classes
		*/
		ODESystemObjectVirtualClass()
		{
			number_of_equations_ = 0;
		}

		//! Equations of the ODE system
		/*!
		  This function must return the equations of the ODE system
		  \param x the current abscissa
		  \param y the current solution (independent variable)
		  \param dy_over_dx the current equation
		*/
		virtual void GetFunctions(const double x, const Vector& y, Vector& dy_over_dx) = 0;

		//! Jacobian matrix of the ODE system
		/*!
		  This function must return the Jacobian matrix of the ODE system
		  \param x the current abscissa
		  \param y the current solution (independent variable)
		  \param J the Jacobian matrix
		*/
		virtual void GetNumericalJacobian(const double x, const Vector& y, Matrix& J) = 0;

		//! Jacobian matrix of the ODE system
		/*!
		  Given the matrix A and the vector b, this function must return the solution z of the linear system Az=b
		  \param A the matrix
		  \param b the vector of known terms
		  \param z the solution of the linear system
		*/
		virtual void GetLinearSystemSolution(const Matrix& A, const Vector& b, Vector& z) = 0;

		//! Print information at each step
		/*!
		  This function is called at every step
		  \param x the current abscissa
		  \param y the current solution (independent variable)
		  \param dy_over_dx the current derivatives
		*/
		virtual void PrintStep(const double x, const Vector& y, const Vector& dy_over_dx) = 0;

		//! Returns the number of equations of the system
		/*!
		  Returns the number of equations of the system
		  \return the number of equations
		*/
		unsigned int number_of_equations() const { return number_of_equations_; }

	protected: 

		unsigned int number_of_equations_ ;		/*!< number of equations */
	};

}

#endif	//ODESMOKE_ODESYSTEMOBJECTCLASS_H