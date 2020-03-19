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

namespace OdeSMOKE
{
	template<typename Vector, typename Method>
	OdeRungeKuttaFamily<Vector,Method>::OdeRungeKuttaFamily()
	{
		this->SetDefaultConditions();

		this->number_steps_ = 0;
		this->number_function_evaluations_ = 0;
		this->number_accepted_steps_ = 0;
		this->number_rejected_steps_ = 0;

		this->user_defined_max_step_ = false;
		this->user_defined_min_step_ = false;
		this->user_defined_first_step_ = false;
	}

	template<typename Vector, typename Method>
	void OdeRungeKuttaFamily<Vector,Method>::SetSafetyCoefficient(const double safety_coefficient)
	{
		if (safety_coefficient >1. || safety_coefficient <= 0.)
			this->FatalError("The safety coefficient must be strictly larger than 0. and smaller than 1.");
		this->safety_coefficient_ = safety_coefficient;
	}

	template<typename Vector, typename Method>
	void OdeRungeKuttaFamily<Vector,Method>::SetMaximumEnlargementRatio(const double max_enlargement_ratio)
	{
		if (max_enlargement_ratio < 1.)
			this->FatalError("The maximum enlargement ratio must be larger than 1.");
		this->max_enlargement_ratio_ = max_enlargement_ratio;
	}

	template<typename Vector, typename Method>
	void OdeRungeKuttaFamily<Vector,Method>::SetStabilizationFactor(const double stabilization_factor)
	{
		if (stabilization_factor < 0. || stabilization_factor > 0.1)
			this->FatalError("The stabilization factor must be larger (or equal) than 0. and smaller than 0.1");
		this->stabilization_factor_ = stabilization_factor;
	}

	template<typename Vector, typename Method>
	void OdeRungeKuttaFamily<Vector,Method>::Solve(const double x0, const double xF)
	{
		this->x0_ = x0;
		this->xF_ = xF;
		this->x_ = x0;
		this->y_ = this->y0_;

		// First step
		if (this->user_defined_first_step_ == false)
			this->h_ = (this->xF_ - this->x0_)/20.;
		//else
		//	h_ = EstimateFirstStep();

		// Maximum step
		if (this->user_defined_max_step_ == false)
			this->max_step_ = fabs(this->xF_ - this->x0_);

		// Minimum step
		if (this->user_defined_min_step_ == false)
			this->min_step_ = 0.;

		// Advance
		this->AdvanceAdaptiveStepSize(this->xF_);
	}

	template<typename Vector, typename Method>
	void OdeRungeKuttaFamily<Vector,Method>::Solve(const double x0, const double dx, const unsigned int nsteps)
	{
		this->x0_ = x0;
		this->xF_ = this->x0_+nsteps*dx;
		this->x_ = this->x0_;
		this->y_ = this->y0_;
		this->h_ = dx;

		this->AdvanceFixedStepSize(nsteps);
	}

	template<typename Vector, typename Method>
	unsigned int OdeRungeKuttaFamily<Vector,Method>::number_steps() const
	{
		return this->number_steps_;
	}

	template<typename Vector, typename Method>
	unsigned int OdeRungeKuttaFamily<Vector,Method>::number_accepted_steps() const
	{
		return this->number_accepted_steps_;
	}

	template<typename Vector, typename Method>
	unsigned int OdeRungeKuttaFamily<Vector,Method>::number_rejected_steps() const
	{
		return this->number_rejected_steps_;
	}

	template<typename Vector, typename Method>
	unsigned int OdeRungeKuttaFamily<Vector,Method>::number_function_evaluations() const
	{
		return this->number_function_evaluations_;
	}
}