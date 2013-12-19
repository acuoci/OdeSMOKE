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
	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKutta4th<ODESystemObject, Vector>::a21 = 0.5;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKutta4th<ODESystemObject, Vector>::a32 = 0.5;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKutta4th<ODESystemObject, Vector>::c2 = 0.5;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKutta4th<ODESystemObject, Vector>::c3 = 0.5;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKutta4th<ODESystemObject, Vector>::c4 = 1.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKutta4th<ODESystemObject, Vector>::b1 = 1./6.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKutta4th<ODESystemObject, Vector>::b2 = 1./3.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKutta4th<ODESystemObject, Vector>::b3 = 1./3.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKutta4th<ODESystemObject, Vector>::b4 = 1./6.;


	template<typename ODESystemObject, typename Vector>
	OdeRungeKutta4th<ODESystemObject, Vector>::OdeRungeKutta4th()
	{
		safety_coefficient_ = 0.70;
		max_enlargement_ratio_ = 4.;
		stabilization_factor_ = 0.;

		OdeSMOKE::resize<Vector>(number_of_equations_,vb);
		OdeSMOKE::resize<Vector>(number_of_equations_,k2_);
		OdeSMOKE::resize<Vector>(number_of_equations_,k3_);
		OdeSMOKE::resize<Vector>(number_of_equations_,k4_);

		OdeSMOKE::resize<Vector>(number_of_equations_,dy_over_dx_);
		OdeSMOKE::resize<Vector>(number_of_equations_,dy_over_dx_old_);
		OdeSMOKE::resize<Vector>(number_of_equations_,y_2h_);
		OdeSMOKE::resize<Vector>(number_of_equations_,y_h_);
	}

	template<typename ODESystemObject, typename Vector>
	void OdeRungeKutta4th<ODESystemObject, Vector>::SetDefaultConditions()
	{
		OdeSMOKE::resize<Vector>(number_of_equations_,y_);
		OdeSMOKE::resize<Vector>(number_of_equations_,y0_);
		OdeSMOKE::resize<Vector>(number_of_equations_,abs_tolerances_);
		OdeSMOKE::resize<Vector>(number_of_equations_,rel_tolerances_);

		for(unsigned int i=0; i<number_of_equations_;i++)
			rel_tolerances_[i]=1.e-7;
		
		for(unsigned int i=0; i<number_of_equations_;i++)
			abs_tolerances_[i]=1.e-12;

		first_step_   = 0.;
		max_step_ = 1.e16;
		min_step_ = 1.e-16;
		x0_ = 0.;
		xF_ = 0.;
		max_number_steps_ = 500000;
		analytical_jacobian_ = false;

		user_defined_max_step_ = false;
		user_defined_min_step_ = false;
		user_defined_first_step_ = false;
	}

	template<typename ODESystemObject, typename Vector>
	void OdeRungeKutta4th<ODESystemObject, Vector>::AdvanceAdaptiveStepSize(const double tF)
	{
		GetFunctions(x_, y_, dy_over_dx_);
		number_function_evaluations_++;
		PrintStep(x_, y_, dy_over_dx_);

		for(;;)
		{
		//	if (number_steps_ > maximum_number_of_steps_)
		//		ErrorMessage();

			// Check if the next step is over the final abscissa
			// If yes automatic reduction of the step is done to match excatly the final abscissa
			if (h_ > 0. && x_ + 2.*h_ > xF_) h_ = 0.5*(xF_ - x_);
			if (h_ < 0. && x_ + 2.*h_ < xF_) h_ = 0.5*(xF_ - x_);

			// Advance one step
			AdvanceOverSingleStep();
			PrintStep(x_, y_, dy_over_dx_);
			number_steps_++;
			
			// If the current step is the final abscissa the required integration is done
			if ( fabs(x_-xF_) <= fabs(xF_)*10.*MACHINE_EPSILON_FLOAT)
				break; 

			// Update the functions
			GetFunctions(x_, y_, dy_over_dx_);
			number_function_evaluations_++;
		}
	}

	template<typename ODESystemObject, typename Vector>
	void OdeRungeKutta4th<ODESystemObject, Vector>::AdvanceFixedStepSize(const unsigned int nsteps)
	{
		PrintStep(x_, y_, dy_over_dx_);

		for(unsigned int i=1;i<=nsteps;i++)
		{
			number_function_evaluations_ += 4;
			GetFunctions(x_, y_, dy_over_dx_);
			Step(x_, y_, dy_over_dx_, y_h_);
			x_ += h_;
			y_ = y_h_;
			PrintStep(x_, y_, dy_over_dx_);
		}
	}

	template<typename ODESystemObject, typename Vector>
	void OdeRungeKutta4th<ODESystemObject, Vector>::AdvanceOverSingleStep()
	{
		const double x_old = x_;
		double h_old = h_;
	
		dy_over_dx_old_ = dy_over_dx_;

		// Searching for the most appropriate step
		// The Richardson extrapolation is used
		for(;;)
		{
			number_function_evaluations_ += 10;

			// Perform 2 integrations with step equl to h
			// from y_ to y_h_
			Step(x_old,y_,dy_over_dx_old_,y_2h_);
			x_ = x_old + h_;
			GetFunctions(x_,y_2h_,dy_over_dx_);
			Step(x_, y_2h_,dy_over_dx_,y_h_);

			// Perform a single integration with step equal to 2h
			// from y_ to y_2h_
			h_ *= 2.;
			Step(x_old, y_, dy_over_dx_old_,y_2h_);

			// Searching for the maximum error between the 2 solutions
			// This is to estimate the local error
			double maximum = 1.e-32;
			for(unsigned int i=0; i<number_of_equations_;i++)
			{
				const double aux = fabs(y_h_[i]-y_2h_[i])/(abs_tolerances_[i] + rel_tolerances_[i]*fabs(y_h_[i]));
				if(aux > maximum)
					maximum = aux;
			}
			maximum /= 30.; // The number 30 is from [2^(p+1)-2]
		
			// Estimation of the local error
			const double E = fabs(maximum);
			double h_new = safety_coefficient_*h_old/pow(E,0.20-stabilization_factor_*0.75);

			// The new step cannot be too much larger than the previous one 
			// by a factor equal to max_enlargement_ratio (default 4)
			if(fabs(h_new) > max_enlargement_ratio_*fabs(h_old))
				h_new = max_enlargement_ratio_*h_old;

			// Check with respect to the user-defined min and max steps
			h_new = std::min(h_new, max_step_);
			h_new = std::max(h_new, min_step_);
			
			// In case the local error is large, the new step cannot be accepted and 
			// a new estimation is needed
			if(E > 1.)
			{
				number_rejected_steps_++;

				h_old = h_ = h_new;
				x_ = x_old;
			}
			// If the local error is sufficiently small, the new estimation can be accepted
			// Moreover, according to the Richardson extrapolation, the new estimation can be improved
			else
			{
				number_accepted_steps_++;
				
				for(unsigned int i=0; i<number_of_equations_;++i)
					y_[i] = y_h_[i]+(y_h_[i]-y_2h_[i])/15.;

				x_ = x_old + h_;
				h_ = h_new;
				break;
			}
		}
	}

	template<typename ODESystemObject, typename Vector>
	void OdeRungeKutta4th<ODESystemObject, Vector>::Step(const double x_in, const Vector& y_in, const Vector& dy_over_dx_in, Vector& y_out)
	{
		const double a21star = a21*h_;
		const double a32star = a32*h_;
		const double b1star  = b1*h_;
		const double b2star  = b2*h_;
		const double b3star  = b3*h_;
		const double b4star  = b4*h_;

		//k1_ = dy_over_dx_in;
		//vb = y_in+a21*k1_*h_;
		sum_plus_scalar_multiplication(&vb, y_in, a21star, dy_over_dx_in);
		GetFunctions(x_in+h_*c2, vb, k2_);

		//vb = y_in+a32*k2_*h_;
		sum_plus_scalar_multiplication(&vb, y_in, a32star, k2_);
		GetFunctions(x_in+h_*c3, vb, k3_);

		//vb = y_in+k3_*h_;
		sum_plus_scalar_multiplication(&vb, y_in, h_, k3_);
		GetFunctions(x_in+h_*c4, vb, k4_);

		//y_out = y_in + (b1*k1_+b2*k2_+b3*k3_+b4*k4_)*h_;
		sum_plus_scalar_multiplication(&y_out, y_in, b1star, dy_over_dx_in, b2star, k2_, b3star, k3_, b4star, k4_);
	}

	template<typename ODESystemObject, typename Vector>
	void OdeRungeKutta4th<ODESystemObject, Vector>::Info() const
	{
		std::cout << "Number of steps:                " << number_steps_ << std::endl;
		std::cout << "Number of accepted steps:       " << number_accepted_steps_ << std::endl;
		std::cout << "Number of rejected steps:       " << number_rejected_steps_ << std::endl;
		std::cout << "Number of function evaluations: " << number_function_evaluations_ << std::endl;
	}
}