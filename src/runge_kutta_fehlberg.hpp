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
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a21 = 0.25;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a31 = 3./32.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a32 = 9./32.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a41 = 1932./2197.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a42 = -7200./2197.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a43 = 7296./2197.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a51 = 439./216.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a52 = -8.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a53 = 3680./513.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a54 =-845./4104.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a61 = -8./27.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a62 = 2.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a63 = -3544./2565.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a64 = 1859./4104.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::a65 = -11./40.;


	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::c2 = 0.25;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::c3 = 3./8.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::c4 = 12./13.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::c5 = 1.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::c6 = 0.5;


	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::b1Tilde = 25./216.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::b3Tilde = 1408./2565.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::b4Tilde = 2197./4104.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::b5Tilde = -0.2;


	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::b1 = 16./135.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::b3 = 6656./12825.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::b4 = 28561./56430.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::b5 = -9./50.;

	template<typename ODESystemObject, typename Vector>
	const double OdeRungeKuttaFehlberg<ODESystemObject, Vector>::b6 = 2./55.;

	template<typename ODESystemObject, typename Vector>
	OdeRungeKuttaFehlberg<ODESystemObject, Vector>::OdeRungeKuttaFehlberg()
	{
		safety_coefficient_ = 0.70;
		max_enlargement_ratio_ = 4.;
		stabilization_factor_ = 0;

		OdeSMOKE::resize<Vector>(this->number_of_equations_,k2_);
		OdeSMOKE::resize<Vector>(this->number_of_equations_,k3_);
		OdeSMOKE::resize<Vector>(this->number_of_equations_,k4_);
		OdeSMOKE::resize<Vector>(this->number_of_equations_,k5_);
		OdeSMOKE::resize<Vector>(this->number_of_equations_,k6_);
		OdeSMOKE::resize<Vector>(this->number_of_equations_,y4_);

		OdeSMOKE::resize<Vector>(this->number_of_equations_,dy_over_dx_);
		OdeSMOKE::resize<Vector>(this->number_of_equations_,y_h_);
		OdeSMOKE::resize<Vector>(this->number_of_equations_,vb);
	}

	template<typename ODESystemObject, typename Vector>
	void OdeRungeKuttaFehlberg<ODESystemObject, Vector>::SetDefaultConditions()
	{
		OdeSMOKE::resize<Vector>(this->number_of_equations_,y_);
		OdeSMOKE::resize<Vector>(this->number_of_equations_,y0_);
		OdeSMOKE::resize<Vector>(this->number_of_equations_,abs_tolerances_);
		OdeSMOKE::resize<Vector>(this->number_of_equations_,rel_tolerances_);

		for(unsigned int i=0; i< this->number_of_equations_;i++)
			rel_tolerances_[i]=1.e-7;
		
		for(unsigned int i=0; i< this->number_of_equations_;i++)
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
	void OdeRungeKuttaFehlberg<ODESystemObject, Vector>::AdvanceAdaptiveStepSize(const double tF)
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
			if (h_ > 0. && x_ + h_ > xF_) h_ = (xF_ - x_);
			if (h_ < 0. && x_ + h_ < xF_) h_ = (xF_ - x_);

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
	void OdeRungeKuttaFehlberg<ODESystemObject, Vector>::AdvanceFixedStepSize(const unsigned int nsteps)
	{
		for(unsigned int i=1;i<=nsteps;i++)
		{
			number_function_evaluations_ += 6;
			GetFunctions(x_, y_, dy_over_dx_);
			Step(x_, y_, dy_over_dx_, y_h_);
			x_ += h_;
			y_ = y_h_;
		}
	}

	template<typename ODESystemObject, typename Vector>
	void OdeRungeKuttaFehlberg<ODESystemObject, Vector>::AdvanceOverSingleStep()
	{
		// Searching for the most appropriate step
		for(;;)
		{
			number_function_evaluations_ += 5;

			Step(x_,y_,dy_over_dx_,y_h_);

			// Searching for the maximum error between the 2 solutions
			// This is to estimate the local error
			double maximum = 1.e-64;
			for(unsigned int i=0; i< this->number_of_equations_;i++)
			{
				const double aux = fabs(y_h_[i]-y4_[i])/(abs_tolerances_[i] + rel_tolerances_[i]*fabs(y_h_[i]));
				if(aux > maximum)
					maximum = aux;
			}
		
			// Estimation of the local error
			const double E = fabs(maximum);
			double h_new = safety_coefficient_*h_/pow(E,0.20-stabilization_factor_*0.75);

			// The new step cannot be too much larger than the previous one 
			// by a factor equal to max_enlargement_ratio (default 4)
			if(fabs(h_new) > max_enlargement_ratio_*fabs(h_))
				h_new = max_enlargement_ratio_*h_;

			// Check with respect to the user-defined min and max steps
			h_new = std::min(h_new, max_step_);
			h_new = std::max(h_new, min_step_);
			
			// In case the local error is large, the new step cannot be accepted and 
			// a new estimation is needed
			if(E > 1.)
			{
				number_rejected_steps_++;
				h_ = h_new;
			}
			// If the local error is sufficiently small, the new estimation can be accepted
			else
			{
				number_accepted_steps_++;
				
				y_ = y_h_;
				x_ += h_;
				h_ = h_new;
				break;
			}
		}
	}

	template<typename ODESystemObject, typename Vector>
	void OdeRungeKuttaFehlberg<ODESystemObject, Vector>::Step(const double x_in, const Vector& y_in, const Vector& dy_over_dx_in, Vector& y_out)
	{
		const double a21star = a21*h_;
		const double a31star = a31*h_;
		const double a32star = a32*h_;
		const double a41star = a41*h_;
		const double a42star = a42*h_;
		const double a43star = a43*h_;
		const double a44star = a44*h_;
		const double a51star = a51*h_;
		const double a52star = a52*h_;
		const double a53star = a53*h_;
		const double a54star = a54*h_;
		const double a61star = a61*h_;
		const double a62star = a62*h_;
		const double a63star = a63*h_;
		const double a64star = a64*h_;
		const double a65star = a65*h_;

		const double b1TildeStar = b1Tilde*h_;
		const double b3TildeStar = b3Tilde*h_;
		const double b4TildeStar = b4Tilde*h_;
		const double b5TildeStar = b5Tilde*h_;

		const double b1Star = b1*h_;
		const double b3Star = b3*h_;
		const double b4Star = b4*h_;
		const double b5Star = b5*h_;
		const double b6Star = b6*h_;

		//k1_ = dy_over_dx_in;
		//vb = y_in+a21*k1_;
		sum_plus_scalar_multiplication(&vb, y_in, a21star, dy_over_dx_in);
		GetFunctions(x_in+h_*c2, vb, k2_);

		//vb = y_in+a31*k1_+a32*k2_;
		sum_plus_scalar_multiplication(&vb, y_in, a31star, dy_over_dx_in, a32star, k2_);
		GetFunctions(x_in+h_*c3, vb, k3_);

		//vb = y_in+a41*k1_+a42*k2_+a43*k3_;
		sum_plus_scalar_multiplication(&vb, y_in, a41star, dy_over_dx_in, a42star, k2_, a43star, k3_);
		GetFunctions(x_in+h_*c4, vb, k4_);

		//vb = y_in+a51*k1_+a52*k2_+a53*k3_+a54*k4_;
		sum_plus_scalar_multiplication(&vb, y_in, a51star, dy_over_dx_in, a52star, k2_, a53star, k3_, a54star, k4_);
		GetFunctions(x_in+h_*c5, vb, k5_);

		//vb = y_in+a61*k1_+a62*k2_+a63*k3_+a64*k4_+a65*k5_;
		sum_plus_scalar_multiplication(&vb, y_in, a61star, dy_over_dx_in, a62star, k2_, a63star, k3_, a64star, k4_, a65star, k5_);
		GetFunctions(x_in+h_*c6, vb, k6_);

		sum_plus_scalar_multiplication(&y4_, y_in, b1TildeStar, dy_over_dx_in,	b3TildeStar, k3_, 
												   b4TildeStar, k4_,			b5TildeStar, k5_);

		sum_plus_scalar_multiplication(&y_out, y_in, b1Star, dy_over_dx_in, b3Star, k3_, 
													 b4Star, k4_, b5Star, k5_, b6Star, k6_);

		//y4_   = y_in + b1Star*k1_+b2Star*k2_+b3Star*k3_+b4Star*k4_+b5Star*k5_+b6Star*k6_;

		//y_out = y_in + b1*k1_+b2*k2_+b3*k3_+b4*k4_+b5*k5_+b6*k6_;


	/*	k1_ = h_*dy_over_dx_in;
		
		vb = y_in+a21*k1_;
		GetFunctions(x_in+h_*c2, vb, k2_);
		k2_ *= h_;

		vb = y_in+a31*k1_+a32*k2_;
		GetFunctions(x_in+h_*c3, vb, k3_);
		k3_ *= h_;

		vb = y_in+a41*k1_+a42*k2_+a43*k3_;
		GetFunctions(x_in+h_*c4, vb, k4_);
		k4_ *= h_;

		vb = y_in+a51*k1_+a52*k2_+a53*k3_+a54*k4_;
		GetFunctions(x_in+h_*c5, vb, k5_);
		k5_ *= h_;

		vb = y_in+a61*k1_+a62*k2_+a63*k3_+a64*k4_+a65*k5_;
		GetFunctions(x_in+h_*c6, vb, k6_);
		k6_ *= h_;

		y4_   = y_in + b1Star*k1_+b2Star*k2_+b3Star*k3_+b4Star*k4_+b5Star*k5_+b6Star*k6_;

		y_out = y_in + b1*k1_+b2*k2_+b3*k3_+b4*k4_+b5*k5_+b6*k6_;
		*/
	}

	template<typename ODESystemObject, typename Vector>
	void OdeRungeKuttaFehlberg<ODESystemObject, Vector>::Info() const
	{
		std::cout << " * Number of steps:                " << number_steps_ << std::endl;
		std::cout << " * Number of accepted steps:       " << number_accepted_steps_ << std::endl;
		std::cout << " * Number of rejected steps:       " << number_rejected_steps_ << std::endl;
		std::cout << " * Number of function evaluations: " << number_function_evaluations_ << std::endl;
	}

}