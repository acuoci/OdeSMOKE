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
#define EIGEN_RUNTIME_NO_MALLOC // Define this symbol to enable runtime tests for allocations
#define MKL_SUPPORT 1
#include <iostream>
#include <stdint.h>
#include <iomanip>

#include <Eigen/Dense>

#include "ode_solver_virtual_class.h"
#include "runge_kutta_4th.h"
#include "runge_kutta_fehlberg.h"

#include "ode_system_test_01.h"
#include "ode_system_test_02.h"
#include "ode_system_test_03.h"
#include "ode_system_test_04.h"
#include "ode_system_test_05.h"
#include "ode_system_test_06.h"
#include "ode_system_test_07.h"
#include "ode_system_test_08.h"
#include "ode_system_test_09.h"
#include "ode_system_test_10.h"
#include "ode_system_test_11.h"
#include "ode_system_test_12.h"
#include "ode_system_test_13.h"
#include "ode_system_test_14.h"
#include "ode_system_test_15.h"
#include "ode_system_lorenz.h"
#include "ode_system_test_diffusion.h"

#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;
#define n_equations 3000
#define n_steps 20000

typedef boost::array< double , n_equations > boost_array;
//typedef boost_array state_type;
typedef std::vector<double> state_type;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
	for(unsigned int i=0;i<n_equations/3;i++)
	{
		unsigned int index = 3*i;
		dxdt[0+index] = sigma * ( x[1+index] - x[0+index] );
		dxdt[1+index] = R * x[0+index] - x[1+index] - x[0+index] * x[2+index];
		dxdt[2+index] = -b * x[2+index] + x[0+index] * x[1+index];
	}
}

void write_lorenz( const state_type &x , const double t )
{
   // cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
}

/*
#include "BzzMath.hpp"
void OdeSystem_01(BzzVector &y,double x,BzzVector &dy_over_dx);
void OdeSystem_05(BzzVector &y,double x,BzzVector &dy_over_dx);
void OdeSystem_Diffusion(BzzVector &y,double x,BzzVector &dy_over_dx);

// Solve BzzMath
void SolveBzzMath(const Eigen::VectorXd& initial_condition, const double xF, BzzVector& y, void pt(BzzVector &y,double x,BzzVector &dy_over_dx))
{
		BzzVector y0(initial_condition.size());
		ChangeDimensions(initial_condition.size(), &y); 
		for(unsigned int i=1;i<=initial_condition.size();i++)
			y0[i]=initial_condition(i-1);
		BzzOdeRK o(y0,0.,*pt);
		o.SetTollAbs(1.e-12);
		o.SetTollRel(1.e-7);
		y = o(xF);

		o.BzzPrint();
}

// Estimation of the error
void EstimationOfError(const Eigen::VectorXd& odesmoke_solution, const BzzVector& bzzmath_solution, const Eigen::VectorXd& analytical_solution)
{
		Eigen::VectorXd e1 = odesmoke_solution - analytical_solution;
		std::cout << "OdeSMOKE Norm2: " << e1.lpNorm<2>() << std::endl;
		Eigen::VectorXd e2(odesmoke_solution.size());
		for(unsigned int i=1;i<=odesmoke_solution.size();i++)
			e2(i-1) = bzzmath_solution[i]-analytical_solution(i-1);
		std::cout << "BzzMath  Norm2: " << e2.lpNorm<2>() << std::endl;
}
*/
template<typename ode_solver_type, typename ode_system_type>
void PrintSummary(const std::string tag, const ode_solver_type& ode_solver_rk, const ode_system_type& ode_system)
{
	std::cout << "--------------------------------------------------------------------------------------" << std::endl;
	std::cout << "                             " << tag                                                 << std::endl;
	std::cout << "--------------------------------------------------------------------------------------" << std::endl;

	Eigen::VectorXd e(ode_solver_rk.y().size());
	Eigen::VectorXd y(ode_solver_rk.y().size());

	for(unsigned int i=0;i<ode_solver_rk.y().size();i++)
		y[i] = ode_solver_rk.y()[i];

	for(unsigned int i=0;i<ode_solver_rk.y().size();i++)
		e[i] = ode_solver_rk.y()[i] - ode_system.analytical_solution(ode_system.final_abscissa())[i];

	std::cout << " * Error norm2:              " << e.lpNorm<2>() << std::endl;
    std::cout << " * Solution norm2:           " << y.lpNorm<2>() << std::endl;
    std::cout << " * Ratio error/solution:     " << e.lpNorm<2>()/y.lpNorm<2>() << std::endl;
	std::cout << " * Number of steps:          " << ode_solver_rk.number_steps() << std::endl;
	std::cout << " * Number of rejected steps: " << ode_solver_rk.number_rejected_steps() << std::endl;
	std::cout << " * Number of functions:      " << ode_solver_rk.number_function_evaluations() << std::endl;

	std::cout << " * Comparison between solutions" << std::endl;
	for(unsigned int i=0;i<e.size();i++)
        std::cout << "   " << y[i] << " " <<  ode_system.analytical_solution(ode_system.final_abscissa())[i] << std::endl;

	std::cout << std::endl;
}

int main()
{
	std::cout.setf(std::ios::scientific);

	// Test 01
	{
		typedef Eigen::VectorXd vector_type;

		typedef ODESystemObject_Test_01<vector_type, vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 01", ode_solver_rk, ode_solver_rk);
	}
	
	// Test 02
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 1 > vector_type;

		typedef ODESystemObject_Test_02<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 02", ode_solver_rk, ode_solver_rk);
	}

	// Test 03
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 2 > vector_type;

		typedef ODESystemObject_Test_03<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 03", ode_solver_rk, ode_solver_rk);
	}

	// Test 04
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 2 > vector_type;

		typedef ODESystemObject_Test_04<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());

		PrintSummary("Test 04", ode_solver_rk, ode_solver_rk);
	}

	// Test 05
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 3 > vector_type;

		typedef ODESystemObject_Test_05<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());

		PrintSummary("Test 05", ode_solver_rk, ode_solver_rk);
	}

	// Test 06
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 1 > vector_type;

		typedef ODESystemObject_Test_06<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());

		PrintSummary("Test 06", ode_solver_rk, ode_solver_rk);
	}

	// Test 07
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 1 > vector_type;

		typedef ODESystemObject_Test_07<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 07", ode_solver_rk, ode_solver_rk);
	}

	// Test 08
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 1 > vector_type;

		typedef ODESystemObject_Test_08<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 08", ode_solver_rk, ode_solver_rk);
	}

	// Test 09
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 1 > vector_type;

		typedef ODESystemObject_Test_09<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 09", ode_solver_rk, ode_solver_rk);
	}

	// Test 10
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 1 > vector_type;

		typedef ODESystemObject_Test_10<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 10", ode_solver_rk, ode_solver_rk);
	}

	// Test 11
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 1 > vector_type;

		typedef ODESystemObject_Test_11<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 11", ode_solver_rk, ode_solver_rk);
	}

	// Test 12
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 1 > vector_type;

		typedef ODESystemObject_Test_12<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 12", ode_solver_rk, ode_solver_rk);
	}

	// Test 13
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 3 > vector_type;

		typedef ODESystemObject_Test_13<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 13", ode_solver_rk, ode_solver_rk);
	}

	// Test 14
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 4 > vector_type;

		typedef ODESystemObject_Test_14<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 14", ode_solver_rk, ode_solver_rk);
	}

	// Test 15
	{
		//typedef Eigen::VectorXd vector_type;
		typedef boost::array< double , 1 > vector_type;

		typedef ODESystemObject_Test_15<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method_RungeKutta4thOrder;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 15", ode_solver_rk, ode_solver_rk);
	}

	
	// Test 14/Fehlberg
	{
		//typedef Eigen::VectorXd vector_type;
		//typedef Eigen::Vector4d vector_type;
		typedef std::vector<double> vector_type;
		//typedef boost::array< double , 4 > vector_type;

		typedef ODESystemObject_Test_14<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKuttaFehlberg< ode_system, vector_type > Method_RungeKuttaFehlberg;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method_RungeKuttaFehlberg> Solver_RungeKutta4thOrder;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver_RungeKutta4thOrder > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
		ode_solver_rk.Solve(0, ode_solver_rk.final_abscissa());
		
		PrintSummary("Test 14 (Fehlberg)", ode_solver_rk, ode_solver_rk);
	}
	
	
	// Lorenz attractor
	{
		typedef boost_array vector_type;
	//	typedef Eigen::VectorXd vector_type;

		const double t1 = OdeSMOKE::GetCpuTime();
		typedef ODESystemObject_Test_Lorenz<vector_type> ode_system;
		typedef OdeSMOKE::OdeRungeKuttaFehlberg< ode_system, vector_type > Method;
//		typedef OdeSMOKE::OdeRungeKutta4th< ode_system, vector_type > Method;
		typedef OdeSMOKE::OdeRungeKuttaFamily< vector_type, Method> Solver;
		OdeSMOKE::ODESolverVirtualClass< vector_type, Solver > ode_solver_rk;
		ode_solver_rk.SetInitialConditions(ode_solver_rk.initial_condition());
//		ode_solver_rk.Solve(0, 24.9);
		ode_solver_rk.Solve(0, 0.1, n_steps);
        const double t2 = OdeSMOKE::GetCpuTime();
//		PrintSummary("Test Lorenz attractor", ode_solver_rk, ode_system);

		// state_type = double
        const double t3 = OdeSMOKE::GetCpuTime();
        runge_kutta4< state_type > stepper;

		{
            state_type x(n_equations);
			
			for(unsigned int i=0;i<n_equations/3;i++)
			{
				x[0+3*i] = 10.;
				x[1+3*i] = 1.;
				x[2+3*i] = 1.;
			}
			integrate_const(	 stepper, lorenz , x , 0.0 , n_steps*0.1 , 0.1 , write_lorenz );
		}
        const double t4 = OdeSMOKE::GetCpuTime();

		std::cout << "OdeSMOKE " << t2-t1 << std::endl;
		std::cout << "OdeInt " << t4-t3 << std::endl;
	}
	
	/*
	// Test Predator-Prey
	{
		ODESystemObject_Test_Diffusion<Eigen::VectorXd, Eigen::MatrixXd> ode_system;
		OdeSMOKE::OdeRungeKutta4th<ODESystemObject_Test_Diffusion<Eigen::VectorXd, Eigen::MatrixXd>, Eigen::VectorXd > ode_solver_rk(&ode_system);
		ode_solver_rk.SetInitialConditions(ode_system.initial_condition());
		ode_solver_rk.Solve(0, ode_system.final_abscissa());
	}
	*/
	/*
	// Test 01
	{
		ODESystemObject_Test_01<Eigen::VectorXd, Eigen::MatrixXd> ode_system;
		OdeSMOKE::OdeRungeKutta4th<ODESystemObject_Test_01<Eigen::VectorXd, Eigen::MatrixXd>, Eigen::VectorXd > ode_solver_rk(&ode_system);
		ode_solver_rk.SetInitialConditions(ode_system.initial_condition());

		const double tStart = OdeSMOKE::GetCpuTime();
		ode_solver_rk.Solve(0, ode_system.final_abscissa());
		const double tEnd = OdeSMOKE::GetCpuTime();
		std::cout << "Test 01 (OdeSMOKE), Cpu Time: " << (tEnd - tStart)*1000. << " ms" << std::endl; 


		// BzzMath
		BzzVector y;
		SolveBzzMath(ode_system.initial_condition(), ode_system.final_abscissa(), y,OdeSystem_01);
		EstimationOfError(ode_solver_rk.y(), y, ode_system.analytical_solution(ode_system.final_abscissa()) );
	}
	
	{
		ODESystemObject_Test_05<Eigen::VectorXd, Eigen::MatrixXd> ode_system;
		OdeSMOKE::OdeRungeKutta4th<ODESystemObject_Test_05<Eigen::VectorXd, Eigen::MatrixXd>, Eigen::VectorXd > ode_solver_rk(&ode_system);
		ode_solver_rk.SetInitialConditions(ode_system.initial_condition());

		const double tStart = OdeSMOKE::GetCpuTime();
		ode_solver_rk.Solve(0, ode_system.final_abscissa());
		const double tEnd = OdeSMOKE::GetCpuTime();
		std::cout << "Test 05 (OdeSMOKE), Cpu Time: " << (tEnd - tStart)*1000. << " ms" << std::endl; 
		ode_solver_rk.Info();

		// BzzMath
		BzzVector y;
		SolveBzzMath(ode_system.initial_condition(), ode_system.final_abscissa(), y,OdeSystem_05);
		EstimationOfError(ode_solver_rk.y(), y, ode_system.analytical_solution(ode_system.final_abscissa()) );
	}
	/*
	{
		ODESystemObject_Test_05<Eigen::Vector3d, Eigen::Matrix3d> ode_system;
		OdeSMOKE::OdeRungeKutta4th<ODESystemObject_Test_05<Eigen::Vector3d, Eigen::Matrix3d>, Eigen::Vector3d > ode_solver_rk(&ode_system);
		ode_solver_rk.SetInitialConditions(ode_system.initial_condition());

		const double tStart = OdeSMOKE::GetCpuTime();
		ode_solver_rk.Solve(0, ode_system.final_abscissa());
		const double tEnd = OdeSMOKE::GetCpuTime();
		std::cout << "Test 05 (OdeSMOKE, 3d), Cpu Time: " << (tEnd - tStart)*1000. << " ms" << std::endl;

		// BzzMath
		BzzVector y;
		SolveBzzMath(ode_system.initial_condition(), ode_system.final_abscissa(), y,OdeSystem_05);
		EstimationOfError(ode_solver_rk.y(), y, ode_system.analytical_solution(ode_system.final_abscissa()) );
	}
	*/
	/*
	{
		const double tStart = OdeSMOKE::GetCpuTime();
		for(unsigned int i=1;i<=1000;i++)
		{
			ODESystemObject_Test_05<Eigen::Vector3d, Eigen::Matrix3d> ode_system;
			OdeSMOKE::OdeRungeKutta4th<ODESystemObject_Test_05<Eigen::Vector3d, Eigen::Matrix3d>, Eigen::Vector3d > ode_solver_rk(&ode_system);
			ode_solver_rk.SetInitialConditions(ode_system.initial_condition());
			ode_solver_rk.Solve(0, ode_system.final_abscissa());
		}
		const double tEnd = OdeSMOKE::GetCpuTime();
		std::cout << "Test 05 (OdeSMOKE, 3d), Cpu Time: " << (tEnd - tStart)*1000. << " ms" << std::endl;

		{
			ODESystemObject_Test_05<Eigen::Vector3d, Eigen::Matrix3d> ode_system;
			BzzVector y;
			BzzVector y0(ode_system.initial_condition().size());
			ChangeDimensions(ode_system.initial_condition().size(), &y); 
			for(unsigned int i=1;i<=ode_system.initial_condition().size();i++)
				y0[i]=ode_system.initial_condition()(i-1);

			const double tStart = OdeSMOKE::GetCpuTime();
			for(unsigned int i=1;i<=1000;i++)
			{
				BzzOdeRK o(y0,0.,OdeSystem_05);
				o.SetTollAbs(1.e-12);
				o.SetTollRel(1.e-7);
				y = o(ode_system.final_abscissa());
			}
			const double tEnd = OdeSMOKE::GetCpuTime();
			std::cout << "Test 05 (BzzMath, 3d), Cpu Time: " << (tEnd - tStart)*1000. << " ms" << std::endl;
		}
	}
	*/

	return EXIT_SUCCESS;
}
/*
void OdeSystem_01(BzzVector &y,double x,BzzVector &dy_over_dx)
{
	dy_over_dx[1] = -y[1];
}

void OdeSystem_05(BzzVector &y,double x,BzzVector &dy_over_dx)
{
	dy_over_dx[1] = -0.1*y[1]-49.9*y[2];
	dy_over_dx[2] = -50.*y[2];
	dy_over_dx[3] =  70.*y[2]-120.*y[3];
}

void OdeSystem_Diffusion(BzzVector &y,double x,BzzVector &dy_over_dx)
{
	const double dx2 = 1./double(y.Size()-1)/double(y.Size()-1);

	dy_over_dx[1] = 0.;
	for(unsigned int i=2;i<y.Size();i++)
		dy_over_dx[i] = 1e-10*(y[i+1]-2.*y[i]+y[i-1])/dx2;
	dy_over_dx[y.Size()] = 0.;
}
*/
