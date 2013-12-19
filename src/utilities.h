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

#ifndef ODESMOKE_UTILITIES_H
#define	ODESMOKE_UTILITIES_H

#include <cstdint>
#include <time.h>

#if MKL_SUPPORT == 1
	#include <mkl.h>
#endif

#include <Eigen/Dense>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

namespace OdeSMOKE
{
	/*
	 * by default any type is not resizable
	 */
	template< class Container>
	struct is_resizeable
	{
		typedef boost::false_type type;
		const static bool value = type::value;
	};

	/*
	 * specialization for std::vector
	 */
	template< class V>
	struct is_resizeable< std::vector<V> >
	{
		//struct type : public boost::true_type { };
		typedef boost::true_type type;
		const static bool value = type::value;
	};

	template<>
	struct is_resizeable< Eigen::VectorXd >
	{
		//struct type : public boost::true_type { };
		typedef boost::true_type type;
		const static bool value = type::value;
	};

	template< class Vector>
	void resize_implementation(const unsigned int n, Vector& v, const boost::true_type )
	{
		v.resize(n);
	};

	template< class Vector>
	void resize_implementation(const unsigned int n, Vector& v, const boost::false_type )
	{
	};

	template< class Vector>
	void resize(const unsigned int n, Vector& v )
	{
		resize_implementation(n, v, typename OdeSMOKE::is_resizeable< Vector >::type() );
	}

	template< typename Vector>
	void sum_plus_scalar_multiplication(Vector* result, const Vector& v, const double alpha1, const Vector& v1)
	{
		#if MKL_SUPPORT == 1
			cblas_dcopy(result->size(), &v[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha1, &v1[0], 1, &(*result)[0], 1);
		#else
			double* pt_result = &(*result)[0];
			const double* pt_v = &v[0];
			const double* pt_v1 = &v1[0];
			for(unsigned int i=0;i<result->size();++i)
				*pt_result++ = *pt_v++ + alpha1*(*pt_v1++);
		#endif
	}

	template< typename Vector>
	void sum_plus_scalar_multiplication(Vector* result, const Vector& v, const double alpha1, const Vector& v1, 
		                                                                 const double alpha2, const Vector& v2 )
	{
		#if MKL_SUPPORT == 1
			cblas_dcopy(result->size(), &v[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha1, &v1[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha2, &v2[0], 1, &(*result)[0], 1);
		#else
			double* pt_result = &(*result)[0];
			const double* pt_v = &v[0];
			const double* pt_v1 = &v1[0];
			const double* pt_v2 = &v2[0];
			for(unsigned int i=0;i<result->size();++i)
				*pt_result++ = *pt_v++ + alpha1*(*pt_v1++) + alpha2*(*pt_v2++);
		#endif
	}

	template< typename Vector>
	void sum_plus_scalar_multiplication(Vector* result, const Vector& v, const double alpha1, const Vector& v1, 
		                                                                 const double alpha2, const Vector& v2,
		                                                                 const double alpha3, const Vector& v3 )	
	{
		#if MKL_SUPPORT == 1
			cblas_dcopy(result->size(), &v[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha1, &v1[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha2, &v2[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha3, &v3[0], 1, &(*result)[0], 1);
		#else
			double* pt_result = &(*result)[0];
			const double* pt_v = &v[0];
			const double* pt_v1 = &v1[0];
			const double* pt_v2 = &v2[0];
			const double* pt_v3 = &v3[0];
			for(unsigned int i=0;i<result->size();++i)
				*pt_result++ = *pt_v++ + alpha1*(*pt_v1++) + alpha2*(*pt_v2++) + alpha3*(*pt_v3++);
		#endif
	}

	template< typename Vector>
	void sum_plus_scalar_multiplication(Vector* result, const Vector& v, const double alpha1, const Vector& v1, 
		                                                                 const double alpha2, const Vector& v2,
		                                                                 const double alpha3, const Vector& v3,
		                                                                 const double alpha4, const Vector& v4 )	
	{
		#if MKL_SUPPORT == 1
			cblas_dcopy(result->size(), &v[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha1, &v1[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha2, &v2[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha3, &v3[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha4, &v4[0], 1, &(*result)[0], 1);
		#else
			double* pt_result = &(*result)[0];
			const double* pt_v = &v[0];
			const double* pt_v1 = &v1[0];
			const double* pt_v2 = &v2[0];
			const double* pt_v3 = &v3[0];
			const double* pt_v4 = &v4[0];
			for(unsigned int i=0;i<result->size();++i)
				*pt_result++ = *pt_v++ + alpha1*(*pt_v1++) + alpha2*(*pt_v2++) + alpha3*(*pt_v3++) + alpha4*(*pt_v4++);
		#endif
	}

	template< typename Vector>
	void sum_plus_scalar_multiplication(Vector* result, const Vector& v, const double alpha1, const Vector& v1, 
		                                                                 const double alpha2, const Vector& v2,
		                                                                 const double alpha3, const Vector& v3,
		                                                                 const double alpha4, const Vector& v4,
		                                                                 const double alpha5, const Vector& v5 )	
	{
		#if MKL_SUPPORT == 1
			cblas_dcopy(result->size(), &v[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha1, &v1[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha2, &v2[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha3, &v3[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha4, &v4[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha5, &v5[0], 1, &(*result)[0], 1);
		#else
			double* pt_result = &(*result)[0];
			const double* pt_v = &v[0];
			const double* pt_v1 = &v1[0];
			const double* pt_v2 = &v2[0];
			const double* pt_v3 = &v3[0];
			const double* pt_v4 = &v4[0];
			const double* pt_v5 = &v5[0];
			for(unsigned int i=0;i<result->size();++i)
				*pt_result++ = *pt_v++ + alpha1*(*pt_v1++) + alpha2*(*pt_v2++) + alpha3*(*pt_v3++) + alpha4*(*pt_v4++) + alpha5*(*pt_v5++);
		#endif
	}

	template< typename Vector>
	void sum_plus_scalar_multiplication(Vector* result, const Vector& v, const double alpha1, const Vector& v1, 
		                                                                 const double alpha2, const Vector& v2,
		                                                                 const double alpha3, const Vector& v3,
		                                                                 const double alpha4, const Vector& v4,
		                                                                 const double alpha5, const Vector& v5,
		                                                                 const double alpha6, const Vector& v6 )	
	{
		#if MKL_SUPPORT == 1
			cblas_dcopy(result->size(), &v[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha1, &v1[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha2, &v2[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha3, &v3[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha4, &v4[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha5, &v5[0], 1, &(*result)[0], 1);
			cblas_daxpy(result->size(), alpha6, &v6[0], 1, &(*result)[0], 1);
		#else
			double* pt_result = &(*result)[0];
			const double* pt_v = &v[0];
			const double* pt_v1 = &v1[0];
			const double* pt_v2 = &v2[0];
			const double* pt_v3 = &v3[0];
			const double* pt_v4 = &v4[0];
			const double* pt_v5 = &v5[0];
			for(unsigned int i=0;i<result->size();++i)
				*pt_result++ = *pt_v++ + alpha1*(*pt_v1++) + alpha2*(*pt_v2++) + alpha3*(*pt_v3++) + alpha4*(*pt_v4++) + alpha5*(*pt_v5++) + alpha6*(*pt_v6++);
		#endif
	}
	
	
	template<typename float_t, typename int_t>
	float_t machine_eps()
	{
		union
		{
			float_t f;
			int_t   i;
		} 
		one, one_plus, little, last_little;
 
		one.f    = 1.0;
		little.f = 1.0;
		last_little.f = little.f;
 
		while(true)
		{
			one_plus.f = one.f;
			one_plus.f += little.f;
 
			if( one.i != one_plus.i )
			{
					last_little.f = little.f;
					little.f /= 2.0;
			}
			else
			{
					return last_little.f;
			}
		}
	}

	static const float  MACHINE_EPSILON_FLOAT = machine_eps<float, uint32_t>();
	static const double MACHINE_EPSILON_DOUBLE = machine_eps<double, uint64_t>();

	double Clock()
	{
		return (double)(clock())/CLOCKS_PER_SEC;
	}

    double GetCpuTime()
    {
		#if MKL_SUPPORT == 1
			return dsecnd();
		#else
            return Clock();
		#endif
    }
}

#endif	//ODESMOKE_UTILITIES_H