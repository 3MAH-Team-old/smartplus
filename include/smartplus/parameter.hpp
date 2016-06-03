///@file parameter.hpp
///@brief parameters of SMART+
///@author Chemisky & Despringre
///@version 1.0
///@date 01/27/2014

#include <boost/math/constants/constants.hpp>
#define UNUSED(x) [&x]{}()

#ifndef version_full
#define version_full 1
#endif

namespace smart {

#ifndef pi
#define pi boost::math::constants::pi<double>()
#endif
    
#ifndef axis_psi
#define axis_psi 3
#endif

#ifndef axis_theta
#define axis_theta 1
#endif

#ifndef axis_phi
#define axis_phi 3
#endif

#ifndef limit
#define limit 1E-9
#endif

#ifndef iota
#define iota 1E-12
#endif

#ifndef miniter_umat
#define miniter_umat 10
#endif

#ifndef maxiter_umat
#define maxiter_umat 100
#endif

#ifndef precision_umat
#define precision_umat 1E-9
#endif

#ifndef div_tnew_dt_umat
#define div_tnew_dt_umat 0.2
#endif

#ifndef mul_tnew_dt_umat
#define mul_tnew_dt_umat 2
#endif

#ifndef lambda_solver
#define lambda_solver 10000
#endif
    
#ifndef miniter_solver
#define miniter_solver 10
#endif

#ifndef maxiter_solver
#define maxiter_solver 100
#endif

#ifndef precision_solver
#define precision_solver 1E-6
#endif

#ifndef inforce_solver
#define inforce_solver 0
#endif

#ifndef div_tnew_dt_solver
#define div_tnew_dt_solver 0.5
#endif

#ifndef mul_tnew_dt_solver
#define mul_tnew_dt_solver 2
#endif


#ifndef maxiter_micro
#define maxiter_micro 100
#endif

#ifndef precision_micro
#define precision_micro 1E-6
#endif

} //end of namespace smart
