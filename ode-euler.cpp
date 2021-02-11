/*                  _____ _     _     _____ ____    ____  _____  _____ ____  ____  _____ ____ 
 *                 /  __// \ /\/ \   /  __//  __\  / ___\/__ __\/  __//  __\/  __\/  __//  __\
 *                 |  \  | | ||| |   |  \  |  \/|  |    \  / \  |  \  |  \/||  \/||  \  |  \/|
 *       oo   k    |  /_ | \_/|| |_/\|  /_ |    /  \___ |  | |  |  /_ |  __/|  __/|  /_ |    /
 *   x   __  x     \____\\____/\____/\____\\_/\_\  \____/  \_/  \____\\_/   \_/   \____\\_/\_\
 *  e  = \  ---     _  _      ____  ____  _  ____  _____ ____    ____ ___  _   ____  ____  _____ _  _      _____ 
 *       /_  k!    / \/ \  /|/ ___\/  __\/ \/  __\/  __//  _ \  /  __\\  \//  /  _ \/  _ \/  __// \/ \  /|/__ __\
 *       k=1       | || |\ |||    \|  \/|| ||  \/||  \  | | \|  | | // \  /   | / \|| | \||  \  | || |\ ||  / \  
 *                 | || | \||\___ ||  __/| ||    /|  /_ | |_/|  | |_\\ / /    | \_/|| |_/||  /_ | || | \||  | |  
 * Alexander       \_/\_/  \|\____/\_/   \_/\_/\_\\____\\____/  \____//_/     \____/\____/\____\\_/\_/  \|  \_/  
 * Toepfer 2020
 * This work was inspired by headmyshoulder's odeint v2.
 * WORK IN PROGRESS
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#define tab '\t'
#define newl '\n'

template< class InOutIter, class InIter, class T >
void increment( InOutIter first1, InOutIter last1, InIter first2, T alpha ) {
    while( first1 != last1 ) (*first1++) += alpha * (*first2++);
} /* computes y += alpha * x1 with alpha = dt */

template< class Container > struct container_traits;
template< class Container, class Time = double, class Traits = container_traits< Container > >
class euler_stepper {
	public:
		typedef Time time_type;
		typedef Traits traits_type;
		
		typedef typename traits_type::container_type container_type;
		typedef typename traits_type::value_type value_type;
    
		/* performs one step with knowledge of dxdt */
		template< class DynamicSystem >
		void doStep( DynamicSystem& system, container_type& x,
			     const container_type& dxdt, time_type t, time_type dt ) {
			increment( traits_type::begin( x ), traits_type::end( x ), traits_type::begin( dxdt ), dt );
		}
		template< class DynamicSystem >
		void doStep( DynamicSystem& system, container_type& x, time_type t, time_type dt ) {
			/* get dxdt from system to perform step */
			system( x, this->m_dxdt, t );
			doStep( system, x, this->m_dxdt, t, dt );
		}
		euler_stepper() = default;
		euler_stepper( const container_type& x ) {
			this->adjustSize( x );
		}
		void adjustSize( const container_type& x ) {
			/* adjust dxdt's size to same as x */
			traits_type::adjustSize( x, m_dxdt );
		}
		
	private:
		container_type m_dxdt;
};

template< class Container >
struct container_traits {
    typedef Container container_type;
	
    typedef typename container_type::value_type value_type;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

    static void resize( const container_type& x, container_type& dxdt ) {
        dxdt.resize( x.size() );
    }
    static bool sameSize( const container_type& x1, const container_type& x2 ) {
        return( x1.size() == x2.size() );
    }
    static void adjustSize( const container_type& x1, container_type& x2 ) {
        if( !sameSize( x1, x2 ) ) resize( x1, x2 );
    }
    static iterator begin( container_type& x ) {
        return x.begin();
    }
	/* const is needed here for dxdt */
    static const_iterator begin( const container_type& x ) {
        return x.begin();
    }
    static iterator end( container_type& x ) {
        return x.end();
    }
    static const_iterator end( const container_type& x ) {
        return x.begin();
    }
};

/* directly pass stepper onto integrate function */
template< class T >
std::remove_reference<T>::type&& (*consumeStepper)(T&& t) = std::move;

/* iterate the state of the ordinary differential equation */
template< class Stepper, class DynamicSystem, class Observer >
void integrate( Stepper&& stepper, DynamicSystem& system, 
                typename Stepper::container_type& x0,
                typename Stepper::time_type t0, 
                typename Stepper::time_type t1,
                typename Stepper::time_type dt, 
                Observer& observer ) {
    stepper.adjustSize( x0 );
    for( size_t i = 0; t0 < t1; ++i, t0 += dt ) {
        stepper.doStep( system, x0, t0, dt );
        observer( x0, t0 );
    }
}

/* factors for lorenz attractor
 * https://en.wikipedia.org/wiki/Lorenz_system */
const double A = 10.0;
const double B = 28.0;
const double C = 8.0 / 3.0;

/* n-equations to define a system of ordinary differential equations */
typedef std::vector< double > state_type;
typedef euler_stepper< state_type > stepper_type;

/* lorenz system see wiki link */
void lorenzSystem( state_type& x, state_type& dxdt, double t ) {
    dxdt[0] = A * ( x[1] - x[0] );
    dxdt[1] = ( x[0] * ( B - x[2] ) ) - x[1];
    dxdt[2] = x[0] * x[1] - C * x[2];
}

/* observer function for printing */
void observeState( state_type& x, double t ) {
    std::cout << tab << std::setw( 15 ) << std::setprecision( 10 ) << x[0]
    << tab << std::setw( 15 ) << std::setprecision( 10 ) << x[1]
    << tab << std::setw( 15 ) << std::setprecision( 10 ) << x[2]
    << newl;
}

/* example for numerical solutions to the lorenz attractor
 * TODO : implement rk4 stepper additionally */
int main() {
    std::cout << tab << std::right << std::setw( 42 ) << std::setfill( ' ' )
    << "Euler stepper to solve lorenz system" << newl
    << tab << "================================================" << newl;
    state_type x0 = { 1.0, 0.0, 0.0 }; /* starting conditions */
    double dt = 0.01;
    integrate( consumeStepper< stepper_type >( stepper_type() ),
	       lorenzSystem, x0, 0.0, 1, dt, observeState );
}
