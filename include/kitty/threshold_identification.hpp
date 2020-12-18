
#pragma once
#include <iostream>
#include <vector>
#include <lpsolve/lp_lib.h>
#include "constructors.hpp"
#include "dynamic_truth_table.hpp"
#include "static_truth_table.hpp"
#include "isop.hpp"
#include "cube.hpp"



namespace kitty
{

/*! \brief Threshold logic function identification

*/
template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool is_threshold( const TT& tt, std::vector<int64_t>* plf = nullptr )
{
	TT tt_star = tt;
	std::vector<int64_t> inv;
	
	if ( unate( tt_star, &inv ) )
	{
		lp_solver( tt_star, plf, &inv ) ;
	}
	else 
	{
		 return false;
	}

}

template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool lp_solver( TT& tt,  std::vector<int64_t>* plf, std::vector<int64_t>* inv )
{
	lprec* lp;
	int Ncol, *variables = NULL, j;
	REAL* weights = NULL;
	Ncol = tt.num_vars() + 1; // we have num_vars variables + 1 for the thresold value
	

	/* We will build the model row by row, so we star with 0 rows and Ncol columns */

	lp = make_lp( 0, Ncol );

	if ( lp == NULL )
		return false; /* couldn't construct a new model... */

	/* create space large enough for one row */
	variables = (int*)malloc( Ncol * sizeof( *variables ) );
	weights = (REAL*)malloc( Ncol * sizeof( *weights ) );

	if ( ( variables == NULL ) || ( weights == NULL ) )
		return false;  
	
	/* We compute the OFFSET and ONSET */
	auto ONSET = isop( tt );
	auto OFFSET = isop( unary_not( tt ) );
	//we set the variables to integer
	for ( auto i = 1; i <= Ncol; i++ )
	{
	    set_int( lp, i, TRUE ); 
	}
	  
	set_add_rowmode( lp, TRUE ); /*building the model is faster if done row by row */
	//we begin with th ONSET cube 

	for ( auto& cube : ONSET )
	{
		j = 0;
		for ( unsigned int i = 0; i < tt.num_vars(); i++ )
		{
			if ( cube.get_mask( i ) ) //Check if the variable is in the cube
			{
				variables[j] = i + 1;
				weights[j] = 1;
				j++;
			}
		}
		variables[j] = Ncol; 
		weights[j] = -1; 
		j++;
		
		add_constraintex( lp, j, weights, variables, GE, 0 );
	}


	for ( auto& cube : OFFSET )
	{
		j = 0;
		for ( auto i = 0; i < tt.num_vars(); i++ )
		{
			if ( !cube.get_mask( i ) ) //Check if the variable is NOT in the cube
			{
				variables[j] = i + 1;
				weights[j] = 1;
				j++;
			}
		}
		variables[j] = Ncol; //this is the right part of the constraint : weight put to -1
		weights[j] = -1;
		j++;
	
		add_constraintex( lp, j, weights, variables, LE, -1 );
	}
	/* we add a constraint on the sign of the variables */
	 for (auto i = 0; i < tt.num_vars() + 1; i++) 
	 {
            for (auto j = 0; j < tt.num_vars() + 1; j++) 
            {
                variables[j] = j + 1;
                weights[j] = (i == j) ? 1 : 0;
            }
            
            add_constraintex(lp, tt.num_vars() + 1, weights, variables, GE, 0);
        }

	set_add_rowmode( lp, FALSE ); /* turned off because we are done building the model */

	
	for ( auto i = 0; i < Ncol; i++ )
	{
		variables[i] = i + 1;
		weights[i] = 1;
	}
	set_obj_fnex( lp, Ncol, weights, variables ); /* set the objective in lpsolve */
	
	//we need to minimize the objective function
	set_minim( lp );
	set_verbose(lp,0);
	/* Do the trick */
	if ( solve( lp ) != OPTIMAL ) //This is not a Threshold function if LP doesn't give an optimal solution
		return false;

	
	get_variables( lp, weights );
	int sum = 0;
	
	for ( j = 0; j < Ncol; j++ )
	{
		bool not_empty = ( inv->size() != 0 );
		if ( ( not_empty ) && ( inv->at( 0 ) == j ) )
		{
			sum += weights[j];
			if ( j == Ncol - 1 )
			{
				plf->push_back( -weights[j] - sum );
			}
			else
			{
				plf->push_back( -weights[j] );
			}
			inv->erase( inv->begin() );
		}
		else
		{
			if ( j == Ncol - 1 )
			{
				plf->push_back( weights[j] - sum );
			}
			else
			{
				plf->push_back( weights[j] );
			}
		}
	}
    
	/* free allocated memory */
	if ( weights != NULL )
	{
		free( weights );
	}
	if ( variables != NULL )
	{
		free( variables );
	}

	if ( lp != NULL )
	{
		delete_lp( lp );
	}
	
	return true;
}

template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool unate( TT& tt, std::vector<int64_t>* inv )
{
  
	bool pos_unate = false;
	bool neg_unate = false;
	auto numvars = tt.num_vars();

	for ( auto i = 0u; i < numvars; i++ )
	{
		pos_unate = false; //reinitialization because we check the unateness for each variable
		neg_unate = false;

		auto const tt_neg = cofactor0( tt, i );
		auto const tt_pos = cofactor1( tt, i );
		for ( auto bit = 0; bit < ( 1 <<  numvars ); bit++ )
		{
			if ( get_bit( tt_neg, bit ) < get_bit( tt_pos, bit ) )
			{
				pos_unate = true;
			}

			else if ( get_bit( tt_neg, bit ) > get_bit( tt_pos, bit ) )
			{
				neg_unate = true;
			}
		}


		if ( neg_unate && pos_unate )
		{
			return false;
		}
		
		if ( neg_unate )
		{
			inv->push_back( i ); 			//we save the inv variables for the linear transfo
			
			for ( int position  = ( 1 << ( i + 1 ) ); position  < ( 1 << numvars ); position  = position  +( 1 << ( i + 2 ) ) )
			{
				for ( int j = 0; j < ( 1 << ( i + 1 ) ); j++ )
				{
					
					if ( get_bit( tt, position  + j ) != get_bit( tt,  position  + j - ( 1 << ( i + 1 ) ) ) )
					{
						flip_bit( tt, position  + j );
						flip_bit( tt, position  + j - ( 1 << ( i + 1 ) ) );
					}
				}
			}
		}
	}
	return true;
}
}
