
#pragma once
#include <iostream>
#include "isop.hpp"
#include "cube.hpp"
#include <lpsolve/lp_lib.h>
#include "constructors.hpp"
#include "dynamic_truth_table.hpp"
#include "static_truth_table.hpp"
#include <vector>


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
	Ncol = tt.num_vars() + 1;
	

	// Start lp problem configurations
	lp = make_lp( 0, Ncol );

	if ( lp == NULL )
		return false;

	variables = (int*)malloc( Ncol * sizeof( *variables ) );
	weights = (REAL*)malloc( Ncol * sizeof( *weights ) );

	if ( ( variables == NULL ) || ( weights == NULL ) )
		return false;
		
	std::vector<cube> ONSET = isop( tt );
	std::vector<cube> OFFSET = isop( unary_not( tt ) );
	//we set the variables to integer
	for ( uint8_t i = 1; i <= Ncol; i++ )
	{
	    set_int( lp, i, TRUE ); 
	}
	  
	set_add_rowmode( lp, TRUE );

	for ( cube cube : ONSET )
	{
		j = 0;
		for ( uint8_t = 0; i < tt.num_vars(); i++ )
		{
			if ( cube.get_mask( i ) ) //Check if the variable is in the cube
			{
				variables[j] = i + 1;
				weights[j] = 1;
				j++;
			}
		}
		variables[j] = Ncol; /* this is the right part of the constraint : weight put to -1 */
		weights[j] = -1; 
		j++;
		
		add_constraintex( lp, j, weights, variables, GE, 0 );
	}


	for ( cube cube : OFFSET )
	{
		j = 0;
		for ( uint8_t i = 0; i < tt.num_vars(); i++ )
		{
			if ( !cube.get_mask( i ) ) 
			{
				variables[j] = i + 1;
				weights[j] = 1;
				j++;
			}
		}
		variables[j] = Ncol; /* this is the right part of the constraint : weight put to -1 */
		weights[j] = -1;
		j++;
	
		add_constraintex( lp, j, weights, variables, LE, -1 );
	}
	
	 for (uint8_t i = 0; i < tt.num_vars() + 1; i++) 
	 {
            for (uint8_t j = 0; j < tt.num_vars() + 1; j++) 
            {
                variables[j] = j + 1;
                weights[j] = (i == j) ? 1 : 0;
            }
            
            add_constraintex(lp, tt.num_vars() + 1, weights, variables, GE, 0);
        }

	set_add_rowmode( lp, FALSE ); /* turned off because we are done building the model */

	
	for ( uint8_t i = 0; i < Ncol; i++ )
	{
		variables[i] = i + 1;
		weights[i] = 1;
	}
	set_obj_fnex( lp, Ncol, weights, variables );
	/* set the objective in lpsolve */
	
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
	

	for ( uint8_t i = 0u; i < tt.num_vars(); i++ )
	{
		pos_unate = false; //reinitialization because we check the unateness for each variable
		neg_unate = false;

		TT const tt_neg = cofactor0( tt, i );
		TT const tt_pos = cofactor1( tt, i );
		for ( TT position = 0; position < ( 1 <<  tt.num_vars() ); position++ )
		{
			if ( get_bit( tt_neg, position ) < get_bit( tt_pos, position ) )
			{
				pos_unate = true;
			}

			else if ( get_bit( tt_neg, position ) > get_bit( tt_pos, position ) )
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
			
			for ( int position  = ( 1 << ( i + 1 ) ); position  < ( 1 << tt.num_vars() ); position  = position  +( 1 << ( i + 2 ) ) )
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
