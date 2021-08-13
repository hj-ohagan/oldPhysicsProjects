//
//  microstate.h
//  ising
//
//  Created by Henry O'Hagan on 05/12/2013.
//  Copyright (c) 2013 Henry O'Hagan. All rights reserved.
//
/*-------------------------------------------------------NOTES-------------------------------------------------------*/
/*The class is basically for a microstate of an ising model, essentially a square array of spins which can be parallel
or anti-parallel wrt each other. 
 
 The system uses periodic boundary conditions so that the RHS is connected to the LHS and the top is connected to the 
bottom, making the system (topologically) toroidal. This is the reason why the indices used to reference neighbouring
spins will often contain the modulo operation.*/

#ifndef ising_spinstate_h
#define ising_spinstate_h
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <gsl/gsl_rng.h>

gsl_rng * r;
const gsl_rng_type * T;

class microstate
{
private:
    std::vector<std::vector<double>> state;
    double en; //energy total spin
    double ts; //total spin
    unsigned long spin_number; //number of spins
    double bfield; //external magnetic field
    
/*N.B, the reason why these functions are private is so that values are saved in the object, rather than needed to
be recalculated, and to make the constructor more compact when calculating the intial values.*/
    
//CALCULATE_TOTAL_ENERGY_AND_TOTAL_SPIN
    void calc_tspin_energy()
    {
        double e=0.0, sp=0.0;
        unsigned long s=state.size();
        for(unsigned long i=0; i<s; i++)
        {
            for(unsigned long j=0; j<s; j++)
            {
                e+=state[i][j]*(state[i][(j+1)%s]+state[(i+1)%s][j]);
                sp+=state[i][j];
            }
        }
        ts=sp;
        en= -e-bfield*ts;
    }
/*Calculates total energy by adding the interaction with the spin at the bottom and to the right for every spin where
there are such neighbours.(Starts at top left, ends at bottom right). Periodic boundary conditions are used, the
modulo operation is so that the index of vector<vector<>> is always in range.
 Also calculates total spin.*/
    
    std::vector<std::vector<long>> reference_neighbours(std::vector<long> reference_centre)
    {
        long s=state.size();
        std::vector<std::vector<long>> neighbours;
        std::vector<long> ref(2);
        neighbours.assign(4, ref);
        neighbours[0][0]=reference_centre[0];
        neighbours[0][1]=(reference_centre[1]+1)%s;
        neighbours[1][0]=(reference_centre[0]+1)%s;
        neighbours[1][1]=reference_centre[1];
        neighbours[2][0]=reference_centre[0];
        neighbours[2][1]=((reference_centre[1]-1)%s+s)%s;
        neighbours[3][0]=((reference_centre[0]-1)%s+s)%s;
        neighbours[3][1]=reference_centre[1];
        return neighbours;
    }
/*returns a set of references for the neighbours of given spin reference in the following order:
               neighbours[0]
 
neighbours[3] reference_centre  neighbours[1]
              
               neighbours [2]
 */
    
    bool element_of(std::vector<long> a , std::vector<std::vector<long>> b)
    {
        bool c=false;
        for(int i=0; i<b.size() ; i++)
        {
            c=(c||a==b[i]);
        }
        return c;
    }
/*checks if "a" is an element of "b"*/
public:
/*-------------------------------------------------CONSTRUCTORS-----------------------------------------------------*/
    microstate()
    {}
    microstate(unsigned long n, bool alignment, double b, bool matchfield)
    {
        spin_number=n*n;
        if (alignment)//case for aligned spins
        {
            std::vector<double> inter;
            if(matchfield)
            {
                if(b>=0)
                {
                inter.assign(n, 1.0);
                }else
                {
                inter.assign(n,-1.0);
                }
            }else
            {
                if(b>=0)
                {
                inter.assign(n,-1.0);
                }else
                {
                    inter.assign(n, 1.0);
                }
            }
            state.assign(n,inter);
        }else//case for randomised spins
        {state.resize(n);
            for(unsigned long i=0; i<n; i++)
            {   state[i].resize(n);
                for(unsigned long j=0; j<n; j++)
                {   //if statement discretises the uniform distribution
                    if (gsl_rng_uniform(r)<0.5){state[i][j]=-1.0;}
                    else{state[i][j]= 1.0;}
                }
            }
        }
        bfield=b;
        calc_tspin_energy();
    }
/*----------------------------------------------------ACCESSORS------------------------------------------------------*/
//ENERGY_DIFFERENCE_IF_(i,j)_IS_FLIPPED
    double ediff( long i,  long j)
    {
        long s=state.size();
        return 2.0*state[i][j]*
        ((state[(i+1)%s][j]+state[((i-1)%s+s)%s][j]+state[i][(j+1)%s]+state[i][((j-1)%s+s)%s])+bfield);
    }
/*Calculates the difference in energy if spin at position n,m is flipped (in units of j)
 N.B. the "%" function returns the REMAINDER as opposed to modulo, so -A%N = -A for A<N , not N-1 as desired, so in 
the cases where the cell being referenced is possible to return a negative number, we use (-A % N + N) % N, which 
returns N-A for A<N, but still returns A for -A positive.
 N.B. signage of the reference is important when referencing a negative, doing n-1 where n is unsigned causes 
 errors due to underflow.*/
    
//ENERGY
    double get_energy()
    {return en;}
    
//MAGNETISATION
    double get_magnetisation()
    {return ts/double(spin_number);}
        
//SYSTEM_LENGTH_(OR_WIDTH)
    unsigned long get_sys_len()
    {return state.size();}
    
//INDIVIDUAL_SPIN
    double get_spin(unsigned long i, unsigned long j)
    {return state[i][j];}
    
//NUMBER_OF_SPINS
    unsigned long get_spin_number()
    {return spin_number;}
    
//TOTAL_SPIN
    double get_total_spin()
    {return ts;}
    
//EXTERNAL_FIELD
    double get_bfield()
    {return bfield;}
/*----------------------------------------------------MODIFIERS------------------------------------------------------*/
//EVOLVE_THE_SYSTEM_VIA_METROPOLIS
    void evolve(double beta)
    {
        long m,n;
        double ed;
        bool condition;
        do
        {
        m=gsl_rng_uniform_int(r,get_sys_len());
        n=gsl_rng_uniform_int(r,get_sys_len());
        ed=ediff(n,m);
        condition=(ed<=0||gsl_rng_uniform(r)<exp(-ed*beta));
            //it is not necessary to fret over inefficiency, the RHS of the || is only evaluated if ed<=0 returns false
            //this is known as 'Lazy evaluation'
        if (condition)
        {
            state[n][m]= -state[n][m];
            en+=ed;
            ts+=2*state[n][m];
        }
        }while(!condition);
    }
    
//EVOLVE_VIA_WOLFF_ALGORITHM
    void wolff(double beta)
    {
        bool condition=false;
        do
        {
        std::vector<long> active_index(2);
        active_index[0]=gsl_rng_uniform_int(r,get_sys_len());
        active_index[1]=gsl_rng_uniform_int(r,get_sys_len());
        std::vector<std::vector<long>> cluster, boundary_old(0), neighbours;
        cluster.push_back(active_index);
        boundary_old.push_back(active_index);
        do
        {
            std::vector<std::vector<long>> boundary_new(0);
            for (int i=0; i<boundary_old.size(); i++)
            {
                neighbours=reference_neighbours(boundary_old[i]);
                for ( int j=0; j<neighbours.size() ; j++)
                {
                    if( state[ neighbours[j][0] ][ neighbours[j][1] ] == state[ boundary_old[i][0] ][ boundary_old[i][1] ]
                       && !element_of(neighbours[j], cluster)
                       && gsl_rng_uniform(r)<(1-exp(-2.0*beta)) ){
                        //and if the spin of the neighbour is equal to the spin of the centre
                        //if the neighbour j is not an element of the cluster,
                        //and if the cluster addition conditions are met
                        boundary_new.push_back(neighbours[j]);
                        cluster.push_back(neighbours[j]);
                    }
                }
            }
            boundary_old=boundary_new;
        }while(boundary_old.size()!=0);
        double ed_bfield=bfield*2.0*state[active_index[0]][active_index[1]]*double(cluster.size());
        condition=(ed_bfield<=0||gsl_rng_uniform(r)<exp(-beta*ed_bfield));
        if(condition)
        {
            for(int i=0; i<cluster.size() ; i++)
            {
                //en+=ediff(cluster[i][0],cluster[i][1]);
                state[ cluster[i][0] ][ cluster[i][1] ] = - state[ cluster[i][0] ][ cluster[i][1] ];
            }//flip the entire selected cluster
            //ts+=2.0*state[active_index[0]][active_index[1]]*double(cluster.size());
            calc_tspin_energy();
        }
        }while(!condition); //keep iterating if the condition is not met
    }
/*This algorithm was originally detailed in the following paper:
Wolff, U. (1989). Collective Monte Carlo Updating for Spin Systems. Phys. Rev. Lett., 62:361.
This update algorithm is only for temperatures close to Tc. This algorithm has a 100% success rate for no external
 field.*/

//INCREMENT_BFIELD
    void inc_bfield(double x)
    {
        bfield+=x;
        en-=x*ts;
    }
/*------------------------------------------------------MISC---------------------------------------------------------*/
//PRINT_TO_FILE
    void print_microstate(std::ofstream& fout)
    {
        for(int i=0 ; i<state.size() ; i++)
        {
            for(int j=0 ; j<state[i].size() ; j++)
                {fout<<state[i][j]<<'\t';
            }fout<<std::endl;
        }fout<<std::endl;
    }
//ERASE_VECTOR
    void clearspins()
    {
        state.erase(state.begin(),state.end());
    }
/*This method is used to free memory*/
};
#endif
