//
//  macrostate.h
//  ising
//
//  Created by Henry O'Hagan on 05/12/2013.
//  Copyright (c) 2013 Henry O'Hagan. All rights reserved.
//
/*-------------------------------------------------------NOTES-------------------------------------------------------*/
/*A class for describing macroscopic properties, using many microstates. In the main program the user will only be 
allowed to access equilibrium states.*/
#ifndef ising_macrostate_h
#define ising_macrostate_h
#include "microstate.h"

class macrostate
{
private:
    std::vector<microstate> state;
    double beta, magsusc, sphc, menergy, mspin, msenergy, msspin;
    //beta, magnetic susceptibility, specific heat capacity, means (which must be stored to be updatable)
    unsigned long active, meansize;
    //used to reference the latest microstate to be updated.
    
/*The following functions are made private to ensure they are only ever called once (during construction).*/    
//COMPUTE_MEANS
    void compute_means()
    {
        double e=0.0, se=0.0, s=0.0, ss=0.0;
        //energy, square energy, total spin, square total spin
        for(unsigned long i=0; i<meansize; i++)
        {
            e+=state[i].get_energy();
            se+=(state[i].get_energy()*state[i].get_energy());
            s+=state[i].get_total_spin();
            ss+=state[i].get_total_spin()*state[i].get_total_spin();
        }
        menergy=e/double(meansize);
        msenergy=se/double(meansize);
        mspin=s/double(meansize);
        msspin=ss/double(meansize);
    }
    
    
    //EVOLVE_TO_STEADY_STATE
    void equilibrate(bool aligned)
    {
        double mean_energy_difference, init_med, bf=state[active].get_bfield();
        macrostate reference;
        long slength=state[active].get_sys_len();
        for( int i=0 ; i<(meansize+slength*slength) ; i++)
        {
            evolve();
        }
        if((state[active].get_magnetisation()>0)==(bf>0))
        {
            macrostate ref(meansize,slength,beta,true,bf,true,true);
            reference=ref;
            /*if the sign of the bfield is equal to the sign of the magnetisation then the reference state needs
             to also match the alignment, and vice versa, because the magnetisation and magnetic field can be of
             opposite sign for sufficiently high beta, otherwise the reference state and current macrostate will
             settle into different equilibrium configurations*/
        }else
        {
            macrostate ref(meansize,slength,beta,true,bf,true,false);
            reference=ref;
        }
        if(!aligned)
            //if the state is aligned, then we need to translate the reference state back, so that the conditions are not
        {
            for( int i=0 ; i<(meansize+slength*slength) ; i++)
            {
                reference.evolve();//use to update the reference state to the same time as the state of the system
            }
        }
        init_med=reference.get_mean_energy()-menergy;//initial mean energy difference
        do
        {
            evolve();
            reference.evolve();
            mean_energy_difference=reference.get_mean_energy()-menergy;
        }while(mean_energy_difference<0.0==init_med<0.0);
        for( int i=0 ; i<(meansize) ; i++)
        {
            evolve();
        }//after it reaches equilibrium update all microstates use in the running mean
    }
    /*Evolves the macrostate until a steady state is reached.
     This creates a macrostate that is not in equilibrium, then allows both systems to evolve until the energy difference
     changes sign (which indicates that the energies have intersected)*/

public:
/*-------------------------------------------------CONSTRUCTORS-----------------------------------------------------*/
    macrostate()
    {}
    macrostate(unsigned long ms,unsigned long size, double b, bool alignment_initial, double bf, bool equilibrating, bool matchfield)
    {
        meansize=ms;
        active=meansize-1;
        beta=b;
        state.resize(meansize);
        microstate intermediate(size,alignment_initial,bf,matchfield);
        for(unsigned long i=0; i<meansize; i++ )
        {
            state[i]=intermediate;
            intermediate.evolve(b);
            if(i!=active)
            {state[i].clearspins();}
            //Since the actual array of spins aren't used in any calculations, this information can be forgotten
            //which reduces memory usage by a factor of ~30.
        }
        compute_means();
        magsusc=beta*(msspin-mspin*mspin)/double(state[0].get_spin_number());
        sphc=beta*beta*(msenergy-menergy*menergy)/double(state[0].get_spin_number());
        if(!equilibrating)
        {
            equilibrate(alignment_initial);
        }
    }
/*----------------------------------------------------ACCESSORS------------------------------------------------------*/
//GET_NUMBER_OF_MICROSTATES
    unsigned long get_microstate_number()
    {return state.size();}
    
//RETURN_BETA
    double get_beta()
    {return beta;}
    
//SPECIFIC_HEAT_CAPACITY
    double get_sphc()
    {return sphc;}

//MAGNETIC_SUSCEPTIBILITY
    double get_magsusc()
    {return magsusc;}
    
//RETURN_PARTICULAR_MICROSTATE
    microstate get_microstate(unsigned long i)
    {return state[i];}
    
//RETURN_ACTIVE_MICROSTATE
    microstate get_active_state()
    {return state[active];}
    
//RETURN_ACTIVE_INDEX/REFERENCE
    unsigned long get_active()
    {return active;}
    
//RETURN_MEAN_ENERGY
    double get_mean_energy()
    {return menergy;}
/*----------------------------------------------------MODIFIERS------------------------------------------------------*/
//INCREMENT_BETA
    void incbeta(double x)
    {
        beta+=x;
        equilibrate(false);
    }
//INCREMENT_BFIELD
    void incbfield(double x)
    {
        menergy-=state[active].get_energy()/double(meansize);
        msenergy-=state[active].get_energy()*state[active].get_energy()/double(meansize);
        state[active].inc_bfield(x);
        menergy+=state[active].get_energy()/double(meansize);
        msenergy+=state[active].get_energy()*state[active].get_energy()/double(meansize);
        sphc=beta*beta*(msenergy-menergy*menergy)/double(state[0].get_spin_number());
        equilibrate(false);
    }
    

    
    
    //EVOLVE_ACTIVE_MICROSTATE
    void evolve()
    {
        microstate intermediate=state[active];

        intermediate.evolve(beta);
        state[active].clearspins();
        active=(active+1)%meansize;
        
        double en_old=state[active].get_energy(),       en_new=intermediate.get_energy();
        double sp_old=state[active].get_total_spin(),   sp_new=intermediate.get_total_spin();
        double msize=double(meansize);
        
        menergy+=(en_new-en_old)/msize;
        msenergy+=(en_new*en_new-en_old*en_old)/msize;
        mspin+=(sp_new-sp_old)/msize;
        msspin+=(sp_new*sp_new-sp_old*sp_old)/msize;
        
        magsusc=beta*(msspin-mspin*mspin)/double(state[active].get_spin_number());
        sphc=beta*beta*(msenergy-menergy*menergy)/double(state[active].get_spin_number());
        
        state[active]=intermediate;
    }
    /*Evolves the latest microstate by a single step, replaces new microstate with the oldest microstate,
     and updates all quantities of interest.*/
    
    //EVOLVE_ACTIVE_MICROSTATE_via wolff algorithm
    void wolff()
    {
        microstate intermediate=state[active];
        intermediate.wolff(beta);
        
        state[active].clearspins();
        active=(active+1)%meansize;
        
        double en_old=state[active].get_energy(),       en_new=intermediate.get_energy();
        double sp_old=state[active].get_total_spin(),   sp_new=intermediate.get_total_spin();
        double msize=double(meansize);
        
        menergy+=(en_new-en_old)/msize;
        msenergy+=(en_new*en_new-en_old*en_old)/msize;
        mspin+=(sp_new-sp_old)/msize;
        msspin+=(sp_new*sp_new-sp_old*sp_old)/msize;
        
        magsusc=beta*(msspin-mspin*mspin)/double(state[active].get_spin_number());
        sphc=beta*beta*(msenergy-menergy*menergy)/double(state[active].get_spin_number());
        
        state[active]=intermediate;
    }
    /*N.B. This is an extension, has been tested to converge to the correct values but not used for */
    
};
#endif
