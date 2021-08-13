//
//  ising.cpp
//  ising
//
//  Created by Henry O'Hagan on 01/12/2013.
//  Copyright (c) 2013 Henry O'Hagan. All rights reserved.
//
/*-------------------------------------------------------NOTES-------------------------------------------------------*/
/*
Use of gsl will require custom header/library search paths, use the command gsl-config --cflags --libs-without-cblas in
 the shell in order to find the necessary information.
GSL REFERENCE MANUAL:
http://www.gnu.org/software/gsl/manual/gsl-ref.html
*/

#include "ising.h"
using namespace std;

/*-------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------2D_ISING_MODEL---------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
int main(void)
{
o    ofstream outfile10("fixbf2.dat");
    ofstream outfile11("fixbf4.dat");
    ofstream outfile12("fixbf6.dat");
    ofstream outfile13("fixbf8.dat");
    ofstream outfile14("fixbfm2.dat");

    ofstream outfile2000("fixedbeta00for.dat");
    ofstream outfile2010("fixedbeta00back.dat");
    ofstream outfile200("fixedbeta0for.dat");
    ofstream outfile201("fixedbeta0back.dat");
    ofstream outfile210("fixedbeta1for.dat");
    ofstream outfile211("fixedbeta1back.dat");
    ofstream outfile220("fixedbeta2for.dat");
    ofstream outfile221("fixedbeta2back.dat");
    
    ofstream outfile310("spin_heat_maps_with_bfield2.dat");
    ofstream outfile311("spin_heat_maps_with_bfield4.dat");
    ofstream outfile312("spin_heat_maps_with_bfield6.dat");
    ofstream outfile313("spin_heat_maps_with_bfield8.dat");
    ofstream outfile314("spin_heat_maps_with_bfieldm2.dat");
    ofstream outfile40("spin_heat_maps_with_bfield0.dat");

    ofstream  outfile00("test0.dat");
    
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T); //select rng, the default is "Mersenne twister" (mt19937)

    cout<<"RNG :"<<'\t'<<gsl_rng_name(r)<<endl; //this is just the RNG algorithm used, FYI for the user
    
    simulate_constantbfield(outfile0 , 0.0, outfile40 );
    simulate_constantbfield(outfile10, 0.2, outfile310);
    simulate_constantbfield(outfile11, 0.4, outfile311);
    simulate_constantbfield(outfile12, 0.6, outfile312);
    simulate_constantbfield(outfile13, 0.8, outfile313);
    simulate_constantbfield(outfile14,-0.2, outfile314);

    simulate_constanttemp(outfile2010, outfile2000, 0.0);
    simulate_constanttemp(outfile201,outfile200, 0.35);
    simulate_constanttemp(outfile211,outfile210, 0.7);
    simulate_constanttemp(outfile221,outfile220, 1.0);
    
    test(outfile00);
    
    gsl_rng_free (r); //frees all the memory associated with the generator r.
    return 0;
}//END OF MAIN PROGRAM
/*-------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------SIMULATION-----------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/

void simulate_constantbfield(ofstream& fout,double fixed_variable, ofstream& spindomains)
{
    unsigned long reseed;
    time_t timer;
    time_t now;
    int reseed_number=100;
    vector<double> inres;
    inres.assign(100, 0.0);
    vector<vector<double>> results;
    results.assign(4, inres);
    for(int j=1;j<=reseed_number;j++)//RESEED LOOP
    {
        macrostate s(1000,25,0.0,false,fixed_variable,false,true);
        time(&timer);
        for(int i=0;i<100;i++)//VARY BETA LOOP
        {
            if(j==1){s.get_active_state().print_microstate(spindomains);}
            addresults(i, results, s, true);
            printresults(reseed_number, i,j, results, s, fout, false);
            s.incbeta(0.01);
            cout<<i<<endl;
        }
        reseed=(unsigned long) time(NULL);
        time(&now);
        printtime(now);
        cout<<double(100*j)/double(reseed_number)<<"% done , "<<difftime(now, timer)<<" s/seed"<<endl;
        cout<<"Fixed Bfield: reseeding with "<<reseed<<endl;
        gsl_rng_set(r,reseed);
    }
}

void simulate_constanttemp(ofstream& fout1,ofstream& fout2,double fixed_variable)
{
    unsigned long reseed;
    time_t timer;
    time_t now;
    int reseed_number=100;
    vector<double> inres;
    inres.assign(400, 0.0);
    vector<vector<double>> results;
    results.assign(4, inres);
    for(int j=1;j<=reseed_number;j++)//RESEED LOOP
    {
        macrostate sb(1000,25,fixed_variable,false,1.0,false,true);
        time(&timer);
        for(int i=0;i<200;i++)//DECREASE BFIELD LOOP
        {
            addresults(i, results, sb,false);
            printresults(reseed_number, i,j, results, sb, fout1, true);
            sb.incbfield(-0.01);
        }
        for(int i=200;i<400;i++)//INCREASE BFIELD LOOP
        {
            addresults(i, results, sb,false);
            printresults(reseed_number, i,j, results, sb, fout2, true);
            sb.incbfield(0.01);
        }
        reseed=(unsigned long) time(NULL);
        time(&now);
        printtime(now);
        cout<<double(100*j)/double(reseed_number)<<"% done, "<<difftime(now, timer)<<" s/seed"<<endl;
        cout<<"Fixed Temperature: reseeding with "<<reseed<<'\t'<<endl;
        gsl_rng_set(r,reseed);
    }
}

/*------------------------------------------------------MISC---------------------------------------------------------*/
void printtime(time_t now)
{
    tm *ltm=localtime(&now);
    cout<<'['<< 1 + ltm->tm_hour << ":"<< 1 + ltm->tm_min <<":"<< 1 + ltm->tm_sec<<']'<<'\t';
}//prints time given by "time_t now"

void addresults(int i, vector<vector<double>> &results , macrostate s, bool absmagnetisation)
{
    results[0][i]+=s.get_active_state().get_energy();
    if(absmagnetisation)
    {
        results[1][i]+=fabs(s.get_active_state().get_magnetisation());
    }else
    {
        results[1][i]+=s.get_active_state().get_magnetisation();
    }
    results[2][i]+=s.get_magsusc();
    results[3][i]+=s.get_sphc();
}
//add results for different seeds

void printresults(int reseed_number,int i,int j,vector<vector<double>> results , macrostate s ,ofstream& fout, bool bfield)
{
    if(j==reseed_number)
    {
        if(bfield)
        {
            fout<<s.get_active_state().get_bfield()<<'\t';
        }else
        {
        fout<<s.get_beta()<<'\t';
        }
        for(int k=0 ; k<4 ; k++)//WRITE RESULTS LOOP
        {
            fout<<results[k][i]/double(reseed_number)<<'\t';
        }
        fout<<endl;
    }
}

void test(ofstream& fout)
{
    //non-equilibrium states
    macrostate saw(1000, 25, 0.4,false, 0.0, true, false);
    macrostate sae(1000, 25, 0.4,false, 0.0, true, false);

    macrostate smw(1000, 25, 0.4, true, 0.0, true, false);
    macrostate sme(1000, 25, 0.4, true, 0.0, true, false);


    for(int i=0 ; i<5000 ; i++)
    {
        fout<<saw.get_mean_energy()<<'\t'<<smw.get_mean_energy()<<'\t'<<sae.get_mean_energy()<<'\t'<<sme.get_mean_energy()<<endl;
        sae.evolve();
        sme.evolve();
        saw.wolff();
        smw.wolff();
    }
}
/*
Test for whether the wolff algorithm converges to the same energies as the metropolis algorithm

*/