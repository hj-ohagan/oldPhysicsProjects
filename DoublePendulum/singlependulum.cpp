//
//  singlependulum.cpp
//  CompPhysSinglePendulum
//
//  Created by Henry O'Hagan on 04/11/2013.
//  Copyright (c) 2013 Henry O'Hagan. All rights reserved.
//

#include "singleheader.h" //includes function declarations and classes used
using namespace std;

/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------SINGLE_PENDULUM---------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/

int main()
{
    
    ofstream outfile1("damped_eulersingle.dat");
    ofstream outfile2("damped_leapfrogsingle.dat");
    ofstream outfile3("damped_rungekuttasingle.dat");
    ofstream outfile4("undamped_eulersingle.dat");
    ofstream outfile5("undamped_leapfrogsingle.dat");
    ofstream outfile6("undamped_rungekuttasingle.dat");
    ofstream outfile7("critical_eulersingle.dat");
    ofstream outfile8("critical_rungekuttasingle.dat");
    complex<double> check(1.0,2.0);
    
    double h;
    double gamma;
    double thetainitial;
    vect theta(2); //vector theta contains theta in the top, theta-dot on the bottom
    
    cout<<"Step size:"<<endl; //choose step size, initial theta, and gamma from within program
    cin>>h;
    cout<<endl<<"Theta initial:"<<endl;
    cin>>thetainitial;
    cout<<endl<<"Gamma:"<<endl;
    cin>>gamma;
    theta.seti(0,thetainitial);
    theta.seti(1,0.0);

    //damped solutions:
    single_step_updater(h,gamma,theta,outfile1, 'e');
    leapfrog(h, gamma, theta, outfile2);
    single_step_updater(h,gamma,theta,outfile3, 'r');
    //undamped solutions:
    single_step_updater(h, 0.0, theta, outfile4, 'e');
    leapfrog(h, 0.0, theta, outfile5);
    single_step_updater(h, 0.0, theta, outfile6, 'r');
    //critical step sizes
    criticalstep(outfile7, 'e');
    criticalstep(outfile8, 'r');
    return 0;
}//END OF MAIN PROGRAM

/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------INTEGRATION-------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/


/*-------------------------------------FUNCTION_FOR_ANY_SINGLE_STEP_METHOD-------------------------------------------*/

void single_step_updater(double h, double gamma, vect theta, ofstream &out1, char u)//perform iterations with update matrix
{
    matrix updater(2,2);
    switch(u)
    {
        case 'e':
            updater=euler_updater_initialise(h,gamma);
        break;
        case 'r':
            updater=rk_updater_initialise(h,gamma);
        break;
        default:
            cout<<"choose 'r' or 'e' for RK4 or euler method"<<endl;
            exit(99);
    }
    int imax=int(100/h);
    double t, thzero=theta.geti(0);
    vect thetaa(2);
    for(int i=0; i<=imax ; i++)
    {
        t=h*i;
        thetaa=analytic(t, gamma, thzero);
        out1<<t<<'\t'<<theta.geti(0)<<'\t'<<theta.geti(1)<<'\t'<<hdifference(thetaa, theta )<<'\t'<<energy(thetaa)<<'\t'<<energy(theta)<<endl;
        theta=updater*theta;
    }
    
}

/*------------------------------------------------EULER_METHOD-------------------------------------------------------*/

matrix euler_updater_initialise(double h, double gamma)//initialise update matrix for euler method
{
    matrix updater(2,2);
    updater.setij(0,0,1.0); updater.setij(0,1,h);
    updater.setij(1,0,-h);  updater.setij(1,1,1.0-h*gamma);
    return updater;
}


/*-----------------------------------------------LEAPFROG_METHOD-----------------------------------------------------*/

//initialise update matrix for leapfrog method
matrix leapfrog_updater_initialise(double h, double gamma)
{
    matrix updater(2,2);
    updater.setij(0,0,0.0);
    updater.setij(0,1,2*h);
    updater.setij(1,0,-2*h);
    updater.setij(1,1,-2*h*gamma);
    return updater;
}

//double step updater
void leapfrog(double h, double gamma, vect theta, ofstream &out1)
{
    double thzero=theta.geti(0), t;
    matrix updater(2,2);
    vect theta2(2), theta3(2), thetaa(2); //extra variables to keep track of theta_n-1
    theta2=theta;
    thetaa=analytic(0.0, gamma, thzero);
    //one step using euler method
    updater=euler_updater_initialise(h,gamma);
    out1<<0.0<<'\t'<<theta.geti(0)<<'\t'<<theta.geti(1)<<'\t'<<hdifference(thetaa, theta )<<'\t'<<fabs(thetaa.geti(0)-theta.geti(0))<<'\t'<<fabs(thetaa.geti(1)-theta.geti(1))<<endl;
    theta=updater*theta;
    
    //leapfrog:
    updater=leapfrog_updater_initialise(h, gamma);
    int imax=(5/h);
    for( int i=0; i<imax; i++)
    {
        t=h*(i+1);
        theta3=theta;
        thetaa=analytic(t, gamma, thzero);
        out1<<t<<'\t'<<theta.geti(0)<<'\t'<<theta.geti(1)<<'\t'<<hdifference(thetaa, theta )<<'\t'<<fabs(thetaa.geti(0)-theta.geti(0))<<'\t'<<fabs(thetaa.geti(1)-theta.geti(1))<<endl;
        theta=updater*theta3+theta2;
        theta2=theta3;
    }
    
}


/*-------------------------------------------------RK4_METHOD--------------------------------------------------------*/

matrix rk_updater_initialise(double h, double gamma) //update matrix for RK4 method
{
    matrix updater(2,2), ddt(2,2), identity(2,2); //a is the matrix for the derivative, i is identity
    ddt.setij(0,0,0.0);  ddt.setij(0,1,1.0);     identity.setij(0,0,1.0); identity.setij(0,1,0.0);
    ddt.setij(1,0,-1.0); ddt.setij(1,1,-gamma);   identity.setij(1,0,0.0); identity.setij(1,1,1.0);
    matrix saver0(2,2), saver1(2,2), saver2(2,2), saver3(2,2); //used for saving progress so less operations are computed
    saver0=ddt*h;
    saver1=(saver0^2);
    saver2=saver1*saver0;
    saver3=saver2*saver0;
    updater=identity+saver0+saver1*0.5+saver2*(1.0/6.0)+saver3*(1.0/24.0);
    return updater;
}

/*-------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------STABILITY--------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------ENERGY_CHECKS------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/

//ENERGY DIFFERENCE OF 2 STATES
double hdifference(vect thetanumerical, vect thetaanalytic)
{
    double hanalytic = energy(thetanumerical);
    double hnumerical= energy(thetaanalytic);
    return fabs((hanalytic-hnumerical)/hanalytic);
}
//ENERGY OF A STATE
double energy(vect theta)
{
    return ( theta.geti(1)*theta.geti(1)/2.0-cos(theta.geti(0)) );
}

/*----------------------------------ANALYTIC_STABILITY_FOR_SINGLE_STEP_METHODS---------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/


//eigenvalues of derivative matrix
vector<complex<double>> deigen(double h, double gamma)
{
    vector<complex<double>> eigen(2);
    if (gamma*gamma-4.0<0)
    {
        complex<double> e1(-h*gamma/2.0 , (h/2.0)*sqrt(4.0-gamma*gamma));
        complex<double> e2(-h*gamma/2.0 , -(h/2)*sqrt(4.0-gamma*gamma));
        eigen[0]=e1; eigen[1]=e2;
    }
    else
    {
        eigen[0]=(h/2.0)*(-gamma + sqrt(gamma*gamma -4.0) );
        eigen[1]=(h/2.0)*(-gamma - sqrt(gamma*gamma -4.0) );
    }
    return eigen;
}

/*------------------------------------------------EULER_METHOD-------------------------------------------------------*/
//eigenvalues of update matrix
vector<complex<double>> eeigen(double h, double gamma)
{
    vector<complex<double>> eigen(2);
    vector<complex<double>> d=deigen(h, gamma);
    eigen[0]=1.0 + d[0];
    eigen[1]=1.0 + d[1];
    return eigen;
}

/*-------------------------------------------------RK4_METHOD--------------------------------------------------------*/
//eigenvalues of update matrix
vector<complex<double>> rkeigen(double h, double gamma)
{
    vector<complex<double>> d=deigen(h, gamma);
    vector<complex<double>> eigen(2); eigen[0]=0.0; eigen[1]=0.0;
    vector<complex<double>> e1(2), e2(2), e3(2) ;
    for (int i=0 ; i<2 ; i++) {e1[i] = d[i]*d[i]; e2[i] = e1[i]*d[i] ; e3[i]=e2[i]*d[i] ; //intermediate values for slightly extra efficiency
        eigen[i]=1.0+d[i]+ e1[i]/2.0 + e2[i]/6.0 + e3[i]/24.0;}
    return eigen;
}


/*---------------------------------------------CRITICAL_STEP_FINDER--------------------------------------------------*/

void criticalstep(ofstream& out1, char u)
{
    double hh=0.00001;
    double h=hh;
    double gamma=0;
    double gg=0.005;
    vector<complex<double>> eigenvalues;
    switch(u)
    {
        case 'e':
            do
            {
            do
            {
            eigenvalues=eeigen(h, gamma);
            h+=hh;
            }while(abs(eigenvalues[0])<1.0 && abs(eigenvalues[1])<1.0);
            out1<<h<<'\t'<<gamma<<endl;;
            h=hh;
            gamma+=gg;
                cout<<gamma<<endl;
            }while(gamma<5.0);
            break;
        case 'r':
            do
            {
                do
                {
                    eigenvalues=rkeigen(h, gamma);
                    h+=hh;
                }while(abs(eigenvalues[0])<1.0 && abs(eigenvalues[1])<1.0);
                out1<<h<<'\t'<<gamma<<endl;
                h=hh;
                gamma+=gg;
                cout<<gamma<<endl;
            }while(gamma<5.0);
            break;
        default:
            cout<<"choose 'r' or 'e' for RK4 or euler method"<<endl;
            exit(99);
    }
}


/*----------------------------------------------ANALYTIC_SOLUTION----------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/

//for calculating energy differences and also solution differences
vect analytic(double t, double gamma, double thetazero)
{
    vect theta(2);
    if (gamma*gamma-4<0)//complex solutions of the characteristic equation
    {
        complex<double> l(-gamma/2.0 , sqrt(4.0-gamma*gamma)/2.0);
        theta.seti(0,real((thetazero/(1.0-l/conj(l)))*(exp(l*t)-(l/conj(l))*exp(conj(l)*t))));
        theta.seti(1,real((thetazero/(1.0-l/conj(l)))*l*(exp(l*t)-exp(conj(l)*t))));
    }
    else if (gamma*gamma-4==0)//repeated solutions of the characteristic equation
    {
        double l=-gamma/2.0;
        theta.seti(0,thetazero*exp(l*t)*(1.0-l*t));
        theta.seti(1,-l*l*thetazero*exp(l*t));
    }
    else//real solutions of the characteristic equation
    {
        double lplus=(-gamma/2.0 + sqrt(gamma*gamma-4.0)/2.0), lminus=(-gamma/2.0 + sqrt(gamma*gamma-4.0)/2.0);
        theta.seti(0,(thetazero/(1.0-lplus/lminus))*(exp(lplus*t)-(lplus/lminus)*exp(lminus*t)));
        theta.seti(1,(thetazero/(1.0-lplus/lminus))*lplus*(exp(lplus*t)-exp(lminus*t)));
    }
    return theta;
}







