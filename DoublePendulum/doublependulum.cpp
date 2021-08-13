//
//  doublependulum.cpp
//  doublependulum
//
//  Created by Henry O'Hagan on 06/11/2013.
//  Copyright (c) 2013 Henry O'Hagan. All rights reserved.
//

#include "doubleheader.h"
using namespace std;
/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------DOUBLE_PENDULUM---------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
int main()
{
    ofstream outfile1("double.dat");
    ofstream outfile2("doublewor.dat");
    double  g;
    vect theta(4);
    theta.seti(0,0.1); theta.seti(1, 0.0); theta.seti(2, 0.0); theta.seti(3,0.0);// theta initial, phi initial, omega, nu

    
    cout<<endl<<"G:"<<endl;
    cin>>g;
    
    rungekutta(g,theta,outfile1,outfile2);

    return 0;
}//END OF MAIN PROGRAM

/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------INTEGRATION-------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------FUNCTION_FOR_ANY_SINGLE_STEP_METHOD-------------------------------------------*/


void rungekutta(double g, vect theta, ofstream &out1, ofstream &out2)//perform iterations with update matrix
{
    double r=0.01, h=0.005;
    matrix updaterwithg(4,4), updaterzerog(4,4);

    int imax=int(100/h);
    for( ; r<=100 ; r*=pow(10.0,0.01) )
    {
        vect thetazerog=theta;
        vect thetawithg=theta;
        updaterwithg=rk_updater_initialise(h, g, r);
        updaterzerog=rk_updater_initialise(h, 0.0, r);
    for(int i=0; i<=imax ; i++)
    {
        if(r<1.01 && r>0.99 )
        {
        out2<<h*i<<'\t'<<thetawithg.geti(0)<<'\t'<<thetawithg.geti(1)<<'\t'<<thetawithg.geti(2)<<'\t'<<thetawithg.geti(3)<<'\t'<<
                   '\t'<<thetazerog.geti(0)<<'\t'<<thetazerog.geti(1)<<'\t'<<thetazerog.geti(2)<<'\t'<<thetazerog.geti(3)<<endl;
        }
        out1<<h*i<<'\t'<<thetawithg.geti(0)<<'\t'<<thetawithg.geti(1)<<'\t'<<thetawithg.geti(2)<<'\t'<<thetawithg.geti(3)<<'\t'<<
                   '\t'<<thetazerog.geti(0)<<'\t'<<thetazerog.geti(1)<<'\t'<<thetazerog.geti(2)<<'\t'<<thetazerog.geti(3)<<'\t'<<log10(r)<<endl<<endl;
        thetawithg=updaterwithg*thetawithg;
        thetazerog=updaterzerog*thetazerog;
    }
        cout<<r<<endl; //used to monitor progress
    }
}

/*----------------------------------------------RK4_UPDATE_MATRIX----------------------------------------------------*/

matrix rk_updater_initialise(double h, double g, double r) //update matrix for RK4 method
{
    matrix updater(4,4), ddt(4,4), identity(4,4); //a is the matrix for the derivative, i is identity
    ddt.setij(0,0,0.0);      ddt.setij(0,1,0.0);    ddt.setij(0,2,1.0);             ddt.setij(0,3,0.0);
    ddt.setij(1,0,0.0);      ddt.setij(1,1,0.0);    ddt.setij(1,2,0.0);             ddt.setij(1,3,1.0);
    ddt.setij(2,0,-(r+1.0)); ddt.setij(2,1,r);      ddt.setij(2,2,-g);              ddt.setij(2,3,0.0);
    ddt.setij(3,0,r+1.0);    ddt.setij(3,1,-(r+1)); ddt.setij(3,2,g*(1.0-(1.0/r))); ddt.setij(3,3,-(g/r));
    for(int j=0; j<4 ; j++){
    for(int i=0; i<4 ; i++){
        if (i==j) { identity.setij(i,j,1.0); }
        else      { identity.setij(i,j,0.0); }
    }} //initialise identity matrix
    matrix saver0(4,4), saver1(4,4), saver2(4,4), saver3(4,4); //used for saving progress so less operations are computed
    saver0=ddt*h;
    saver1=(saver0^2);
    saver2=saver1*saver0;
    saver3=saver2*saver0;
    updater=identity+saver0+saver1*0.5+saver2*(1.0/6.0)+saver3*(1.0/24.0);
    return updater;
}




