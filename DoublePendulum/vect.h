//
//  vect.h
//  CompPhysSinglePendulum
//
//  Created by Henry O'Hagan on 05/11/2013.
//  Copyright (c) 2013 Henry O'Hagan. All rights reserved.
//

//vector class with additional operations

#ifndef CompPhysSinglePendulum_vect_h
#define CompPhysSinglePendulum_vect_h
#include <vector>
#include <iostream>

class vect
{
private:
    std::vector<double> vec;
    
public:
    
/*--------------------------------------------------CONSTRUCTORS-----------------------------------------------------*/
    
    vect(unsigned int m)
    {
        vec.resize(m);
    }
    
/*----------------------------------------------------ACCESSORS------------------------------------------------------*/
    double geti(unsigned int i)
    {
        return vec[i];
    }
    
    unsigned int getsize()
    {
        return int(vec.size());
    }
    
/*----------------------------------------------------MODIFIERS------------------------------------------------------*/
    void seti(unsigned int i, double x)
    {
        vec[i]=x;
    }
    
    void vresize(unsigned int i)
    {
        vec.resize(i);
    }
    

/*-----------------------------------------------OVERLOADED_OPERATORS------------------------------------------------*/
	// + OPERATOR
	vect operator+(vect s)
	{
        vect v_sum(int(vec.size()));
        if (s.getsize()==int(vec.size()))
        {
                for( unsigned int i=0 ; i<int(vec.size()) ; i++)
                {
                    v_sum.seti(i,s.geti(i)+vec[i]);
                }
        }
        else
        {
                std::cout<<"Vector std::size do not match."<<std::endl;
                exit(99);
        }
		
		return v_sum;
	}
    // - OPERATOR
	vect operator-(vect s)
	{
        vect v_new(int(vec.size()));
        if (s.getsize()==int(vec.size()))
        {
                for( unsigned int i=0 ; i< int(vec.size()) ; i++)
                {
                    v_new.seti(i,vec[i]-s.geti(i));
                }
        }
        else
        {
        std::cout<<"Vector size do not match."<<std::endl;
                exit(99);
        }
		
		return v_new;
	}
    // SCALAR MULTIPLICATION
    vect operator*(double s)
    {
        vect v_new(int(vec.size()));
        
        for(unsigned int i=0 ; i< int(vec.size()); i++)
        {
        v_new.seti(i,vec[i]*s);
        }
        return v_new;
    }
    // INNER PRODUCT
    double operator*(vect s)
    {
        double d_new=0.0;
        if(s.getsize()==int(vec.size()))
        {
                for (unsigned int i=0 ; i<int(vec.size()) ; i++)
                {
                    d_new=d_new+(vec[i]*s.geti(i));
                }
        }
        else
        {
            std::cout<<"vector size do not match"<<std::endl;
             exit(99);
        }
        return d_new;
    }
    
/*-------------------------------------------------------MISC--------------------------------------------------------*/
    void print() //prints as column vector
    {
        for(unsigned int i=0; i<vec.size() ; i++)
        {
            std::cout<<vec[i]<<std::endl;
        }
        std::cout<<std::endl;
    }
    
};


#endif
