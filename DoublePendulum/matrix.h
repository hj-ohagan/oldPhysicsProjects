//
//  matrix.h
//  CompPhysSinglePendulum
//
//  Created by Henry O'Hagan on 04/11/2013.
//  Copyright (c) 2013 Henry O'Hagan. All rights reserved.
//

//class for basic linear algebra
//important operations this class will acheive:
//matrix*matrix, matrix*vector,scalar*matrix, transpose, matrix+matrix, matrix-matrix, eigenvalues


#ifndef CompPhysSinglePendulum_matrix_h
#define CompPhysSinglePendulum_matrix_h
#include "vect.h"


class matrix
{
private:
    std::vector<std::vector<double>> mat;
    
public:
    
/*--------------------------------------------------CONSTRUCTORS-----------------------------------------------------*/
    
    //m by n matrix:
    matrix(unsigned int m,unsigned int n)
    {
        mat.resize(m);
        for(unsigned int i=0 ; i<m ; i++)
        {
            mat[i].resize(n);
        }
    }
    
/*----------------------------------------------------ACCESSORS------------------------------------------------------*/
    
    double getij(unsigned int i,unsigned int j)//get value of row i column j
    {
        return mat[i][j];
    }
    unsigned int getisize()//get number of rows
    {
        return int(mat.size());
    }
    unsigned int getjsize()//get number of columns
    {
        return int(mat[0].size());
    }
    
/*----------------------------------------------------MODIFIERS------------------------------------------------------*/
    void setij(unsigned int i,unsigned int j, double x)//set value of row i column j
    {
        mat[i][j]=x;
    }
    
    void resizei(unsigned int i)//resize rows
    {
        mat.resize(i);
    }
    void resizej(unsigned int j)//resize columns
    {
        for(unsigned int i=0 ; i<mat.size() ; i++ )
        {
            mat[i].resize(j);
        }
    }
    
    //TRANSPOSE
    void transpose()//needed to left-multiply
    {
        unsigned int rows=int(mat.size());
        unsigned int columns=int(mat[0].size());
        matrix m1(int(mat.size()), int(mat[0].size()));
        for(unsigned int i=0; i<rows; i++)
        {
            for(unsigned int j=0 ; i<columns; j++)
            {
                m1.setij(i,j,mat[i][j]);
            }
        }
        

        mat.resize(columns);
        for(unsigned int i=0; i<columns ; i++)
        {
            mat[i].resize(rows);
        }
        
        for(unsigned int i=0 ; i < columns ; i++)
        {
            for (unsigned int j=0 ; j < rows ; j++ )
            {
                mat[i][j]=m1.getij(j,i);
            }
        }
    }
/*-----------------------------------------------OVERLOADED_OPERATORS------------------------------------------------*/
    // MATRIX * MATRIX
    matrix operator*(matrix s)
    {
        //rows of s must = columns of m_before
        if( s.getisize() == int(mat[0].size()) )
        {
            matrix m_new(int(mat.size()),s.getjsize());
            matrix m_before(int(mat.size()),int(mat[0].size()));
            for(unsigned int i=0 ; i<int(mat.size());i++)
            {
                for(unsigned int j=0 ; j<int(mat[0].size()) ; j++)
                {
                    m_before.setij(i,j,mat[i][j]); //using matrix as opposed to vector<vector> to use matrix*vector
                }
            }
            vect v_intermediate1(s.getisize());
            vect v_intermediate2(m_before.getisize());
            for (int j=0; j<s.getjsize() ;j++)
            {
            for (int i=0; i<s.getisize();i++)
            {
                v_intermediate1.seti(i,s.getij(i,j));
            }
                v_intermediate2=(m_before*v_intermediate1);
                for(int i=0; i<v_intermediate2.getsize() ; i++)
                {
                    m_new.setij( i, j, v_intermediate2.geti(i));
                    
                }
            }
            

        return m_new;
        }
        else
        {
            std::cout<<"#columns of m_left != #rows of m_right"<<std::endl;
            exit(99);
        }
    }
    //POWERS OF MATRICES
    matrix operator^(int n)
    {
        matrix m_new(int(mat.size()), int(mat[0].size())), m_1(int(mat.size()), int(mat[0].size()));
        for (unsigned int i=0; i<int(mat.size()); i++)
        {
            for(unsigned int j=0; j<int(mat[0].size());j++)
            {
                m_new.setij(i,j,mat[i][j]);
            }
        }
        m_1=m_new;
        for (unsigned int i=1 ; i<n; i++)
        {
            m_new= (m_new*m_1);
            
        }
        
        return m_new;
    }
    // MATRIX * VECTOR
    vect operator*(vect s)
    {
        vect v_new(int(mat.size()));
        
        for(unsigned int i=0 ; i<int(mat.size()); i++)
        {
            vect v_intermediate(int(mat[i].size()));
        if (mat[i].size()==s.getsize()) //does the number of columns of row i, match the number of entries of vector s?
        {
                for(unsigned int j=0 ; j<int(mat[i].size()) ; j++)
                {
                    v_intermediate.seti(j,mat[i][j]); //make new vector that inherits dot product
                }
                v_new.seti(i,v_intermediate*s); //use inherited properties of vect
        }
        else
        {
                std::cout<<"size of vector != #columns of matrix"<<std::endl;
                exit(99);
        }
        }
        return v_new;
    }
    //SCALAR MULTIPLICATION
    matrix operator*(double s)
    {
        matrix m_new(int(mat.size()),int(mat[0].size()));
        for(unsigned int i=0 ; i< mat.size(); i++)
        {
            for(unsigned int j=0 ; j<int(mat[i].size()) ; j++)
            {
                m_new.setij(i,j,mat[i][j]*s);
            }
        }
        return m_new;
    }
    // + OPERATOR
    matrix operator+(matrix s)
    {
        matrix m_new(int(mat.size()),int(mat[0].size()));
        if(int(mat.size())==s.getisize())
        {
            for (unsigned int i=0; i<int(mat.size()) ; i++)
            {
                if (int(mat[i].size())==s.getjsize())
                {
                    for( unsigned int j=0 ; j<int(mat[i].size()); j++)
                    {
                        m_new.setij(i,j,mat[i][j]+s.getij(i,j));
                    }
                }
                else
                {
                    std::cout<<"#columns of matrix_1 != #columns of matrix_2"<<std::endl;
                    exit(99);
                }
                        
            }
        }
        else
        {
            std::cout<<"#rows of matrix_1 != #rows of matrix_2"<<std::endl;
            exit(99);
        }
        return m_new;
    }
    // - OPERATOR
    matrix operator-(matrix s)
    {
        matrix m_new(int(mat.size()),int(mat[0].size()));
        if(int(mat.size())==s.getisize())
        {
            for (unsigned int i=0; i<int(mat.size()) ; i++)
            {
                if (int(mat[i].size())==s.getjsize())
                {
                    for( unsigned int j=0 ; j<int(mat[i].size()); j++)
                    {
                        m_new.setij(i,j,mat[i][j]-s.getij(i,j));
                    }
                }
                else
                {
                    std::cout<<"#columns of matrix_1 != #columns of matrix_2"<<std::endl;
                    exit(99);
                }
                
            }
        }
        else
        {
            std::cout<<"#rows of matrix_1 != #rows of matrix_2"<<std::endl;
            exit(99);
        }
        return m_new;
    }

/*-------------------------------------------------------MISC--------------------------------------------------------*/
void print()//print entire matrix
{
    for(unsigned int i=0; i<int(mat.size()) ; i++ )
    {
        for(unsigned int j=0; j<int(mat[0].size()) ; j++)
        {
            std::cout<<mat[i][j]<<'\t';
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}
void printcolj(int j) //print column vector j
{
    for(unsigned int i=0; i<mat.size(); i++)
    {
        std::cout<<mat[i][j]<<std::endl;
    }
    std::cout<<std::endl;
    
}
void printrowi(int i) //print row vector i
{
    for(unsigned int j=0; j<mat[i].size(); j++)
    {
        std::cout<<mat[i][j]<<'\t';
    }
    std::cout<<std::endl;
    
}

};

#endif
