//
//  main.cpp
//  Banerjee-Wellner
//
//  Created by Piet on 22/05/15.
//  Copyright (c) 2015 Piet. All rights reserved.
//


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <chrono>
#include <random>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

typedef struct
{
    double t;
    int delta;
}
SampleTime;

double  F0(double x);
double  bias(double t, double h);
int     CompareTime(const void *a, const void *b);
void    cumsum(int m, double cs[], int delta[]);
void    convexmin(int n, double cumw[], double cs[], double y[]);
void    curstatgen(int n, double tt[], int delta[], int seed);

// [[Rcpp::export]]

List ComputeIntervals(int N=1000,int NumB=50)
{
    double          B,*data,*data2,*grid,step;
    double          *data0,*lowbound,*upbound;
    double          *cumw,*cs,*F,*y,*y1,*y2;
    double          *MLE;
    double          a,b,*cumw1,*cs1,L,L0;
    int             i,j,k,m,m1,n,*delta,*delta2;
    int             iter,ngrid,NumIt,npoints,seed;
    int             *percentage;
    clock_t         StartTime, StopTime;
    double          Time_bootstrap;
    
    Rcout << "Piet Groeneboom 2015" << std::endl << "For more information please see:" << std::endl;
    Rcout << "Nonparametric Estimation under Shape Constraints, Section 9.5" << std::endl;
    Rcout << "Piet Groeneboom & Geurt Jongbloed, Cambridge University Press, 2014." << std::endl;
    Rcout << "Simulation of Banerjee-Wellner confidence intervals." << std::endl << std::endl;
    
    NumIt=(int)NumB;
    ngrid=1000;
    n=(int)N;
    npoints=100;
    
    seed=100;
    
    delta = new int[n+1];
    delta2 = new int[n+1];
    data = new double[n+1];
    data0 = new double[n+1];
    data2 = new double[n+1];
    
    delta[0]=delta2[0]=0;
    data[0]=data0[0]=data2[0]=0;
    
    B = 2.0;
    
    step = B/ngrid;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]= i*step;
    
    F= new double[n+1];
    cumw= new double[n+1];
    cs= new double[n+1];
    cumw1= new double[n+1];
    cs1= new double[n+1];
    y= new double[n+1];
    y1= new double[n+1];
    y2= new double[n+1];
    
    
    MLE= new double[ngrid+1];
    
    lowbound=new double[npoints];
    upbound=new double[npoints];
    percentage = new int[npoints];
    
    for (i=0;i<npoints;i++)
        percentage[i]=0;
    
    
    F[0]=0;
    cumw[0]=cs[0]=cumw1[0]=cs1[0]=0;
    for (i=1;i<=n;i++)
        cumw[i]=cumw1[i]=i*1.0;
    
    y[0]=y1[0]=y2[0]=0;
    
    
    StartTime = clock();
    
    for (iter=1;iter<=NumIt;iter++)
    {
        seed++;
        curstatgen(n,data,delta,seed);
        
        cumsum(n,cs,delta);
        convexmin(n,cumw,cs,y);
        
        j=0;
        
        for (i=1;i<=n;i++)
        {
            if (y[i]>y[i-1])
            {
                j++;
                F[j]=y[i];
                data0[j]=data[i];
            }
        }
        
        m=j;
        
        for (i=0;i<=ngrid;i++)
        {
            if (grid[i]<data0[1])
                MLE[i]=0;
            for (j=1;j<m;j++)
            {
                if (data0[j]<= grid[i] && grid[i]<data0[j+1])
                    MLE[i] = F[j];
            }
            if (data0[m]<= grid[i])
                MLE[i]=F[m];
        }
        
        L=0;
        for (i=1;i<=n;i++)
        {
            if (delta[i]==1)
                L += log(y[i]);
            else
                L += log(1-y[i]);
        }
        
        
        for (k=1;k<npoints;k++)
        {
            j=1;
            while (data[j]<grid[10*k] && j<=n)
                j++;
            m1=j-1;
            
            a=y[m1];
            
            for (i=1;i<=n-m1;i++)
                delta2[i]=1-delta[n-i+1];
            
            for (i=1;i<=n-m1;i++)
                cs1[i]=cs1[i-1]+delta2[i];
            
            convexmin(m1,cumw,cs,y1);
            convexmin(n-m1,cumw1,cs1,y2);
            
            b=0;
            j=0;
            
            while (b<2.26916 && a+j*0.0001<=0.9999)
            {
                L0=0;
                
                j++;
                
                for (i=1;i<=m1;i++)
                {
                    if (delta[i]==1)
                        L0 += log(fmin(a+j*0.0001,y1[i]));
                    else
                        L0 += log(1-fmin(a+j*0.0001,y1[i]));
                }
                
                for (i=1;i<=n-m1;i++)
                {
                    if (delta2[i]==1)
                        L0 += log(fmin(1-a-j*0.0001,y2[i]));
                    else
                        L0 += log(1-fmin(1-a-j*0.0001,y2[i]));
                }
                
                b=2*(L-L0);
            }
            
            upbound[k]=a+j*0.0001;
            
            j=0;
            b=0;
            
            while (b<2.26916 && a-j*0.0001>=0.0001)
            {
                L0=0;
                
                j++;
                
                for (i=1;i<=m1;i++)
                {
                    if (delta[i]==1)
                        L0 += log(fmin(a-j*0.0001,y1[i]));
                    else
                        L0 += log(1-fmin(a-j*0.0001,y1[i]));
                }
                
                for (i=1;i<=n-m1;i++)
                {
                    if (delta2[i]==1)
                        L0 += log(fmin(1-a+j*0.0001,y2[i]));
                    else
                        L0 += log(1-fmin(1-a+j*0.0001,y2[i]));
                }
                
                b=2*(L-L0);
                
                //printf(" value LR statistics at %5.3f = %14.10f\n",a+j*0.001,b);
            }
            
            lowbound[k]=a-j*0.0001;
            
            if (F0(grid[k*10])<lowbound[k] || F0(grid[k*10])>upbound[k])
                percentage[k]++;
        }
        
        Rcout  << setw(10) << iter << setprecision(6) <<  setw(15) << F0(grid[500]) << setprecision(6) <<  setw(15) << lowbound[50] << setprecision(6) <<  setw(15) << upbound[50] << setw(10) << percentage[50] << std::endl;
        
    }
    
    StopTime  = clock();
    Time_bootstrap   = (double)(StopTime - StartTime)/(double)CLOCKS_PER_SEC;
    
    Rcout << std::endl << std::endl;
    Rcout << "The computations took    " << setprecision(10) << Time_bootstrap << "   seconds"  << std::endl;

    NumericMatrix out1 = NumericMatrix(m+1,2);
    
    for (i=0;i<=m;i++)
    {
      out1(i,0)=data0[i];
      out1(i,1) = F[i];
    }
    
    NumericMatrix out2 = NumericMatrix(npoints-1,3);
    
    for (i=0;i<npoints-1;i++)
    {
        out2(i,0)=grid[10*(i+1)];
        out2(i,1)=lowbound[i+1];
        out2(i,2)=upbound[i+1];
    }
    
    NumericMatrix out3 = NumericMatrix(npoints-1,2);
    
    for (i=0;i<npoints-1;i++)
    {
        out3(i,0)=grid[10*(i+1)];
        out3(i,1)=(double)percentage[i+1]/NumIt;
    }

    
    ofstream file0_("MLE.txt");
    
    if (file0_.is_open())
    {
        for (i=0;i<=m;i++)
        {
            file0_ << setprecision(11) << setw(20) << data0[i];
            file0_ << setprecision(11) <<  setw(20) << F[i];
            file0_ << "\n";
        }
        file0_ << setprecision(10) << setw(20) << B;
        file0_ << setprecision(11) <<  setw(20) << F[m];
        file0_ << "\n";
        file0_.close();
    }
    
    ofstream file_("CI_MLE.txt");
    
    if (file_.is_open())
    {
        for (i=1;i<npoints;i++)
        {
            file_ << setprecision(10) << setw(20) << grid[10*i];
            file_ << setprecision(11) <<  setw(20) << lowbound[i];
            file_ << setprecision(11) <<  setw(20) << upbound[i];
            file_ << "\n";
        }
        file_.close();
    }
    
    ofstream file2_("percentages.txt");
    
    if (file2_.is_open())
    {
        for (i=1;i<npoints;i++)
        {
            file2_ << setprecision(10) << setw(20) << grid[10*i];
            file2_ << setprecision(11) <<  setw(20) << (double)percentage[i]/NumIt;
            file2_ << "\n";
        }
        file2_.close();
    }

    Rcout << "Making output list" << std::endl;
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    
    List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("CI_Banerjee_Wellner")=out2,Rcpp::Named("percentages")=out3);
    
    
    Rcout << "Freeing memory" << std::endl;
    
    // free memory
    
    delete[] data, delete[] data0, delete[] data2, delete[] delta, delete[] delta2,
    delete[] F, delete[] cumw, delete[] cs, delete[] cumw1, delete[] cs1, delete[] MLE,
    delete[] y, delete[] y2, delete[] lowbound, delete[] upbound;
    
    return out;
    
}

void cumsum(int m, double cs[], int delta[])
{
    int	i;
    
    cs[1]= delta[1];
    for (i=2;i<=m;i++)
        cs[i] = cs[i-1] + delta[i];
}



void convexmin(int n, double cumw[], double cs[], double y[])
{
    int	i, j, m;
    
    y[1] = cs[1]/cumw[1];
    for (i=2;i<=n;i++)
    {
        y[i] = (cs[i]-cs[i-1])/(cumw[i]-cumw[i-1]);
        if (y[i-1]>y[i])
        {
            j = i;
            while (y[j-1] > y[i] && j>1)
            {
                j--;
                if (j>1)
                    y[i] = (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1]);
                else
                    y[i] = cs[i]/cumw[i];
                for (m=j;m<i;m++)	y[m] = y[i];
            }
        }
    }
}


int CompareTime(const void *a, const void *b)
{
    if ((*(SampleTime *) a).t < (*(SampleTime *) b).t)
        return -1;
    if ((*(SampleTime *) a).t > (*(SampleTime *) b).t)
        return 1;
    return 0;
}

void curstatgen(int n, double tt[], int delta[], int seed)
{
    int    i;
    SampleTime *obs;
    double x;
    
    obs = new SampleTime[n];
    
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<> dis(0, 1);
    
    for (i = 1; i <= n;i++)
    {
        x= dis(gen);
        x=-log(1-(1-exp(-2))*x);
        tt[i]=2*dis(gen);
        if (x<=tt[i]) delta[i]=1;
        else delta[i]=0;
    }
    
    for (i=0;i<n;i++)
    {
        obs[i].t=tt[i+1];
        obs[i].delta=delta[i+1];
    }
    
    qsort(obs,n,sizeof(SampleTime),CompareTime);
    
    for (i=1;i<=n;i++)
    {
        tt[i]=obs[i-1].t;
        delta[i]=obs[i-1].delta;
    }
    
    delete[] obs;
}

double F0(double x)
{
    return (1-exp(-x))/(1-exp(-2.0));
}

