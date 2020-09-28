//  curstat.cpp
//  LR test based confidence intervals for MLE
//
//  Created by Piet Groeneboom on 28-05-15.
//  Copyright (c) 2015 Piet Groeneboom. All rights reserved.
//

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
//#include <chrono>
//#include <random>
#include <Rcpp.h>



using namespace std;
using namespace Rcpp;

void    convexmin(int n, double cumw[], double cs[], double y[]);
double  bdf(double A, double B,  int njumps, double *jumploc, double *p, double h, double u);
double  KK(double x);


// [[Rcpp::export]]

List ComputeIntervals(DataFrame input)
{
    double  A,B,c,*grid,*p,step,h;
    double  *Fsmooth,*lowbound,*upbound;
    double  a,b,*cumw,*cs,*F,*jumploc,*y,*y1,*y2;
    int     i,j,k,m,N,n,njumps;
    double  *cumw1,*cs1,L0,L;
    int     *tt,*uu,*tt1,*uu1;
    int     ngrid=1000,npoints=100;
    clock_t StartTime, StopTime;
    double  Time_ConfInt;

    
    
    Rcpp::DataFrame DF = Rcpp::DataFrame(input);
    Rcpp::NumericVector data = DF["V1"];
    Rcpp::IntegerVector freq1 = DF["V2"];
    Rcpp::IntegerVector freq2 = DF["V3"];
    
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "Piet Groeneboom 2015" << std::endl << "For further information see:" << std::endl;
    Rcpp::Rcout << "Nonparametric Estimation under Shape Constraints." << std::endl;
    Rcpp::Rcout << "Piet Groeneboom & Geurt Jongbloed, Cambridge University Press, 2014." << std::endl << std::endl;
    Rcpp::Rcout << "The program produces LR test based (pointwise) 95% confidence intervals for the cdf," << std::endl;
    Rcpp::Rcout << "using the procedure proposed by Banerjee and Wellner." << std::endl << std::endl;
    Rcpp::Rcout << "See Banerjee, M. and Wellner, J.A, Confidence Intervals for current status data," << std::endl;
    Rcpp::Rcout << "Scandinavian Journal of Statistics, 32, 405-424." << std::endl << std::endl;
    Rcpp::Rcout << "The program also computes the SMLE (blue curve, MLE is red)." << std::endl << std::endl;
    
    
    // determine the number of rows of the data frame

    n = (int)data.size();
    
    Rcpp::Rcout << "Number of unique observations:" << std::setw(10) << n << std::endl;
    
    A = 0.0;
    B = (double)data[n-1];
    
    c=B-A;
    step = (B-A)/ngrid;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]= A + i*step;
    
    N=0;
    for (i=0;i<n;i++)
        N += (int)freq2[i];
    
    Rcpp::Rcout << "Total sample size N =" << std::setw(10) << N << std::endl;
    
    Rcpp::Rcout << "Left point estimation interval A = " <<  std::setw(5) << A << std::endl;
    
    Rcpp::Rcout << "Right point estimation interval B = " <<  std::setw(5) << B << std::endl;
    
    // form the data vector for the bootstrap
    
    
    F= new double[n+1];
    cumw= new double[n+1];
    cs= new double[n+1];
    cumw1= new double[n+1];
    cs1= new double[n+1];
    y= new double[n+1];
    y1= new double[n+1];
    y2= new double[n+1];
    jumploc= new double[n+1];
    
    tt1= new int[n+1];
    uu1= new int[n+1];
    tt= new int[n+1];
    uu= new int[n+1];
    
    Fsmooth= new double[ngrid+1];
    p= new double[n+1];
    
    
    F[0]=0;
    cumw[0]=cs[0]=cumw1[0]=cs1[0]=0;
    
    y[0]=y1[0]=y2[0]=0;
    
    StartTime = clock();
    
    for (i=1;i<=n;i++)
    {
        cs[i]=cs[i-1]+(double)freq1[i-1];
        cumw[i]=cumw[i-1]+(double)freq2[i-1];
    }
    
    convexmin(n,cumw,cs,y);
    
    j=0;
    
    jumploc[0]=0;
    
    for (i=1;i<=n;i++)
    {
        if (y[i]>y[i-1])
        {
            j++;
            p[j]=y[i]-y[i-1];
            F[j]=y[i];
            jumploc[j]=data[i-1];
        }
    }
    
    njumps=j;
    
    // bandwidth for SMLE
    h = c*pow(N,-1.0/5);
    
    NumericMatrix out1 = NumericMatrix(njumps+1,2);

    for (i=0;i<=njumps;i++)
    {
        out1(i,0)=jumploc[i];
        out1(i,1) = F[i];
    }
    
     // computation of the SMLE
    
    for (i=0;i<=ngrid;i++)
        Fsmooth[i]=bdf(A,B,njumps,jumploc,p,grid[i],h);

    
    NumericMatrix out2 = NumericMatrix(ngrid+1,2);
    
    for (i=0;i<=ngrid;i++)
    {
        out2(i,0)=grid[i];
        out2(i,1) = Fsmooth[i];
    }
    
    lowbound=new double[npoints];
    upbound=new double[npoints];

    
    Rcpp::Rcout << "Computing Banerjee-Wellner LR test based confidence intervals" << std::endl;
    
    tt[0]=uu[0]=0;
    
    for (i=1;i<=n;i++)
    {
        tt[i]=(double)freq1[i-1];
        uu[i]=(double)freq2[i-1];
    }
    
    L=0;
    for (i=1;i<=n;i++)
    {
        if (tt[i]>0)
            L += log(y[i])*tt[i];
        if (uu[i]-tt[i]>0)
            L += log(1-y[i])*(uu[i]-tt[i]);
    }
    
    for (k=1;k<npoints;k++)
    {
        j=1;
        while ((double)data[j-1]<grid[10*k] && j<=n)
            j++;
        m=j-1;
        
        a=y[m];
      
        for (i=1;i<=n-m;i++)
        {
            uu1[i]=uu[n-i+1];
            tt1[i]=uu[n-i+1]-tt[n-i+1];
        }
        
        for (i=1;i<=n-m;i++)
        {
            cumw1[i]=cumw1[i-1]+(double)uu1[i];
            cs1[i]=cs1[i-1]+(double)tt1[i];
        }
        
        b=0;
        j=0;
        
        while (b<2.26916 && a+j*0.0001<=0.9999)
        {
            L0=0;
            
            j++;
            
            convexmin(m,cumw,cs,y1);
            
            for (i=1;i<=m;i++)
            {
                if (tt[i]>0)
                    L0 += log(fmin(a+j*0.0001,y1[i]))*tt[i];
                if (uu[i]-tt[i]>0)
                    L0 += log(1-fmin(a+j*0.0001,y1[i]))*(uu[i]-tt[i]);
            }
            
            convexmin(n-m,cumw1,cs1,y2);
            
            for (i=1;i<=n-m;i++)
            {
                if (tt1[i]>0)
                    L0 += log(fmin(1-a-j*0.0001,y2[i]))*(tt1[i]);
                if (uu1[i]-tt1[i]>0)
                    L0 += log(1-fmin(1-a-j*0.0001,y2[i]))*(uu1[i]-tt1[i]);
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
            
            convexmin(m,cumw,cs,y1);
            
            for (i=1;i<=m;i++)
            {
                if (tt[i]>0)
                    L0 += log(fmin(a-j*0.0001,y1[i]))*tt[i];
                if (uu[i]-tt[i]>0)
                    L0 += log(1-fmin(a-j*0.0001,y1[i]))*(uu[i]-tt[i]);
            }
            
            convexmin(n-m,cumw1,cs1,y2);
            
            for (i=1;i<=n-m;i++)
            {
                if (tt1[i]>0)
                    L0 += log(fmin(1-a+j*0.0001,y2[i]))*(tt1[i]);
                if (uu1[i]-tt1[i]>0)
                    L0 += log(1-fmin(1-a+j*0.0001,y2[i]))*(uu1[i]-tt1[i]);
            }
            
            b=2*(L-L0);
        }
        
        lowbound[k]=a-j*0.0001;
    }
    
    StopTime  = clock();
    Time_ConfInt   = (double)(StopTime - StartTime)/(double)CLOCKS_PER_SEC;
    
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "The computations took    " << setprecision(10) << Time_ConfInt << "   seconds"  << std::endl;
    
    ofstream file0_("MLE.txt");
    
    if (file0_.is_open())
    {
        for (i=0;i<=njumps;i++)
        {
            file0_ << setprecision(11) << setw(20) << jumploc[i];
            file0_ << setprecision(11) <<  setw(20) << F[i];
            file0_ << "\n";
        }
        file0_ << setprecision(10) << setw(20) << grid[ngrid];
        file0_ << setprecision(11) <<  setw(20) << F[njumps];
        file0_ << "\n";
        file0_.close();
    }


    ofstream file_("CI_Banerjee_Wellner.txt");
    
    if (file_.is_open())
    {
        for (i=0;i<npoints-1;i++)
        {
            file_ << setprecision(10) << setw(20) << grid[10*(i+1)];
            file_ << setprecision(11) <<  setw(20) << lowbound[i+1]
            << setprecision(11) <<  setw(20) << upbound[i+1];
            file_ << "\n";
         }
         file_.close();
    }
    
    NumericMatrix out3 = NumericMatrix(npoints-1,3);
    
    for (i=0;i<npoints-1;i++)
    {
        out3(i,0)=grid[10*(i+1)];
        out3(i,1)=lowbound[i+1];
        out3(i,2)=upbound[i+1];
    }

    
    List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("SMLE")=out2,Rcpp::Named("CI_Banerjee_Wellner")=out3);

    
    Rcpp::Rcout << "Freeing memory" << std::endl;
    
    // free memory

    delete[] F, delete[] cumw, delete[] cs, delete[] cumw1, delete[] cs1,
    delete[] y, delete[] y1, delete[] y2, delete[] jumploc,
    delete[] Fsmooth, delete[] p, delete[] lowbound, delete[] upbound,
    delete[] tt, delete[] uu, delete[] tt1, delete[] uu1;
    
    return out;
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


double KK(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y = (16.0 + 35*x - 35*pow(x,3) + 21*pow(x,5) - 5*pow(x,7))/32.0;
    else
    {
        if (x>1)
            y=1;
        else
            y=0;
        
    }
    
    return y;
}

double bdf(double A, double B, int m, double t[], double p[], double u, double h)
{
    int			k;
    double		t1,t2,t3,sum;
    
    
    sum=0;
    for (k=1;k<=m;k++)
    {
        t1=(u-t[k])/h;
        t2=(u+t[k]-2*A)/h;
        t3=(2*B-u-t[k])/h;
        sum+= (KK(t1)+KK(t2)-KK(t3))*p[k];
    }
    return fmax(0,sum);
}
