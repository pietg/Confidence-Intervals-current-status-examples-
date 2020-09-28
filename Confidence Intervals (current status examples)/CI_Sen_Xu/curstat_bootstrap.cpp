//
//  main.cpp
//  trial
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
//#include <chrono>
//#include <random>
#include <Rcpp.h>


using namespace std;
using namespace Rcpp;


typedef struct
{
    double t;
    int freq1;
    int freq2;
}
SampleTime;

typedef struct
{
    double t;
    int delta;
}
SampleTime2;


int     CompareTime(const void *a, const void *b);
int     compare(const void *a, const void *b);
void    convexmin(int n, double cumw[], double cs[], double y[]);
double  bdf(double A, double B,  int njumps, double *jumploc, double *p, double h, double u);
double  KK(double x);
void    data_binom(int N, int n, double data[], int delta2[], int **freq, double F[]);


// [[Rcpp::export]]

List ComputeIntervals(DataFrame input)
{
    double          A,B,c,*data,*data1,*data2,*grid,*p,*p2,step,h;
    double          *Fsmooth,*Fsmooth2,*tt,**f3,*lowbound,*upbound,*f4;
    double          *cumw,*cs,*F,*F1,*F2,*jumploc,*y,*y2,*SMLE,*MLE,*MLE2;
    int             i,j,k,m,N,n,*delta,**freq,njumps,*delta2,*freq1,*freq2;
    int             percentile1,percentile2,iter,ngrid=1000,NumIt=1000,npoints=100;
    clock_t         StartTime, StopTime;
    double          Time_bootstrap;

    
    DataFrame DF = Rcpp::DataFrame(input);
    NumericVector data0 = DF["V1"];
    IntegerVector freq01 = DF["V2"];
    IntegerVector freq02 = DF["V3"];
    
    Rcout << std::endl;
    Rcout << "Piet Groeneboom 2015" << std::endl << "For further information see:" << std::endl;
    Rcout << "Piet Groeneboom & Geurt Jongbloed, Cambridge University Press, 2014." << std::endl << std::endl;
    Rcout << "The program produces (pointwise) 95% confidence intervals for the cdf," << std::endl;
    Rcout << "using the SMLE in the smooth bootstrap procedure proposed by Sen and Xu." << std::endl << std::endl;
    Rcout << "See Sen and Xu, Model based bootstrap methods for interval censored data" << std::endl;
    Rcout << "Computational Statistics & Data Analysis, 2015, Pages 121–129" << std::endl;
    Rcout << "The program also computes the SMLE (blue curve, MLE is red)." << std::endl << std::endl;
    
    
    // determine the number of rows of the data frame
    
    percentile1=(int)(0.025*NumIt);
    percentile2=(int)(0.975*NumIt);
    
    n = (int)data0.size();
    
    Rcout << "Number of unique observations:" << std::setw(7) << n << std::endl << std::endl;
    Rcout << "Number of bootstrap samples:" << std::setw(10) << NumIt << std::endl;
    
    Rcout << "percentile1: " <<  std::setw(7) << percentile1 << std::endl;
    Rcout << "percentile2: " <<  std::setw(7) << percentile2 << std::endl;
    
    
    
    data= new double[n+1];
    freq1 = new int[n+1];
    freq2 = new int[n+1];
    
    for (i=1;i<=n;i++)
    {
        data[i]=(double)data0[i-1];
        freq1[i]=(int)freq01[i-1];
        freq2[i]=(int)freq02[i-1];
    }
 
    
    A = 0.0;
    B = data[n];
    c=B-A;
    
    data[0]=0;
    
    step = (B-A)/ngrid;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]= A + i*step;
    
    N=0;
    for (i=1;i<=n;i++)
        N += freq2[i];
    
    Rcout << "Sample size N =" << std::setw(5) << N << std::endl;
    Rcout << "Left point of observation interval A = " <<  std::setw(5) << A << std::endl;
    Rcout << "Right point of observation interval B = " <<  std::setw(5) << B << std::endl;
    
    // form the data vector for the bootstrap
    
    data = new double[N+1];
    delta = new int[N+1];
    tt=new double[N+1];
    delta2 = new int[N+1];
    
    
    data[0]=0;
    
    j=0;
    
    for (i=1;i<=n;i++)
    {
        for (k=1;k<=freq1[i];k++)
        {
            j++;
            data[j]=data0[i];
            delta[j]=1;

        }
        for (k=1;k<=freq2[i]-freq1[i];k++)
        {
            j++;
            data[j]=data0[i];
            delta[j]=0;
        }
    }
    
    
    F= new double[n+1];
    F1= new double[N+1];
    F2= new double[n+1];
    cumw= new double[n+1];
    cs= new double[n+1];
    y= new double[n+1];
    y2= new double[n+1];
    jumploc= new double[n+1];
    data1= new double[n+1];
    data2= new double[n+1];
    
    Fsmooth= new double[ngrid+1];
    Fsmooth2= new double[ngrid+1];
    SMLE= new double[ngrid+1];
    p= new double[n+1];
    p2= new double[n+1];
    MLE= new double[ngrid+1];
    MLE2= new double[ngrid+1];
    
    freq = new int*[2];
    for (i=0;i<2;i++)
        freq[i] = new int[n+1];
    
    f3  = new double*[NumIt+1];
    
    for (iter=0;iter<NumIt+1;iter++)
        f3[iter] = new double[npoints];
    
    
    F[0]=F1[0]=F2[0]=0;
    cumw[0]=cs[0]=0;
    
    y[0]=y2[0]=0;
    
    StartTime = clock();
    
    for (i=1;i<=n;i++)
    {
        cs[i]=cs[i-1]+(double)freq1[i];
        cumw[i]=cumw[i-1]+(double)freq2[i];
    }

    
    convexmin(n,cumw,cs,y);
    
    jumploc[0]=0;
    F[0]=0;
    
    j=0;
    
    for (i=1;i<=n;i++)
    {
        if (y[i]>y[i-1])
        {
            j++;
            p[j]=y[i]-y[i-1];
            F[j]=y[i];
            jumploc[j]=data0[i];
        }
    }
    
    njumps=j;
    
    NumericMatrix out1 = NumericMatrix(njumps+1,2);
    
    for (i=0;i<=njumps;i++)
    {
      out1(i,0)=jumploc[i];
      out1(i,1) = F[i];
    }
    
    
    // bandwidth for SMLE
    h = c*pow(N,-1.0/5);

    
    for (i=1;i<=N;i++)
        F1[i]=bdf(A,B,njumps,jumploc,p,data[i],h);
    
    SMLE[0]=0;
    
    for (i=1;i<=ngrid;i++)
        SMLE[i]=bdf(A,B,njumps,jumploc,p,grid[i],h);
    
    for (i=0;i<=ngrid;i++)
    {
        if (0<= grid[i] && grid[i]<jumploc[1])
            MLE[i]=0;
        for (j=1;j<njumps;j++)
        {
            if (jumploc[j]<= grid[i] && grid[i]<jumploc[j+1])
                MLE[i] = F[j];
        }
        if (jumploc[njumps]<= grid[i])
            MLE[i]=F[njumps];
    }
    
    
    for (iter=0;iter<NumIt;iter++)
    {
        data_binom(N,n,data,delta2,freq,F1);
        
        cumw[0]=cs[0]=0;
        
        for (i=1;i<=n;i++)
        {
            cs[i]=cs[i-1]+(double)freq[1][i];
            cumw[i]=cumw[i-1]+(double)(freq[0][i]+freq[1][i]);
        }
        
        convexmin(n,cumw,cs,y2);
        
        j=0;
        
        for (i=1;i<=n;i++)
        {
            if (y2[i]>y2[i-1])
            {
                j++;
                data2[j]=data0[i];
                p2[j]=y2[i]-y2[i-1];
                F2[j]=y2[i];
            }
        }
        
        m=j;
    
        
        for (i=0;i<=ngrid;i++)
        {
            if (grid[i]<data2[1])
                MLE2[i]=0;
            for (j=1;j<m;j++)
            {
                if (data2[j]<= grid[i] && grid[i]<data2[j+1])
                    MLE2[i] = F2[j];
            }
            if (data2[m]<= grid[i])
                MLE2[i]=F2[m];
        }
        
        for (i=1;i<npoints;i++)
            f3[iter][i]=MLE2[10*i]-SMLE[10*i];
    }
    
    StopTime  = clock();
    Time_bootstrap   = (double)(StopTime - StartTime)/(double)CLOCKS_PER_SEC;
    
    f4= new double[NumIt+1];
    
    lowbound=new double[npoints];
    upbound=new double[npoints];
    
    for (i=1;i<npoints;i++)
    {
        for (iter=0;iter<NumIt;iter++)
            f4[iter]=f3[iter][i];
        
        qsort(f4,NumIt,sizeof(double),compare);
        
        lowbound[i]= fmax(0,MLE[10*i]-f4[percentile2-1]);
        upbound[i]= fmin(1,MLE[10*i]-f4[percentile1-1]);
    }
    
    
    
    Rcout << std::endl << std::endl;
    Rcout << "The computations took    " << setprecision(10) << Time_bootstrap << "   seconds"  << std::endl;
    

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
    
    // determination of the bandwidth of the SMLE and density using a simple rule of thumb
    
    
    Rcout << std::endl << std::endl;
    
    Rcout << "Sorting results of bootstrap" << std::endl << std::endl;
    
    Rcout << "Forming output with SMLE" << std::endl;
    
    NumericMatrix out2 = NumericMatrix(ngrid+1,2);
    
    // computation of the SMLE
    
    ofstream file_("CI_Sen_Xu.txt");
    
    if (file_.is_open())
    {
        for (i=1;i<npoints;i++)
        {
            file_ << setprecision(10) << setw(20) << grid[10*i];
            file_ << setprecision(11) <<  setw(20) << lowbound[i]
            << setprecision(11) <<  setw(20) << upbound[i];
            file_ << "\n";
        }
        file_.close();
    }
    
    
    ofstream file1_("SMLE.txt");
    
    if (file1_.is_open())
    {
        for (i=0;i<=ngrid;i++)
        {
            file1_ << setprecision(10) << setw(20) << grid[i];
            file1_ << setprecision(11) <<  setw(20) << SMLE[i];
            file1_ << "\n";
        }
        file1_.close();
    }
    
    
    for (i=0;i<=ngrid;i++)
    {
        out2(i,0)=grid[i];
        out2(i,1) = SMLE[i];
    }
    
    Rcout << "Forming matrix with confidence intervals" << std::endl;
    
    NumericMatrix out3 = NumericMatrix(npoints-1,3);
    
    for (i=0;i<npoints-1;i++)
    {
        out3(i,0)=grid[10*(i+1)];
        out3(i,1)=lowbound[i+1];
        out3(i,2)=upbound[i+1];
    }
    
    Rcout << "Making output list" << std::endl;
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    
    List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("SMLE")=out2,Rcpp::Named("CI_Sen_Xu")=out3);
    
    
    Rcout << "Freeing memory" << std::endl;
    
    // free memory
    
    delete[] data, delete[] delta, delete[] freq1, delete[] freq2, delete[] tt, delete[] delta2,
    delete[] F, delete[] F1, delete[] F2, delete[] cumw,
    delete[] cs, delete[] y, delete[] y2, delete[] jumploc,  delete[] data1, delete[] data2,
    delete[] Fsmooth,  delete[] p, delete[] p2, delete[] lowbound, delete[] upbound;
    
    for (i = 0;i<2;i++)
        delete[] freq[i];
    delete[] freq;
    
    for (iter = 0;iter < NumIt;iter++)
        delete[] f3[iter];
    delete[] f3;
    
    return out;
    
    return 0;
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
        //sum+= KK(t1)*p[k];
    }
    return fmax(0,sum);
}


int compare(const void *a, const void *b)
{
    double x = *(double*)a;
    double y = *(double*)b;
    
    if (x < y)
        return -1;
    if (x > y)
        return 1;
    return 0;
}

int CompareTime(const void *a, const void *b)
{
    if ((*(SampleTime *) a).t < (*(SampleTime *) b).t)
        return -1;
    if ((*(SampleTime *) a).t > (*(SampleTime *) b).t)
        return 1;
    return 0;
}


void data_binom(int N, int n, double data[], int delta2[], int **freq, double F[])
{
    int	i,j;
    double x;
    
    for (i=0;i<2;i++)
    {
        for (j=0;j<=n;j++)
            freq[i][j]=0;
    }
  

    for (i=1;i<=N;i++)
    {
        j = rand();
        x= ((double) j)/(double)RAND_MAX;
        
        if (x<=F[i])
            delta2[i]=1;
        else
            delta2[i]=0;
    }
    
    j=0;
    
    for (i=1;i<=N;i++)
    {
        if (data[i]>data[i-1])
        {
            j++;
            freq[delta2[i]][j]=1;
        }
        else
            freq[delta2[i]][j]++;
    }
}



/* void data_binom(int N, int n, double data[], int delta2[], int **freq, double F[])
{
    int	i,j;
    double x;
    
    for (i=0;i<2;i++)
    {
        for (j=0;j<=n;j++)
            freq[i][j]=0;
    }
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis(0,1);
    
    for (i=1;i<=N;i++)
    {
        x= dis(gen);
        
        if (x<=F[i])
            delta2[i]=1;
        else
            delta2[i]=0;
    }
    
    j=0;
    
    for (i=1;i<=N;i++)
    {
        if (data[i]>data[i-1])
        {
            j++;
            freq[delta2[i]][j]=1;
        }
        else
            freq[delta2[i]][j]++;
    }
}*/

