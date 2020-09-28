//
//  curstat_smooth_bootstrapSMLE.cpp
//  CI_SMLE_smoothbootstrap
//
//  Created by Piet Groeneboom on 22/05/15.
//  Copyright (c) 2015 Piet Groeneboom. All rights reserved.
//


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <Rcpp.h>


#define SQR(x) ((x)*(x))

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
double  bdf_conv(double B, int m, double data[], double p[], double u, double h);
double  K(double x);
double  KK(double x);
double  KK2(double x);
double  varF(int N, int n, int **freq, double *F, double A, double B, double t[], double h, double u);
void    data_binom(int N, int n, double data[], int delta2[], int **freq, double F[]);

// [[Rcpp::export]]

List ComputeIntervals(DataFrame input)
{
    double          A,B,c,*data0,*data,*data1,*data2,*grid,*p,*p2,step,h;
    double          *tt,**f3,*lowbound,*upbound,*f4;
    double          *cumw,*cs,*F,*F2,*jumploc,*y,*y2,*SMLE,*SMLE1,*SMLE2,*Fsmooth,*cc,*cc1;
    int             i,j,k,m,N,n,*delta,njumps,*delta2,*freq1,*freq2,**frequence1,**frequence2;
    int             percentile1,percentile2,iter,ngrid=1000,NumIt=1000,npoints=100;
    clock_t         StartTime, StopTime;
    double          Time_bootstrap;

    
    DataFrame DF = Rcpp::DataFrame(input);
    NumericVector xx = DF["V1"];
    IntegerVector freq01 = DF["V2"];
    IntegerVector freq02 = DF["V3"];
    
    Rcout << std::endl;
    Rcout << "For further information see:" << std::endl;
    Rcout << "Piet Groeneboom and Kim Hendrickx (2016) Confidence intervals for the current status model." << std::endl;
    Rcout << "https://arxiv.org/abs/1611.08299" << std::endl << std::endl;
    Rcout << "The program produces smooth bootstrap based (pointwise) 95% confidence intervals for the cdf," << std::endl;
    Rcout << "using the resampling from the SMLE." << std::endl << std::endl;
    Rcout << "The data set is the hepatitis A data set, provided to us by Niels Keiding." << std::endl;
    Rcout << "The program also computes the MLE (red curve, SMLE is blue)." << std::endl << std::endl;
    
    
    // determine the number of rows of the data frame
    
    percentile1=(int)(0.025*NumIt);
    percentile2=(int)(0.975*NumIt);
    
    n = (int)xx.size();
    
    Rcout << "Number of unique observations:" << std::setw(7) << n << std::endl << std::endl;
    Rcout << "Number of bootstrap samples:" << std::setw(10) << NumIt << std::endl;
    
    Rcout << "percentile1: " <<  std::setw(7) << percentile1 << std::endl;
    Rcout << "percentile2: " <<  std::setw(7) << percentile2 << std::endl;
    
    data0= new double[n+1];
    freq1= new int[n+1];
    freq2= new int[n+1];
    
    data0[0]=0;
    
    for (i=1;i<=n;i++)
    {
      data0[i]=(double)xx[i-1];
      freq1[i]=(int)freq01[i-1];
      freq2[i]=(int)freq02[i-1];
    }
    
    A = 0;
    B = data0[n];
    c= B-A;
    
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
    cc = new double[ngrid+1];
    
    for (i=1;i<=ngrid;i++)
        cc[i]=0.5+grid[i]/grid[ngrid];
    
  
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
    F2= new double[n+1];
    cumw= new double[n+1];
    cs= new double[n+1];
    y= new double[n+1];
    y2= new double[n+1];
    jumploc= new double[n+1];
    data1= new double[n+1];
    data2= new double[n+1];
    cc1 = new double[N+1];
 
    Fsmooth = new double[N+1];
    SMLE= new double[ngrid+1];
    SMLE1= new double[ngrid+1];
    SMLE2= new double[ngrid+1];
    p= new double[n+1];
    p2= new double[n+1];
    
    frequence1 = new int*[2];
    for (i=0;i<2;i++)
        frequence1[i] = new int[n+1];
    
    frequence2 = new int*[2];
    for (i=0;i<2;i++)
        frequence2[i] = new int[n+1];
    
    f3  = new double*[NumIt+1];
    
    for (iter=0;iter<NumIt+1;iter++)
        f3[iter] = new double[npoints];
    
    
    F[0]=F2[0]=0;
    cumw[0]=cs[0]=0;
    
    y[0]=y2[0]=0;
    
    StartTime = clock();
    
    for (i=1;i<=n;i++)
    {
        cs[i]=cs[i-1]+(double)freq1[i];
        cumw[i]=cumw[i-1]+(double)freq2[i];
        frequence1[1][i]=freq1[i];
        frequence1[0][i]=freq2[i]-freq1[i];
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
    
    SMLE[0]=0;
    
    for (i=1;i<=ngrid;i++)
        SMLE[i]=bdf(A,B,njumps,jumploc,p,grid[i],cc[i]*h);
    
    // SMLE1 is the estimate which uses the convolution kernel
    
    for (i=1;i<=ngrid;i++)
        SMLE1[i]=bdf_conv(B,njumps,jumploc,p,grid[i],cc[i]*h);
    
    
    NumericMatrix out2 = NumericMatrix(ngrid+1,2);
    
    for (i=0;i<=ngrid;i++)
    {
        out2(i,0)=grid[i];
        out2(i,1) = SMLE[i];
    }
    
    j=0;
    
    for (i=1;i<=n;i++)
    {
        for (k=1;k<=freq2[i];k++)
        {
            j++;
            cc1[j]=0.5+data0[i]/grid[ngrid];
        }
    }

    
    for (i=1;i<=N;i++)
        Fsmooth[i]=bdf(A,B,njumps,jumploc,p,data[i],cc1[i]*h);

    
    for (iter=0;iter<NumIt;iter++)
    {
        data_binom(N,n,data,delta2,frequence2,Fsmooth);
        
        cumw[0]=cs[0]=0;
        
        for (i=1;i<=n;i++)
        {
            cs[i]=cs[i-1]+(double)frequence2[1][i];
            cumw[i]=cumw[i-1]+(double)(frequence2[0][i]+frequence2[1][i]);
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
        
        SMLE2[0]=0;
        
        for (i=1;i<=ngrid;i++)
            SMLE2[i]=bdf(A,B,m,data2,p2,grid[i],cc[i]*h);
        
        
        for (i=1;i<npoints;i++)
            f3[iter][i]=(SMLE2[10*i]-SMLE1[10*i])/sqrt(varF(N,n,frequence2,y2,0,B,data0,cc[10*i]*h,grid[10*i]));
        
        // alternatively, one can use non-Studentized intervals:
        //for (i=1;i<npoints;i++)
            //f3[iter][i]=SMLE2[10*i]-SMLE1[10*i];
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
        
        lowbound[i]= SMLE[10*i]-f4[percentile2-1]*sqrt(varF(N,n,frequence1,y,0,B,data0,cc[10*i]*h,grid[10*i]));
        upbound[i] = SMLE[10*i]-f4[percentile1-1]*sqrt(varF(N,n,frequence1,y,0,B,data0,cc[10*i]*h,grid[10*i]));
        
        // alternatively, one can use non-Studentized intervals:
        //lowbound[i]= SMLE[10*i]-f4[percentile2-1];
        //upbound[i] = SMLE[10*i]-f4[percentile1-1];
        
    }
    
    
    Rcout << std::endl << std::endl;
    Rcout << "The computations took    " << setprecision(10) << Time_bootstrap << "   seconds"  << std::endl;
    
    NumericMatrix out3 = NumericMatrix(npoints-1,3);
    
    for (i=0;i<npoints-1;i++)
    {
        out3(i,0)=grid[10*(i+1)];
        out3(i,1)=lowbound[i+1];
        out3(i,2)=upbound[i+1];
    }

    
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
    
    
    ofstream file_("CI_SMLE.txt");
    
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
        for (i=1;i<=ngrid;i++)
        {
            file1_ << setprecision(10) << setw(20) << grid[i];
            file1_ << setprecision(11) <<  setw(20) << SMLE[i];
            file1_ << "\n";
        }
        file1_.close();
    }
    
    Rcout << std::endl;
    
    Rcout << "Making output list" << std::endl;
    
    // make the list for the output, containing the MLE, SMLE and the bootstrap confidence intervals
    
    List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("SMLE")=out2,Rcpp::Named("CI_SMLE")=out3);

    
    // free memory
    
    delete[] data0, delete[] freq1, delete[] freq2,
    delete[] data, delete[] delta, delete[] tt, delete[] delta2,
    delete[] F, delete[] F2, delete[] cumw, delete[] cs, delete[] y, delete[] y2,
    delete[] jumploc, delete[] data1, delete[] data2, delete[] Fsmooth,
    delete[] SMLE, delete[] SMLE1, delete[] SMLE2, delete[] cc, delete[] cc1, delete[] p, delete[] p2;
    
    for (i = 0;i<2;i++)
        delete[] frequence1[i];
    delete[] frequence1;
    
    for (i = 0;i<2;i++)
        delete[] frequence2[i];
    delete[] frequence2;
    
    delete[] lowbound, delete[] upbound, delete[] f4;
    
    for (iter = 0;iter < NumIt;iter++)
        delete[] f3[iter];
    delete[] f3;
    
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

double K(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y=(35.0/32)*pow(1-u,3);
    else
        y=0.0;
    
    return y;
}


double varF(int N, int n, int **freq, double *F, double A, double B, double t[], double h, double u)
{
    int			i;
    double		t1,t2,t3,sum;
    
    
    sum=0;
    
    for (i=1;i<=n;i++)
    {
        t1=(u-t[i])/h;
        t2=(u+t[i]-2*A)/h;
        t3=(2*B-u-t[i])/h;
        
        sum += SQR(K(t1)-K(t2)-K(t3))*(SQR(F[i]-1)*freq[1][i]+SQR(F[i])*freq[0][i]);
    }
    
    sum = sum/N;
    
    return sum;
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
    return sum;
}

double KK2(double x)
{
    double y;
    
    y=0;
    
    if (x<=-2)
        y=0;
    
    if (x>=2)
        y=1;
    
    if (x>-2 && x<0)
        y=pow(2.0 + x,8)*(6864 - 16256*x + 16976*SQR(x) - 9440*pow(x,3) + 2690*pow(x,4) - 400*pow(x,5) + 25*pow(x,6))/3514368.0;
    
    
    if (x>=0 && x<2)
        y = 0.5 + 350*x/429.0 - 35*pow(x,3)/66.0 + 7*pow(x,5)/24.0 - 5*pow(x,7)/32.0 + 35*pow(x,8)/512.0 - 7*pow(x,10)/1536.0 + 35*pow(x,12)/135168.0 - 25*pow(x,14)/3514368.0;
    
    if (x==2)
        y = (1.0/32)*(70*pow(x,4) - 84*pow(x,5) + 35*pow(x,6) - 5*pow(x,7));
    
    return y;
}

double bdf_conv(double B, int m, double data[], double p[], double u, double h)
{
    int			i;
    double		t1,t2,t3,sum;
    
    
    sum=0;
    
    for (i=1;i<=m;i++)
    {
        t1=(u-data[i])/h;
        t2=(u+data[i])/h;
        t3=(2*B-u-data[i])/h;
        
        sum+= (KK2(t1)+KK2(t2)-KK2(t3))*p[i];
        
    }
    
    return sum;
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

