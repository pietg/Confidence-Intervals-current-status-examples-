
//  curstat_bootstrapSMLE.cpp
//  CI_SMLE
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
#include <chrono>
#include <random>
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

typedef struct{double alpha, beta;}
weight_t;

double  F0(double x);
double  f1(double x);
double  bias(double t, double h);
double  bias(double t, int m, double tt[], double pp[], double h, double h0);
int     CompareTime(const void *a, const void *b);
int     compare(const void *a, const void *b);
void    cumsum(int m, double cs[], int delta[]);
void    convexmin(int n, double cumw[], double cs[], double y[]);
double  bdf(double A, double B,  int m, double *tt, double *p, double h, double u);
//double  bdf_conv(double B, int m, double data[], double p[], double u, double h);
double  bdf_conv(double B, int ngrid, double grid[], double u, double SMLE0[], double h);
double  K(double x);
double  KK(double x);
double  KK2(double x);
double  Kprime(double x);
void    data_binom(int N, int n, double data[], int delta2[], int **freq, double F[], int seed);
double  varF(int N, int n, int **freq, double *F, double A, double B, double t[], double h, double u);
void curstatgen(int n, double tt[], int delta[], int seed);
weight_t weight(double x);
weight_t weight2(double x);
double dens(double A, double B,  int m, double t[], double p[], double u, double h);
double densprime(double B, int m, double t[], double p[], double u, double h);
double Compute_bandwidth(int n,int m, double tt[], double pp[], double h, double h0, double u);
void    Solve_boundary(double h, double function[], double a[]);
double bdf_conv(double A, double B, int m, double data[], double p[], double u, double h1, double h0);
double KK2(double x, double c);
double Compute_bandwidth(double h0, double u);


// [[Rcpp::export]]

List ComputeIntervals(int N=1000, int NumB=1000)
{
    double          A,B,c,*data,*tt,*data1,*data2,*grid,*pp,*pp2,step,h,h0;
    double          *lowbound,*upbound;
    double          *cumw,*cs,*F,*F2,*y,*y2,*SMLE,*SMLE0,*SMLE1,*SMLE2,*Fsmooth,*cc,*cc1,*cc2;
    double          **f3,*f4,*diff;
    int             i,j,n,m,m_bootstrap,*delta,*delta2,**freq,**freq2;
    int             iter,iter2,ngrid,NumIt,NumIt2,npoints,seed;
    int             below,above,*percentage;
    clock_t         StartTime, StopTime;
    double          Time_bootstrap;

    
    Rcout << "For further information see:" << std::endl;
    Rcout << "Piet Groeneboom and Kim Hendrickx (2016) Confidence intervals for the current status model." << std::endl;
    Rcout << "https://arxiv.org/abs/1611.08299" << std::endl << std::endl;
    Rcout << "The program produces the smooth bootstrap (pointwise) 95% confidence intervals for the cdf," << std::endl;
    Rcout << "using the SMLE." << std::endl << std::endl;
    
    NumIt=(int)NumB;
    NumIt2=1000;
    ngrid=1000;
    n=(int)N;
    npoints=100;
    
    // bandwidth for SMLE
    h  = pow(n,-1.0/5);
    h0 = pow(n,-1.0/9);
    
    seed=1;
    
    Rcout << "Number of observations:" << std::setw(7) << n << std::endl << std::endl;
    Rcout << "Number of samples:" << std::setw(10) << NumIt << std::endl;
    Rcout << "Number of bootstrap samples from each sample:" << std::setw(10) << NumIt2 << std::endl;
 
    below=(int)(0.025*NumIt2-1);
    above=(int)(0.975*NumIt2-1);
    
    m=0;
    
    delta = new int[n+1];
    data = new double[n+1];
    tt = new double[n+1];
    data1 = new double[n+1];
    delta2 = new int[n+1];
    data2 = new double[n+1];
    
    delta[0]=0;
    data[0]=data1[0]=data2[0]=0;
    tt[0]=0;
    
    A = 0.0;
    B = 2.0;
    c=B;
    
    step = B/ngrid;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]= i*step;
    
    F= new double[n+1];
    F2= new double[n+1];
    cumw= new double[n+1];
    cs= new double[n+1];
    y= new double[n+1];
    y2= new double[n+1];
    
    pp= new double[n+1];
    pp2= new double[n+1];
    
    SMLE= new double[ngrid+1];
    SMLE0= new double[ngrid+1];
    SMLE1= new double[ngrid+1];
    SMLE2= new double[ngrid+1];
    Fsmooth= new double[n+1];
    cc = new double[ngrid+1];
    cc2 = new double[ngrid+1];
    cc1 = new double[n+1];
    
    freq = new int*[2];
    for (i=0;i<=1;i++)
        freq[i] = new int[n+1];
    
    freq2 = new int*[2];
    for (i=0;i<=1;i++)
        freq2[i] = new int[n+1];
    
    f4= new double[NumIt2+1];
    
    f3  = new double*[NumIt2+1];
    
    for (iter2=0;iter2<NumIt2+1;iter2++)
        f3[iter2] = new double[npoints];
    
    lowbound=new double[npoints];
    upbound=new double[npoints];
    percentage = new int[npoints];
    
    diff = new double[npoints+1];
    
    for (i=0;i<npoints;i++)
        percentage[i]=0;
    
    
    F[0]=F2[0]=0;
    cumw[0]=cs[0]=0;
    
    SMLE[0]=SMLE1[0]=SMLE2[0]=0;
    
    y[0]=y2[0]=0;
    
    for (i=1;i<=ngrid;i++)
        cc[i]=2;
    
    for (i=1;i<=n;i++)
        cc1[i]=2;
    
    StartTime = clock();
    
    Rcout << "F_0 is standard truncated exponential on [0,2]" << std::endl;
    Rcout << "Observation distribution is Uniform on [0,2]" << std::endl << std::endl;
    
    Rcout << "1000 samples with 1000 bootstrap samples from each sample:" << std::endl << std::endl;
    
    Rcout << "     Iteration  " << "  F_0(1)  " << "     lower bound  "<< "  upper bound  " << "#{F_0(1) not in interval}  " << std::endl << std::endl;
    
    //for (i=1;i<=ngrid;i++)
        //cc[i]=Compute_bandwidth(n,m,tt,pp,h0,h,grid[i]);
            
    //for (i=1;i<=n;i++)
        //cc1[i]=Compute_bandwidth(n,m,tt,pp,h,h0,data[i]);
    
    for (iter=1;iter<=NumIt;iter++)
    {
        seed++;
        
        for (i=1;i<npoints;i++)
            diff[i]=0;
        
        curstatgen(n,data,delta,seed);
        
        for (i=1;i<=n;i++)
        {
            freq[1][i]=delta[i];
            freq[0][i]=1-delta[i];
        }
        
        for (i=1;i<=n;i++)
            cumw[i]=i*1.0;
        
        
        cumsum(n,cs,delta);
        convexmin(n,cumw,cs,y);
        
        j=0;
        
        for (i=1;i<=n;i++)
        {
            if (y[i]>y[i-1])
            {
                j++;
                pp[j]=y[i]-y[i-1];
                F[j]=y[i];
                tt[j]=data[i];
            }
        }
        
        m=j;
        
        //for (i=1;i<=ngrid;i++)
            //cc[i]=Compute_bandwidth(2*h,grid[i]);
        
        for (i=1;i<=ngrid;i++)
            SMLE[i]=bdf(A,B,m,tt,pp,grid[i],cc[i]*h);
        
        for (i=1;i<=ngrid;i++)
            SMLE0[i]=bdf(A,B,m,tt,pp,grid[i],cc[i]*h0);
        
        for (i=1;i<=ngrid;i++)
            SMLE1[i]=bdf_conv(A,B,m,tt,pp,grid[i],cc[i]*h,2*h0);
            //SMLE1[i]=bdf_conv(B,ngrid,grid,grid[i],SMLE0,cc[i]*h);
                    
        for (i=1;i<=n;i++)
            Fsmooth[i]=bdf(A,B,m,tt,pp,data[i],2*h0);
        
        for (iter2=1;iter2<=NumIt2;iter2++)
        {
            seed++;
            data_binom(n,n,data,delta2,freq2,Fsmooth,seed);
            
            for (i=1;i<=n;i++)
            {
                cs[i]=cs[i-1]+(double)freq2[1][i];
                cumw[i]=cumw[i-1]+(double)(freq2[0][i]+freq2[1][i]);
            }
            
            convexmin(n,cumw,cs,y2);
            
            j=0;
            
            for (i=1;i<=n;i++)
            {
                if (y2[i]>y2[i-1])
                {
                    j++;
                    data2[j]=data[i];
                    pp2[j]=y2[i]-y2[i-1];
                    F2[j]=y2[i];
                }
            }
            
            m_bootstrap=j;
            
            //for (i=1;i<=ngrid;i++)
                //cc2[i]=Compute_bandwidth(n,m,tt,pp,2*h,2*h0,grid[i]);
            
            for (i=1;i<=ngrid;i++)
                SMLE2[i]=bdf(A,B,m_bootstrap,data2,pp2,grid[i],cc[i]*h);
            
            
            for (i=1;i<npoints;i++)
                f3[iter2][i]=SMLE2[10*i]-SMLE1[10*i];
                //f3[iter2][i]=(SMLE2[10*i]-SMLE1[10*i])/sqrt(varF(n,n,freq2,y2,0.0,B,data,cc[10*i]*h,grid[10*i]));
            
            for (i=1;i<npoints;i++)
                //diff[i] += SMLE2[10*i]-SMLE0[10*i];
                diff[i] = bias(grid[10*i],cc[10*i]*h);
                //diff[i] = bias(grid[10*i],m,data,pp,2*h,2*h0);
        }
    
        
        for (i=1;i<npoints;i++)
        {
            //diff[i]/=NumIt2;
            
            //diff[i]=fmin(0,diff[i]);
            
            for (iter2=0;iter2<NumIt2;iter2++)
                f4[iter2]=f3[iter2+1][i];
            
            qsort(f4,NumIt2,sizeof(double),compare);
                    
            //lowbound[i] = SMLE[10*i]-f4[above]-bias(grid[10*i],cc[10*i]*h);
            lowbound[i] = SMLE[10*i]-f4[above]-diff[i];
            //lowbound[i] = SMLE[10*i]-f4[above];
            //lowbound[i] = SMLE[10*i]-f4[above]-bias(grid[10*i],m,tt,pp,2*h,2*h0);
            //lowbound[i] = SMLE[10*i]-diff[i]-f4[above]*sqrt(varF(n,n,freq,y,0.0,B,data,cc[10*i]*h,grid[10*i]));
            //lowbound[i] = SMLE[10*i]-bias(grid[10*i],cc[10*i]*h)-f4[above]*sqrt(varF(n,n,freq,y,0.0,B,data,cc[10*i]*pow(n,-1.0/5),grid[10*i]));
            
            upbound[i] = SMLE[10*i]-f4[below]-diff[i];
            //upbound[i] = SMLE[10*i]-f4[below];
            
            //upbound[i] = SMLE[10*i]-f4[below]-bias(grid[10*i],m,tt,pp,2*h,2*h0);

            //upbound[i]= SMLE[10*i]-diff[i]-f4[below]*sqrt(varF(n,n,freq,y,0.0,B,data,cc[10*i]*h,grid[10*i]));
            
            if (F0(grid[i*10])<lowbound[i] || F0(grid[i*10])>upbound[i])
                percentage[i]++;
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
      out1(i,0)=tt[i];
      out1(i,1) = F[i];
    }
    
    NumericMatrix out2 = NumericMatrix(ngrid+1,2);
    
    for (i=0;i<=ngrid;i++)
    {
      out2(i,0)=grid[i];
      out2(i,1) = SMLE[i];
    }
    
    
    NumericMatrix out3 = NumericMatrix(npoints-1,3);
    
    for (i=0;i<npoints-1;i++)
    {
        out3(i,0)=grid[10*(i+1)];
        out3(i,1)=lowbound[i+1];
        out3(i,2)=upbound[i+1];
    }
    
    NumericMatrix out4 = NumericMatrix(npoints-1,2);
    
    for (i=0;i<npoints-1;i++)
    {
        out4(i,0)=grid[10*(i+1)];
        out4(i,1)=(double)percentage[i+1]/NumIt2;
    }


    
    ofstream file0_("MLE.txt");
    
    if (file0_.is_open())
    {
        for (i=0;i<=m;i++)
        {
            file0_ << setprecision(11) << setw(20) << tt[i];
            file0_ << setprecision(11) <<  setw(20) << F[i];
            file0_ << "\n";
        }
        file0_ << setprecision(10) << setw(20) << grid[ngrid];
        file0_ << setprecision(11) <<  setw(20) << F[m];
        file0_ << "\n";
        file0_.close();
    }
    
    
    ofstream file_("CI_SMLE.txt");
    
    if (file_.is_open())
    {
        for (i=1;i<npoints;i++)
        {
            file_ << setprecision(10) << setw(20) << grid[10*i];
            file_ << setprecision(11) <<  setw(20) << lowbound[i] << setprecision(11) <<  setw(20) << upbound[i];
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
    
    ofstream file3_("diff.txt");
    
    if (file3_.is_open())
    {
        for (i=1;i<npoints;i++)
        {
            file3_ << setprecision(10) << setw(20) << grid[10*i];
            file3_ << setprecision(11) <<  setw(20) << diff[i];
            file3_ << "\n";
        }
        file3_.close();
    }

    
    Rcout << std::endl;
    
    Rcout << "Making output list" << std::endl;
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    
    List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("SMLE")=out2,Rcpp::Named("CI_SMLE")=out3,Rcpp::Named("percentages")=out4);

    
    // free memory
    
    delete[] data; delete[] delta; delete[] tt; delete[] delta2; delete[] SMLE;  delete[] SMLE1; delete[] SMLE2; delete[] SMLE0;
    delete[] F; delete[] F2; delete[] Fsmooth; delete[] cumw;
    delete[] cs; delete[] y; delete[] y2;  delete[] data1; delete[] data2;
    delete[] pp; delete[] pp2; delete[] lowbound; delete[] upbound; delete[] cc; delete[] cc1; delete[] cc2;
    
    for (i = 0;i<2;i++)
    {
        delete[] freq[i]; delete[] freq2[i];
    }
    
    delete[] freq; delete[] freq2;
    
    for (iter2 = 0;iter2 < NumIt2;iter2++)
        delete[] f3[iter2];
    delete[] f3;
    
    return out;
}

/*double Compute_bandwidth(double h, double u)
{
    double v,w;
    double numerator, denominator;
    // c is the length of the support of the observation distribution
    
    if (u<h)
        numerator = pow(((70*pow(u,2)*(39*pow(h,11) - 143*pow(h,9)*pow(u,2) + 429*pow(h,7)*pow(u,4) - 429*pow(h,6)*pow(u,5)+143*pow(h,4)*pow(u,7) - 39*pow(h,2)*pow(u,9) + 5*pow(u,11)))/(429.*pow(h,13)))*F0(u)*(1-F0(u)),1.0/5);
    else
    {
        if (u<=2-h)
            numerator = pow((350.0/429)*F0(u)*(1-F0(u)),1.0/5);
        else
            numerator = pow(((70*pow(2-u,2)*(39*pow(h,11) - 143*pow(h,9)*pow(2-u,2) + 429*pow(h,7)*pow(2-u,4) - 429*pow(h,6)*pow(2-u,5)+143*pow(h,4)*pow(2-u,7) - 39*pow(h,2)*pow(2-u,9) + 5*pow(2 - u,11)))/(429.*pow(h,13)))*F0(u)*(1-F0(u)),1.0/5);
    }
    
    
    if (u<h)
    {
        v=u/h;
        w=(pow(1-v,6)*(64 + v*(69 + 5*v*(6 + v))))/5576;
        denominator = pow(fabs(f1(u))*(1.0/9-w),2.0/5);
    }
    else
    {
        if (u<=2-h)
            denominator = pow(fabs(f1(u))/9,2.0/5);
        else
        {
            v=(2-u)/h;
            w=(pow(1-v,6)*(64 + v*(69 + 5*v*(6 + v))))/5576;
            denominator = pow(fabs(f1(u))*(1.0/9-w),2.0/5);
        }
    }
    
    return numerator/denominator;
}*/

double bias(double t, double h)
{
    double v;
    if (t>= h && t<=2-h)
        return -SQR(h)*exp(-t)/(18*(1-exp(-2)));
    else
    {
        if (t<h)
        {
            v=t/h;
            return -SQR(h)*exp(-t)*(35*v/64 - v*v + 35*pow(v,3)/48 - 7*pow(v,5)/32 + pow(v,7)/16 - 5*pow(v,9)/576)/(2*(1-exp(-2)));
        }
        else
        {
            v=(2.0-t)/h;
            return -SQR(h)*exp(-t)*(35*v/64 - v*v + 35*pow(v,3)/48 - 7*pow(v,5)/32 + pow(v,7)/16 - 5*pow(v,9)/576)/(2*(1-exp(-2)));
        }
    }
    
}

double bias(double t, int m, double tt[], double pp[], double h,  double h0)
{
    double v;
    if (t>= h && t<=2-h)
        return (1.0/18)*SQR(h)*densprime(2,m,tt,pp,t,h0);
    else
    {
        if (t<h)
        {
            v=t/h;
            return 0.5*SQR(h)*densprime(2,m,tt,pp,t,h0)*(35*v/64 - v*v + 35*pow(v,3)/48 - 7*pow(v,5)/32 + pow(v,7)/16 - 5*pow(v,9)/576);
        }
        else
        {
            v=(2.0-t)/h;
            return 0.5*SQR(h)*densprime(2,m,tt,pp,t,h0)*(35*v/64 - v*v + 35*pow(v,3)/48 - 7*pow(v,5)/32 + pow(v,7)/16 - 5*pow(v,9)/576);
        }
    }
    
}


double F0(double t)
{
    if (t>=0 && t<=2)
        return (1-exp(-t))/(1-exp(-2));
    else
    {
        if (t<0)
            return 0;
        else
            return 1;
    }
}

double f1(double t)
{
    if (t>=0 && t<=2)
        return -exp(-t)/(1-exp(-2));
    else
    {
        return 0;
    }
}

double Compute_bandwidth(int n,int m, double tt[], double pp[], double h0, double h, double u)
{
    double numerator,denominator;
        
    //numerator = 2*(350.0/429)*bdf(0,2,m,tt,pp,u,h)*(1-bdf(0,2,m,tt,pp,u,2*h));
    numerator = 2*(350.0/429)*F0(u)*(1-F0(u));
    numerator = pow(numerator,1.0/5);
    
    //denominator = (1.0/9)*fmax(0.1,fabs(densprime(2,m,tt,pp,u,2*h0)));
    denominator = (1.0/9)*exp(2-u)/(exp(2)-1);
    denominator = pow(denominator,-2.0/5);
    
    return numerator*denominator;
}

weight_t weight(double x)
{
    short i;
    double y1, y2, y3, y4;
    double *xx;
    weight_t temp;
    
    xx = new double[10];;
    
    xx[1] = x;
    for (i=2;i<=9;i++)
        xx[i] = x * xx[i - 1];
    
    y1 = 0.5 + (35.0*xx[1])/32.0-(35.0*xx[3])/32.0+(21.0*xx[5])/32.0-(5.0*xx[7])/32.0;
    y2 = -35.0/256.0+(35.0 * xx[2])/64.0-(105.0*xx[4])/128.0+(35.0*xx[6])/64.0
    -(35.0*xx[8])/256.0;
    y3 = 1.0/18 + 25.0*xx[3]/96 - 21.0*xx[5]/32 + 15.0*xx[7]/32 - 35.0*xx[9]/288;
    
    
    y4 = y1 * y3 - (y2*y2);
    temp.alpha = y3 / y4;
    temp.beta = -y2 / y4;
    
    delete[] xx;
    return temp;
}


weight_t weight2(double x)
{
    double y1, y2;
    weight_t temp;
    
    y1 = -71680*pow(1-x,3)*SQR(128 + 5*x*(-65 + x*(69 + 7*(-5 + x)*x)))/
    ((pow(1 + x,7)*SQR(5359 + 5*x*(-3550 + x*(4909 + x*(-3620 + x*(1517 + 35*(-10 + x)*x)))))));
    
    y2 = 645120*pow(1 - x,3)*(-35 + x*(47 + 5*(-5 + x)*x))*
    (128 + 5*x*(-65 + x*(69 + 7*(-5 + x)*x)))/
    (pow(1 + x,7)*SQR(5359 + 5*x*(-3550 + x*(4909 + x*(-3620 + x*(1517 + 35*(-10 + x)*x))))));
    
    temp.alpha = y1;
    temp.beta = y2;
    
    return temp;
}

double Kprime(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y=-(105.0/16)*x*SQR(1-u);
    else
        y=0.0;
    
    return y;
}

/*double densprime(double B, int m, double t[], double p[], double u, double h)
{
    int k;
    double        rho,x,sum;
    weight_t    weights,weights2;
    
    sum=0;
    
    if (u>=h && u<=B-h)
    {
        for (k=1;k<=m;k++)
        {
            x=(u-t[k])/h;
            sum += Kprime(x)*p[k]/SQR(h);
        }
    }
    else
    {
        if (u<h)
        {
            rho = u/h;
            weights=weight(rho);
            weights2=weight2(rho);
            for (k=1;k<=m;k++)
            {
                x=(u-t[k])/h;
                sum += (Kprime(x)*weights.alpha + (x*Kprime(x)+K(x))*weights.beta +K(x)*weights2.alpha + x*K(x)*weights2.beta)*p[k]/SQR(h);
            }
        }
        else
        {
            if (u>B-h)
            {
                rho = (B-u)/h;
                weights=weight(rho);
                weights2=weight2(rho);
                for (k=1;k<=m;k++)
                {
                    x=(u-t[k])/h;
                    sum += (Kprime(x)*weights.alpha - (x*Kprime(x)+K(x))*weights.beta
                            - K(x)*weights2.alpha + x*K(x)*weights2.beta)*p[k]/SQR(h);
                }
            }
            
        }
    }
    return sum;
}*/


double densprime(double B, int m, double t[], double p[], double u, double h)
{
    int            k;
    double        t1,t2,t3,sum;
    
    
    sum=0;
    for (k=1;k<=m;k++)
    {
        t1=(u-t[k])/h;
        t2=(u+t[k])/h;
        t3=(2*B-u-t[k])/h;
        sum+= (Kprime(t1)+Kprime(t2)-Kprime(t3))*p[k]/SQR(h);
    }
    return fmax(0,sum);
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

/*double bdf(double A, double B, int m, double tt[], double pp[], double u, double h)
{
    int     k;
    double  t,sum,sum2;
    double *a,*function;
    
    a = new double[4];
    function = new double[4];
    
    sum=sum2=0;
    
    if (u>=h && u<=B-h)
    {
        for (k=1;k<=m;k++)
        {
            t=(u-tt[k])/h;
            sum+= KK(t)*pp[k];
        }
    }
    else
    {
        if (u<h)
        {
            for (k=1;k<=m;k++)
            {
                t=(h-tt[k])/h;
                sum2 += KK(t)*pp[k];
            }

            function[1] = sum2;
            function[2] = dens(A,B,m,tt,pp,h,h);
            function[3] = densprime(B,m,tt,pp,h,h);
            
            Solve_boundary(h,function,a);
            sum = a[1]*u+a[2]*SQR(u)+a[3]*pow(u,3);
        }
        else
        {
            for (k=1;k<=m;k++)
            {
                t=(B-h-tt[k])/h;
                sum2 += KK(t)*pp[k];
            }
            sum2=1-sum2;
            
            function[1] = sum2;
            function[2] = dens(A,B,m,tt,pp,B-h,h);
            function[3] = -densprime(B,m,tt,pp,B-h,h);
            
            Solve_boundary(h,function,a);
            sum = 1-(a[1]*(B-u)+a[2]*SQR(B-u)+a[3]*pow(B-u,3));
        }
    }
    delete[] a; delete[] function;
    return sum;
}*/

void Solve_boundary(double h, double function[], double a[])
{
    double y1,y2,y3;
    
    y1= function[1];
    y2= function[2];
    y3= function[3];
    
    a[1] = -0.5*(-6*y1 + 4*h*y2 - SQR(h)*y3)/h;
    a[2] = -((3*y1 - 3*h*y2 + SQR(h)*y3)/SQR(h));
    a[3] = -0.5*(-2*y1 + 2*h*y2 - SQR(h)*y3)/pow(h,3);
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

double  bdf_conv(double B, int ngrid, double grid[], double u, double SMLE0[], double h)
{
    int            i;
    double        t1,t2,t3,sum;
    
    
    sum=0;
    
    for (i=1;i<=ngrid;i++)
    {
        t1=(u-grid[i])/h;
        t2=(u+grid[i])/h;
        t3=(2*B-u-grid[i])/h;
        
        sum+= (KK(t1)+KK(t2)-KK(t3))*(SMLE0[i]-SMLE0[i-1]);
        
    }
    
    return sum;
}

/*double bdf_conv(double B, int m, double data[], double p[], double u, double h)
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
}*/


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




void curstatgen(int n, double tt[], int delta[], int seed)
{
    int	i;
    SampleTime2 *obs;
    double x;
    
    obs = new SampleTime2[n];
    
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
    
    qsort(obs,n,sizeof(SampleTime2),CompareTime);
    
    for (i=1;i<=n;i++)
    {
        tt[i]=obs[i-1].t;
        delta[i]=obs[i-1].delta;
    }
    
    delete[] obs;
}


void data_binom(int N, int n, double data[], int delta2[], int **freq, double F[], int seed)
{
    int	i,j;
    double x;
    
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<> dis(0, 1);
    
    for (i=0;i<2;i++)
    {
        for (j=0;j<=n;j++)
            freq[i][j]=0;
    }
    
    
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
    
    sum = sum/(N*N);
    
    return sum;
}

double dens(double A, double B,  int m, double t[], double p[], double u, double h)
{
    int k;
    double      t1,sum;
    
    sum=0;
    
    for (k=1;k<=m;k++)
    {
        t1=(u-t[k])/h;
        sum += K(t1)*p[k]/h;
    }
    
    return sum;
}

double bdf_conv(double A, double B, int m, double data[], double p[], double u, double h1, double h0)
{
    int            i;
    double        t1,t2,t3,sum;
    
    
    sum=0;
    
    for (i=1;i<=m;i++)
    {
        t1=(u-data[i])/h1;
        t2=(u+data[i])/h1;
        t3=(2*B-u-data[i])/h1;
        
        sum+= (KK2(t1,h0/h1)+KK2(t2,h0/h1)-KK2(t3,h0/h1))*p[i];
        
    }
    
    return sum;
}

double KK2(double x, double c)
{
    double y=0;
    
    if (0<c && c<1)
    {
        if (1 + c < x)
            y=1;
        else
        {
            if (1 + c == x && c + x > 1)
                y = pow(-1 + x,7)/pow(c,7);
            else
            {
                if (c >= 1 + x && 1 + c + x > 0)
                    y = -2.84546182983683e-7*(pow(1 + c + x,8)*(2145*pow(c,6) + 40*pow(c,5)*(-429 + 131*x) +pow(c,4)*(56199 + 5*x*(-4952 + 625*x)) + 40*c*pow(1 + x,2)*(-429 + x*(239 + 5*(-11 + x)*x)) -5*pow(1 + x,3)*(-429 + x*(239 + 5*(-11 + x)*x)) - 16*pow(c,3)*(5577 + x*(-1728 + 5*x*(3 + 8*x))) -pow(c,2)*(1 + x)*(-56199 + x*(28551 + 5*x*(-1161 + 89*x)))))/pow(c,7);
                else
                {
                    if (c < 1 + x && c + x <= 1)
                        y = (6864 - 35*(-429 + 143*pow(c,2) - 39*pow(c,4) + 5*pow(c,6))*x - 455*(33 - 22*pow(c,2) + 5*pow(c,4))*pow(x,3) + 1001*(9 - 5*pow(c,2))*pow(x,5) - 2145*pow(x,7))/13728;
                     else
                     {
                         if (1 + c > x && c + x > 1)
                         y = (2145*pow(c,14) - 22400*pow(c,13)*x - 640640*pow(c,9)*x*pow(-1 + pow(x,2),2) +
                              525525*pow(c,8)*pow(-1 + pow(x,2),3) - 58240*pow(c,11)*x*(-3 + 5*pow(x,2)) +
                              21021*pow(c,12)*(-1 + 5*pow(x,2)) + 105105*pow(c,10)*(1 - 6*pow(x,2) + 5*pow(x,4)) -
                              54912*pow(c,7)*(-48 - 35*x + 35*pow(x,3) - 21*pow(x,5) + 5*pow(x,7)) -
                              5005*pow(c,4)*pow(-1 + x,7)*(3 + x)*(7 + x*(4 + x)) + 15015*pow(c,6)*pow(-1 + x,5)*(35 + x*(47 + 5*x*(5 + x))) + 91*pow(c,2)*pow(-1 + x,9)*(231 + x*(159 + 5*x*(9 + x))) - 5*pow(-1 + x,11)*(429 + x*(239 + 5*x*(11 + x))))/(3.514368e6*pow(c,7));
                           else
                               y=0;
                     }
                }
            }
        }
    }
    if (c>=1)
    {
        if ((c == 1 && x >= 2) || (c > 1 && 1 + c <= x))
            y=1;
        else
        {
            if (c > 1 && -1 < c + x && c+x <= 1)
                y = -2.84546182983683e-7*(pow(1 + c + x,8)*(2145*pow(c,6) + 40*pow(c,5)*(-429 + 131*x) +
                    pow(c,4)*(56199 + 5*x*(-4952 + 625*x)) + 40*c*pow(1 + x,2)*(-429 + x*(239 + 5*(-11 + x)*x)) -
                    5*pow(1 + x,3)*(-429 + x*(239 + 5*(-11 + x)*x)) - 16*pow(c,3)*(5577 + x*(-1728 + 5*x*(3 + 8*x))) -
                    pow(c,2)*(1 + x)*(-56199 + x*(28551 + 5*x*(-1161 + 89*x)))))/pow(c,7);
            else
            {
                if (c == 1 && -2 < x && x<= 0)
                    y=(pow(2 + x,8)*(6864 + x*(-16256 + x*(16976 + 5*x*(-1888 + x*(538 + 5*(-16 + x)*x))))))/3.514368e6;
                else
                {
                    if (c > 1 + x && c + x > 1)
                        y = (-175*x + 13*(528*pow(c,7) + 35*pow(c,2)*(3 - 11*pow(c,2) + 33*pow(c,4))*x -
                            35*(5 - 22*pow(c,2) + 33*pow(c,4))*pow(x,3) + 77*(-5 + 9*pow(c,2))*pow(x,5) - 165*pow(x,7)))/
                       (13728.*pow(c,7));
                    else
                    {
                        if (c == 1 && 0 < x && x < 2)
                            y = 0.5 + (350*x)/429. - (35*pow(x,3))/66. + (7*pow(x,5))/24. - (5*pow(x,7))/32. + (35*pow(x,8))/512. -
                           (7*pow(x,10))/1536. + (35*pow(x,12))/135168. - (25*pow(x,14))/3.514368e6;
                        else
                        {
                            if (c > 1 && x < 1 + c && c <= 1 + x)
                                y = (2145*pow(c,14) - 22400*pow(c,13)*x - 640640*pow(c,9)*x*pow(-1 + pow(x,2),2) +
                                 525525*pow(c,8)*pow(-1 + pow(x,2),3) - 58240*pow(c,11)*x*(-3 + 5*pow(x,2)) +
                                 21021*pow(c,12)*(-1 + 5*pow(x,2)) + 105105*pow(c,10)*(1 - 6*pow(x,2) + 5*pow(x,4)) -
                                 54912*pow(c,7)*(-48 - 35*x + 35*pow(x,3) - 21*pow(x,5) + 5*pow(x,7)) -
                                 5005*pow(c,4)*pow(-1 + x,7)*(3 + x)*(7 + x*(4 + x)) + 15015*pow(c,6)*pow(-1 + x,5)*(35 + x*(47 + 5*x*(5 + x))) +
                                 91*pow(c,2)*pow(-1 + x,9)*(231 + x*(159 + 5*x*(9 + x))) - 5*pow(-1 + x,11)*(429 + x*(239 + 5*x*(11 + x))))/(3.514368e6*pow(c,7));
                        }
                    }
                }
            }
        }
    }
    return y;
}
