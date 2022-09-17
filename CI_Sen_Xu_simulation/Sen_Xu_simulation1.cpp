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
#include <chrono>
#include <random>
#include <Rcpp.h>

#define SQR(x) ((x)*(x))

using namespace std;
using namespace Rcpp;

typedef struct
{
    double t;
    int delta;
}
SampleTime;


typedef struct{double alpha, beta;}
weight_t;

double  F0(double x);
double  bias(double t, double h);
int     CompareTime(const void *a, const void *b);
int     compare(const void *a, const void *b);
void    cumsum(int m, double cs[], int delta[]);
void    convexmin(int n, double cumw[], double cs[], double y[]);
double  bdf(double A, double B, int m, double tt[], double pp[], double u, double h);
double  bdf(double A, double B, int m, double tt[], double pp[], double u, double h, double h1, double h0);
double  Kprime(double x);
double  K(double x);
double  KK(double x);
void    data_binom(int n,  int delta2[], double F[]);
double  varF(int N, int n, int **freq, double *F, double A, double B, double t[], double h, double u);
void    curstatgen(int n, double tt[], int delta[]);
weight_t weight(double x);
weight_t weight2(double x);
double  dens(double A, double B,  int m, double t[], double p[], double u, double h);
double  densprime(double B, int m, double t[], double p[], double u, double h);
double  Compute_bandwidth(int n,int m, double tt[], double pp[], double u, double h, double h1, double h0);
void    Solve_boundary(double h, double function[], double a[]);

// [[Rcpp::export]]

List ComputeIntervals(int N=1000)
{
    double          A,B,c,*data0,*data,*data2,*grid,*p,*p2,step,h,h0,h1;
    double          *lowbound,*upbound,*cc,*cc1;
    double          *cumw,*cs,*F,*F1,*F2,*y,*y2,*SMLE;
    double          **f3,*f4,*MLE,*MLE2;
    int             i,j,m,m2,n,*delta,*delta2;
    int             iter,iter2,ngrid,NumIt,NumIt2,npoints;
    int             below,above,*percentage;
    clock_t         StartTime, StopTime;
    double          Time_bootstrap;
    
    Rcout << "Piet Groeneboom 2015" << std::endl << "For more information please see:" << std::endl;
    Rcout << "Nonparametric Estimation under Shape Constraints, pp. 8-10 and 276-279," << std::endl;
    Rcout << "Piet Groeneboom & Geurt Jongbloed, Cambridge University Press, 2014." << std::endl << std::endl;
    Rcout << "See also Sen and Xu, Model based bootstrap methods for interval censored data" << std::endl;
    Rcout << "Computational Statistics & Data Analysis, 2015, Pages 121–129" << std::endl << std::endl;

    NumIt=1000;
    NumIt2=1000;
    ngrid=1000;
    n=(int)N;
    npoints=100;
    
    // bandwidth for SMLE
    h   = 2*pow(n,-1.0/5);
    h1  = 2*pow(n,-1.0/7);
    h0  = 2*pow(n,-1.0/9);
    
    Rcout << "Number of observations:" << std::setw(7) << n << std::endl << std::endl;
    Rcout << "Number of bootstrap samples:" << std::setw(10) << NumIt2 << std::endl;
    
    below=(int)(0.025*NumIt2-1);
    above=(int)(0.975*NumIt2-1);
    
    delta = new int[n+1];
    delta2 = new int[n+1];
    data = new double[n+1];
    data0 = new double[n+1];
    data2 = new double[n+1];
    
    delta[0]=delta2[0]=0;
    data[0]=data0[0]=data2[0]=0;
    
    A = 0.0;
    B = 2.0;
    c=B;
    
    step = B/ngrid;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]= i*step;
    
    F= new double[n+1];
    F1= new double[n+1];
    F2= new double[n+1];
    cumw= new double[n+1];
    cs= new double[n+1];
    y= new double[n+1];
    y2= new double[n+1];
    
    p= new double[n+1];
    p2= new double[n+1];

    SMLE  = new double[ngrid+1];
    MLE   = new double[ngrid+1];
    MLE2  = new double[ngrid+1];
    cc    = new double[ngrid+1];
    cc1   = new double[n+1];
    

    f4= new double[NumIt2+1];
    
    f3  = new double*[NumIt2+1];
    
    for (iter2=0;iter2<NumIt2+1;iter2++)
        f3[iter2] = new double[npoints];

    lowbound=new double[npoints];
    upbound=new double[npoints];
    percentage = new int[npoints];
    
    for (i=0;i<npoints;i++)
        percentage[i]=0;
    
    
    F[0]=F1[0]=F2[0]=0;
    cumw[0]=cs[0]=0;
    for (i=1;i<=n;i++)
        cumw[i]=i*1.0;
    
    y[0]=y2[0]=0;
    
    StartTime = clock();
    
    Rcout << "F_0 is standard truncated exponential on [0,2]" << std::endl;
    Rcout << "Observation distribution is Uniform on [0,2]" << std::endl << std::endl;
    
    Rcout << "1000 samples with 1000 bootstrap samples from each sample:" << std::endl << std::endl;
    
    Rcout << "     Iteration  " << "  F_0(1)  " << "     lower bound  "<< "  upper bound  " << "#{F_0(1) not in interval}  " << std::endl << std::endl;
    

    
    for (iter=1;iter<=NumIt;iter++)
    {
        curstatgen(n,data,delta);
        
        cumsum(n,cs,delta);
        convexmin(n,cumw,cs,y);
        
        j=0;
        
        for (i=1;i<=n;i++)
        {
            if (y[i]>y[i-1])
            {
                j++;
                p[j]=y[i]-y[i-1];
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

        for (i=1;i<=ngrid;i++)
            cc[i]=Compute_bandwidth(n,m,data0,p,grid[i],c*h0,c*h0,c*h0);
        
        for (i=1;i<=n;i++)
            cc1[i]=Compute_bandwidth(n,m,data0,p,data[i],c*h0,c*h0,c*h0);
        
        SMLE[0]=0;
        
        for (i=1;i<=ngrid;i++)
            SMLE[i]=bdf(A,B,m,data0,p,grid[i],cc[i]*h0,cc[i]*h0,cc[i]*h0);
        
        for (i=1;i<=n;i++)
            F1[i]= bdf(A,B,m,data0,p,data[i],cc1[i]*h0,cc1[i]*h0,cc1[i]*h0);
        
        
        for (iter2=1;iter2<=NumIt2;iter2++)
        {
            data_binom(n,delta2,F1);
            
            cumsum(n,cs,delta2);
            convexmin(n,cumw,cs,y2);
            
            j=0;
            
            for (i=1;i<=n;i++)
            {
                if (y2[i]>y2[i-1])
                {
                    j++;
                    data2[j]=data[i];
                    p2[j]=y2[i]-y2[i-1];
                    F2[j]=y2[i];
                }
            }
            
            m2=j;
            
            for (i=0;i<=ngrid;i++)
            {
                if (grid[i]<data2[1])
                    MLE2[i]=0;
                for (j=1;j<m2;j++)
                {
                    if (data2[j]<= grid[i] && grid[i]<data2[j+1])
                        MLE2[i] = F2[j];
                }
                if (data2[m2]<= grid[i])
                    MLE2[i]=F2[m2];
            }
            
            for (i=1;i<npoints;i++)
                f3[iter2][i]=MLE2[10*i]-SMLE[10*i];
        }
        
        for (i=1;i<npoints;i++)
        {
            for (iter2=0;iter2<NumIt2;iter2++)
                f4[iter2]=f3[iter2+1][i];
            
            qsort(f4,NumIt2,sizeof(double),compare);
    
            lowbound[i]= MLE[10*i]-f4[above];
            upbound[i]= MLE[10*i]-f4[below];
            
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
      out1(i,0)=data0[i];
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
        file0_ << setprecision(10) << setw(20) << A;
        file0_ << setprecision(11) <<  setw(20) << F[0];
        file0_ << "\n";
        for (i=1;i<=m;i++)
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
            file_ << setprecision(11) <<  setw(20) << lowbound[i]
            << setprecision(11) <<  setw(20) << upbound[i];
            file_ << "\n";
        }
        file_.close();
    }
    
    
   ofstream file1_("MLE.txt");
    
    if (file1_.is_open())
    {
        file1_ << setprecision(10) << setw(20) << A;
        file1_ << setprecision(11) <<  setw(20) << F[0];
        file1_ << "\n";
        for (i=1;i<=m;i++)
        {
            file1_ << setprecision(10) << setw(20) << data0[i];
            file1_ << setprecision(11) <<  setw(20) << F[i];
            file1_ << "\n";
        }
        file1_ << setprecision(10) << setw(20) << B;
        file1_ << setprecision(11) <<  setw(20) << F[m];
        file1_ << "\n";

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
    
    ofstream file3_("SMLE.txt");
    
    if (file3_.is_open())
    {
        file3_ << setprecision(10) << setw(20) << A;
        file3_ << setprecision(11) <<  setw(20) << 0.0;
        file3_ << "\n";
        for (i=1;i<=ngrid;i++)
        {
            file3_ << setprecision(10) << setw(20) << grid[i];
            file3_ << setprecision(11) <<  setw(20) << SMLE[i];
            file3_ << "\n";
        }
        file3_.close();
    }

    Rcpp::Rcout << "Making output list" << std::endl;
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    
    List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("SMLE")=out2,Rcpp::Named("CI_Sen_Xu")=out3,Rcpp::Named("percentages")=out4);
    
    
    Rcpp::Rcout << "Freeing memory" << std::endl;
    
    // free memory
    
    delete[] data; delete[] delta; delete[] delta2; delete[] SMLE;
    delete[] F;  delete[] F1; delete[] F2; delete[] cumw; delete[] cs; delete[] MLE; delete[] MLE2;
    delete[] y; delete[] y2; delete[] data0;  delete[] data2;
    delete[] p; delete[] p2; delete[] lowbound; delete[] upbound; delete[] f4;
    delete[] cc; delete[] cc1;
    
    for (iter2 = 0;iter2 < NumIt2;iter2++)
        delete[] f3[iter2];
    delete[] f3;
    
    return out;
    
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

double Compute_bandwidth(int n,int m, double tt[], double pp[], double u, double h, double h1, double h0)
{
    double numerator,denominator;
            
    numerator = 2*(350.0/429)*bdf(0,2,m,tt,pp,u,h,h1,h0)*(1-bdf(0,2,m,tt,pp,u,h,h1,h0));
    //numerator = 2*(350.0/429)*F0(u)*(1-F0(u));
    numerator = pow(numerator,1.0/5);
    
    denominator = (1.0/9)*fmax(0.01,fabs(densprime(2,m,tt,pp,u,h0)));
    //denominator = (1.0/9)*exp(2-u)/(exp(2)-1);
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

double densprime(double B, int m, double t[], double p[], double u, double h)
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
}

void cumsum(int m, double cs[], int delta[])
{
    int    i;
    
    cs[1]= delta[1];
    for (i=2;i<=m;i++)
        cs[i] = cs[i-1] + delta[i];
}



void convexmin(int n, double cumw[], double cs[], double y[])
{
    int    i, j, m;
    
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
                for (m=j;m<i;m++)    y[m] = y[i];
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

/*double bdf(double A, double B, int m, double t[], double p[], double u, double h)
{
    int            k;
    double        t1,t2,t3,sum;
    
    
    sum=0;
    for (k=1;k<=m;k++)
    {
        t1=(u-t[k])/h;
        t2=(u+t[k]-2*A)/h;
        t3=(2*B-u-t[k])/h;
        sum+= (KK(t1)+KK(t2)-KK(t3))*p[k];
        //sum+= KK(t1)*p[k];
    }
    return sum;
}*/

double bdf(double A, double B, int m, double tt[], double pp[], double u, double h, double h2, double h3)
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
            
            //function[2] = dens(A,B,m,tt,pp,h,h2);
            //function[3] = densprime(B,m,tt,pp,h,h3);
            
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
            
            //function[2] = dens(A,B,m,tt,pp,B-h,h2);
            //function[3] = -densprime(B,m,tt,pp,B-h,h3);
            
            Solve_boundary(h,function,a);
            sum = 1-(a[1]*(B-u)+a[2]*SQR(B-u)+a[3]*pow(B-u,3));
        }
    }
    delete[] a; delete[] function;
    return sum;
}

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




void curstatgen(int n, double tt[], int delta[])
{
    int    i,j;
    SampleTime *obs;
    double x;
    
    obs = new SampleTime[n];
    
    //std::random_device rd;
    //std::mt19937 generator(rd());
    //std::uniform_real_distribution<> dis(0, 1);
    
    for (i = 1; i <= n;i++)
    {
        j = rand();
        x= ((double) j)/(double)RAND_MAX;
        x=-log(1-(1-exp(-2))*x);
        j = rand();
        tt[i]=2*((double) j)/(double)RAND_MAX;
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

void data_binom(int n,  int delta2[], double F[])
{
    int    i,j;
    double x;
    
    for (i=1;i<=n;i++)
    {
        j = rand();
        x= ((double) j)/(double)RAND_MAX;
        
        if (x<=F[i])
            delta2[i]=1;
        else
            delta2[i]=0;
    }
}

/*double dens(double A, double B,  int m, double t[], double p[], double u, double h)
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
}*/


double dens(double A, double B,  int njumps, double jumploc[], double p[], double u, double h)
{
    int i;
    double        rho,sum,x;
    weight_t    weights;
    
    sum=0;
    
    if (u>=A+h && u<=B-h)
    {
        for (i=1;i<=njumps;i++)
        {
            x=(u-jumploc[i])/h;
            sum+= K(x)*p[i]/h;
        }
    }
    else
    {
        if (u<A+h)
        {
            rho=(u-A)/h;
            weights=weight(rho);
            for (i=1;i<=njumps;i++)
            {
                x=(u-jumploc[i])/h;
                sum += (K(x)*weights.alpha + x*K(x)*weights.beta)*p[i]/h;
            }
        }
        else
        {
            if (u>B-h)
            {
                rho = (B-u)/h;
                weights=weight(rho);
                for (i=1;i<=njumps;i++)
                {
                    x=(u-jumploc[i])/h;
                    sum += (K(x)*weights.alpha - x*K(x)*weights.beta)*p[i]/h;
                }
            }
        }
    }
    return fmax(0,sum);
}

