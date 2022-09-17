//
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
double  bias(double t, double h);
int     CompareTime(const void *a, const void *b);
int     compare(const void *a, const void *b);
void    cumsum(int m, double cs[], int delta[]);
void    convexmin(int n, double cumw[], double cs[], double y[]);
double  bdf(double A, double B,  int njumps, double *jumploc, double *p, double h, double u);
double  K(double x);
double  KK(double x);
double  Kprime(double x);
void    data_binom(int N, int n, double data[], int delta2[], int **freq, double F[]);
void    data_bootstrap(int n, int *m, double x[], double x2[], double data2[],
                       int **freq, int delta[], int delta2[], int seed);
double  varF(int N, int n, int **freq, double *F, double A, double B, double t[], double h, double u);
void    curstatgen(int n, double tt[], int delta[], int seed);
double  c0(double t);
weight_t weight(double x);
weight_t weight2(double x);
double  dens(double A, double B,  int m, double t[], double p[], double u, double h);
double  densprime(double B, int m, double t[], double p[], double u, double h);
double  Compute_bandwidth(int n,int m, double tt[], double pp[], double h, double h0, double u);
void    Solve_boundary(double h, double function[], double a[]);

// [[Rcpp::export]]

List ComputeIntervals(int N=1000, int NumB =1000)
{
    double          A,B,c,*data,*data1,*data2,*grid,*p,*p2,step,h;
    double          *data0,*tt,*lowbound,*upbound;
    double          *cumw,*cs,*F,*F2,*jumploc,*y,*y2,*SMLE,*SMLE2,*cc;
    double          **f3,*f4,*diff;
    int             i,j,m,m2,n,*delta,*delta2,**freq,**freq2,njumps,seed;
    int             iter,iter2,ngrid,NumIt,NumIt2,npoints;
    int             below,above,*percentage;
    clock_t         StartTime, StopTime;
    double          Time_bootstrap;

    
    Rcout << "Piet Groeneboom 2015" << std::endl << "For further information see:" << std::endl;
    Rcout << "Nonparametric Estimation under Shape Constraints, pp. 8-10 and Section 9.5," << std::endl;
    Rcout << "Piet Groeneboom & Geurt Jongbloed, Cambridge University Press, 2014." << std::endl << std::endl;
    Rcout << "The program produces the (naive) bootstrap (pointwise) 95% confidence intervals for the cdf," << std::endl;
    Rcout << "using the SMLE." << std::endl << std::endl;
      
    // determine the number of rows of the data frame
    
    NumIt= (int)NumB;
    NumIt2=1000;
    ngrid=1000;
    n= (int)N;
    npoints=100;
    
    seed=100;
    
    Rcout << "Number of observations:" << std::setw(7) << n << std::endl << std::endl;
    Rcout << "Number of bootstrap samples:" << std::setw(10) << NumIt2 << std::endl;
 
    below=(int)(0.025*NumIt2-1);
    above=(int)(0.975*NumIt2-1);
    
    njumps=0;
    
    delta = new int[n+1];
    data = new double[n+1];
    data0 = new double[n+1];
    data1 = new double[n+1];
    delta2 = new int[n+1];
    data2 = new double[n+1];
    
    delta[0]=0;
    data[0]=data1[0]=data2[0]=0;
    
    A = 0.0;
    B = 2.0;
    c=B;
    
    step = B/ngrid;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]= i*step;
    
    F= new double[n+1];
    F2= new double[n+1];
    tt=new double[n+1];
    cumw= new double[n+1];
    cs= new double[n+1];
    y= new double[n+1];
    y2= new double[n+1];
    
    p= new double[n+1];
    p2= new double[n+1];
    jumploc= new double[n+1];
    
    SMLE= new double[ngrid+1];
    SMLE2= new double[ngrid+1];
    cc = new double[ngrid+1];
    
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
    
    for (i=1;i<=ngrid;i++)
        cc[i]=c;
        //cc[i]=c0(grid[i]);
        //cc[i]=0.5+0.5*grid[i]/grid[ngrid];
    
    lowbound=new double[npoints];
    upbound=new double[npoints];
    percentage = new int[npoints];
    
    diff = new double[npoints];
    
    for (i=0;i<npoints;i++)
        percentage[i]=0;
    
    
    F[0]=F2[0]=0;
    cumw[0]=cs[0]=0;
    
    y[0]=y2[0]=0;
    
    // bandwidth for SMLE
    h = pow(n,-1.0/5);
    
    StartTime = clock();
    
    Rcout << "F_0 is standard truncated exponential on [0,2]" << std::endl;
    Rcout << "Observation distribution is Uniform on [0,2]" << std::endl << std::endl;
    
    Rcout << "Number of observations:" << std::setw(7) << n << std::endl << std::endl;
    Rcout << "Number of bootstrap samples for eac sample:" << std::setw(10) << NumIt2 << std::endl;
    
    Rcout << "     Iteration  " << "  F_0(1)  " << "     lower bound  "<< "  upper bound  " << "#{F_0(1) not in interval}  " << std::endl << std::endl;
    
    for (i=1;i<npoints;i++)
        diff[i] = bias(grid[10*i],2*h);
    
    for (iter=1;iter<=NumIt;iter++)
    {
        seed++;
        curstatgen(n,data,delta,seed);
        
        for (i=1;i<=n;i++)
        {
            freq[1][i]=delta[i];y
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
                p[j]=y[i]-y[i-1];
                F[j]=y[i];
                data0[j]=data[i];
            }
        }
        
        njumps=j;
        
        
        SMLE[0]=0;
        
        for (i=1;i<=ngrid;i++)
            SMLE[i]=bdf(A,B,njumps,data0,p,grid[i],cc[i]*h);
        
        
        for (iter2=1;iter2<=NumIt2;iter2++)
        {
            seed++;
            data_bootstrap(n,&m,data,tt,data1,freq2,delta,delta2,seed);
            
            for (i=1;i<=m;i++)
            {
                cs[i]=cs[i-1]+(double)freq2[1][i];
                cumw[i]=cumw[i-1]+(double)(freq2[0][i]+freq2[1][i]);
            }
            
            convexmin(m,cumw,cs,y2);
            
            j=0;
            
            for (i=1;i<=m;i++)
            {
                if (y2[i]>y2[i-1])
                {
                    j++;
                    data2[j]=data1[i];
                    p2[j]=y2[i]-y2[i-1];
                    F2[j]=y2[i];
                }
            }
            
            m2=j;
            
            SMLE2[0]=0;
            
            for (i=1;i<=ngrid;i++)
                SMLE2[i]=bdf(A,B,m2,data2,p2,grid[i],cc[i]*h);
            
            
            for (i=1;i<npoints;i++)
                //f3[iter2][i]=SMLE2[10*i]-SMLE[10*i];
                f3[iter2][i]=(SMLE2[10*i]-SMLE[10*i])/sqrt(varF(n,m,freq2,y2,0.0,B,data1,cc[10*i]*h,grid[10*i]));
            
            //for (i=1;i<npoints;i++)
                //diff[i] = bias(grid[10*i],cc[10*i]*h);
                //diff[i] += SMLE2[10*i]-SMLE[10*i];
        }
        
        for (i=1;i<npoints;i++)
        {
            //diff[i]/=NumIt2;
            
            for (iter2=0;iter2<NumIt2;iter2++)
                f4[iter2]=f3[iter2+1][i];
            
            qsort(f4,NumIt2,sizeof(double),compare);
            
            //lowbound[i] = SMLE[10*i]-diff[i]-f4[above];
            lowbound[i] = SMLE[10*i]-diff[i]-f4[above]*sqrt(varF(n,n,freq,y,0.0,B,data,cc[10*i]*h,grid[10*i]));
            //upbound[i]= SMLE[10*i]-diff[i]-f4[below];
            upbound[i]= SMLE[10*i]-diff[i]-f4[below]*sqrt(varF(n,n,freq,y,0.0,B,data,cc[10*i]*h,grid[10*i]));
            
            if (F0(grid[i*10])<lowbound[i] || F0(grid[i*10])>upbound[i])
                percentage[i]++;
        }
        
        
        
        Rcout  << setw(10) << iter << setprecision(6) <<  setw(15) << F0(grid[500]) << setprecision(6) <<  setw(15) << lowbound[50] << setprecision(6) <<  setw(15) << upbound[50] << setw(10) << percentage[50] << std::endl;
    }
    
    StopTime  = clock();
    Time_bootstrap   = (double)(StopTime - StartTime)/(double)CLOCKS_PER_SEC;
    
    Rcout << std::endl << std::endl;
    Rcout << "The computations took    " << setprecision(10) << Time_bootstrap << "   seconds"  << std::endl;
    
    NumericMatrix out1 = NumericMatrix(njumps+1,2);
    
    for (i=0;i<=njumps;i++)
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
        for (i=0;i<=njumps;i++)
        {
            file0_ << setprecision(11) << setw(20) << data0[i];
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
    
    delete[] data, delete[] delta, delete[] tt, delete[] delta2, delete[] SMLE, delete[] SMLE2,
    delete[] F,  delete[] F2, delete[] cumw,
    delete[] cs, delete[] y, delete[] y2, delete[] jumploc,  delete[] data1, delete[] data2,
    delete[] p, delete[] p2, delete[] lowbound, delete[] upbound, delete[] cc; delete[] diff;
    
    for (i = 0;i<2;i++)
        delete[] freq[i], delete[] freq2[i];
    
    delete[] freq, delete[] freq2;
    
    for (iter2 = 0;iter2 < NumIt2;iter2++)
        delete[] f3[iter2];
    delete[] f3;
    
    return out;
}

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


double c0(double t)
{
    return pow(7.0/143,0.2)*pow(3,0.6)*pow(5/exp(1),0.4)*
    pow((exp(2) - exp(t))*(exp(t)-1),0.2);
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

double bdf(double A, double B, int m, double tt[], double pp[], double u, double h)
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

/*double bdf(double A, double B, int m, double t[], double p[], double u, double h)
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
    int    i;
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


void data_bootstrap(int n, int *m, double x[], double x2[], double data2[], int **freq, int delta[], int delta2[], int seed)
{
    int	i,j,k;
    SampleTime2 *obs;
    
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<int> dis(1,n);
    
    obs = new SampleTime2[n+1];
    
    for (k=1;k<=n;k++)
        freq[0][k]=freq[1][k]=0;
    
    for (k=1;k<=n;k++)
    {
        j=dis(gen);
        
        x2[k]=x[j];
        delta2[k]=delta[j];
    }
    
    for (i=0;i<n;i++)
    {
        obs[i].t=x2[i+1];
        obs[i].delta=delta2[i+1];
    }
    
    qsort(obs,n,sizeof(SampleTime2),CompareTime);
    
    for (i=1;i<=n;i++)
    {
        x2[i]=obs[i-1].t;
        delta2[i]=obs[i-1].delta;
    }
    
    
    x2[0]=0;
    
    j=0;
    
    for (i=1;i<=n;i++)
    {
        if (x2[i]>x2[i-1])
        {
            j++;
            data2[j]=x2[i];
            freq[delta2[i]][j]=1;
        }
        else
        {
            data2[j]=x2[i];
            freq[delta2[i]][j]++;
        }
    }
    
    *m=j;
    
    
    delete[] obs;
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

double F0(double x)
{
    return (1-exp(-x))/(1-exp(-2.0));
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
