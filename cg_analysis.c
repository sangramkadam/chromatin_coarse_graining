#include<stdio.h>
#include<stdlib.h>
#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))	//Nearest integer
#include<math.h>
#define Sqr(x)     ((x) * (x))

int main()
{
int i,j,k;
int ii;
int n;
int N;
N=400;
int t,nn,clust_id[N];
int clust_size[N];
int size_dist[N];
int count,ncount;
int count_bond;
int count_angle;
int n_sys;
int n_ens;
int n_flip,n_bind,n_slide,ncoh;
int J,mu;
int nupdate;
int flag;
int seg;
int pp;
double dist_bond[500];
double dist_angle[500];
double bin_bond;
double bin_angle;
double sum_bond;
double sum_spring;
double xr,yr,zr,rr;
double sigma,rc,rc2,theta;
FILE *ip,*op;
char file[100];
sigma=1.0;
rc=1.5*sigma;
rc2=rc*rc;
seg=10;
int dot_count[N/seg];
double dot_bond[N/seg];
double x[N/seg],y[N/seg],z[N/seg];
double bond[N/seg];
double bond2[N/seg];
double tx[N/seg];
double ty[N/seg];
double tz[N/seg];
int nframes=61;
n_ens=0;
count=0;
count_bond=0;
count_angle=0;


  for(i=0;i<N/seg;i++)
  {
  bond[i]=0.0;
  bond2[i]=0.0;
  dot_bond[i]=0.0;
  dot_count[i]=0;
  }
  for(i=0;i<500;i++)
  {
      dist_bond[i]=0.0;
      dist_angle[i]=0.0;
  }
  bin_bond=0.1;
  sum_bond=0.0;
  sum_spring=0.0;
  bin_angle=0.05;

   //Bond length and bond angle
    sprintf(file,"com%d.xyz",seg);
    ip=fopen(file,"r");
    if(ip==NULL)printf("ERROR: File not found\n");
    ncount=0;
    while(ncount!=(N/seg+2)*nframes) 
    {
    ncount=ncount+N/seg+2;
    fscanf(ip,"%d\n",&nn);
    fscanf(ip,"%d\n",&t);
        for(i=0;i<N/seg;i++)
        {
        fscanf(ip,"%d %lf %lf %lf\n",&pp,&x[i],&y[i],&z[i]);
        }
        if(ncount>(N/seg+2)*0) 
        {
            n_ens++;
            for(i=0;i<N/seg-1;i++)
            {
            xr=x[i+1]-x[i];
            yr=y[i+1]-y[i];
            zr=z[i+1]-z[i];
            rr=Sqr(xr)+Sqr(yr)+Sqr(zr);
            rr=sqrt(rr);
            bond[i]+=rr;
            count_bond++;
            sum_bond+=rr;
            dist_bond[(int)(rr/bin_bond)]+=1.0;
            tx[i]=xr/rr;
            ty[i]=yr/rr;
            tz[i]=zr/rr;
            }
            for(i=0;i<N/seg-2;i++)
            {
            theta=-(tx[i]*tx[i+1]+ty[i]*ty[i+1]+tz[i]*tz[i+1]);
            theta=acos(theta);
            dist_angle[(int)(theta/bin_angle)]+=1.0;
            count_angle++;
            }
            for(i=0;i<N/seg-1;i++)
            {
                for(j=i;j<N/seg-1;j++)
                {
                dot_bond[j-i]+=tx[i]*tx[j]+ty[i]*ty[j]+tz[i]*tz[j];
                dot_count[j-i]++;
                }
            }
        }
    }
  fclose(ip);
  sprintf(file,"bond_length_seg%d.dat",seg);
  op=fopen(file,"w");
  for(i=0;i<N/seg-1;i++)
  {
  bond[i]=bond[i]/(double)n_ens;
  fprintf(op,"%d %lf\n",i*(int)seg+(int)seg,bond[i]);
  }
  fclose(op);

  sprintf(file,"dist_bond_length_seg%d.dat",seg);
  op=fopen(file,"w");
  for(i=0;i<500;i++)
  {
  dist_bond[i]=dist_bond[i]/(double)(count_bond);
  fprintf(op,"%lf %lf\n",(double)i*bin_bond+bin_bond/2.0,dist_bond[i]/bin_bond);
  }
  fclose(op);
  sprintf(file,"avg_bond_length.dat");
  op=fopen(file,"a");
  fprintf(op,"%d %lf\n",seg,sum_bond/(double)(count_bond));
  fclose(op);

  sprintf(file,"dist_angle_seg%d.dat",seg);
  op=fopen(file,"w");
  for(i=0;i<500;i++)
  fprintf(op,"%lf %lf\n",(double)i*bin_angle+bin_angle/2.0,dist_angle[i]/((double)(count_angle)*bin_angle));
  fclose(op);

// Spring constant calculation
n_ens=0;
    sprintf(file,"com%d.xyz",seg);
    ip=fopen(file,"r");
    if(ip==NULL)printf("ERROR: File not found\n");
    ncount=0;
    while(ncount!=(N/seg+2)*nframes) 
    {
    ncount=ncount+N/seg+2;
    fscanf(ip,"%d\n%d\n",&nn,&t);
        for(i=0;i<N/seg;i++)
        {
        fscanf(ip,"%d %lf %lf %lf\n",&pp,&x[i],&y[i],&z[i]);
        }
        if(ncount>(N/seg+2)*0) 
        {
        n_ens++;
            for(i=0;i<N/seg-1;i++)
            {
            xr=x[i+1]-x[i];
            yr=y[i+1]-y[i];
            zr=z[i+1]-z[i];
            rr=Sqr(xr)+Sqr(yr)+Sqr(zr);
            rr=sqrt(rr);
            bond2[i]+=Sqr(rr-bond[i]);
            }
        }
    }
  fclose(ip);

  sprintf(file,"spring_const_seg%d.dat",seg);
  op=fopen(file,"w");
  for(i=0;i<N/seg-1;i++)
  {
  bond2[i]=bond2[i]/(double)n_ens;
  fprintf(op,"%d %lf\n",i*(int)seg+(int)seg,1.0/bond2[i]);
  sum_spring+=1.0/bond2[i];
  }
  fclose(op);

  sprintf(file,"avg_spring.dat");
  op=fopen(file,"a");
  fprintf(op,"%d %lf\n",seg,sum_spring/(double)(N/seg));
  fclose(op);
return 0;
}
