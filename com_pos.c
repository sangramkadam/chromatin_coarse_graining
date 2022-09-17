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
int s[N];
int t,nn,clust_id[N];
int clust_size[N];
int size_dist[N];
int count,ncount;
int n_sys;
int n_ens;
int n_flip,n_bind,n_slide,ncoh;
int J,mu;
int nupdate;
int flag;
int seg;
int nframes;
double bin;
double sum;
double pp;
double xr,yr,zr,rr;
double cx,cy,cz,rg[1000];
double sigma,rc,rc2;
FILE *ip,*op;
char file[100];
sigma=1.0;
rc=1.1*sigma;
rc2=rc*rc;
n_sys=1000; //Number of trajectories
bin=0.05;
seg=10;     //segment size for coarse-graining (n_b)
nframes=61;

double x[N],y[N],z[N];
double rg_avg[N/seg];
int count_avg[N/seg];
n_ens=0;
count=0;


  for(i=0;i<N/seg;i++)
  {
  rg_avg[i]=0;
  count_avg[i]=0;
  }

  for(i=0;i<1000;i++) 
  rg[i]=0;

    
    sprintf(file,"pos.xyz");
    ip=fopen(file,"r");
        if(ip==NULL)
            printf("ERROR: File not found\n");
    ncount=0;
    while(ncount!=(N+2)*nframes) 
    {
    ncount=ncount+N+2;
    fscanf(ip,"%d\n",&nn);
    fscanf(ip,"%d %d\n",&nn,&t);
        for(i=0;i<N;i++)
        {
        fscanf(ip,"%d %lf %lf %lf\n",&s[i],&x[i],&y[i],&z[i]);
        }
        if(ncount>(N+2)*0) 
        {
        n_ens++;
        //Radius of gyration
        j=0;
        sprintf(file,"com%d.xyz",seg);
        op=fopen(file,"a");
        fprintf(op,"%d\n%d\n",N/seg,t);
        fclose(op);
        while(j<N-seg+1)
        {
            cx=0.0;
            cy=0.0;
            cz=0.0;
            for(i=j;i<j+seg;i++)
            {
            cx=cx+x[i];
            cy=cy+y[i];
            cz=cz+z[i];
            }
            cx=cx/(double)(seg);
            cy=cy/(double)(seg);
            cz=cz/(double)(seg);

            sprintf(file,"com%d.xyz",seg);
            op=fopen(file,"a");
            fprintf(op,"%d %lf %lf %lf\n",0,cx,cy,cz);
            fclose(op);
            sum=0.0;
            for(i=j;i<j+seg;i++)
            {
            sum=sum+Sqr(x[i]-cx)+Sqr(y[i]-cy)+Sqr(z[i]-cz);
            }
            sum=sum/(double)(seg);
            sum=sqrt(sum);
            rg[(int)(sum/bin)]++;
            rg_avg[j/seg]+=sum;
            count_avg[j/seg]++;
            count++;
        j=j+seg;
        }
        }
    }
  fclose(ip);

    //Average value of radius of gyration
    sprintf(file,"rg_avg_seg%d.dat",seg);
    op=fopen(file,"w");
    for(i=0;i<N/seg;i++)
    {
    fprintf(op,"%lf %lf\n",seg*(double)i+seg/2.0,rg_avg[i]/(double)count_avg[i]);
    }
    fclose(op);

    //Distribution of radius of gyration values
    sprintf(file,"rg_dist_seg%d.dat",seg);
    op=fopen(file,"w");
    for(i=0;i<1000;i++)
    {
    fprintf(op,"%lf %lf\n",bin/2.0+bin*(double)i,rg[i]/(double)count);
    }
    fclose(op);

return 0;
}
