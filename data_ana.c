/************************************************************/
/*********  Interval Censoring data analysis   ***********/
/*********  by Hao Liu on 6/7/07                  **********/
/************************************************************/

/*** accepted data file: first p colmuns are covariates, last 2 columns are observed time intervals, 9999=censored ***/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#define MAX_ITER 10000
#define TOL 1.0e-6
#define sim_num 300

/** sim_num is the bootstrap replicates **/

#define machine_p 1.0e-15
#define myPHO 0.5
#define myPI 0.01
#define Sigma 20.0

struct fun_params {
    double *dij; 
    double *zij; 
    double *pj;
    int nn;
    int np;
    int mk;
};

double fdf_eplike(int, const double *, double *, void *);
void myxij(double *, const double *, const double *, void *);
void deveplike(double *, double *, const double *, double *, double *, void *);

double maxp(double *, double *, double, const double *, double *, void *);
double normf(double *, double *, double, double, const double *, double *, void *);

extern int conmin(int, double *, double *, double *, int *, int *, double, int, double *, int, double, int,
                  double (*)(int, const double *, double *, void *), void *);

int main(int argc, char *argv[])
{
    FILE *fpi, *fpo1, *fpo2;
    if(argc<=4) {
        printf("need: 1) data file name, 2) sample size n, 3) number of covariate p, and 4) the output file1 5) files2...\n");
        return 0;
    }
    fpi=fopen(argv[1], "r");
    fpo1=fopen(argv[4], "w");
    fpo2=fopen(argv[5], "w");
    
    int nn=0, np=0, i, j;
    char *cs=malloc(1000*sizeof(char));
    
    sscanf(argv[2], "%i", &nn);
    sscanf(argv[3], "%i", &np);
    
    if (np>nn) {
        printf("sample size should be larger than the number of covariates.");
        return 0;
    }
    
    double *li= malloc(nn*sizeof(double)), *lia= malloc(nn*sizeof(double)), *li0=malloc(nn*sizeof(double));
    double *ri= malloc(nn*sizeof(double)), *ria= malloc(nn*sizeof(double)), *ri0=malloc(nn*sizeof(double));
    double *zij= malloc(nn*np*sizeof(double)), *zij0= malloc(nn*np*sizeof(double)), *dij=malloc(nn*nn*sizeof(double));
    double *q=malloc(nn*sizeof(double)), *p=malloc(nn*sizeof(double)), *gam= malloc(nn*nn*sizeof(double));
    
    double *pj=calloc(nn, sizeof(double)), *mypj=calloc(nn, sizeof(double));
        
    double alpha0=1.0, theta0=1.0, nvold; 
    register double *tau= calloc(nn, sizeof(double));
    double past_time, now_time, mysize, tmp00;

    double curt, wt, u, likval1, likval2, likval;
    
	int j1, k, mk, ss, iter, status, iter2, cure;
    int *bsi=malloc(nn*sizeof(int)), *bsi2=malloc(nn*sizeof(int));
    
    double *ww=calloc(2*(np+1)*(np+8), sizeof(double)), 
        *xx=calloc(2*(np+1), sizeof(double)), *gg=calloc(2*(np+1), sizeof(double)), fval;
    double acc=1e-11, eps=1e-5;;
    int nflag=-1, ifun=-1, ifunter=-1;
    
    struct fun_params *mydata=malloc(sizeof(struct fun_params));

    const gsl_rng_type *Ty;
    gsl_rng *rn;
    
    gsl_rng_env_setup();
    Ty = gsl_rng_default;
    rn = gsl_rng_alloc (Ty);
    gsl_rng_set(rn, 100);
    
    /** read the data, first line is the hearder and thus removed ***/
    fgets(cs, 1000, fpi);
    
    for (i = 0; i<nn; i++) {
        for (j=0;j<np;j++) fscanf(fpi, "%lf,", zij0+np*i+j);
        fscanf(fpi, "%lf,%lf", li0+i, ri0+i);
        bsi[i]=i;
    }
    fclose(fpi);
    mydata->nn = nn;
    mydata->np = np;
    
    /** sim_num is the bootstrap replicates **/
            
    for (ss=0;ss<sim_num;ss++){
        
        /** select the bootstrap sample ***/        
        /*** random sample with replacement **/

        gsl_ran_sample (rn, bsi2, nn, bsi, nn, sizeof(int));
        for(i=0;i<nn;i++) {
            for (j=0;j<np;j++) zij[np*i+j] = zij0[np*bsi2[i]+j];
            lia[i]=li0[bsi2[i]];
            li[i]=li0[bsi2[i]];
            ria[i]=ri0[bsi2[i]];
            ri[i]=ri0[bsi2[i]];
        }
        
        mydata->zij = zij;        
        gsl_sort(ria,1,nn);
        gsl_sort(lia,1,nn);

        /*** generate p, q data ***/
        p[0]=ria[0];
        i=j=k=0;
        do {
            while (lia[i]<=p[k] && i<nn) i++;
            q[k]=lia[i-1];
            k++;
            while (ria[j]<lia[i] && j<nn) j++;
            if (j<nn) p[k]=ria[j];
            else p[k]=ria[j-1];
        } while (i<nn && j<nn);
                
        mk=k-1;
        mydata->mk=mk;
        for (i=0; i<nn; i++)
            for (j = 0; j <= mk; j++)
                dij[i*(mk+1)+j]=(li[i]<=q[j] && p[j]<=ri[i]) ? 1.0: 0.0;
        mydata->dij= dij;
        
        /** Initial values **/
        
        for (j=0; j<mk; j++) pj[j] = 1.0/mk;
        (mydata->pj)=pj;
        
        for (j=0; j<=np; j++) xx[j] = 0.0;
        
        nflag=2;
        fval=9999.0;
        iter2=0;
        do {
            iter2++;
            likval2=fval;
            
            for (j = 0; j < mk; j++) tau[j]=1.0/(mydata->pj)[j];
            nvold = 0.1;

            /*** calculate expected Xij ***/
            myxij(gam, xx, pj, mydata);

            /*** Primal-dual interior-point method ***/
                        
            nvold = maxp(pj, tau, nvold, xx, gam, mydata);
            
            /** minimize to get theta **/
            
            if (nflag!=0) for (j=0; j<=np; j++) xx[j] = 0.1;            
            nflag= conmin(np+1, xx, &fval, gg, &ifun, &ifunter, 0.001, 300, ww, 400, 1e-18, 1, fdf_eplike, mydata);

            likval1=fval;

        } while(fabs(likval1-likval2)>0.0001 & iter2 < MAX_ITER);
        
        now_time = (double) clock ();

        /*** outputs ****/
         
        for(j=0; j<=np; j++) fprintf(fpo1, "%.10f, ", xx[j]); 
        fprintf(fpo1, "%.10f, %.2f \n", fval, now_time/(double) CLOCKS_PER_SEC);
        for(j=0; j<mk; j++) fprintf(fpo2, "%.10f,", (mydata->pj)[j]);
        fprintf(fpo2, "%i\n", mk);
    }
    
    fclose(fpo1);
    fclose(fpo2);
    
    free(li);
    free(ri);    
    free(zij);

    free(li0);
    free(ri0);
    free(zij0);
    
    free(lia);
    free(ria);
    free(p);
    free(q);
    free(mypj);
    
    free(gam);    
    free(tau);
    
    free(xx);
    free(ww);
    free(gg);
    
    free(mydata);
    
    return 0;
}

/*** the expected likelihood, 0=alpha, 1=theta and the derivative ***/

double fdf_eplike(int np1, const double *theta, double *df, void * par)
{
    struct fun_params * param = (struct fun_params *) par;
    int i,j,j1, m=(param->mk), nn=(param->nn), tmpi;
    
    double alp0, sz, likval, msumval, tmpz, tmpz1, tmpa, tmpc, tmpd, tmp0, dfz1, dfz2;
    double *xij, *cij, *pps, pp;
    
    xij=calloc(m+1, sizeof(double));
    pps=calloc(m+1, sizeof(double));
    cij=calloc(m+1, sizeof(double));

    alp0=theta[0];
        
    for(j=0;j<m;j++) {
        pps[j]=0.0;
        for(j1=0;j1<j;j1++) pps[j] += ((param->pj)[j1]) ;
    }

    for(j=0;j<np1;j++) df[j]=0.0;
    
    likval=0.0;
    for(i=0;i<nn;i++) {
        tmpz1=0.0;
        for(j=1;j<np1;j++) tmpz1 += theta[j] * (param->zij)[i*(np1-1)+j-1];
        tmpz=exp(alp0 + tmpz1);
                
        tmpc=0.0;
        for(j=0;j<m;j++) {
            pp = (param->pj)[j];
            xij[j] = ((param->dij)[i*(m+1)+j]) * exp( -tmpz* pps[j]) *(1.0- exp(-pp*tmpz));
            tmpc += xij[j];
        }
        xij[m]= (param -> dij)[i*(m+1)+m] * exp(-tmpz);
        tmpc += xij[m];
        
        if (tmpc == 0.0) {
            for(j=0;j<np1;j++) df[j]=-9999999.0;
            return 9999999.0;
        }
        
        for(j=0;j<=m;j++) {
            cij[j]=0.0;
            for(j1=j+1; j1<m; j1++) cij[j] += xij[j1];
            xij[j] /=tmpc;
        }
        
        msumval=0.0;
        for(j=0;j<m;j++) {
            cij[j] /=tmpc;
            pp= (param -> pj)[j];
            likval += ( xij[j]* log(1.0- exp(- pp* tmpz )) - cij[j] * pp * tmpz);
            msumval += (( xij[j] *exp(-pp* tmpz) /(1.0- exp(-pp* tmpz)) - cij[j] )*pp);
        }
        df[0] += ((msumval-xij[m])*tmpz);
        
        for(j=1;j<np1;j++) df[j] += ((msumval-xij[m])*(param->zij)[i*(np1-1)+j-1]*tmpz);
        likval -= xij[m]*tmpz;
    }
    
    for(j=0;j<np1;j++) df[j] = - df[j]/((double) nn);

    free(cij);
    free(xij);
    free(pps);

    return -likval/((double) nn);
}


/*** Calculate the updated Xij ***/

void myxij(double *lamij, const double *the, const double *mypj, void * parm)
{   
    struct fun_params * para = (struct fun_params *) parm;

    double alp0, tmp0, tmpc, pp, tmpz, tmpz1;    
    int i,j,j1, tmpi, m=(para->mk), nn=(para->nn), np=(para->np);
    double *gj=calloc(m+1, sizeof(double)), *xij=calloc(m+1, sizeof(double));
    double *dij=(para->dij), *zij=(para->zij);
    
    alp0=the[0];
    
    for(j=0;j<m;j++) {
        gj[j]=0.0;
        for(j1=0;j1<j;j1++) gj[j] += mypj[j1];
    }    
    for(i=0;i<nn;i++) {
        tmpz1=0.0;
        for(j=0;j<np;j++) tmpz1 += the[j+1] * zij[j+i*np];
        tmpz=exp(alp0+ tmpz1);
        
        tmpc=0.0;
        for(j=0;j<m;j++) {
            pp = mypj[j];
            xij[j] = dij[i*(m+1)+j] * exp( -tmpz* gj[j]) *(1.0- exp(-pp*tmpz));
            tmpc += xij[j];
        }
        xij[m]= dij[i*(m+1)+m] * exp(-tmpz);
        tmpc += xij[m];
        
        for(j=0;j<m;j++)  lamij[i*m+j] = xij[j] / tmpc;
    }
    
    free(xij);
    free(gj);
}

/*** Calculate the first derivative and second derivative wrt p of the expected liklihood ***/
void deveplike(double *df, double *ddf, const double *the, double *gam, double * myp, void *parm)
{
    struct fun_params * para = (struct fun_params *) parm;
    double alp0, tmp1, tmpz1, tmp2, tmpa, pp;
    int i,j,k,j1, m=(para->mk), nn=(para->nn), np=(para->np);
    double *zij=(para->zij);
    
    alp0=the[0];
    for(j=0;j<m;j++) {
        
        df[j] = ddf[j] = 0.0;
        pp = myp[j];
        
        for(i=0;i<nn;i++) {      
 
            tmpz1=0.0;
            for(k=0;k<np;k++) tmpz1 += the[k+1] * zij[k+i*np];
            tmp1=exp(alp0+ tmpz1);
            tmp2=exp(tmp1 * pp);
            
            tmpa =0.0;
            for(j1=j+1; j1<m; j1++) tmpa += gam[i*m+j1];
            df[j] += tmp1*( gam[i*m+j] /( tmp2- 1.0) - tmpa);
            ddf[j] -= (gam[i*m+j]*tmp1*tmp1*tmp2 / (tmp2 - 1.0) / (tmp2 - 1.0));
        }
    }
}


/*** NormF ***/

double normf(double *p, double *tau, double nv, double eta, const double *the, double *gam, void *parm)
{
    struct fun_params * para = (struct fun_params *) parm;
	double *fj, *ffj, sum1, tmp1, tmp2;
	int j, m=(para->mk);
    
    fj=calloc(m, sizeof(double));
    ffj=calloc(m, sizeof(double));
    
    deveplike(fj, ffj, the, gam, p, parm);
    
    sum1=0.0;
    for(j=0;j<m;j++) {
        tmp1 = (fj[j] + tau[j] - nv);
        tmp2 = p[j] * tau[j] - 1.0/eta;
        sum1 += (tmp1*tmp1 + tmp2*tmp2);
    }

    free(fj);
    free(ffj);
    return sqrt(sum1);
}

/**** Maximize for p ****/

double maxp(double *pj, double *tau, double nvold, const double *the, double *gam, void *parm)
{
    struct fun_params * para = (struct fun_params *) parm;
    
	int j, iter, m=(para->mk), nn=(para->nn), np=(para->np);
    double eta, nv, curnv, steps, sum1, sum2, tmp1, tmp2;
    double *fj, *ffj, *delp, *curp, *taunew, *curtau;

    fj = calloc(m, sizeof(double)); 
    ffj = calloc(m, sizeof(double)); 
    delp = calloc(m, sizeof(double)); 
    curp= calloc(m, sizeof(double));
    taunew = calloc(m, sizeof(double));
    curtau = calloc(m, sizeof(double));
    
    iter=0;
    do {
        iter++;
        
        // calculate the df=fj and diag(ddf)=ffj
        
        for(j=0;j<m;j++) eta += pj[j]*tau[j];
        eta = m*Sigma/eta;
        deveplike(fj, ffj, the, gam, pj, parm);
        
        sum1 = sum2 = 0.0;
        for(j=0;j<m;j++) {
            sum1 += ((1.0/eta + pj[j] * fj[j]) / (tau[j] - pj[j] * ffj[j]));
            sum2 += (pj[j] / (tau[j] - pj[j] * ffj[j]));
        }
        nv=sum1/sum2;
        
        // backtracking: find the largest positive step = steps
        
        steps=1.0;
        for(j=0;j<m;j++) { 
            delp[j] =  (pj[j] * (fj[j]-nv) + 1.0/eta) / (tau[j] - pj[j] * ffj[j]);
            taunew[j] =   nv - ffj[j] * delp[j] - fj[j];
            tmp1 = taunew[j] - tau[j];
            if (tmp1 <0) steps = fmin(steps, -tau[j]/tmp1);
        }
        steps *= 0.99;
        
        // backtracking: to keep the minimal p > 0;
        tmp2 = 1.0;
        for(j=0;j<m;j++) tmp2=fmin(tmp2,  pj[j] + steps * delp[j]);
        while(tmp2 <= 0) {
            steps *= myPHO;
            tmp2 = 1.0;
            for(j=0;j<m;j++) tmp2=fmin(tmp2, pj[j] + steps * delp[j]);
        }

        // backtracking: 
        
        for(j=0;j<m;j++) {
            curp[j] = pj[j] + steps * delp[j];
            curtau[j] = tau[j] + steps * (taunew[j]-tau[j]);
        }
        curnv=nvold + steps * (nv-nvold);
        tmp2 = normf(curp, curtau, curnv, eta, the, gam, parm);
        tmp1 = normf(pj, tau, nvold, eta, the, gam, parm);
        
        while(tmp2 > (1.0-myPI*steps)*tmp1) {
            steps *= myPHO;
            for(j=0;j<m;j++) {
                curp[j] = pj[j] + steps * delp[j];
                curtau[j] = tau[j] + steps * (taunew[j]-tau[j]);
            }
            curnv=nvold + steps * (nv-nvold);
            tmp2=normf(curp, curtau, curnv, eta, the, gam, parm);
        }
        
        deveplike(fj, ffj, the, gam, curp, parm);

        sum1=sum2=0.0;
        for(j=0;j<m;j++) {
            tmp1 = (fj[j] + curtau[j] - curnv);
            sum1 += (tmp1*tmp1);
            sum2 += curp[j] * curtau[j];
            pj[j] = curp[j];
            tau[j] = curtau[j];
        }
        
        nvold=curnv;
//         printf("iter=%i: sum1=%.20f sum2=%f \n", iter, sqrt(sum1), sum2/m);
        
    } while ( (sqrt(sum1)> TOL)& iter < MAX_ITER*10 );
    
    free(fj);
    free(ffj);
    free(delp);
    free(curtau);
    free(taunew);
    free(curp);
    
    return nvold;
}

