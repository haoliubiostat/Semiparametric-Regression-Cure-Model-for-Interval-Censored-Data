#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* PURPOSE:    SUBROUTINE CONMIN MINIMIZES AN UNCONSTRAINED NONLINEAR */
/*             SCALAR VALUED FUNCTION OF A VECTOR VARIABLE X */
/*             EITHER BY THE BFGS VARIABLE METRIC ALGORITHM OR BY A */
/*             BEALE RESTARTED CONJUGATE GRADIENT ALGORITHM. */

/* USAGE:     
   int conmin(int n, double *x, double *f, double *g, int *ifun, int *iter,
           double eps, int mxfun, double *w, int mdim, double acc, int nmeth);
**/

/* PARAMETERS: N      THE NUMBER OF VARIABLES IN THE FUNCTION TO */
/*                    BE MINIMIZED. */
/*             X      THE VECTOR CONTAINING THE CURRENT ESTIMATE TO */
/*                    THE MINIMIZER. ON ENTRY TO CONMIN,X MUST CONTAIN */
/*                    AN INITIAL ESTIMATE SUPPLIED BY THE USER. */
/*                    ON EXITING,X WILL HOLD THE BEST ESTIMATE TO THE */
/*                    MINIMIZER OBTAINED BY CONMIN. X MUST BE DOUBLE */
/*                    PRECISIONED AND DIMENSIONED N. */
/*             F      ON EXITING FROM CONMIN,F WILL CONTAIN THE LOWEST */
/*                    VALUE OF THE OBJECT FUNCTION OBTAINED. */
/*                    F IS DOUBLE PRECISIONED. */
/*             G      ON EXITING FROM CONMIN,G WILL CONTAIN THE */
/*                    ELEMENTS OF THE GRADIENT OF F EVALUATED AT THE */
/*                    POINT CONTAINED IN X. G MUST BE DOUBLE */
/*                    PRECISIONED AND DIMENSIONED N. */
/*             IFUN   UPON EXITING FROM CONMIN,IFUN CONTAINS THE */
/*                    NUMBER OF TIMES THE FUNCTION AND GRADIENT */
/*                    HAVE BEEN EVALUATED. */
/*             ITER   UPON EXITING FROM CONMIN,ITER CONTAINS THE */
/*                    TOTAL NUMBER OF SEARCH DIRECTIONS CALCULATED */
/*                    TO OBTAIN THE CURRENT ESTIMATE TO THE MINIZER. */
/*             EPS    EPS IS THE USER SUPPLIED CONVERGENCE PARAMETER. */
/*                    CONVERGENCE OCCURS WHEN THE NORM OF THE GRADIENT */
/*                    IS LESS THAN OR EQUAL TO EPS TIMES THE MAXIMUM */
/*                    OF ONE AND THE NORM OF THE VECTOR X. EPS */
/*                    MUST BE DOUBLE PRECISIONED. */
/*             NFLAG  UPON EXITING FROM CONMIN,NFLAG STATES WHICH */
/*                    CONDITION CAUSED THE EXIT. */
/*                    IF NFLAG=0, THE ALGORITHM HAS CONVERGED. */
/*                    IF NFLAG=1, THE MAXIMUM NUMBER OF FUNCTION */
/*                       EVALUATIONS HAVE BEEN USED. */
/*                    IF NFLAG=2, THE LINEAR SEARCH HAS FAILED TO */
/*                       IMPROVE THE FUNCTION VALUE. THIS IS THE */
/*                       USUAL EXIT IF EITHER THE FUNCTION OR THE */
/*                       GRADIENT IS INCORRECTLY CODED. */
/*                    IF NFLAG=3, THE SEARCH VECTOR WAS NOT */
/*                       A DESCENT DIRECTION. THIS CAN ONLY BE CAUSED */
/*                       BY ROUNDOFF,AND MAY SUGGEST THAT THE */
/*                       CONVERGENCE CRITERION IS TOO STRICT. */
/*             MXFUN  MXFUN IS THE USER SUPPLIED MAXIMUM NUMBER OF */
/*                    FUNCTION AND GRADIENT CALLS THAT CONMIN WILL */
/*                    BE ALLOWED TO MAKE. */
/*             W      W IS A VECTOR OF WORKING STORAGE.IF NMETH=0, */
/*                    W MUST BE DIMENSIONED 5*N+2. IF NMETH=1, */
/*                    W MUST BE DIMENSIONED N*(N+7)/2. IN BOTH CASES, */
/*                    W MUST BE DOUBLE PRECISIONED. */
/*             IOUT   IOUT IS A USER  SUPPLIED OUTPUT PARAMETER. */
/*                    IF IOUT = 0, THERE IS NO PRINTED OUTPUT FROM */
/*                    CONMIN. IF IOUT > 0,THE VALUE OF F AND THE */
/*                    NORM OF THE GRADIENT SQUARED,AS WELL AS ITER */
/*                    AND IFUN,ARE WRITTEN EVERY IOUT ITERATIONS. */
/*             MDIM   MDIM IS THE USER SUPPLIED DIMENSION OF THE */
/*                    VECTOR W. IF NMETH=0,MDIM=5*N+2. IF NMETH=1, */
/*                    MDIM=N*(N+7)/2. */
/*             IDEV   IDEV IS THE USER SUPPLIED NUMBER OF THE OUTPUT */
/*                    DEVICE ON WHICH OUTPUT IS TO BE WRITTEN WHEN */
/*                    IOUT>0. */
/*             ACC    ACC IS A USER SUPPLIED ESTIMATE OF MACHINE */
/*                    ACCURACY. A LINEAR SEARCH IS UNSUCCESSFULLY */
/*                    TERMINATED WHEN THE NORM OF THE STEP SIZE */
/*                    BECOMES SMALLER THAN ACC. IN PRACTICE, */
/*                    ACC=10.D-20 HAS PROVED SATISFACTORY. ACC IS */
/*                    DOUBLE PRECISIONED. */
/*             NMETH  NMETH IS THE USER SUPPLIED VARIABLE WHICH */
/*                    CHOOSES THE METHOD OF OPTIMIZATION. IF */
/*                    NMETH=0,A CONJUGATE GRADIENT METHOD IS */
/*                    USED. IF NMETH=1, THE BFGS METHOD IS USED. */

/* REMARKS:    IN ADDITION TO THE SPECIFIED VALUES IN THE ABOVE */
/*             ARGUMENT LIST, THE USER MUST SUPPLY A SUBROUTINE */
/*             CALCFG WHICH CALCULATES THE FUNCTION AND GRADIENT AT */
/*             X AND PLACES THEM IN F AND G(1),...,G(N) RESPECTIVELY. */
/*             THE SUBROUTINE MUST HAVE THE FORM: */
/*                    SUBROUTINE CALCFG(N,X,F,G) */
/*                    DOUBLE PRECISION X(N),G(N),F */

/*** return nflag ***/

///extern void calcfg(int, double *, double *, double *);

int conmin(int n, double *x, double *fval, double *g, int *ifun, int *iter,
           double eps, int mxfun, double *w, int mdim, double acc, int nmeth,
           double (*calcfg)(int, const double *, double *, void *), void *param)
{
    int i, j,  ii, ij, ng, nx, nrd, nry, rsw, ngpi, ngpj, nxpi, nrst;
    double u1, u2, u3, u4, dg, ap, dp, at, fp, dg1, dal, gsq, xsq, funmin, step, rtst, alpha;
    int nrdpi, ncons, nrypi, ncons1, ncons2, ncalls;
    
    /* INITIALIZE ITER, IFUN. */
    
    *iter = *ifun = 0;
    
    /* SET PARAMETERS TO EXTRACT VECTORS FROM W. */
    /* W(I) HOLDS THE SEARCH VECTOR, W(NX+I) HOLDS THE BEST CURRENT ESTIMATE TO THE MINIMIZER,
        AND W(NG+I) HOLDS THE GRADIENT AT THE BEST CURRENT ESTIMATE. */
    
    nx = n;
    ng = nx + n;
    
    /* TEST WHICH METHOD IS BEING USED. */
    /* IF NMETH=0, W(NRY+I) HOLDS THE RESTART Y VECTOR AND */
    /* W(NRD+I) HOLDS THE RESTART SEARCH VECTOR. */
    
    if (nmeth == 0) {
        nry = ng + n;
        nrd = nry + n;
        ncons = n * 5;
        ncons1 = ncons + 0;
        ncons2 = ncons + 1;
    }
    else ncons = n * 3;   /* IF NMETH=1,W(NCONS+I) HOLDS THE APPROXIMATE INVERSE HESSIAN. */        
    
    /* CALCULATETHE FUNCTION AND GRADIENT AT THE INITIAL POINT AND INITIALIZE 
        nrst, WHICH IS USED TO DETERMINE WHETHER A BEALE RESTART IS BEING DONE. 
        NRST=N MEANS THAT THIS ITERATION IS A RESTART ITERATION. INITIALIZE 
        rsw, WHICH INDICATES THAT THE CURRENT SEARCH DIRECTION IS A GRADIENT DIRECTION. */
L20:
        
    fval[0]=calcfg(n, x, g, param);
    (*ifun)++;
    nrst = n;
    rsw = 1;
    
    /* CALCULATE THE INITIAL SEARCH DIRECTION , THE NORM OF X SQUARED, */
    /* AND THE NORM OF G SQUARED. DG1 IS THE CURRENT DIRECTIONAL */
    /* DERIVATIVE, WHILE XSQ AND GSQ ARE THE SQUARED NORMS. */
    
    dg1 = xsq = 0.0;
    for (i=0; i<n; i++) {
        w[i] = -g[i];
        xsq += x[i] * x[i];        
        dg1 -= g[i] * g[i];
    }
    gsq = -dg1;
    dg=dg1;
    
    /* TEST IF THE INITIAL POINT IS THE MINIMIZER. */
    
    if (gsq <= eps * eps * fmax(1.0, xsq)) return 0;
    
    /* BEGIN THE MAJOR ITERATION LOOP. NCALLS IS USED TO GUARANTEE THAT 
        AT LEAST TWO POINTS HAVE BEEN TRIED WHEN NMETH=0. FunMIN IS THE CURRENT FUNCTION VALUE. */
    
    for(;;) {
        
        funmin = fval[0];
        ncalls = *ifun;
                
        /* BEGIN LINEAR SEARCH. ALPHA IS THE STEPLENGTH. */
        /* SET ALPHA TO THE NONRESTART CONJUGATE GRADIENT ALPHA. */
        
        alpha = alpha * dg / dg1;
        
        /* IF NMETH=1 OR A RESTART HAS BEEN PERFORMED, SET ALPHA=1.0. */
        
        if (nrst == 1 || nmeth == 1) alpha = 1.0;
        
        /* IF A GRADIENT DIRECTION IS USED, SET ALPHA=1.0/DSQRT(GSQ), */
        /* WHICH SCALES THE INITIAL SEARCH VECTOR TO UNITY. */
        
        if (rsw==1) alpha = 1.0 / sqrt(gsq);
        
        /* THE LINEAR SEARCH FITS A CUBIC TO F AND DAL, THE FUNCTION AND ITS */
        /* DERIVATIVE AT ALPHA, AND TO FP AND DP,THE FUNCTION */
        /* AND DERIVATIVE AT THE PREVIOUS TRIAL POINT AP. */
        
        /* INITIALIZE AP ,FP,AND DP. */
        
        ap = 0.0;
        fp = funmin;
        dp = dg1;
        
        /* SAVE THE CURRENT DERIVATIVE TO SCALE THE NEXT SEARCH VECTOR. */
        
        dg = dg1;
        
        /* UPDATE THE ITERATION. */
        
        *iter++;
        
        /* CALCULATE THE CURRENT STEPLENGTH AND STORE THE CURRENT X AND G. */
        
        step = 0.0;
        for (i=0; i<n; i++) {
            step += w[i] * w[i];
            nxpi = nx + i;
            ngpi = ng + i;
            w[nxpi] = x[i];
            w[ngpi] = g[i];
        }
        step = sqrt(step);
        
        /* BEGIN THE LINEAR SEARCH ITERATION. */
        /* TEST FOR FAILURE OF THE LINEAR SEARCH. */
        
        for(;;) {
            if (alpha * step <= acc) {
                
                /* TEST IF DIRECTION IS A GRADIENT DIRECTION. */
                
                if (rsw==0) goto L20;
                return 2;
            }
            
            /* CALCULATE THE TRIAL POINT. */
            
            for (i=0; i<n; i++) {
                nxpi = nx + i;
                x[i] = w[nxpi] + alpha * w[i];
            }
            
            /* EVALUATE THE FUNCTION AT THE TRIAL POINT. */
            
            fval[0]=calcfg(n, x, g, param);
            
//            printf ("line search %5d %.5f %.5f %10.5f\n", ifun[0], x[0], x[1], fval[0]);
            
            /* TEST IF THE MAXIMUM NUMBER OF FUNCTION CALLS HAVE BEEN USED. */
            
            (*ifun)++;
            if (*ifun > mxfun) return 1;
                
            /* COMPUTE THE DERIVATIVE OF F AT ALPHA. */                
            
            dal = 0.0;
            for (i=0; i<n; i++) dal += g[i] * w[i];
            
            /* TEST WHETHER THE NEW POINT HAS A NEGATIVE SLOPE BUT A HIGHER */
            /* FUNCTION VALUE THAN ALPHA=0. IF THIS IS THE CASE,THE SEARCH */
            /* HAS PASSED THROUGH A LOCAL MAX AND IS HEADING FOR A DISTANT LOCAL MINIMUM. */
            
            if (fval[0] > funmin && dal < 0.0) {
                /* A RELATIVE MAX HAS BEEN PASSED. REDUCE ALPHA AND RESTART THE SEARCH. */
                
                alpha /= 3.0;
                ap = 0.0;
                fp = funmin;
                dp = dg;
            }
            else {
                
                /* IF NOT, TEST WHETHER THE STEPLENGTH CRITERIA HAVE BEEN MET. */
                
                if (!( fval[0] > funmin + alpha * 0.0001 * dg || fabs(dal/dg) > 0.9)) {
                    /* IF THEY HAVE BEEN MET, TEST IF TWO POINTS HAVE BEEN TRIED */
                    /* IF NMETH=0 AND IF THE TRUE LINE MINIMUM HAS NOT BEEN FOUND. */
                    
                    if (!(*ifun - ncalls <= 1 && fabs(dal / dg) > eps && nmeth == 0)) break;
                }
                
                /* A NEW POINT MUST BE TRIED. USE CUBIC INTERPOLATION TO FIND THE TRIAL POINT AT. */
                
                u1 = dp + dal - (fp - fval[0]) * 3.0 / (ap - alpha);
                u2 = u1 * u1 - dp * dal;
                if (u2 < 0.0) u2 = 0.0;
                u2 = sqrt(u2);
                at = alpha - (alpha - ap) * (dal + u2 - u1) / (dal - dp + u2 * 2.0);
                
                /* TEST WHETHER THE LINE MINIMUM HAS BEEN BRACKETED. */
                
                if (dal/dp > 0.0) {
                    
                    if (!(dal > 0.0 && 0.0 < at && at < fmin(ap,alpha) * .99)){
                        
                        /* TEST IF BOTH POINTS ARE LESS THAN THE MINIMUM AND THE TRIAL POINT */
                        /* IS SUFFICIENTLY LARGE. */
                        
                        if (!(dal <= 0.0 && at > fmax(ap,alpha) * 1.01)) {
                            /* IF THE TRIAL POINT IS TOO SMALL,DOUBLE THE LARGEST PRIOR POINT. */
                            
                            if (dal <= 0.0) at = fmax(ap,alpha) * 2.0;
                            
                            /* IF THE TRIAL POINT IS TOO LARGE, HALVE THE SMALLEST PRIOR POINT. */
                            
                            if (dal > 0.0) at = fmin(ap,alpha) / 2.0;                                
                        }
                    }
                }
                else{
                    /* THE MINIMUM HAS BEEN BRACKETED. TEST WHETHER THE TRIAL POINT LIES */
                    /* SUFFICIENTLY WITHIN THE BRACKETED INTERVAL. */
                    /* IF IT DOES NOT, CHOOSE AT AS THE MIDPOINT OF THE INTERVAL. */
                    
                    if (at < fmin(alpha,ap) * 1.01 || at > fmax(alpha,ap) * 0.99) at = (alpha + ap) / 2.0;                        
                }
                
                /* THE MINIMUM HAS NOT BEEN BRACKETED. TEST IF BOTH POINTS ARE GREATER THAN THE MINIMUM 
                    AND THE TRIAL POINT IS SUFFICIENTLY SMALLER THAN EITHER. */
                
                /* SET AP=ALPHA, ALPHA=AT,AND CONTINUE SEARCH. */            
                ap = alpha;
                fp = fval[0];
                dp = dal;
                alpha = at;
            }
        }
        
        /* THE LINE SEARCH HAS CONVERGED. TEST FOR CONVERGENCE OF THE ALGORITHM. */
        
        gsq = xsq = 0.0;
        for (i=0; i<n; i++) {
            gsq += g[i] * g[i];
            xsq += x[i] * x[i];
        }
        if (gsq <= eps * eps * fmax(1.0, xsq)) return 0;
                
        /* SEARCH CONTINUES. SET W(I)=ALPHA*W(I),THE FULL STEP VECTOR. */
        
        for (i=0; i<n; i++) w[i] = alpha * w[i];
        
        /* COMPUTE THE NEW SEARCH VECTOR. FIRST TEST WHETHER A */
        /* CONJUGATE GRADIENT OR A VARIABLE METRIC VECTOR IS USED. */
        
        if (nmeth == 0) {
            
            /* CONJUGATE GRADIENT UPDATE SECTION. */
            /* TEST IF A POWELL RESTART IS INDICATED. */
            
            rtst = 0.0;
            for (i=0; i<n; i++) {
                ngpi = ng + i;
                rtst += g[i] * w[ngpi];
            }
            if (fabs(rtst/gsq) > 0.2) nrst = n;
            
            /* IF A RESTART IS INDICATED, SAVE THE CURRENT D AND Y */
            /* AS THE BEALE RESTART VECTORS AND SAVE D'Y AND Y'Y */
            /* IN W(NCONS+1) AND W(NCONS+2). */
            
            if (nrst == n) {
                w[ncons1] = 0.0;
                w[ncons2] = 0.0;
                for (i=0; i<n; i++) {
                    nrdpi = nrd + i;
                    nrypi = nry + i;
                    ngpi = ng + i;
                    w[nrypi] = g[i] - w[ngpi];
                    w[nrdpi] = w[i];
                    w[ncons1] += w[nrypi] * w[nrypi];
                    w[ncons2] += w[i] * w[nrypi];
                }
            }
            
            /* CALCULATE THE RESTART HESSIAN TIMES THE CURRENT GRADIENT. */
            
            u1 = u2 = 0.0;
            for (i=0; i<n; i++) {
                nrdpi = nrd + i;
                nrypi = nry + i;
                u1 -= w[nrdpi] * g[i] / w[ncons1];
                u2 += w[nrdpi] * g[i] * 2.0 / w[ncons2] - w[nrypi] * g[i] / w[ncons1];
            }
            u3 = w[ncons2] / w[ncons1];
            for (i=0; i<n; i++) {
                nxpi = nx + i;
                nrdpi = nrd + i;
                nrypi = nry + i;
                w[nxpi] = -u3 * g[i] - u1 * w[nrypi] - u2 * w[nrdpi];
            }
            
            /* IF THIS IS A RESTART ITERATION,W(NX+I) CONTAINS THE NEW SEARCH VECTOR. */
            
            if (nrst != n) {
                /* NOT A RESTART ITERATION. CALCULATE THE RESTART HESSIAN */
                /* TIMES THE CURRENT Y. */
                
                u1 = u2 = u3 = u4 = 0.0;
                for (i=0; i<n; i++) {
                    ngpi = ng + i;
                    nrdpi = nrd + i;
                    nrypi = nry + i;
                    u1 -= (g[i] - w[ngpi]) * w[nrdpi] / w[ncons1];
                    u2 -= (g[i] - w[ngpi]) * w[nrypi] / w[ncons1] + w[nrdpi] * 2.0 * (g[i] - w[ngpi]) / w[ncons2];
                    u3 += w[i] * (g[i] - w[ngpi]);
                }
                step = 0.0;
                for (i=0; i<n; i++) {
                    ngpi = ng + i;
                    nrdpi = nrd + i;
                    nrypi = nry + i;
                    step = w[ncons2] / w[ncons1] * (g[i] - w[ngpi]) + u1 * w[nrypi] + u2 * w[nrdpi];
                    u4 += step * (g[i] - w[ngpi]);
                    w[ngpi] = step;
                }
                
                /* CALCULATE THE DOUBLY UPDATED HESSIAN TIMES THE CURRENT */
                /* GRADIENT TO OBTAIN THE SEARCH VECTOR. */
                
                u1 = u2 = 0.0;
                for (i=0; i<n; i++) {
                    u1 -= w[i] * g[i] / u3;
                    ngpi = ng + i;
                    u2 += (u4 / u3 + 1.0) * w[i] * g[i] / u3 - w[ngpi] * g[i] / u3;
                }
                for (i=0; i<n; i++) {
                    ngpi = ng + i;
                    nxpi = nx + i;
                    w[nxpi] -= (u1 * w[ngpi] + u2 * w[i]);
                } 
            }
            
            /* CALCULATE THE DERIVATIVE ALONG THE NEW SEARCH VECTOR. */
            
            dg1 = 0.0;
            for (i=0; i<n; i++) {
                nxpi = nx + i;
                w[i] = w[nxpi];
                dg1 += w[i] * g[i];
            }
            
//            printf ("after line search %5d %.5f %.5f %10.5f\n", ifun[0], x[0], x[1], fval[0]);
//            printf("  1st dg1=%.10f\n", dg1);
            
            /* IF THE NEW DIRECTION IS NOT A DESCENT DIRECTION,STOP. */
            /* ROUNDOFF HAS PRODUCED A BAD DIRECTION. */
            
            if (dg1 > 0.0) return 3;
            
            /* UPDATE NRST TO ASSURE AT LEAST ONE RESTART EVERY N ITERATIONS. */
            
            if (nrst == n) nrst = 0;
            ++nrst;
        }
        else {
            /* A VARIABLE METRIC ALGORITM IS BEING USED. CALCULATE Y AND D'Y. */
            
            u1 = 0.0;
            for (i=0; i<n; i++) {
                ngpi = ng + i;
                w[ngpi] = g[i] - w[ngpi];
                u1 += w[i] * w[ngpi];
            }
            
            /* IF RSW=.TRUE.,SET UP THE INITIAL SCALED APPROXIMATE HESSIAN. */
            
            if (rsw==1) {
                /* CALCULATE Y'Y. */
                
                u2 = 0.0;
                for (i=0; i<n; i++) {
                    ngpi = ng + i;
                    u2 += w[ngpi] * w[ngpi];
                }
                
                /* CALCULATE THE INITIAL HESSIAN AS H=(P'Y/Y'Y)*I */
                /* AND THE INITIAL U2=Y'HY AND W(NX+I)=HY. */
                
                ij = 0;
                u3 = u1 / u2;
                for (i=0; i<n; i++) {
                    for (j=i; j<n; j++) {
                        ncons1 = ncons + ij;
                        w[ncons1] = 0.0;
                        if (i == j) w[ncons1] = u3;
                        ++ij;
                    }
                    nxpi = nx + i;
                    ngpi = ng + i;
                    w[nxpi] = u3 * w[ngpi];
                }
                u2 *= u3 ;
            }
            else {
                /* CALCULATE W(NX+I)=HY AND U2=Y'HY. */
                
                u2 = 0.0;
                for (i=0; i<n; i++) {
                    u3 = 0.0;
                    ij = i;
                    if (i != 0) {
                        for (j = 0; j < i; j++) {
                            ngpj = ng + j;
                            ncons1 = ncons + ij;
                            u3 += w[ncons1] * w[ngpj];
                            ij += (n-j-1);
                        }
                    }
                    for (j=i; j<n; j++) {
                        ncons1 = ncons + ij;
                        ngpj = ng + j;
                        u3 += w[ncons1] * w[ngpj];
                        ++ij;
                    }
                    ngpi = ng + i;
                    u2 += (u3 * w[ngpi]);
                    nxpi = nx + i;
                    w[nxpi] = u3;
                }            
            }
            
            /* CALCULATE THE UPDATED APPROXIMATE HESSIAN. */
            
            u4 = u2 / u1 + 1.0;
            for (i=0; i<n; i++) {
                nxpi = nx + i;
                ngpi = ng + i;
                w[ngpi] = u4 * w[i] - w[nxpi];
            }
            
            ij = 0;
            for (i=0; i<n; i++) {
                nxpi = nx + i;
                u3 = w[i] / u1;
                u4 = w[nxpi] / u1;
                for (j=i; j<n; j++) {
                    ncons1 = ncons + ij;
                    ngpj = ng + j;
                    w[ncons1] = w[ncons1] + u3 * w[ngpj] - u4 * w[j];
                    ++ij;
                }
            }
            
            /* CALCULATE THE NEW SEARCH DIRECTION W(I)=-HG AND ITS DERIVATIVE. */
            
            dg1 = 0.0;
            for (i=0; i<n; i++) {
                u3 = 0.0;
                ij = i;
                if (i != 0){
                    for (j = 0; j < i; j++) {
                        ncons1 = ncons + ij;
                        u3 -= w[ncons1] * g[j];
                        ij += (n-j-1);
                    } 
                }
                for (j = i; j <n; j++)  {
                    ncons1 = ncons + ij;
                    u3 -= w[ncons1] * g[j];
                    ++ij;
                }
                dg1 += u3 * g[i];
                w[i] = u3;
            }

            /* TEST FOR A DOWNHILL DIRECTION. */
            /* ROUNDOFF HAS PRODUCED A BAD DIRECTION. */            
            if (dg1 > 0.0) return 3;
        }
        rsw =0;
    }
} 

