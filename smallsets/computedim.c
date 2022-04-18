#include "get_polyfuns.h"

void readmats(slong *multmatrix, slong *mapmat, bool_mat_t temp, slong *multi, slong *comulti, slong total, slong uni, slong couni, slong maxmul, int num) {
    
    char numers[30], filename[30]; 
    FILE *fp;
    slong k0, k1, nmax;  

    sprintf(filename, "multmatrix%d.dat", num);
    printf("multiplicity data = %s \n", filename);
  
    //fp = fopen("multmatrix.dat","r");
    fp = fopen(filename,"r");
    for (k0=0; k0< uni; k0++) {
        fscanf(fp, "%s", numers); 
        nmax = atol(numers);  
        multmatrix[k0*(maxmul+1)] = nmax; 
        for (k1=1; k1 < nmax-1; k1 ++) {
            fscanf(fp, "%s", numers); 
            multmatrix[k0*(maxmul+1)+k1] = atol(numers); 
            multi[atol(numers)-1] = k0;
        }
    }
    fclose(fp);

    sprintf(filename, "comultmatrix%d.dat", num);
    printf("co-multiplicity data = %s \n", filename);

    //fp = fopen("comultmatrix.dat","r");
    fp = fopen(filename,"r");
    for (k0=0; k0 < couni; k0++) {
        fscanf(fp, "%s", numers); 
        nmax = atol(numers);  
        for (k1=1; k1 < nmax-1; k1 ++) {
            fscanf(fp, "%s", numers); 
            comulti[atol(numers)-1] = k0;
        }
    }
    fclose(fp);
   
    fp = fopen("multi.c.dat","w"); 
    for (k0=0; k0<total; k0++) {  flint_fprintf(fp,"%wd \t",multi[k0]); }
    fclose(fp);
 
    fp = fopen("comulti.c.dat","w"); 
    for (k0=0; k0<total; k0++) {  flint_fprintf(fp,"%wd \t",comulti[k0]); }
    fclose(fp);

    sprintf(filename, "smallmatrices%d.dat", num);
    printf("small matrices file = %s \n", filename);

    fp = fopen(filename,"r");
    //fp = fopen("smallmatrices10.dat","r");
    for (k0=0; k0 < total; k0++) {
        for (k1=0; k1 < 4; k1 ++) {
            fscanf(fp, "%s", numers); 
            mapmat[k0*4+k1] = atol(numers); 
        }
    }
    fclose(fp);

    sprintf(filename, "minimatrix%d.dat", num);
    printf("minimatrix file = %s \n", filename);

    fp = fopen(filename,"r");
//    fp = fopen("minimatrix.dat","r"); 
    for (k0 = 0; k0<uni; k0++) {
        for (k1 = 0; k1<couni; k1++) {
            fscanf(fp, "%s", numers); 
            bool_mat_set_entry(temp,k0,k1,atoi(numers));
        }
    }
    fclose(fp);
}

int main(int argc, char *argv[]) {

    slong print_dig = 20;
    slong prec, uni, couni, total, maxmul, digits;
    arb_t beta; arb_init(beta); 
    int num;

    FILE *fp;
    char numers[20];
    fp = fopen("init.dat","r");

    fscanf(fp, "%s", numers); 
    num = atoi(numers); 
    fscanf(fp, "%s", numers); 
    prec = atol(numers); 
    fscanf(fp, "%s", numers); 
    uni = atol(numers);
    fscanf(fp, "%s", numers); 
    couni = atol(numers);
    fscanf(fp, "%s", numers); 
    maxmul = atol(numers);
    fscanf(fp, "%s", numers); 
    total = atol(numers);
    fscanf(fp, "%s", numers); 
    digits = atol(numers);

    fclose(fp);

//    flint_printf("Checking the value beta = \t");
//    arb_printd(beta,print_dig);
//    flint_printf("\n");

    slong *multmatrix = (slong *)malloc(uni*(maxmul+1)*sizeof(slong)); 
    slong *mapmats = (slong *)malloc(total*4*sizeof(slong));
    slong *comulti = (slong *)malloc(total*sizeof(slong));
    slong *multi = (slong *)malloc(total*sizeof(slong));
    
    slong *params = (slong *)malloc(2*sizeof(slong));

    bool_mat_t temp;
    bool_mat_init(temp,uni,couni);
 
    readmats(multmatrix,mapmats,temp,multi,comulti,total,uni,couni,maxmul,num);

    arb_t eps; arb_init(eps); arb_one(eps); arb_mul_2exp_si(eps,eps,-4*(slong)digits);
    flint_printf("Computing Hausdorff dimension with an error up to \t");
    arb_printd(eps,20); 

    flint_printf(" Matrices loaded \n");

    slong k0, k1;     

    int sw0 = 1, maxc0 = 10, act0 = 0; //success switch --- dimension guess confirmed with m given; attempts counters; 
    int sw1 = 0, sw3 = 0; // switches
    slong m = 6, m0, logpart = -8;    

    
    arb_ptr nuls;     
    arb_t Lamq; arb_init(Lamq); 
    arb_poly_t testfs[uni]; // test functions

    arf_ptr bb; bb = _arf_vec_init(2); 
    arf_t une, a, b; 
    arf_init(une); arf_one(une); 
    arf_init(a); arf_zero(a);
    arf_init(b); arf_one(b);
    
    arb_t t0, t1, ll; arb_init(t0); arb_init(t1); arb_zero(t0); arb_one(t1);
    arb_add(beta,t0,t1,prec);
    arb_div_si(beta,beta,2,prec);   

    while ( (intv_length(t0,t1,eps,prec) > 0) && (act0 < maxc0) && (logpart > -20)) {
        nuls = _arb_vec_init(m); 
        cheb_zeros(nuls,m,prec);

        for (k0=0; k0<uni; k0++) {
                    arb_poly_init(*(testfs+k0)); 
                    arb_poly_fit_length(*(testfs+k0),m);
        }

        params[0] = prec;  
        params[1] = logpart;

        sw1 = leading_vector(testfs,Lamq,params,multmatrix,mapmats,comulti,temp,beta,m,nuls,uni,total,maxmul,a,b);

        if (sw1!=0) {
            while ((sw1!=0) && (prec < 2048) && (m*uni < 3500) ) {
                act0++;
                if (sw1 == 1) { //no eigenvector --- perhaps we should try the usual method  
                    prec = prec*2 + 128;
                    params[0] = prec;
                    cheb_zeros(nuls,m,prec);
                    flint_printf("sw1= %wd, \t prec = %wd, \t m = %wd, \t logpart = %wd \n",sw1,prec,m,logpart);
                    sw1 = leading_vector(testfs,Lamq,params,multmatrix,mapmats,comulti,temp,beta,m,nuls,uni,total,maxmul,a,b);
                } 
                if (sw1 == 2) { //the resulting function is not positive
                    prec = prec*2 + 128;
                    params[0] = prec;
                    _arb_vec_clear(nuls,m); 
                    m = m + 2;
                    nuls = _arb_vec_init(m); 
                    cheb_zeros(nuls,m,prec);
                    flint_printf("sw1= %wd, \t prec = %wd, \t m = %wd, \t logpart = %wd \n",sw1,prec,m,logpart);
                    sw1 = leading_vector(testfs,Lamq,params,multmatrix,mapmats,comulti,temp,beta,m,nuls,uni,total,maxmul,a,b);
                }
                if (sw1 == 3) { //m is too small for the eigenvalue
                    _arb_vec_clear(nuls,m); 
                    m = m + 2;
                    nuls = _arb_vec_init(m); 
                    cheb_zeros(nuls,m,prec);
                    flint_printf("sw1= %wd, \t prec = %wd, \t m = %wd, \t logpart = %wd \n",sw1,prec,m,logpart);
                    sw1 = leading_vector(testfs,Lamq,params,multmatrix,mapmats,comulti,temp,beta,m,nuls,uni,total,maxmul,a,b);
                }
            }
        } 

        prec = params[0];
        logpart = params[1];

        if (sw1==0) {    
            ratbounds(bb, logpart, testfs, beta, nuls, a, b, m-1, uni, temp, mapmats, comulti, multi, total, prec);

                if ( arf_cmp(bb,une)>0 )  {
                    sw0 = 0; m0 = m; prec = prec*2; 
                    arf_printd(bb,print_dig);
                    flint_printf("--- lower bound \n"); 
                    arf_printd(bb+1,print_dig);
                    flint_printf("--- upper bound \n");  

                    arb_printd(Lamq,print_dig); 
                    flint_printf("-- eigenvalue; inf > 1. \n ");  
                    flint_printf("Dimension is larger than beta = \t");
                    arb_printd(beta,print_dig);
                    flint_printf("\n");
                    arb_set(t0,beta); arb_add(beta,t0,t1,prec); arb_div_si(beta,beta,2,prec);
                } else if ( arf_cmp(bb+1,une) < 0 ) { 
                    sw0 = 0; m0 = m; prec = prec*2;
                    arf_printd(bb,print_dig);
                    flint_printf("--- lower bound \n"); 
                    arf_printd(bb+1,print_dig);
                    flint_printf("--- upper bound \n");  

                    arb_printd(Lamq,print_dig); 
                    flint_printf("-- eigenvalue; sup < 1. \n ");  
                    flint_printf("Dimension is smaller than beta = \t");
                    arb_printd(beta,print_dig);
                    flint_printf("\n");
                    arb_set(t1,beta); arb_add(beta,t0,t1,prec); arb_div_si(beta,beta,2,prec);
                } else {
                    act0++; 
                    prec = prec*2 + 32; //leading_vector reduced precision by a factor of 2
                    logpart = logpart - 1;

                    if (m*uni < 3500) { sw3 ++; m0 = m + 2;  } else {  m0 = m; }
//                    if (sw3 > 2) { logpart = logpart-2; sw3 = 0;  } // adjust logpart and m 

                    arf_printd(bb,print_dig);
                    flint_printf("--- lower bound \n"); 
                    arf_printd(bb+1,print_dig);
                    flint_printf("--- upper bound \n");  
                    arb_printd(Lamq, print_dig); 
                    flint_printf("-- eigenvalue \n ");  
                    flint_printf("I cannot confirm this estimate yet, sorry. \n");
                } //ratbounds sorted 
                _arb_vec_clear(nuls,m);
                for (k1=0; k1<uni; k1++) { arb_poly_clear(*(testfs+k1)); }
                m = m0;
                flint_printf("prec = %wd, \t m = %wd, \t logpart = %wd \n",prec,m,logpart);
        } else {
            flint_printf("\n Cannot compute the eigenvector. \n");
            sw0 = 0; 
        }  // sw1 = 0 no eigenvector  
    } //  while 

    flint_printf("t0: \t");
    arb_printd(t0,print_dig);
    flint_printf("\n");
    

    flint_printf("t1: \t"); 
    arb_printd(t1,print_dig);
    flint_printf("\n");


    
    free(multmatrix); free(mapmats); 
    free(multi); free(comulti); 
    bool_mat_clear(temp);
    arb_clear(beta);
    arb_clear(t0); 
    arb_clear(t1);
    free(params);

    arf_clear(a);     arf_clear(b);     arf_clear(une);
    _arf_vec_clear(bb,2);
    flint_cleanup();
    return 0;
}
 
