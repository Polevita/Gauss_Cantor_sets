#include "get_polyfuns.h"

void mapmatrices(slong *mapmat, slong An, const slong *albet) {
    
    slong k0;  

    for (k0 = 0; k0 < An; k0++) {
        mapmat[k0*4] = 0;
        mapmat[k0*4+1] = 1;
        mapmat[k0*4+2] = 1;
        mapmat[k0*4+3] = albet[k0];
    
    }

}

int main(int argc, char *argv[]) {

    slong prec, digits, an, *alphbet;
    arb_t beta; arb_init(beta); 
    slong k0;     

    if (argc < 5) {
        flint_printf(" \n Precision, #alphabet, Alphabet, and #digits are required. Nothing to do. Stopping \n");
        return 0;
    }
    else {
        prec = atol(argv[1]); 
        an = atol(argv[2]);
        alphbet = (slong *)malloc(an*sizeof(slong)); 
        for (k0 = 0; k0 < an; k0++) {
            alphbet[k0] = atol(argv[3+k0]);
        }
        digits = atol(argv[3+an]);
    }
   
   flint_printf("Parameters: prec = %wd, an =%wd, digits = %wd \n", prec, an, digits);
   for (k0 = 0; k0 < an; k0++) {
       flint_printf("a_n = %wd, \t n = %wd \t", alphbet[k0], k0);
   }   

//    arb_printd(beta,print_dig);
    flint_printf("\n");

    slong *mapmats = (slong *)malloc(an*4*sizeof(slong));
    
    slong *params = (slong *)malloc(2*sizeof(slong));

    mapmatrices(mapmats,an,alphbet);

    arb_t eps; arb_init(eps); 
    arf_t une, a, b; 
    arf_init(une); 
    arb_set_si(eps,10);
    arb_log_base_ui(eps,eps,2,prec); 
    arb_mul_si(eps,eps,digits,prec);
    arb_add_si(eps,eps,1,prec);

    arb_get_ubound_arf(une,eps,prec); 
    k0 = arf_get_si(une,ARF_RND_CEIL);
    k0 = -k0;
    
    arb_one(eps); arb_mul_2exp_si(eps,eps,k0);
    flint_printf("Computing Hausdorff dimension with an error up to \t");
    arb_printd(eps,20); 

    flint_printf(" Matrices computed \n");


    int sw0 = 1, maxc0 = digits + 10, act0 = 0; //success switch --- dimension guess confirmed with m given; attempts counters; 
    int sw1 = 0, sw3 = 0; // switches
    slong m = 6, m0, logpart = -8;    
    slong print_dig = digits + 5;

    
    arb_ptr nuls;     
    arb_t Lamq; arb_init(Lamq); 
    arb_poly_t testfs; // test functions

    arf_ptr bb; bb = _arf_vec_init(2); 
    arf_one(une); 
    arf_init(a); arf_zero(a);
    arf_init(b); arf_one(b);
    
    arb_t t0, t1; arb_init(t0); arb_init(t1); 
    arb_zero(t0); arb_one(t1);
    arb_add(beta,t0,t1,prec);
    arb_div_si(beta,beta,2,prec);   

    while ( (intv_length(t0,t1,eps,prec) > 0) && (act0 < maxc0) && (logpart > -20)) {
   //     act0++;
        nuls = _arb_vec_init(m); 
        cheb_zeros(nuls,m,prec);

        arb_poly_init(testfs); 
        arb_poly_fit_length(testfs,m);

        params[0] = prec;  
        params[1] = logpart;

        sw1 = leading_vector(testfs,Lamq,params,mapmats,an,beta,m,nuls,a,b);

        if (sw1!=0) {
            while ((sw1!=0) && (prec < 2048)  ) {
                act0++;
                if (sw1 == 1) { //no eigenvector --- perhaps we should try the usual method  
                    prec = prec+128;
                    params[0] = prec;
                    cheb_zeros(nuls,m,prec);
                    flint_printf("sw1= %wd, \t prec = %wd, \t m = %wd, \t logpart = %wd \n",sw1,prec,m,logpart);
                    sw1 = leading_vector(testfs,Lamq,params,mapmats,an,beta,m,nuls,a,b);
                    //sw1 = leading_vector(testfs,Lamq,params,multmatrix,mapmats,comulti,temp,beta,m,nuls,nunique,total,maxmul,a,b);
                } 
                if (sw1 == 2) { //the resulting function is not positive
                    prec = prec+128;
                    params[0] = prec;
                    _arb_vec_clear(nuls,m); 
                    m = m + 2;
                    nuls = _arb_vec_init(m); 
                    cheb_zeros(nuls,m,prec);
                    flint_printf("sw1= %wd, \t prec = %wd, \t m = %wd, \t logpart = %wd \n",sw1,prec,m,logpart);
                    sw1 = leading_vector(testfs,Lamq,params,mapmats,an,beta,m,nuls,a,b);
                }
                if (sw1 == 3) { //m is too small for the eigenvalue
                    _arb_vec_clear(nuls,m); 
                    m = m + 2;
                    nuls = _arb_vec_init(m); 
                    // prec = prec*2;
                    params[0] = prec;
                    cheb_zeros(nuls,m,prec);
                    flint_printf("sw1= %wd, \t prec = %wd, \t m = %wd, \t logpart = %wd \n",sw1,prec,m,logpart);
                    sw1 = leading_vector(testfs,Lamq,params,mapmats,an,beta,m,nuls,a,b);
                }
            }
        } 

        prec = params[0];
        logpart = params[1];

        if (sw1==0) {    
            ratbounds(bb, logpart, testfs, beta, nuls, a, b, m-1, mapmats, an, prec);

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

                    if (m < 3500) { sw3 ++; m0 = m + 2;  } else {  m0 = m; }
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
                arb_poly_clear(testfs);
                m = m0;
                flint_printf("prec = %wd, \t m = %wd, \t logpart = %wd \n",prec,m,logpart);
        } else {
            flint_printf("\n Cannot compute the eigenvector. \n");
            sw0 = 0; 
        }  // sw1 = 0 no eigenvector  
    } //  while 

    flint_printf("Failed attempts = %wd \n t0: \t", act0);
    arb_printd(t0,print_dig);
    flint_printf("\n");
    

    flint_printf("t1: \t"); 
    arb_printd(t1,print_dig);
    flint_printf("\n");

    free(mapmats); 
    arb_clear(beta);
    arb_clear(Lamq);
    arb_clear(t0); 
    arb_clear(t1);
    free(params);
    free(alphbet);

    arf_clear(a);     arf_clear(b);     arf_clear(une);
    _arf_vec_clear(bb,2);
    flint_cleanup();
    return 0;
}
 
