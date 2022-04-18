#include "arb.h"
#include "flint/profiler.h"
#include "arb_mat.h" 
#include "acb_mat.h" 
#include "bool_mat.h" 
#include "math.h" 
#include "string.h" 
#include "arb_fmpz_poly.h"

#include <omp.h>
//#include "immintrin.h"

void test_fun_poly(arb_poly_t res, const arb_ptr Vq, const arb_ptr nuls, slong m, slong prec) {

    // returns test function as a polynomial
    arb_t x0, x1, pm, sum1;
    arb_init(x0); arb_init(x1); arb_init(sum1); arb_init(pm);

    arb_poly_t p1, p2;
    arb_poly_init(p1);
    arb_poly_init(p2);

    slong k, j, k0;

    fmpz *zcoefs;
    arb_ptr tc, pc;
    arb_poly_zero(p2);
    for (k = 1; k<m; k++) {
     //   flint_printf("k = %wd \n", k);
        arb_zero(sum1);
        for (j=0; j<m; j++) {
            cheb_shift_t(x0,k,nuls+j,prec);
            arb_addmul(sum1,x0,Vq+j, prec);
        }
        arb_mul_ui(sum1,sum1,2,prec);
        arb_div_si(sum1,sum1,m,prec);

        tc = _arb_vec_init(k+1);
        tcoefs(tc, k, prec);

        pc = _arb_vec_init(k+1);
        _arb_vec_zero(pc,k+1);
        for (j=0; j<k+1; j++) {
            arb_zero(x0);
            for (k0 = j; k0<k+1; k0++) {
                arb_one(pm); arb_neg(pm,pm); arb_pow_ui(pm,pm,(ulong)k0-j,prec); 
                arb_bin_uiui(x1,(ulong)k0,(ulong)j,prec);
                arb_mul(x1,x1,pm,prec);
                arb_addmul(x0,x1,tc+k0,prec);
            }
            arb_one(x1); arb_mul_ui(x1,x1,2,prec); arb_pow_ui(x1,x1,j,prec);

            arb_mul(pc+j,x0,x1,prec);
        }
        
        zcoefs = _fmpz_vec_init(k+1); 
        _arb_vec_get_unique_fmpz_vec(zcoefs,pc,k+1);
       
        for (j=0; j<k+1; j++) {  
            arb_set_fmpz(x0,zcoefs+j);
            arb_poly_set_coeff_arb(p1,j,x0);   
        }

       // arb_poly_printd(p1,20);
       // flint_printf("\t --- p1 \n");

        arb_poly_scalar_mul(p1,p1,sum1,prec);
        arb_poly_add(p2,p2,p1,prec);
  
        arb_poly_zero(p1);
        _arb_vec_clear(tc,k+1);
        _arb_vec_clear(pc,k+1);
        _fmpz_vec_clear(zcoefs,k+1);
    }
   // arb_poly_printd(p2,20);
   // flint_printf("\t --- p2 \n");

    arb_zero(sum1);
    for (j=0; j<m; j++) {
        arb_add(sum1,sum1,Vq+j,prec);
    }
    arb_div_si(sum1,sum1,m,prec);

    arb_poly_one(p1); 
    arb_poly_scalar_mul(p1,p1,sum1,prec);
    //arb_poly_printd(p1,20);
    //flint_printf("\t --- p1 \n");
   
    arb_poly_add(p2,p2,p1,prec);  // p2 is now denominator polynomial 
    //arb_poly_printd(p2,20);
    //flint_printf("\t --- p2 \n");

    arb_poly_set(res,p2);
    
    arb_clear(x0); arb_clear(x1); arb_clear(sum1); arb_clear(pm);

    arb_poly_clear(p1);
    arb_poly_clear(p2);

}

void function_San0(arb_poly_t res, arb_poly_t reslin, const arb_poly_t q, const slong *an, const slong d, slong prec, arb_poly_t temp1, arb_poly_t temp2) {

    arb_poly_t l1 ; //, m1, m2; 
  //  arb_poly_init(m1); arb_poly_init(m2);
    arb_poly_init(l1);
    arb_poly_fit_length(l1, 2); 
 //   arb_poly_fit_length(m1, d+1);
 //   arb_poly_fit_length(m2, d+1);
    arb_ptr qc; 
    qc = _arb_vec_init(d+1);

    ulong k0; 
    for (k0 = 0; k0 < d+1; k0++) {
        arb_poly_get_coeff_arb(qc+k0,q,k0); 
    }

   // flint_printf("matrix elements: %wd \t %wd \t %wd \t %wd; \t d = %wd \n", an[0],an[1],an[2],an[3], d);
   // arb_poly_printd(q,10);
   // flint_printf("\n");

    arb_poly_set_coeff_si(l1,1,*(an));
    arb_poly_set_coeff_si(l1,0,*(an+1));
    
    arb_poly_set_coeff_si(reslin,1,*(an+2));
    arb_poly_set_coeff_si(reslin,0,*(an+3)); 
    
    arb_poly_zero(res);
    for (k0 = 0; k0 < d+1; k0++) {
        arb_poly_pow_ui(temp1,l1,k0,prec);
        arb_poly_pow_ui(temp2,reslin,d-k0,prec);
        arb_poly_scalar_mul(temp1,temp1,qc+k0,prec); 
        arb_poly_mul(temp1,temp1,temp2,prec);
        arb_poly_add(res,res,temp1,prec); 
    }

    _arb_vec_clear(qc,d+1);
    arb_poly_clear(l1);
//    arb_poly_clear(m1);
//    arb_poly_clear(m2); 

}


void function_San1(arb_poly_t res, const arb_poly_t lin, const arb_poly_t s0, const slong gamma, const arb_t pow, slong prec, arb_poly_t temp1, arb_poly_t temp2) {

    arb_t c0; 
    arb_init(c0); 

    arb_mul_si(c0,pow,gamma,prec); 

    arb_poly_scalar_mul(res,s0,c0,prec);

 //   arb_poly_t difs0; // , lin; 
 //   arb_poly_init(difs0); // arb_poly_init(lin);
 //   arb_poly_fit_length(difs0,d+1); 
       
    arb_poly_derivative(temp1,s0,prec);
    arb_poly_mul(temp1,temp1,lin,prec); 

    arb_poly_add(res,res,temp1,prec);

    arb_clear(c0);
 //   arb_poly_clear(difs0); 
}

void function_Panbn0(arb_poly_t res, const arb_poly_t lin, const arb_poly_t s0, const arb_poly_t s1, const arb_poly_t q, slong prec, arb_poly_t temp1, arb_poly_t temp2) { 


 //   arb_poly_t difq; 
 //   arb_poly_init(difq); 
    
    arb_poly_derivative(temp1,q,prec); 
    arb_poly_mul(temp1,temp1,lin,prec); 
    arb_poly_mul(temp1,temp1,s0,prec); 
    arb_poly_neg(temp1,temp1); 

    arb_poly_mul(res,s1,q,prec); 
    arb_poly_add(res,res,temp1,prec); 

//    arb_poly_clear(difq);
//    arb_poly_clear(lin);  
}

void function_Panbn1(arb_poly_t res, const arb_poly_t lin, const arb_poly_t p0, const arb_t pow, const slong gamma, slong prec, arb_poly_t temp1, arb_poly_t temp2)  {

    arb_t c0; 
    arb_init(c0); 

    arb_mul_si(c0,pow,gamma,prec);

//    arb_poly_t difp0, ptemp; 
//    arb_poly_init(difp0);  arb_poly_init(ptemp); 

    arb_poly_derivative(temp1,p0,prec);
    arb_poly_mul(temp1,temp1,lin,prec); 
    
    arb_poly_scalar_mul(temp2,p0,c0,prec); 

    arb_poly_add(res,temp1,temp2,prec); 

    arb_clear(c0); 
//    arb_poly_clear(difp0); arb_poly_clear(ptemp); 
}

void taylorbound(arb_t res, arb_ptr terms, const arf_t rad, slong n, slong prec) {
// given the values of the first N derivatives at the centre and the last derivative on a ball, return lower bound 
    slong k0; 
    
    arb_t x0;   arb_init(x0);  
    arf_t u0;   arf_init(u0); 

    for (k0=1; k0<n; k0++) {
        arb_get_abs_ubound_arf(u0,terms+n-k0,prec);
        arb_set_arf(x0,u0);
        arb_mul_arf(x0,x0,rad,prec); 
        arb_add_error(terms+n-k0-1,x0);
    }
    arb_set(res,terms); 
    arb_clear(x0); arf_clear(u0);
}

void function_Nbn(arb_mat_t res, const arb_poly_t pqn, const arb_ptr z0, const arb_ptr baz0, const slong totalz0, 
                  const slong ndif, const slong *mapmatrix, const arb_ptr pows,  
                  const slong d, const slong an, const slong prec)   {
    
//    function_Nbn(numdiv, pqn, z0, baz0, npart, ndif, mapmatrix, pows, d, an, prec); 
// take denominators and collection of points produce derivative of the numerator at this points 
// return upper bound on derivative of denominators 
// pows = -2q-d, -2q-d-1, \ldots, -2q - d - (d+1)/2+1
    
    slong k0, k1, k2; 
    arb_poly_t s0, s1, p0, p1, temp1, temp2;
    arb_poly_init(s0);  arb_poly_init(s1); arb_poly_init(p0); arb_poly_init(p1); 
    arb_poly_init(temp1); arb_poly_init(temp2);
    arb_poly_fit_length(s1,(d+1));
    arb_poly_fit_length(s0,(d+1));
    arb_poly_fit_length(p0,2*(d+1));
    arb_poly_fit_length(p1,2*(d+1));

    arb_poly_t lin;
    arb_poly_init(lin);
    
    arb_t x0, x1; 
    arb_init(x0); arb_init(x1);    

    arb_ptr t0; t0 = _arb_vec_init(2); _arb_vec_zero(t0,2);
    arb_t t1, t2; arb_init(t1); arb_init(t2); arb_zero(t1); 
   
    arb_mat_t xv1; arb_mat_init(xv1,totalz0,ndif);  
    arb_ptr xv0; xv0 = _arb_vec_init(totalz0); 
   
    arb_mat_zero(res);

    for (k0 = 0; k0<an; k0++) {
       arb_poly_zero(s1);
       arb_poly_zero(p0);
//       flint_printf("thread #= %wd \n", omp_get_thread_num());

       function_San0(s0,lin,pqn,(mapmatrix+4*k0),d,prec,temp1,temp2);
       
       for (k2 = 0; k2< totalz0; k2++) {
           arb_poly_evaluate(t1,lin,z0+k2,prec); 
           for (k1 = 0; k1 < ndif-1; k1++) {
               arb_pow(arb_mat_entry(xv1,k2,k1),t1,pows+k1,prec);
           }
           
           arb_poly_evaluate(t1,lin,baz0+k2,prec); 
           arb_pow(xv0+k2,t1,pows+k1,prec); //k1 = ndif
       
       }
       
       for (k2 = 0; k2<totalz0; k2++)  {
           arb_poly_evaluate(x1,s0,z0+k2,prec);
           arb_addmul(arb_mat_entry(res,k2,0),x1,arb_mat_entry(xv1,k2,0),prec);
       }

       function_San1(s1,lin,s0,*(mapmatrix+4*k0+2),pows,prec,temp1,temp2);

       function_Panbn0(p0,lin,s0,s1,pqn,prec,temp1,temp2);
       for (k1 = 1; k1< ndif-1; k1++)  {
           for (k2 = 0; k2<totalz0; k2++)  {
               arb_poly_evaluate(x1,p0,z0+k2,prec);
               //                        #pragma omp atomic
               arb_addmul(arb_mat_entry(res,k2,k1),x1,arb_mat_entry(xv1,k2,k1),prec);
           }
           function_Panbn1(p0,lin,p0,pows+k1,*(mapmatrix+4*k0+2),prec,temp1,temp2);
       }  //k1
       for (k2 = 0; k2<totalz0; k2++)  {
           arb_poly_evaluate(x1,p0,baz0+k2,prec);
           arb_addmul(arb_mat_entry(res,k2,k1),x1,xv0+k2,prec);
       } //k2
    } //k0

    _arb_vec_clear(t0,2); arb_clear(t1); arb_clear(t2);
    arb_poly_clear(s0);  arb_poly_clear(s1); arb_poly_clear(p0); arb_poly_clear(p1); 
    arb_poly_clear(lin); 
    arb_clear(x0); arb_clear(x1); 
    arb_mat_clear(xv1);  
    _arb_vec_clear(xv0,totalz0);
    arb_poly_clear(temp1); 
    arb_poly_clear(temp2); 
}

void partition(arb_ptr z0, arb_ptr baz0, const slong logpart, slong npart, const arf_t a, const arf_t b, slong prec) {

arf_t t0, t1, pp, err; 
arf_init (t0); arf_init(t1); arf_init(pp); arf_init(err);
slong k0; 

arf_set_ui_2exp_si(pp,1,logpart);
arf_set_ui_2exp_si(err,1,logpart*3);

for (k0=1; k0<npart+1; k0++) {
      arf_set(t0,a);
      arf_set(t1,a);

      arf_addmul_si(t0,pp,k0-1,prec,ARF_RND_DOWN);
      arf_addmul_si(t1,pp,k0,prec,ARF_RND_UP);

      arf_sub(t0,t0,err,prec,ARF_RND_FLOOR);
      arf_add(t1,t1,err,prec,ARF_RND_CEIL);

      arb_set_interval_arf(baz0+k0-1,t0,t1,prec);
      arb_get_mid_arb(z0+k0-1,baz0+k0-1);
//      arb_printd(baz0+k0-1,16); 
//      flint_printf("\t --- baz0, k0 = %wd \n", k0); 
    }
arf_clear (t0); arf_clear(t1); arf_clear(pp); arf_clear(err);
}

int ratbounds(arf_ptr bb, slong logpart, arb_poly_t pqn, const arb_t q, const arb_ptr nuls, const arf_t a, const arf_t b, const slong d, const slong *mapmatrix, const slong an, slong prec) {
    // compute the ratio of the image and the polynomial 
    arf_set_si(bb,100); 
    arf_set_si(bb+1,-100); 

    arf_t pp, dmax; 
    arf_init(pp); arf_init(dmax);
    slong npart; 

    arf_set_ui_2exp_si(pp,1,-logpart);
    npart = arf_get_si(pp,ARF_RND_CEIL);
    arf_clear(pp); 
    
    arb_ptr z0; z0 = _arb_vec_init(npart);
    arb_ptr baz0; baz0 = _arb_vec_init(npart); 

    partition(z0, baz0, logpart, npart, a, b, prec);
    
    arb_t rad; arb_init(rad); arb_get_rad_arb(rad,baz0); 
    arf_t r0; arf_init(r0); arb_get_abs_ubound_arf(r0,rad,prec);

    arb_t x0, x1, x2;
    arb_init(x0); arb_init(x1); arb_init(x2);
 
    slong k0, k1, k2;  
    
    arb_mat_t numdiv;
    slong ndif = floor((d+1)/2+1);

    arb_ptr nb; nb = _arb_vec_init(ndif-1);
    
    arb_mat_init(numdiv,npart,ndif);  
    arb_ptr pows; pows = _arb_vec_init(ndif);
    arb_mul_si(pows,q,2,prec);
    arb_add_si(pows,pows,d,prec);
    arb_neg(pows,pows); 
       
    
    for (k0 = 1; k0 < ndif; k0++)  {  arb_sub_si(pows+k0, pows+k0-1,1,prec); }
    // pows = -2q-d, -2q-d-1, \ldots, -2q - d - (d+1)/2

    function_Nbn(numdiv, pqn, z0, baz0, npart, ndif, mapmatrix, pows, d, an, prec); 
    
        if (arb_poly_is_zero(pqn) == 0 ) {

            for (k1 = 0; k1 < npart; k1++) { 
                for (k2 = 1; k2 < ndif; k2++) {
                    arb_set(nb+k2-1,arb_mat_entry(numdiv,k1,k2)); 
                }
                arb_poly_evaluate(x1,pqn,z0+k1,prec);  //denominator

                arb_div(x1,arb_mat_entry(numdiv,k1,0),x1,prec); //value of the function [LQ]_k0/Q_k0 at the centre

                //arb_printd(x1,20);
                //flint_printf("\t --- value of the function at the centre \n");

                taylorbound(x0,nb,r0,ndif-1,prec); //upper bound for the numerator of the derivative on the ball
                //arb_printd(x0,20);
                //flint_printf("\t --- numerator of the derivative on a ball \n");

                arb_poly_evaluate(x2,pqn,baz0+k1,prec); //denominator on a ball   
                arb_sqr(x2,x2,prec); //denominator on a ball squared 
                //arb_printd(x2, 10);
                //flint_printf("\t --- denominator on a ball squared \n");

                arb_div(x0,x0,x2,prec); //derivative on a ball
                //arb_printd(x0,20);
                //flint_printf("\t --- derivative on a ball \n");

                arb_mul(x0,x0,rad,prec); // derivative*radius
                arb_get_abs_ubound_arf(dmax,x0,prec); //correction: derivative*radius 

                arb_add_error_arf(x1,dmax); 
                if ( arb_contains_si(x1,1) !=0 ) {
                    flint_printf("nothing good --- return 1 \n");
                    arb_get_interval_arf(bb,bb+1,x1,prec);

                    arf_clear(dmax);  arb_clear(rad); arf_clear(r0);
                    arb_clear(x0); arb_clear(x1); arb_clear(x2);
                    _arb_vec_clear(z0,npart);  _arb_vec_clear(baz0,npart);  _arb_vec_clear(nb,ndif-1);

                    arb_mat_clear(numdiv );    
                    return 1;
                }

                arb_get_abs_ubound_arf(dmax,x1,prec);  
                if (arf_cmp(dmax,bb+1) > 0) {  arf_set(bb+1,dmax); }
                arb_get_abs_lbound_arf(dmax,x1,prec);  
                if (arf_cmp(dmax,bb)< 0) {  arf_set(bb,dmax); }

                arb_set_interval_arf(x0,bb,bb+1,prec); 
                if ( arb_contains_si(x0,1) !=0 ) {
                    flint_printf("nothing good --- return 2 \n");

                    arf_clear(dmax);  arb_clear(rad); arf_clear(r0);
                    arb_clear(x0); arb_clear(x1); arb_clear(x2);
                    _arb_vec_clear(z0,npart);  _arb_vec_clear(baz0,npart);  _arb_vec_clear(nb,ndif-1);
                    arb_mat_clear(numdiv); 

                    return 2;
                }
            } // x0 for 
        } else {
            flint_printf(" function vanished, no estimate is needed");
        }   // functions if 

    arf_clear(dmax);  arb_clear(rad); arf_clear(r0);
    arb_clear(x0); arb_clear(x1); arb_clear(x2);
    _arb_vec_clear(z0,npart);  _arb_vec_clear(baz0,npart);  _arb_vec_clear(nb,floor(d+1)/2);
    _arb_vec_clear(pows,floor(ndif));

    arb_mat_clear(numdiv);     

   
    return 0;
}

