#include "arb.h"
#include "flint/profiler.h"
#include "arb_mat.h" 
#include "acb_mat.h" 
#include "math.h" 
#include "string.h" 
#include "arb_fmpz_poly.h"

void cheb_shift_t(arb_t res, ulong n, const arb_t x, slong prec) {
   
    arb_t t0;
    arb_init(t0);        
    arb_one(t0);
    arb_neg(t0,t0);
    arb_addmul_si(t0,x,2,prec);
   
    arb_chebyshev_t_ui(res,n,t0,prec);
    arb_clear(t0);
}

void cheb_shift_u(arb_t res, ulong n, const arb_t x, slong prec) {

    arb_t t0; 
    arb_init(t0); 
    arb_one(t0);
    arb_neg(t0,t0); 

    arb_addmul_si(t0,x,2,prec);

    arb_chebyshev_u_ui(res,n,t0,prec);

    arb_clear(t0);
}

void cheb_zeros(arb_ptr nuls, slong m, slong prec) {

// returns m zeros of the polynomial T_m shifted to [0,1]
    arb_t x0;
    arb_init(x0);
    slong j; //, condense = 20;
//    flint_printf("Zeros: m = %wd \n",m);
    for (j=1; j<m+1; j++) {
        arb_set_si(x0,2*j-1);
        arb_div_si(x0,x0,2*m,prec);
        arb_cos_pi(x0,x0,prec);
        arb_add_si(x0,x0,1,prec);
        arb_div_si(x0,x0,2,prec);
        arb_set(nuls+j-1,x0); 
    } 
    
    arb_clear(x0);

}

int tcoefs(arb_ptr coefs, ulong deg, slong prec) {
/* returns coefficients of Chebyshev polynomials */
// a_0 , a_1, ... a_{deg}     
    arb_t x1, x2;
    arb_init(x1); arb_init(x2);
    ulong k; 
    slong numc; 
    
    numc = (slong)(deg+1); //number of coefficients

    arb_ptr plusminus;  
    plusminus = _arb_vec_init(numc); 

    for (k=0; k<numc; k++) {
        if (k % 2) {
            arb_one(plusminus+k);
            arb_neg(plusminus+k,plusminus+k);
        } else {
            arb_one(plusminus+k);
        }
    } 
     
    if (deg % 2 != 0) {  
        for (k = 0; k < deg/2+1; k++) { 
             arb_set_ui(x1,deg);  
             arb_div_ui(x1,x1,deg-k,prec);   
             arb_ui_pow_ui(x2,2,deg-2*k-1,prec);
             arb_mul(x1,x1,x2,prec);
             arb_mul(x1,x1,plusminus+k,prec);
             arb_bin_uiui(x2,deg-k,k,prec); 
             arb_mul(coefs+deg-2*k,x1,x2,prec);
        }
    }  else {
            for (k = 0; k < deg/2; k++) { 
                arb_set_ui(x1,deg);  
                arb_div_ui(x1,x1,deg-k,prec);   
                arb_si_pow_ui(x2,2,deg-2*k-1,prec);
                arb_mul(x1,x1,x2,prec);
                arb_mul(x1,x1,plusminus+k,prec);
                arb_bin_uiui(x2,deg-k,k,prec); 
                arb_mul(coefs+deg-2*k,x1,x2,prec);
            }
            arb_set(coefs,plusminus+deg/2);
        }
      
    arb_clear(x1);
    arb_clear(x2); 
    _arb_vec_clear(plusminus,numc);
     
    int sw; 
    fmpz * zcoefs;
    zcoefs = _fmpz_vec_init(numc); 

    sw = _arb_vec_get_unique_fmpz_vec(zcoefs,coefs,numc);
    if (!sw) {
        flint_printf("Error in computing coefficients of the Chebyshev polynomial \n");
    } else {
        for (k = 0; k<numc; k++) {
            arb_set_fmpz(coefs+k,zcoefs+k);
        }
    }
    
    _fmpz_vec_clear(zcoefs,numc);
    if (sw) {
        return 1;
    } else {
        return 0;
    }
}


