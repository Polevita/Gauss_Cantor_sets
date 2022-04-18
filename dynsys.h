#include "arb.h"

void difdynsys(arb_ptr res, const arb_t z, const slong* amat, slong prec) {
// assuming det(amat) = 1   
// computes image of z under z \to (a_1 z + a_2)/(a_3 z + a_4) 
// and the derivative 1/(a_3 z + a_4)^2 
    arb_t x0, x1;
    arb_init(x1); arb_init(x0); 

    arb_mul_si(x0,z,*(amat),prec);
    arb_add_si(x0,x0,*(amat+1),prec);
    arb_mul_si(x1,z,*(amat+2),prec);
    arb_add_si(x1,x1,*(amat+3),prec);
    arb_div(res,x0,x1,prec);

    arb_sqr(x1,x1,prec); 
    arb_inv(res+1,x1,prec); 

    arb_clear(x0); arb_clear(x1); 
}    
