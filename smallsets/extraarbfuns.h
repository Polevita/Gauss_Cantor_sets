#include "arb.h"

arf_ptr _arf_vec_init(slong n)
{
    slong i;
    arf_ptr v = (arf_ptr) flint_malloc(sizeof(arf_struct) * n);

    for (i = 0; i < n; i++)
        arf_init(v + i);

    return v;
}

void _arf_vec_clear(arf_ptr v, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        arf_clear(v + i);
    flint_free(v);
}


int _arb_vec_is_negative(arb_ptr v, slong n) {
    slong j, k=1;
    for (j = 0; j<n; j++) {
        if (arb_is_positive(v+j)) {
           k=0;
        }
    }
    if (k) {
        return 1; // yes 
    } else {
        return 0; // no
    }
}

void mid_point(arb_t z, const arb_t x, const arb_t y, slong prec) {
    arb_union(z,x,y,prec);
    arb_get_mid_arb(z,z);
}

int intv_length(const arb_t t0, const arb_t t1, const arb_t eps, slong prec) {
    // checks if t1-t0>eps 
    arb_t x0; 
    arb_init(x0);
    arb_sub(x0,t1,t0,prec); 
    if (arb_gt(x0,eps)) {
//        flint_printf("The length of the interval > epsilon \n");
        arb_clear(x0); return 1; 
    } else {
//        flint_printf("The length of the interval < epsilon \n");
        arb_clear(x0); return 0;
    }
}


