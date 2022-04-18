#include "bool_mat.h" 
#include "extraarbfuns.h"
#include "mychebysh.h" //chebyshev-related routines
#include "dynsys.h" // dynamical system is here
#include "symroutines.h" 

// #define MIN(x, y) (((x) > (y)) ? (y) : (x))

void arb_mat_window_add(arb_mat_t big, const arb_mat_t small, slong rmin, slong rmax, slong cmin, slong cmax, slong prec) {

    slong j, k;

    for (j=rmin; j<rmax; j++) {
       for (k=cmin; k<cmax; k++) {
         arb_add(arb_mat_entry(big,j,k),arb_mat_entry(big,j,k),arb_mat_entry(small,j-rmin,k-cmin),prec);
       }
    }  
}

void littlematrix(arb_mat_t smallmatrix, const slong *mapmats, const slong an, const arb_t beta, const slong m, const arb_ptr nuls, slong prec)  {

    slong j, k0, k1, k2;

    arb_ptr z, c0;
    z = _arb_vec_init(2);
    c0 = _arb_vec_init(m); 

    for (j=0; j<m; j++) {
        arb_one(c0+j); 
        if (j%2>0) {
            arb_neg(c0+j,c0+j);
        }
    }
   
    arb_mat_t bmatrix, t0;
    arb_mat_init(t0,m,m);
    arb_mat_init(bmatrix,m,m); 

    arb_t r0, r1, r2; 
    arb_init(r0); arb_init(r1); arb_init(r2);

    arb_mat_zero(smallmatrix);

    for (k0 = 0; k0 < an; k0++) {
        arb_mat_zero(bmatrix);
        for (k2 = 0; k2<m; k2++) { 
            difdynsys(z,nuls+k2,(mapmats+k0*4),prec);
            arb_mul_si(r0,z,2,prec);
            arb_sub_si(r0,r0,1,prec);
            for (k1=0; k1<m; k1++) {
                arb_abs(r2,z+1);
                arb_pow(r2,r2,beta,prec); 
                arb_mul(r2,r2,c0+k1,prec);
                arb_chebyshev_t_ui(r1,m,r0,prec);
                arb_mul(r2,r2,r1,prec);
                arb_sub(r1,z,nuls+k1,prec);
                arb_div(r2,r2,r1,prec);

                arb_mul(r1,nuls+k1,nuls+k1,prec);
                arb_sub(r1,nuls+k1,r1,prec);
                arb_sqrt(r1,r1,prec);
                arb_div_si(r1,r1,m,prec); 
                arb_mul(arb_mat_entry(bmatrix,k2,k1),r2,r1,prec);

            } // k1 
        } // k2 
        arb_mat_add(smallmatrix,smallmatrix,bmatrix,prec); 

    } //k0 
            
    arb_mat_transpose(smallmatrix,smallmatrix);

    _arb_vec_clear(z,2);
    _arb_vec_clear(c0,m);

    arb_mat_clear(bmatrix);
    arb_mat_clear(t0);
    arb_clear(r0);
    arb_clear(r1);
    arb_clear(r2);
    
    flint_cleanup();
}

int vector_positive (const arb_poly_t testfuns, slong logpart, const arf_t a, const arf_t b, slong prec) { 
// dim = number of components
// funs = test functions
// logpart = logarithm of the partition intervals

slong k0, npart;

//for (k1 = 0; k1 < dim; k1++)  {
//    arb_poly_printd(*(testfuns+k1),10);
//    flint_printf(" --- polynomial #%wd \n", k1);
//}

arb_t x0, x1; 
arb_init(x0); arb_init(x1);

arf_t ep0, ep1, pp, err, u1; 
arf_init(ep0); arf_init(ep1); arf_init(pp); arf_init(err); arf_init(u1);

arf_set_ui_2exp_si(pp,1,-logpart);
npart = arf_get_si(pp,ARF_RND_NEAR);

arf_set_ui_2exp_si(pp,1,logpart);
arf_set_ui_2exp_si(err,1,logpart*3);
flint_printf("\n taking partition into npart= \t %wd \t intervals \n",npart);

int sw = 1;

for (k0=1; k0<npart+1; k0++) {
    arf_set(ep0,a);
    arf_set(ep1,a);

    arf_addmul_si(ep0,pp,k0-1,prec,ARF_RND_DOWN);
    arf_addmul_si(ep1,pp,k0,prec,ARF_RND_UP);

    arf_sub(ep0,ep0,err,prec,ARF_RND_FLOOR);
    arf_add(ep1,ep1,err,prec,ARF_RND_CEIL);

    arb_set_interval_arf(x0,ep0,ep1,prec);
 
    if (arb_poly_is_zero(testfuns)==0) {
        arb_poly_evaluate(x1,testfuns,x0,prec);
        arb_get_lbound_arf(u1, x1, prec); 
        sw = arf_sgn(u1);
        if ((sw < 0) || sw == 0)  {
            arf_clear(ep1); arf_clear(ep0); arf_clear(pp); arf_clear(u1); arb_clear(x0); arf_clear(err);
            flint_printf("The vector doesn't give a positive function. Refine and recompute \n");
            return 0;
        } // if
    }
} // k0

arf_clear(ep1); arf_clear(ep0); arf_clear(pp); arf_clear(err); arf_clear(u1); arb_clear(x0); arb_clear(x1);
flint_printf("The vector gives a positive function, all good \n ");
return 1;

}

            
int get_eig_poly(arb_poly_t res, slong *params, const arb_ptr Vq, const arb_ptr nuls, const slong m, const arf_t a, const arf_t b) {

  slong  prec, logpart, print_dig=16; 
  int sw;

  // params[0] = prec 
  // params[1] = logpart

  prec = params[0]/2;
  logpart = params[1];

  arb_poly_fit_length(res,m);
  test_fun_poly(res,Vq,nuls,m,prec);
  arb_poly_printd(res,print_dig);
  flint_printf(" --- polynomial \n");
  
  sw = vector_positive(res,logpart,a,b,prec);
 
  while ((sw == 0) && (logpart > -20))   {
      logpart = logpart - 5; 
      sw = vector_positive(res,logpart,a,b,prec);
  }

  params[0] = prec;
  params[1] = logpart;

  if (sw == 1) { return 1; } else { return 0; }

}


int leading_vector(arb_poly_t res, arb_t rlam, slong *params, const slong *mapmats, const slong an, const arb_t beta, const slong m, const arb_ptr nuls, const arf_t a, const arf_t b)  {

    arb_mat_t B, BP, v0, v1;


    arb_ptr Vq;
  //  acb_mat_t BC, L, R;
    arb_t lam0, nv0, t0; 
    mag_t v0norm;
    arf_t anorm, eps, anorm1;
    
    slong k0, k1, count, it1=5, it2=50;
    slong dim = m;

    slong prec = params[0];

    arb_mat_init(B,dim,dim);
    arb_mat_init(BP,dim,dim);
    arb_init(lam0); 

    arb_zero(rlam);

    littlematrix(B, mapmats, an, beta, m, nuls, prec);

    arb_mat_transpose(B,B);
    
//    flint_printf("\n Matrix calculated\n");
//    arb_printd(arb_mat_entry(B,0,0),30);
//    flint_printf("\t -- a_11 \n");
//    arb_printd(arb_mat_entry(B,1,2),30);            
//    flint_printf("\t -- a_12 \n");
    
    arb_mat_set(BP,B);
    arb_init(nv0);
    for (k0 = 0; k0 < it1; k0++) {
        arb_mat_mul(BP,BP,B,prec);
        if (k0 % 4 == 0) {
            arb_mat_frobenius_norm(nv0,BP,prec);
            arb_mat_scalar_div_arb(BP,BP,nv0,prec); 
        }
    }
 //   arb_printd(arb_mat_entry(BP,1,1),20);
 //   flint_printf("--- BP(1,1)\n");

    arb_mat_init(v0,dim,1);  
    arb_mat_ones(v0);
    arf_init(anorm);    
    mag_init(v0norm);
    for (k0 = 0; k0 < it2; k0++)  {
       arb_mat_mul(v0,BP,v0,prec);
       arb_mat_bound_inf_norm(v0norm,v0);
       arf_set_mag(anorm,v0norm);
       arb_set_arf(nv0,anorm); 
       arb_mat_scalar_div_arb(v0,v0,nv0,prec);
    }

    arb_mat_clear(BP);    

    for (k0 = 0; k0 < dim; k0 ++) {
        if (arb_contains_zero(arb_mat_entry(v0,k0,0)) == 0) {  
            break; 
        } else {
            if ( k0 == (dim-1) ) {
                flint_printf("Cannot compute the eigenvalue - increase precision, prec=%wd", prec);

                arb_mat_clear(v0); arb_mat_clear(B);
                arf_clear(anorm); 
                mag_clear(v0norm); 
                arb_clear(nv0); arb_clear(lam0); 
            
            return 1; // increase only prec
            }
        }
    }


    arb_mat_init(v1,dim,1);  
    arb_mat_mul(v1,B,v0,prec);

    //arb_printd(arb_mat_entry(v1,k0,0),20);
    //flint_printf("\t --- v1 \n");
    //arb_printd(arb_mat_entry(v0,k0,0),20);
    //flint_printf("\t --- v0 \n");

    arb_div(lam0,arb_mat_entry(v1,k0,0),arb_mat_entry(v0,k0,0),prec);

    flint_printf("The leading eigenvalue: \n");
    arb_printd(lam0, 50); 
    flint_printf("\n");

    arb_one(nv0);
    arb_sub(nv0,nv0,lam0,prec); 
    arb_abs(nv0,nv0);
    arb_log_base_ui(nv0,nv0,10,prec);
 
    arb_init(t0); 
    arb_set_str(t0,"1.25",prec);
    arb_mul(t0,t0,nv0,prec);
    arb_get_lbound_arf(anorm,t0,prec); 
    arb_clear(t0); 
    k0 = arf_get_si(anorm,ARF_RND_FLOOR);
    k0 = MAX(-k0,6);
    flint_printf(" \n Minimal required m: = %wd \n",k0); 
    if (m < k0) {
        flint_printf(" Increase m \n"); 
        arb_clear(nv0);  arb_clear(lam0);
        arb_mat_clear(v1); arb_mat_clear(v0); 
        arb_mat_clear(B); arf_clear(anorm); 
        mag_clear(v0norm);
        return 3; // better to increase m now
    }

    arb_mul_ui(nv0,nv0,10,prec); // experimental: -1 --> -6 
 //   arb_printd(nv0,20);
 //   flint_printf(" nv0 \n");
    arb_get_lbound_arf(anorm,nv0,prec); 
 //   arf_printd(anorm,20);
 //   flint_printf("anorm \n");
    k0 = arf_get_si(anorm,ARF_RND_FLOOR);
    k0 = MIN(k0,-24);
 //   flint_printf("k0 = %wd \n", k0);
    arf_init(eps); 
    arf_one(eps);
    arf_mul_2exp_si(eps,eps,k0); 
 
    arf_init(anorm1);    
    arb_mat_mul(v1,B,v0,prec); 
    arb_mat_scalar_mul_arb(v0,v0,lam0,prec);
    arb_mat_sub(v1,v0,v1,prec);
    arb_mat_bound_inf_norm(v0norm,v1);
    arf_set_mag(anorm,v0norm);
    
    flint_printf("The leading eigenvector error: \t");
    arf_printd(anorm, 20);
    flint_printf("\n");
    flint_printf("Allowed error: eps=\t ");
    arf_printd(eps, 20);
    flint_printf("\n");

    count = 0;
    int sw = 1;
    if (arf_cmp(anorm,eps)>0) { 
        while((arf_cmp(anorm,eps) > 0) && (count<10) )  {
            count++;     
            arf_set(anorm1,anorm);
            for (k0 = 0; k0<50; k0++)  {  //apply the matrix to the eigenvector 100 times, renormalizing 
                arb_mat_mul(v0,B,v0,prec);
                arb_mat_bound_inf_norm(v0norm,v0);
                arf_set_mag(anorm,v0norm);
                arb_set_arf(nv0,anorm); 
                arb_mat_scalar_div_arb(v0,v0,nv0,prec);
            }   

            arb_mat_mul(v1,B,v0,prec); 
            arb_mat_scalar_mul_arb(v0,v0,lam0,prec);
            arb_mat_sub(v1,v0,v1,prec);
            arb_mat_bound_inf_norm(v0norm,v1);
            arf_set_mag(anorm,v0norm);
            flint_printf("Refinement #=%wd, The leading eigenvector error: \t",count);
            arf_printd(anorm, 20);
            flint_printf("\n");
            flint_printf("Allowed error: eps=\t ");
            arf_printd(eps, 20);
            flint_printf("\n");

            if (arf_cmp(anorm,anorm1) > 0 ) {
        
                flint_printf("\n Error got bigger - increase precision \n");
                
                mag_clear(v0norm); 
                arb_mat_clear(v0);   arb_mat_clear(v1);
                arb_mat_clear(B);    
                arf_clear(anorm);    arf_clear(anorm1);
                arf_clear(eps);      arb_clear(nv0);
                arb_clear(lam0);     
                flint_cleanup();
                return 1; //we are stuck --- increase only prec
            }
        } // while 

        if (arf_cmp(anorm,eps)>0) {  // 12 iterations didn't help -- need more iterations
            flint_printf("Cannot compute the leading eigenvector after 500 iterations. Error:= ");
            arf_printd(anorm, 20);
            flint_printf("\n Allowed error: eps=\t ");
            arf_printd(eps,  20); 
            flint_printf("\n");
            
            arb_clear(lam0);    arb_clear(nv0); 

            arb_mat_clear(v0);    arb_mat_clear(v1);     
            arb_mat_clear(B);

            mag_clear(v0norm); 
            arf_clear(anorm);    arf_clear(eps);
            arf_clear(anorm1);
            flint_cleanup();
            return 1; //we are stuck exit completely
        } else {  // got eigenvector after refinement 
            flint_printf("The leading eigenvector error: \t");
            arf_printd(anorm, 20);
            flint_printf("\n");
            flint_printf("Allowed error: eps=\t ");
            arf_printd(eps, 20); 
            flint_printf("\n");

            arb_mat_bound_inf_norm(v0norm,v0);
            arf_set_mag(anorm,v0norm);
            arb_set_arf(nv0,anorm); 
            arb_mat_scalar_div_arb(v0,v0,nv0,prec);
        
            mag_clear(v0norm);

            flint_printf("The leading eigenvector: \n");
            arb_mat_printd(v0, 30);
            flint_printf("The leading eigenvalue: \n");
            arb_printd(lam0, 20); 
            flint_printf("\n");
            
            Vq = _arb_vec_init(dim);
            
            for (k1 = 0; k1<dim; k1++) { arb_set(Vq+k1,arb_mat_entry(v0,k1,0)); }
            if (_arb_vec_is_negative(Vq,dim)) { _arb_vec_neg(Vq,Vq,dim); }

            arb_set(rlam,lam0);
            
            sw = get_eig_poly(res,params,Vq,nuls,m,a,b);
 
            if (sw) {
                arb_clear(lam0);   arb_mat_clear(B);  
                arb_mat_clear(v0);   arb_mat_clear(v1);

                arf_clear(eps); arf_clear(anorm); arf_clear(anorm1);
                arb_clear(nv0); _arb_vec_clear(Vq,dim); 
                flint_cleanup();

                return 0; // win after an effort
            } else { 
                flint_printf("The eigenvector doesn't give a positive function - increase precision");
                
                arb_clear(lam0);   arb_mat_clear(B);  
                arb_mat_clear(v0);   arb_mat_clear(v1);

                arf_clear(eps); arf_clear(anorm); arf_clear(anorm1);
                arb_clear(nv0); _arb_vec_clear(Vq,dim); 
                flint_cleanup();

                return 2; // increase precision and m 
            }

        } // if check after while

    } else { //immediate win
        flint_printf("The leading eigenvector error: \t");
        arf_printd(anorm, 20);
        flint_printf("\n");
        flint_printf("Allowed error: eps=\t ");
        arf_printd(eps, 20); 
        flint_printf("\n");

        arb_mat_bound_inf_norm(v0norm,v0);
        arf_set_mag(anorm,v0norm);
        arb_set_arf(nv0,anorm); 
        arb_mat_scalar_div_arb(v0,v0,nv0,prec);
        
        flint_printf("The leading eigenvector: \n");
        arb_mat_printd(v0, 30);
        flint_printf("The leading eigenvalue: \n");
        arb_printd(lam0, 20); 
        flint_printf("\n");

     
        Vq = _arb_vec_init(dim);

        for (k1 = 0; k1<dim; k1++) { arb_set(Vq+k1,arb_mat_entry(v0,k1,0)); }
        if (_arb_vec_is_negative(Vq,dim)) { _arb_vec_neg(Vq,Vq,dim); }

        arb_set(rlam,lam0);

        sw = get_eig_poly(res,params,Vq,nuls,m,a,b);

        _arb_vec_clear(Vq,dim);        

        if (sw == 0 ) {
            flint_printf("The eigenvector doesn't give a positive function - increase precision");

            arb_clear(lam0);   arb_mat_clear(B);  
            arb_mat_clear(v0);   arb_mat_clear(v1);
            mag_clear(v0norm); 

            arf_clear(eps); arf_clear(anorm); arf_clear(anorm1);
            arb_clear(nv0); _arb_vec_clear(Vq,dim); 
            flint_cleanup();

            return 2; // first won, then lost = increase precision and m
        }

    }

    arb_mat_clear(B);     arb_mat_clear(v0);      arb_mat_clear(v1);
    arb_clear(nv0);     arb_clear(lam0);
    arf_clear(anorm); arf_clear(anorm1);
    mag_clear(v0norm); 
    arf_clear(eps);   
    flint_cleanup();
    return 0;
}


