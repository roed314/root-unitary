/*
  Low-level code to exhaust over trees of Weil polynomials.
  This code does not implement parallelism; see the Cython wrapper.

  TODO: check for memory leaks.
*/

#include <flint.h>
#include <fmpz_poly.h>
#include <fmpq.h>
#include <fmpq_mat.h>
#include <arith.h>

/*
    Use a subresultant (Sturm-Habicht) sequence to test whether a given
    polynomial has all real roots. Note that this test has an early abort
    mechanism: having all real roots means that the sign sequence has
    the maximal number of sign changes, so the test aborts as soon
    as a sign change is missed.

    This function assumes that:
        - {poly, n} is a normalized vector with n >= 2
        - {w, 2*n+1} is scratch space.

   Based on code by Sebastian Pancratz from the FLINT repository.
   TODO: compare with floating-point interval arithmetic.
*/

int _fmpz_poly_all_real_roots(fmpz *poly, slong n, fmpz *w, int force_squarefree) {
  fmpz *f0     = w + 0*n;
  fmpz *f1     = w + 1*n;
  fmpz *c      = w + 2*n;
  fmpz *d      = w + 2*n+1;
  fmpz *t;
  
  if (n <= 2) return(1);
  _fmpz_vec_set(f0, poly, n);
  _fmpz_poly_derivative(f1, f0, n);
  n--;
  int sgn0_l = fmpz_sgn(f0+n);
  
  for ( ; ; ) {
    /* At this point deg(f0) = n, deg(f1) = n-1.
       We explicitly compute the pseudoremainder of f0 modulo f1:
       f0 := f1[n-1]*f0 - f0[n]*x*f1
       f0 := f0[n-1]*f1 - f1[n-1]*f0
    */
    fmpz_set(c, f0+n);
    _fmpz_vec_scalar_mul_fmpz(f0, f0, n, f1+n-1);
    _fmpz_vec_scalar_submul_fmpz(f0+1, f1, n-1, c);
    n--;
    fmpz_set(c, f0+n);
    fmpz_neg(d, f1+n);
    _fmpz_vec_scalar_mul_fmpz(f0, f0, n, d);
    _fmpz_vec_scalar_addmul_fmpz(f0, f1, n, c);
    
    if (!force_squarefree && _fmpz_vec_is_zero(f0, n)) return(1);
    
    /* If we miss any one sign change, we cannot have enough. */
    if (fmpz_sgn(f0+n-1) != sgn0_l) return(0);

    if (n==1) return(1); /* If f0 is a scalar, it is nonzero and we win. */
    
    /* Extract content from f0; in practice, this seems to do better than
       an explicit subresultant computation. */
    _fmpz_vec_content(c, f0, n);
    _fmpz_vec_scalar_divexact_fmpz(f0, f0, n, c);
    
    /* Swap f0 with f1. */
    t = f0; f0 = f1; f1 = t;
  }
}

/* Primary data structures for tree exhaustion.
 */

typedef struct ps_static_data {
  int d, sign, force_squarefree;
  long node_limit;
  fmpz_t lead, q;
  fmpz_mat_t binom_mat;
  fmpz *cofactor;
  fmpz *modlist;
  fmpq_mat_t *hausdorff_mats;
  fmpq_mat_t *sum_mats;
  fmpq_t *f;
} ps_static_data_t;

typedef struct ps_dynamic_data {
  int d, n, ascend, flag;
  long node_count;
  fmpq_mat_t power_sums, sum_prod, hankel_mat, hankel_dets,
    hausdorff_prod, hausdorff_sums1, hausdorff_sums2;
  fmpz *pol, *sympol, *upper;

  /* Scratch space */
  fmpz *w;
  int wlen; /* = 3*d+8 */
  fmpq *w2;
  int w2len; /* = 5 */
} ps_dynamic_data_t;

/* Set res to floor(a). */
void fmpq_floor(fmpz_t res, const fmpq_t a) {
  fmpz_fdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

/* Set res to ceil(a). */
void fmpq_ceil(fmpz_t res, const fmpq_t a) {
  fmpz_cdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

void fmpz_sqrt_f(fmpz_t res, const fmpz_t a) {
  fmpz_sqrt(res, a);
}

void fmpz_sqrt_c(fmpz_t res, const fmpz_t a) {
  int s = fmpz_is_square(a);
  fmpz_sqrt(res, a);
  if (!s) fmpz_add_ui(res, res, 1);
}

/* Set res to floor(a + b sqrt(q)). 
   For efficiency, we do not assume a and b are canonical;
   we must thus be careful about signs. */
void fmpq_floor_quad(fmpz_t res, const fmpq_t a,
		     const fmpq_t b, const fmpz_t q) {
  if (b==NULL) fmpq_floor(res, a);
  else {
    fmpz *anum = fmpq_numref(a);
    fmpz *aden = fmpq_denref(a);
    int aden_s = fmpz_sgn(aden);
    fmpz *bnum = fmpq_numref(b);
    int bnum_s = fmpz_sgn(bnum);
    fmpz *bden = fmpq_denref(b);
    int bden_s = fmpz_sgn(bden);
    
    fmpz_mul(res, aden, bnum);
    fmpz_mul(res, res, res);
    fmpz_mul(res, res, q);
    if (bnum_s*bden_s >= 0) fmpz_sqrt_f(res, res);
    else {
      fmpz_sqrt_c(res, res);
      fmpz_neg(res, res);
    }
    fmpz_mul_si(res, res, aden_s*bden_s);
    fmpz_addmul(res, anum, bden);
    if (bden_s > 0) fmpz_fdiv_q(res, res, aden);
    else fmpz_cdiv_q(res, res, aden);
    fmpz_fdiv_q(res, res, bden);
  }
}

/* Set res to ceil(a + b sqrt(q)). */
void fmpq_ceil_quad(fmpz_t res, const fmpq_t a,
		     const fmpq_t b, const fmpz_t q) {
  if (b==NULL) fmpq_ceil(res, a);
  else {
    fmpz *anum = fmpq_numref(a);
    fmpz *aden = fmpq_denref(a);
    int aden_s = fmpz_sgn(aden);
    fmpz *bnum = fmpq_numref(b);
    int bnum_s = fmpz_sgn(bnum);
    fmpz *bden = fmpq_denref(b);
    int bden_s = fmpz_sgn(bden);

    fmpz_mul(res, aden, bnum);
    fmpz_mul(res, res, res);
    fmpz_mul(res, res, q);
    if (bnum_s*bden_s >= 0) fmpz_sqrt_c(res, res);
    else {
      fmpz_sqrt_f(res, res);
      fmpz_neg(res, res);
    }
    fmpz_mul_si(res, res, aden_s*bden_s);
    fmpz_addmul(res, anum, bden);
    if (bden_s > 0) fmpz_cdiv_q(res, res, aden);
    else fmpz_fdiv_q(res, res, aden);
    fmpz_cdiv_q(res, res, bden);
  }
}

/* Memory allocation and initialization. */
ps_static_data_t *ps_static_init(int d, fmpz_t q, int coeffsign, fmpz_t lead,
				 int cofactor, fmpz *modlist, long node_limit,
				 int force_squarefree) {
  int i, j, k, l;
  ps_static_data_t *st_data;
  fmpz_poly_t pol;
  fmpz_t m, const1;
  fmpq *k1;

  fmpz_poly_init(pol);
  fmpz_init(m);
  fmpz_init_set_ui(const1, 1);

  st_data = (ps_static_data_t *)malloc(sizeof(ps_static_data_t));

  st_data->d = d;
  st_data->sign = coeffsign;
  fmpz_init(st_data->q);
  fmpz_set(st_data->q, q);
  st_data->node_limit = node_limit;
  st_data->force_squarefree = force_squarefree;

  fmpz_init(st_data->lead);
  fmpz_set(st_data->lead, lead);

  st_data->cofactor = _fmpz_vec_init(3);
  switch (cofactor) {
  case 0: /* Cofactor 1 */
    fmpz_set_si(st_data->cofactor, 1);
    fmpz_set_si(st_data->cofactor+1, 0);
    fmpz_set_si(st_data->cofactor+2, 0);
    break;

  case 1: /* Cofactor x+sqrt(q) */
    fmpz_set(st_data->cofactor, st_data->q);
    fmpz_sqrt(st_data->cofactor, st_data->cofactor);
    fmpz_set_si(st_data->cofactor+1, 1);
    fmpz_set_si(st_data->cofactor+2, 0);
    break;

  case 2:  /* Cofactor x-sqrt(q) */
    fmpz_set(st_data->cofactor, st_data->q);
    fmpz_sqrt(st_data->cofactor, st_data->cofactor);
    fmpz_neg(st_data->cofactor, st_data->cofactor);
    fmpz_set_si(st_data->cofactor+1, 1);
    fmpz_set_si(st_data->cofactor+2, 0);
    break;

  case 3: /* Cofactor x^2-q */
    fmpz_neg(st_data->cofactor, st_data->q);
    fmpz_set_si(st_data->cofactor+1, 0);
    fmpz_set_si(st_data->cofactor+2, 1);
    break;
  }

  st_data->modlist = _fmpz_vec_init(d+1);
  st_data->f = _fmpq_vec_init(d+1);
  for (i=0; i<=d; i++) {
    fmpz_set(st_data->modlist+i, modlist+d-i);
    fmpq_set_si(st_data->f+i, d-i, 1);
    fmpq_div_fmpz(st_data->f+i, st_data->f+i, st_data->lead);
    /* In order to apply power sums and Descartes' rule of signs
       when the modulus is 0, we must pretend that the modulus is 1. */
    if (!fmpz_is_zero(st_data->modlist+i))
      fmpq_mul_fmpz(st_data->f+i, st_data->f+i, st_data->modlist+i);
  }

  fmpz_mat_init(st_data->binom_mat, d+1, d+1);
  for (i=0; i<=d; i++)
    for (j=0; j<=d; j++)
      fmpz_bin_uiui(fmpz_mat_entry(st_data->binom_mat, i, j), i, j);
  
  st_data->hausdorff_mats = (fmpq_mat_t *)malloc((d+1)*sizeof(fmpq_mat_t));
  for (i=0; i<=d; i++) {

    fmpq_mat_init(st_data->hausdorff_mats[i], 2*d+2, d+1);
    fmpq_mat_zero(st_data->hausdorff_mats[i]);

    for (j=0; j<=i; j++)
      for (k=0; k<=i; k++) {
	// The coefficient of t^k in (t-2 sqrt(q))^j (t+2 sqrt(q))^{i-j}, rounding down the exponent of q.
	if ((i-k)%2==0)
	  k1 = fmpq_mat_entry(st_data->hausdorff_mats[i], 2*j, k);
	else
	  k1 = fmpq_mat_entry(st_data->hausdorff_mats[i], 2*j+1, k);
	for (l=0; l<=j; l++) if (k-l >=0 && k-l<=i-j) {
	    fmpz_mul(m, fmpz_mat_entry(st_data->binom_mat, j, l),
		     fmpz_mat_entry(st_data->binom_mat, i-j, k-l));
	    if ((j-l)%2==1) fmpq_neg(m, m);
	    fmpq_add_fmpz(k1, k1, m);
	  }
	fmpz_mul_2exp(k1, k1, i-k);
	for (l=0; l<(i-k)/2; l++) fmpz_mul(k1, k1, q);
      }
  }

  st_data->sum_mats = (fmpq_mat_t *)malloc((d+1)*sizeof(fmpq_mat_t));
  for (i=0; i<=d; i++) {

    fmpq_mat_init(st_data->sum_mats[i], 1, d+1);
    fmpq_mat_zero(st_data->sum_mats[i]);

    arith_chebyshev_t_polynomial(pol, i);
    for (j=0; j<=d; j++) {
      
      /* Coefficients of 2*(i-th Chebyshev polynomial)(x/2). 
         If q != 1, the coeff of x^j is multiplied by q^{floor(i-j)/2}. */
      if (j <= i) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 0, j);
	fmpq_set_fmpz_frac(k1, fmpz_poly_get_coeff_ptr(pol, j), const1);
	fmpz_mul_2exp(m, const1, j);
	fmpq_div_fmpz(k1, k1, m);
	fmpz_set_ui(m, 2);
	fmpq_mul_fmpz(k1, k1, m);
	if (!fmpz_is_one(st_data->q) && i%2==j%2) {
	  fmpz_set(m, st_data->q);
	  fmpz_pow_ui(m, m, (i-j)/2); 
	  fmpq_mul_fmpz(k1, k1, m);
	}
      }

    }
  }
  
  fmpz_poly_clear(pol);
  fmpz_clear(m);
  fmpz_clear(const1);

  return(st_data);
}

ps_dynamic_data_t *ps_dynamic_init(int d, fmpz *coefflist) {
  ps_dynamic_data_t *dy_data;
  int i;

  dy_data = (ps_dynamic_data_t *)malloc(sizeof(ps_dynamic_data_t));
  dy_data->d = d;

  /* Initialize mutable quantities */
  dy_data->n = d;
  dy_data->node_count = 0;
  dy_data->ascend = 0;
  dy_data->pol = _fmpz_vec_init(d+1);
  dy_data->sympol = _fmpz_vec_init(2*d+3);
  if (coefflist != NULL)
    for (i=0; i<=d; i++) 
      fmpz_set(dy_data->pol+i, coefflist+i);
  
  fmpq_mat_init(dy_data->power_sums, d+1, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->power_sums, 0, 0), d, 1);
  fmpq_mat_init(dy_data->hankel_mat, d/2+1, d/2+1);
  fmpq_mat_init(dy_data->hankel_dets, d/2+1, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->hankel_dets, 0, 0), d, 1);
  fmpq_mat_init(dy_data->hausdorff_prod, 2*d+2, 1);
  fmpq_mat_init(dy_data->hausdorff_sums1, d+1, d+1);
  fmpq_mat_init(dy_data->hausdorff_sums2, d+1, d+1);

  dy_data->upper = _fmpz_vec_init(d+1);

  /* Allocate scratch space */
  fmpq_mat_init(dy_data->sum_prod, 1, 1);
  dy_data->wlen = 3*d+8;
  dy_data->w = _fmpz_vec_init(dy_data->wlen);
  dy_data->w2len = 5;
  dy_data->w2 = _fmpq_vec_init(dy_data->w2len);
  return(dy_data);
}

ps_dynamic_data_t *ps_dynamic_clone(ps_dynamic_data_t *dy_data) {
  ps_dynamic_data_t *dy_data2;
  int d = dy_data->d;

  dy_data2 = ps_dynamic_init(d, NULL);
  dy_data2->n = dy_data->n;
  dy_data2->node_count = dy_data->node_count;
  dy_data2->ascend = dy_data->ascend;
  _fmpz_vec_set(dy_data2->pol, dy_data->pol, d+1);
  _fmpz_vec_set(dy_data2->upper, dy_data->upper, d+1);
  fmpq_mat_set(dy_data2->power_sums, dy_data->power_sums);
  fmpq_mat_set(dy_data2->hankel_dets, dy_data->hankel_dets);
  fmpq_mat_set(dy_data2->hausdorff_sums1, dy_data->hausdorff_sums1);
  fmpq_mat_set(dy_data2->hausdorff_sums2, dy_data->hausdorff_sums2);
  return(dy_data2);
}

/* Split off a subtree. */
ps_dynamic_data_t *ps_dynamic_split(ps_dynamic_data_t *dy_data) {
  if (dy_data==NULL) return(NULL);

  ps_dynamic_data_t *dy_data2;
  int i, d = dy_data->d, n = dy_data->n, ascend=dy_data->ascend;

  for (i=d; i>n+ascend; i--)
    if (fmpz_cmp(dy_data->pol+i, dy_data->upper+i) <0) {
      dy_data2 = ps_dynamic_clone(dy_data);
      fmpz_set(dy_data->upper+i, dy_data->pol+i);
      dy_data2->ascend = i-n;
      dy_data2->node_count = 0;
      return(dy_data2);
  }
  return(NULL);
}

/* Memory deallocation. */
void ps_static_clear(ps_static_data_t *st_data) {
  if (st_data == NULL) return(NULL);
  int i, d = st_data->d;
  fmpz_clear(st_data->lead);
  fmpz_clear(st_data->q);
  _fmpz_vec_clear(st_data->cofactor, 3);
  fmpz_mat_clear(st_data->binom_mat);
  _fmpq_vec_clear(st_data->f, d+1);
  _fmpz_vec_clear(st_data->modlist, d+1);
  for (i=0; i<=d; i++) 
    fmpq_mat_clear(st_data->hausdorff_mats[i]);
  free(st_data->hausdorff_mats);
  for (i=0; i<=d; i++) 
    fmpq_mat_clear(st_data->sum_mats[i]);
  free(st_data->sum_mats);
  free(st_data);
}

void ps_dynamic_clear(ps_dynamic_data_t *dy_data) {
  if (dy_data == NULL) return(NULL);
  int d = dy_data->d;
  _fmpz_vec_clear(dy_data->pol, d+1);
  _fmpz_vec_clear(dy_data->sympol, 2*d+3);
  _fmpz_vec_clear(dy_data->upper, d+1);
  fmpq_mat_clear(dy_data->power_sums);
  fmpq_mat_clear(dy_data->sum_prod);
  fmpq_mat_clear(dy_data->hankel_mat);
  fmpq_mat_clear(dy_data->hankel_dets);
  fmpq_mat_clear(dy_data->hausdorff_prod);
  fmpq_mat_clear(dy_data->hausdorff_sums1);
  fmpq_mat_clear(dy_data->hausdorff_sums2);
  _fmpz_vec_clear(dy_data->w, dy_data->wlen);
  _fmpq_vec_clear(dy_data->w2, dy_data->w2len);
  free(dy_data);
}

/* The following is the key subroutine: given some initial coefficients, compute
   a lower and upper bound for the next coefficient, or detect a dead end.

   Return values: 
   -r, r>0: if the n-th truncated polynomial does not have roots in the
       interval, and likewise for all choices of the bottom r-1 coefficients
   1: if lower <= upper
   0: otherwise. 

   The case n=0 is allowed. In this case, we return 1 if the polynomial is
   admissible and 0 otherwise.
   
*/
int set_range_from_power_sums(ps_static_data_t *st_data,
			      ps_dynamic_data_t *dy_data,
			      int test_roots) {
  int i, j, r;
  int d = st_data->d;
  int n = dy_data->n;
  int k = d+1-n;
  fmpz *modulus = st_data->modlist + n-1;
  fmpz *pol = dy_data->pol;
  fmpz *q = st_data->q;
  int q_is_1;
  fmpq *f;
    
  /* Allocate temporary variables from persistent scratch space. */
  fmpz *tpol = dy_data->w;
  fmpz *tpol2 = dy_data->w + d+1;
  fmpz *tpol3 = dy_data->w + 2*d+2;

  fmpz *t0z = dy_data->w+3*d+3;
  fmpz *lower = dy_data->w+3*d+4;
  fmpz *upper = dy_data->w+3*d+5;
  
  fmpq *t0q = dy_data->w2;
  fmpq *t1q = dy_data->w2+1;
  fmpq *t2q = dy_data->w2+2;
  fmpq *t3q = dy_data->w2+3;
  fmpq *t4q = dy_data->w2+4;

  /* Embedded subroutines to adjust lower and upper bounds. 
   These use t0z, t0q, t4q as persistent scratch space.
   The pair (val1, val2) stands for val1 + val2*sqrt(q);
   passing NULL for val2 is a faster variant of passing 0. 

  Usage: if f is a monic linear function of the k-th power sum, then
  set_upper(f) or change_upper(f) imposes the condition f >= 0;
  set_lower(f) or change_lower(f) imposes the condition f <= 0.*/

  void set_lower(const fmpq_t val1, const fmpq_t val2) {
    fmpq_div(t0q, val1, f);
    if (val2==NULL) fmpq_ceil(lower, t0q);
    else {
      fmpq_div(t4q, val2, f);
      fmpq_ceil_quad(lower, t0q, t4q, q);
    }
  }
  
  void set_upper(const fmpq_t val1, const fmpq_t val2) {
    fmpq_div(t0q, val1, f);
    if (val2==NULL) fmpq_floor(upper, t0q);
    else {
      fmpq_div(t4q, val2, f);
      fmpq_floor_quad(upper, t0q, t4q, q);
    }
  }

  void change_lower(const fmpq_t val1, const fmpq_t val2) {
    fmpq_div(t0q, val1, f);
    if (val2==NULL) fmpq_ceil(t0z, t0q);
    else {
      fmpq_div(t4q, val2, f);
      fmpq_ceil_quad(t0z, t0q, t4q, q);
    }
    if (fmpz_cmp(t0z, lower) > 0) fmpz_set(lower, t0z);
  }
  
  void change_upper(const fmpq_t val1, const fmpq_t val2) {
    fmpq_div(t0q, val1, f);
    if (val2==NULL) fmpq_floor(t0z, t0q);
    else {
      fmpq_div(t4q, val2, f);
      fmpq_floor_quad(t0z, t0q, t4q, q);
    }
    if (fmpz_cmp(t0z, upper) < 0) fmpz_set(upper, t0z);
  }

  /* Impose the condition that val1*val3 >= val2, assuming that val1 is a linear monic function
     of the k-th power sum and val2, val3 do not depend on this sum. */
  void impose_quadratic_condition(const fmpq_t val1, const fmpq_t val2, const fmpq_t val3) {
    int s = fmpq_sgn(val3);
    if (s) {
      fmpq_mul(t0q, val2, val2);
      fmpq_div(t0q, t0q, val3);
      fmpq_sub(t0q, val1, t0q);
      if (s>0) change_upper(t0q, NULL);
      else change_lower(t0q, NULL);
    }
  }

  /* End embedded subroutines */
    
  /* Compute the divided n-th derivative of pol, answer in tpol. */
  for (i=0; i<=k-1; i++)
    fmpz_mul(tpol+i, fmpz_mat_entry(st_data->binom_mat, n+i, n), pol+n+i);

  /* Condition: by Rolle's theorem, tpol must have real roots. */
  /* TODO: try using real root isolation instead of Sturm sequences. */
  if (test_roots && !_fmpz_poly_all_real_roots(tpol, k, tpol2, st_data->force_squarefree))
    return(-1);
  
  /* If k>d, no further coefficients to bound. */
  if (k>d) return(1);

  /* Update power_sums[k]. */
  /* TODO: use matrix multiplication here. */
  f = fmpq_mat_entry(dy_data->power_sums, k, 0);
  fmpq_set_si(f, -k, 1);
  fmpq_mul_fmpz(f, f, pol+d-k);
  for (i=1; i<k; i++) {
    fmpq_set_si(t0q, -1, 1);
    fmpq_mul_fmpz(t0q, t0q, pol+d-i);
    fmpq_addmul(f, t0q, fmpq_mat_entry(dy_data->power_sums, k-i, 0));
  }
  fmpq_div_fmpz(f, f, pol+d);
  
  /* Condition: the k-th symmetrized power sum must lie in [-2*sqrt(q), 2*sqrt(q)]. */
  f = st_data->f + n-1;
  fmpq_mat_mul(dy_data->sum_prod, st_data->sum_mats[k], dy_data->power_sums);

  q_is_1 = !fmpz_cmp_ui(q, 1);
  if (k%2==0) {
    fmpq_set_si(t1q, 2*d, 1);
    if (!q_is_1) {
      fmpz_pow_ui(t0z, q, k/2);
      fmpq_mul_fmpz(t1q, t1q, t0z);
    }
    fmpq_sub(t0q, fmpq_mat_entry(dy_data->sum_prod, 0, 0), t1q);
    set_lower(t0q, NULL);
    fmpq_add(t0q, fmpq_mat_entry(dy_data->sum_prod, 0, 0), t1q);
    set_upper(t0q, NULL);
  } else {
    fmpq_zero(t1q); 
    fmpq_set_si(t2q, 2*d, 1);
    if (!q_is_1) {
      fmpz_pow_ui(t0z, q, k/2);
      fmpq_mul_fmpz(t2q, t2q, t0z);
    }
    set_upper(fmpq_mat_entry(dy_data->sum_prod, 0, 0), t2q);
    fmpq_neg(t2q, t2q);
    set_lower(fmpq_mat_entry(dy_data->sum_prod, 0, 0), t2q);
    }
    
  /* Undo one derivative on tpol. */
  for (i=k; i>=1; i--) {
    fmpz_mul_si(tpol+i, tpol+i-1, n);
    fmpz_divexact_si(tpol+i, tpol+i, i);
  }
  fmpz_set(tpol, pol+d-k);
  
  /* Condition: Descartes' rule of signs applies at -2*sqrt(q), +2*sqrt(q). 
   This is only a new condition for the evaluations at these points. */
  
  fmpq_set_si(t3q, -k, 1);
  fmpq_div_fmpz(t3q, t3q, pol+d);

  for (i=0; 2*i <= k; i++) fmpz_mul_2exp(tpol2+i, tpol+2*i, 2*i);
  _fmpz_poly_evaluate_fmpz(t0z, tpol2, (k+2) / 2, q);
  fmpq_mul_fmpz(t1q, t3q, t0z);

  for (i=0; 2*i+1 <= k; i++) fmpz_mul_2exp(tpol2+i, tpol+2*i+1, 2*i+1);
  _fmpz_poly_evaluate_fmpz(t0z, tpol2, (k+1) / 2, q);
  fmpq_mul_fmpz(t2q, t3q, t0z);
  
  change_lower(t1q, t2q);
  
  fmpq_neg(t2q, t2q);
  if (k%2==1) change_upper(t1q, t2q);
  else change_lower(t1q, t2q);

  /* If modulus==0, then return 1 if [lower, upper] contains 0
     and 0 otherwise. After this, we may assume modulus>0.
   */
  if (fmpz_is_zero(modulus)) {
    if ((fmpz_sgn(lower) > 0) || (fmpz_sgn(upper) < 0)) return(0);
      fmpz_zero(lower);
      fmpz_zero(upper);
      return(1);
  }

  /* Condition: nonnegativity of the Hankel determinant. */
  if (k%2==0) {
    fmpq_mat_one(dy_data->hankel_mat);
    for (i=0; i<=k/2; i++)
      for (j=0; j<=k/2; j++)
	fmpq_set(fmpq_mat_entry(dy_data->hankel_mat, i, j),
		 fmpq_mat_entry(dy_data->power_sums, i+j, 0));
    fmpq_mat_det(t0q, dy_data->hankel_mat);
    fmpq_set(fmpq_mat_entry(dy_data->hankel_dets, k/2, 0), t0q);
    fmpq_set(t3q, fmpq_mat_entry(dy_data->hankel_dets, k/2-1, 0));
    if (fmpq_sgn(t3q) > 0) {
      fmpq_div(t0q, t0q, t3q);
      change_upper(t0q, NULL);
      }
    else if (fmpq_sgn(t0q) < 0) return(0);
    else change_upper(fmpq_mat_entry(dy_data->power_sums, k, 0), NULL);
  } 

  if (fmpz_cmp(lower, upper) > 0) return(0);

  /* Condition: the Hausdorff moment criterion for having roots in [-2, 2]. */
  fmpq_mat_mul(dy_data->hausdorff_prod, st_data->hausdorff_mats[k], dy_data->power_sums);
  for (i=0; i<=k; i++) {
    fmpq_set(t1q, fmpq_mat_entry(dy_data->hausdorff_prod, 2*i, 0));
    fmpq_set(t2q, fmpq_mat_entry(dy_data->hausdorff_prod, 2*i+1, 0));
    if (i%2==0) change_upper(t1q, t2q);
    else change_lower(t1q, t2q);
    fmpq_set(fmpq_mat_entry(dy_data->hausdorff_sums1, k, i), t1q);
    fmpq_set(fmpq_mat_entry(dy_data->hausdorff_sums2, k, i), t2q);
  }
  
  if (fmpz_cmp(lower, upper) > 0) return(0);
  
  /* Condition: log convexity based on Cauchy-Schwarz. */
  /* Todo: extend to q != 1 without losing too much efficiency. */
  if (q_is_1) {
    for (i=0; i<=k-2; i++) {
      fmpq_add(t1q, fmpq_mat_entry(dy_data->hausdorff_sums1, k, i),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k, i));
      fmpq_add(t2q, fmpq_mat_entry(dy_data->hausdorff_sums1, k-1, i),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k-1, i));
      fmpq_add(t3q, fmpq_mat_entry(dy_data->hausdorff_sums1, k-2, i),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k-2, i));
      impose_quadratic_condition(t1q, t2q, t3q);
    }
    for (i=2; i<=k; i++) {
      fmpq_add(t1q, fmpq_mat_entry(dy_data->hausdorff_sums1, k, i),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k, i));
      fmpq_add(t2q, fmpq_mat_entry(dy_data->hausdorff_sums1, k-1, i-1),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k-1, i-1));
      fmpq_add(t3q, fmpq_mat_entry(dy_data->hausdorff_sums1, k-2, i-2),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k-2, i-2));
      impose_quadratic_condition(t1q, t2q, t3q);
    }
  }

  if (fmpz_cmp(lower, upper) > 0) return(0);

  fmpz_addmul(tpol, lower, modulus);
  while (1) {
    if (_fmpz_poly_all_real_roots(tpol, k+1, tpol2, st_data->force_squarefree))
      break;
    fmpz_add_ui(lower, lower, 1);
    if (fmpz_cmp(lower, upper) > 0)
      return(0);
    fmpz_add(tpol, tpol, modulus);
    }
  
  /* Set the new upper bound. Note that modulus>0 at this point. */
  fmpz_mul(upper, upper, modulus);
  fmpz_add(dy_data->upper+n-1, pol+n-1, upper);
  
  /* Set the new polynomial value. */
  fmpz_mul(lower, lower, modulus);
  fmpz_add(pol+n-1, pol+n-1, lower);

  /* Correct the k-th power sum and related quantities. */
  t1q = fmpq_mat_entry(dy_data->power_sums, k, 0);
  fmpq_mul_fmpz(t0q, f, lower);

  fmpq_sub(t1q, t1q, t0q);
  for (i=0; i<=k; i++) {
    t1q = fmpq_mat_entry(dy_data->hausdorff_sums1, k, i);
    fmpq_sub(t1q, t1q, t0q);
  }
  if (k%2==0) {
    t1q = fmpq_mat_entry(dy_data->hankel_dets, k/2, 0);
    fmpq_submul(t1q, fmpq_mat_entry(dy_data->hankel_dets, k/2-1, 0), t0q);
  }
  return(1);
}

/* Return value sent back in dy_data->flag:
    1: if a solution has been found
    0: if the tree has been exhausted
   -1: if the maximum number of nodes has been reached
   -2: none of the above
*/

void next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps) {
  if (dy_data==NULL) return(0);

  int d = st_data->d;
  int node_limit = st_data->node_limit;
  fmpz *modlist = st_data->modlist;

  int ascend = dy_data->ascend;
  int n = dy_data->n;
  int count = dy_data->node_count;
  fmpz *upper = dy_data->upper;
  fmpz *pol = dy_data->pol;
  fmpz *sympol = dy_data->sympol;

  int i, j, t, r, count_steps = 0, test_roots = 1;
  fmpq *tq;

  if (n>d) return(0);
  while (1) {
    if (ascend > 0) {
      n += 1;
      if (n>d) {
	/* We have exhausted the entire tree. */
	t=0;
	break;
      }
    } else {
      i = dy_data->n;
      dy_data->n = n;
      r = set_range_from_power_sums(st_data, dy_data, test_roots);
      test_roots = 0;
      if (r > 0) {
	n -= 1;
	if (n<0) { 
	  t=1;
	  /* We have found a solution! Convert it back into symmetric form for output. */
	  _fmpz_vec_zero(sympol, 2*d+3);
	  fmpz *temp = dy_data->w;
	  for (i=0; i<=d; i++) {
	    fmpz_one(temp);
	    for (j=0; j<=i; j++) {
	      fmpz_addmul(sympol+2*d-(d-i+2*j), pol+i, temp);
	      if (j<i) {
		fmpz_mul(temp, temp, st_data->q);
		fmpz_mul_si(temp, temp, i-j);
		fmpz_divexact_si(temp, temp, j+1);
	      }
	    }
	  }
	  _fmpz_vec_scalar_mul_si(sympol, sympol, 2*d+1, st_data->sign);
	  _fmpz_poly_mul_KS(sympol,sympol, 2*d+1, st_data->cofactor, 3);
	  ascend = 1;
	  break; 
	}
	continue;
      } else {
	count += 1;
	if (node_limit != -1 && count >= node_limit) { t= -1; break; }
	if (r<0) {
	  /* Rolle condition failed; it cannot succeed again at this level. */
	  ascend = 1;
	  continue;
	}
      }
    }
    if (ascend>1) ascend -= 1;
    else if (fmpq_is_zero(modlist+n)) ascend = 1;
    else {
      fmpz_add(pol+n, pol+n, modlist+n);
      if (fmpz_cmp(pol+n, upper+n) > 0) ascend = 1;
      else {
	ascend = 0;
	/* Update the (d-n)-th power sum and related quantities. */
	tq = fmpq_mat_entry(dy_data->power_sums, d-n, 0);
	fmpq_sub(tq, tq, st_data->f+n);
	if ((d-n)%2==0)
	  fmpq_submul(fmpq_mat_entry(dy_data->hankel_dets, (d-n)/2, 0),
		      st_data->f+n, fmpq_mat_entry(dy_data->hankel_dets, (d-n)/2-1, 0));
	for (j=0; j<=d-n; j++) {
	  tq = fmpq_mat_entry(dy_data->hausdorff_sums1, d-n, j);
	  fmpq_sub(tq, tq, st_data->f+n);
	}
	test_roots = 1;
      }
    }
    count_steps += 1;
    if (count_steps > max_steps) {
      t = -2;
      break;
      }
  }
  dy_data->ascend = ascend;
  dy_data->n = n;
  dy_data->node_count = count;
  dy_data->flag = t;
}
