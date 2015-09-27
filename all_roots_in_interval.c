#include <fmpz_poly.h>

/* Based on code by Sebastian Pancratz from the FLINT repository.
 */

#define SWAP(f, g)  do { t = f; f = g; g = t; } while (0)

/*
    Assumes that:
        - {poly, n} is a normalized vector with n >= 2
        - {w, 3 * n + 9} is scratch space
 */
int _fmpz_poly_all_roots_in_interval(fmpz *poly, slong n, 
                                     fmpz const * a, fmpz const * b, fmpz *w)
{
    fmpz *f0     = w + 0 * n;
    fmpz *f1     = w + 1 * n;
    fmpz *f2     = w + 2 * n;
    fmpz *val0_a = w + 3 * n;
    fmpz *val0_b = w + 3 * n + 1;
    fmpz *val1_a = w + 3 * n + 2;
    fmpz *val1_b = w + 3 * n + 3;
    fmpz *c      = w + 3 * n + 4;
    fmpz *d      = w + 3 * n + 5;
    fmpz *t1     = w + 3 * n + 6;
    fmpz *t2     = w + 3 * n + 7;
    fmpz *t3     = w + 3 * n + 8;

    fmpz *l0 = f0 + (n - 1);
    fmpz *l1 = f1 + (n - 2);

    fmpz *t;

    int sgn0_a;
    int sgn0_b;

    _fmpz_vec_set(f0, poly, n);
    _fmpz_poly_evaluate_fmpz(val0_a, f0, n, a);

    /* Remove all factors of val0_a */
    while (fmpz_is_zero(val0_a))
    {
        /* {t1, 2} is available */
        fmpz_one(t1 + 1);
        fmpz_neg(t1 + 0, a);

        _fmpz_poly_divrem(f1, f2, f0, n, t1, 2);
        SWAP(f0, f1);
        n--;

        _fmpz_poly_evaluate_fmpz(val0_a, f0, n, a);
    }

    _fmpz_poly_evaluate_fmpz(val0_b, f0, n, b);

    /* Remove all factors of val0_b */
    while (fmpz_is_zero(val0_b))
    {
        /* {t1, 2} is available */
        fmpz_one(t1 + 1);
        fmpz_neg(t1 + 0, b);

        _fmpz_poly_divrem(f1, f2, f0, n, t1, 2);
        SWAP(f0, f1);
        n--;

        _fmpz_poly_evaluate_fmpz(val0_b, f0, n, b);
    }

    if (n == 1)
        return 1;

    _fmpz_poly_derivative(f1, f0, n);
    n--;
    _fmpz_poly_evaluate_fmpz(val1_a, f1, n, a);
    _fmpz_poly_evaluate_fmpz(val1_b, f1, n, b);

    sgn0_a = fmpz_sgn(val0_a);
    sgn0_b = fmpz_sgn(val0_b);

    for ( ; ; )
    {
        /* Invariant:  n = len(f1) */

        /* If we miss any one sign change, we cannot have enough */
        sgn0_a = -sgn0_a;
        if (fmpz_sgn(val1_a) != sgn0_a || fmpz_sgn(val1_b) != sgn0_b)
            return 0;

        /* 
            Explicitly compute the pseudoremainder f2 of f0 modulo f1
            f2 := - l1 * f0 + l0 * x * f1
            c  := - f2[n-1]
            f2 := l1 * f2 + c * f1
         */
        fmpz_zero(f2 + 0);
        _fmpz_vec_scalar_mul_fmpz(f2 + 1, f1, n - 1, l0);
        _fmpz_vec_scalar_submul_fmpz(f2, f0, n, l1);
        fmpz_neg(c, f2 + n - 1);

        _fmpz_vec_scalar_mul_fmpz(f2, f2, n - 1, l1);
        _fmpz_vec_scalar_addmul_fmpz(f2, f1, n - 1, c);

        if (_fmpz_vec_is_zero(f2, n - 1))
            return 1;

        n--;

        /* Cannot have enough sign changes if the degree drops more than 1 */
        if (fmpz_is_zero(f2 + n - 1))
            return 0;

        /* Extract content from f2; in practice, this seems to do better than
	 an explicit subresultant computation. */
        _fmpz_vec_content(d, f2, n);

        /* Evaluate f2 at a and b without an explicit function call */

        /* val2_a = (c*val1_a + lead1*(lead0*val1_a*a - lead1*val0_a)) // d */
        fmpz_mul(t1, val1_a, a);
        fmpz_mul(t2, l0, t1);
        fmpz_submul(t2, l1, val0_a);
        fmpz_mul(t1, c, val1_a);
        fmpz_addmul(t1, l1, t2);
        fmpz_swap(val0_a, val1_a);
        fmpz_divexact(val1_a, t1, d);

        /* val2_b = (c*val1_b + lead1*(lead0*val1_b*b - lead1*val0_b)) // d */
        fmpz_mul(t1, val1_b, b);
        fmpz_mul(t2, l0, t1);
        fmpz_submul(t2, l1, val0_b);
        fmpz_mul(t1, c, val1_b);
        fmpz_addmul(t1, l1, t2);
        fmpz_swap(val0_b, val1_b);
        fmpz_divexact(val1_b, t1, d);

        /* Rotate the polynomials, leading coefficients, and values */
        _fmpz_vec_scalar_divexact_fmpz(f0, f2, n, d);
        SWAP(f0, f1);

        l0 = f0 + n;
        l1 = f1 + n - 1;
    }

    return 1;
}