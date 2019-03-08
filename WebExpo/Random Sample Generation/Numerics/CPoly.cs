namespace Zygotine.Numerics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Numerics;
    using Zygotine.Statistics.Util;
    using Zygotine.Util;

    public class CPoly
    {
        /* 
         * Formerly src/appl/cpoly.c:
         *
         *  Copyright (C) 1997-1998 Ross Ihaka
         *  Copyright (C) 1999-2001 R Core Team
         *
         *	cpoly finds the zeros of a complex polynomial.
         *
         *	On Entry
         *
         *	opr, opi      -	 double precision vectors of real and
         *			 imaginary parts of the coefficients in
         *			 order of decreasing powers.
         *
         *	degree	      -	 int degree of polynomial.
         *
         *
         *	On Return
         *
         *	zeror, zeroi  -	 output double precision vectors of
         *			 real and imaginary parts of the zeros.
         *
         *	fail	      -	 output int parameter,	true  only if
         *			 leading coefficient is zero or if cpoly
         *			 has found fewer than degree zeros.
         *
         *	The program has been written to reduce the chance of overflow
         *	occurring. If it does occur, there is still a possibility that
         *	the zerofinder will work provided the overflowed quantity is
         *	replaced by a large number.
         *
         *	This is a C translation of the following.
         *
         *	TOMS Algorithm 419
         *	Jenkins and Traub.
         *	Comm. ACM 15 (1972) 97-99.
         *
         *	Ross Ihaka
         *	February 1997
         */

        //# include <R_ext/Arith.h> /* for declaration of hypot */
        //# include <R_ext/Memory.h> /* for declaration of R_alloc */

        //# include <float.h> /* for FLT_RADIX */

        //# include <Rmath.h> /* for R_pow_di */

        private int nn;
        private double[] pr, pi, hr, hi, qpr, qpi, qhr, qhi, shr, shi;
        private double sr, si;
        private double tr, ti;
        private double pvr, pvi;

        const double eta = RVaria.DBL_EPSILON;
        const double are = /* eta = */ RVaria.DBL_EPSILON;
        const double mre = 2.0 * RVaria.M_SQRT2 * /* eta, i.e. */ RVaria.DBL_EPSILON;
        const double infin = RVaria.DBL_MAX;

        private CPoly()
        {
        }

        internal static List<Complex> GetRoots(Complex[] op)
        {
            int degree = op.Length - 1;
            double[] opr = op.Real();
            double[] opi = op.Imaginary();
            double[] zeror = new double[degree];
            double[] zeroi = new double[degree];

            CPoly cpoly = new CPoly();
            bool fail = false;
            cpoly.R_cpolyroot(opr, opi, degree, zeror, zeroi, out fail);
            List<Complex> rep = new List<Complex>();
            if(!fail)
            {
                for(int i=0; i<degree; i++)
                {
                    rep.Add( new Complex(zeror[i], zeroi[i]));
                }
            }
            //en cas d'échec la liste sera vide.
            return rep;
        }

        internal static List<Complex> GetRoots(double[] opR)
        {
            return GetRoots(opR.Select(r => new Complex(r, 0)).ToArray());
        }

        private void R_cpolyroot(double[] opr, double[] opi, int degree, double[] zeror, double[] zeroi, out bool fail)
        {
            const double smalno = RVaria.DBL_MIN;
            const double base_ = (double)RVaria.FLT_RADIX;

            // R_cpolyroot variables ...
            int d_n, i, i1, i2;
            double zr = 0, zi = 0, xx, yy;
            double bnd, xxx;

            bool conv;
            int d1;
            const double cosr =/* cos 94 */ -0.06975647374412529990;
            const double sinr =/* sin 94 */  0.99756405025982424767;

            xx = RVaria.M_SQRT1_2;/* 1/Math.Sqrt(2) = 0.707.... */
            yy = -xx;
            fail = false;
            nn = degree;
            d1 = nn - 1;

            /* algorithm fails if the leading coefficient is zero. */

            if ((opr[0] == 0) && (opi[0] == 0))
            {
                fail = true;
                return;
            }

            /* remove the zeros at the origin if any. */

            while ((opr[nn] == 0.0) && (opi[nn] == 0.0))
            {
                d_n = d1 - nn + 1;
                zeror[d_n] = 0.0;
                zeroi[d_n] = 0.0;
                nn--;
            }

            nn++;
            /*-- Now, global var.  nn := #{coefficients} = (relevant degree)+1 */

            if (nn == 1)
            {
                return;
            }

            /* Use a single allocation as these as small */
            //const void* vmax = vmaxget();
            //tmp = new double[nn]; 
            pr = new double[nn];
            pi = new double[nn];

            hr = new double[nn];
            hi = new double[nn];

            qpr = new double[nn];
            qpi = new double[nn];

            qhr = new double[nn];
            qhi = new double[nn];

            shr = new double[nn];
            shi = new double[nn];

            /* make a copy of the coefficients and shr[] = | p[] | */
            for (i = 0; i < nn; i++)
            {
                pr[i] = opr[i];
                pi[i] = opi[i];
                shr[i] = RVaria.hypot(pr[i], pi[i]);
            }

            /* scale the polynomial with factor 'bnd'. */
            bnd = cpoly_scale(nn, shr, eta, infin, smalno, base_);
            if (bnd != 1.0)
            {
                for (i = 0; i < nn; i++)
                {
                    pr[i] *= bnd;
                    pi[i] *= bnd;
                }
            }

            /* start the algorithm for one zero */

            while (nn > 2)
            {
                /* calculate bnd, a lower bound on the modulus of the zeros. */
                for (i = 0; i < nn; i++)
                {
                    shr[i] = RVaria.hypot(pr[i], pi[i]);
                }

                bnd = cpoly_cauchy(nn, shr, shi);

                /* outer loop to control 2 major passes */
                /* with different sequences of shifts */

                for (i1 = 1; i1 <= 2; i1++)
                {
                    /* first stage calculation, no shift */
                    noshft(5);

                    /*	inner loop to select a shift */
                    for (i2 = 1; i2 <= 9; i2++)
                    {
                        /* shift is chosen with modulus bnd */
                        /* and amplitude rotated by 94 degrees */
                        /* from the previous shift */

                        xxx = cosr * xx - sinr * yy;
                        yy = sinr * xx + cosr * yy;
                        xx = xxx;
                        this.sr = bnd * xx;
                        this.si = bnd * yy;

                        /*  second stage calculation, fixed shift */
                        conv = fxshft(i2 * 10, ref zr, ref zi);
                        if (conv)
                        {
                            goto L10;
                        }
                    }
                }

                /* the zerofinder has failed on two major passes */
                /* return empty handed */
                fail = true;
                return;

                /* the second stage jumps directly to the third stage iteration.
                 * if successful, the zero is stored and the polynomial deflated.
                 */
                L10:
                d_n = d1 + 2 - nn;
                zeror[d_n] = zr;
                zeroi[d_n] = zi;
                --nn;

                for (i = 0; i < nn; i++)
                {
                    pr[i] = qpr[i];
                    pi[i] = qpi[i];
                }
            }/*while*/

            /*	calculate the final zero and return */
            cdivid(-pr[1], -pi[1], pr[0], pi[0], out zeror[d1], out zeroi[d1]);
            return;
        }


        /*  Computes the derivative polynomial as the initial
         *  polynomial and computes l1 no-shift h polynomials.	*/
        void noshft(int l1)
        {
            int i, j, jj, n = nn - 1, nm1 = n - 1;
            double t1, t2, xni;
            for (i = 0; i < n; i++)
            {
                xni = (double)(nn - i - 1);
                this.hr[i] = xni * this.pr[i] / n;
                this.hi[i] = xni * this.pi[i] / n;
            }

            for (jj = 1; jj <= l1; jj++)
            {
                if (RVaria.hypot(this.hr[n - 1], this.hi[n - 1]) <= eta * 10.0 * RVaria.hypot(this.pr[n - 1], this.pi[n - 1]))
                {
                    /*	If the constant term is essentially zero, */
                    /*	shift h coefficients. */
                    for (i = 1; i <= nm1; i++)
                    {
                        j = this.nn - i;
                        this.hr[j - 1] = this.hr[j - 2];
                        this.hi[j - 1] = this.hi[j - 2];
                    }

                    this.hr[0] = 0.0;
                    this.hi[0] = 0.0;
                }
                else
                {
                    cdivid(-pr[nn - 1], -pi[nn - 1], hr[n - 1], hi[n - 1], out tr, out ti);
                    for (i = 1; i <= nm1; i++)
                    {
                        j = nn - i;
                        t1 = hr[j - 2];
                        t2 = hi[j - 2];
                        hr[j - 1] = tr * t1 - ti * t2 + pr[j - 1];
                        hi[j - 1] = tr * t2 + ti * t1 + pi[j - 1];
                    }
                    hr[0] = pr[0];
                    hi[0] = pi[0];
                }
            }
        }

        /*  
         *  Computes l2 fixed-shift h polynomials and tests for convergence.
         *  initiates a variable-shift iteration and returns with the
         *  approximate zero if successful.
         *  
         */
        bool fxshft(int l2, ref double zr, ref double zi)
        {
            /*  
             *  l2	  - limit of fixed shift steps
             *  zr,zi - approximate zero if convergence (result TRUE)
             *
             * Return value indicates convergence of stage 3 iteration
             *
             * Uses global (sr,si), nn, pr[], pi[], .. (all args of polyev() !)
             * 
             */
            bool pasd, bool_, test;
            double svsi, svsr;
            int i, j, n;
            double oti, otr;
            n = nn - 1;

            /* evaluate p at s. */
            polyev(nn, sr, si, pr, pi, qpr, qpi, out pvr, out pvi);
            test = true;
            pasd = false;

            /* calculate first t = -p(s)/h(s). */
            calct(out bool_);

            /* main loop for one second stage step. */
            for (j = 1; j <= l2; j++)
            {
                otr = tr;
                oti = ti;

                /* compute next h polynomial and new t. */
                nexth(bool_);
                calct(out bool_);
                zr = sr + tr;
                zi = si + ti;

                /* test for convergence unless stage 3 has */
                /* failed once or this is the last h polynomial. */
                if (!bool_ && test && (j != l2))
                {
                    if (RVaria.hypot(tr - otr, ti - oti) >= RVaria.hypot(zr, zi) * 0.5)
                    {
                        pasd = false;
                    }
                    else if (!pasd)
                    {
                        pasd = true;
                    }
                    else
                    {

                        /* the weak convergence test has been */
                        /* passed twice, start the third stage */
                        /* iteration, after saving the current */
                        /* h polynomial and shift. */
                        for (i = 0; i < n; i++)
                        {
                            shr[i] = hr[i];
                            shi[i] = hi[i];
                        }
                        svsr = sr;
                        svsi = si;
                        if (vrshft(10, ref zr, ref zi))
                        {
                            return true;
                        }

                        /* the iteration failed to converge. */
                        /* turn off testing and restore */
                        /* h, s, pv and t. */

                        test = false;
                        for (i = 1; i <= n; i++)
                        {
                            hr[i - 1] = shr[i - 1];
                            hi[i - 1] = shi[i - 1];
                        }

                        sr = svsr;
                        si = svsi;
                        polyev(nn, sr, si, pr, pi, qpr, qpi, out pvr, out pvi);
                        calct(out bool_);
                    }
                }
            }

            /* attempt an iteration with final h polynomial */
            /* from second stage. */
            return (vrshft(10, ref zr, ref zi));
        }


        /* 
         * carries out the third stage iteration.
         * 
         */
        bool vrshft(int l3, ref double zr, ref double zi)
        {
            /*  l3	    - limit of steps in stage 3.
             *  zr,zi   - on entry contains the initial iterate;
             *	      if the iteration converges it contains
             *	      the final iterate on exit.
             * Returns TRUE if iteration converges
             *
             * Assign and uses  GLOBAL sr, si
             */
            bool bool_, b;
            int i, j;
            double r1, r2, mp, ms, tp, relstp = 0.0;
            double omp = 0.0;

            b = false;
            sr = zr;
            si = zi;

            /* main loop for stage three */
            for (i = 1; i <= l3; i++)
            {
                /* evaluate p at s and test for convergence. */
                polyev(nn, sr, si, pr, pi, qpr, qpi, out pvr, out pvi);
                mp = RVaria.hypot(pvr, pvi);
                ms = RVaria.hypot(sr, si);
                if (mp <= 20.0 * errev(nn, qpr, qpi, ms, mp, /*are=*/eta, mre))
                {
                    goto L_conv;
                }

                /* 
                 * polynomial value is smaller in value than 
                 * a bound on the error in evaluating p, 
                 * terminate the iteration. 
                 * 
                 */
                if (i != 1)
                {
                    if (!b && (mp >= omp) && (relstp < .05))
                    {
                        /* 
                         * iteration has stalled. probably a 
                         * cluster of zeros. do 5 fixed shift 
                         * steps into the cluster to force 
                         * one zero to dominate.
                         * 
                         */
                        tp = relstp;
                        b = true;
                        if (relstp < eta)
                        {
                            tp = eta;
                        }

                        r1 = Math.Sqrt(tp);
                        r2 = sr * (r1 + 1.0) - si * r1;
                        si = sr * r1 + si * (r1 + 1.0);
                        sr = r2;
                        polyev(nn, sr, si, pr, pi, qpr, qpi, out pvr, out pvi);
                        for (j = 1; j <= 5; ++j)
                        {
                            calct(out bool_);
                            nexth(bool_);
                        }
                        omp = infin;
                        goto L10;
                    }
                    else
                    {
                        /* exit if polynomial value */
                        /* increases significantly. */
                        if (mp * .1 > omp)
                        {
                            return false;
                        }
                    }
                }

                omp = mp;

                /* calculate next iterate. */
                L10:
                calct(out bool_);
                nexth(bool_);
                calct(out bool_);
                if (!bool_)
                {
                    relstp = RVaria.hypot(tr, ti) / RVaria.hypot(sr, si);
                    sr += tr;
                    si += ti;
                }
            }

            return false;
            L_conv:
            zr = sr;
            zi = si;
            return true;
        }

        void calct(out bool bool_)
        {
            /* 
             * computes	 t = -p(s)/h(s).
             * bool   - logical, set true if h(s) is essentially zero.
             * 
             */
            int n = nn - 1;
            double hvi, hvr;

            /* evaluate h(s). */
            polyev(n, sr, si, hr, hi, qhr, qhi, out hvr, out hvi);
            bool_ = RVaria.hypot(hvr, hvi) <= are * 10.0 * RVaria.hypot(hr[n - 1], hi[n - 1]);
            if (!bool_)
            {
                cdivid(-pvr, -pvi, hvr, hvi, out tr, out ti);
            }
            else
            {
                tr = 0.0;
                ti = 0.0;
            }
        }

        void nexth(bool bool_)
        {
            /* 
             * calculates the next shifted h polynomial.
             * bool :	if TRUE  h(s) is essentially zero
             * 
             */
            int j, n = nn - 1;
            double t1, t2;

            if (!bool_)
            {
                for (j = 1; j < n; j++)
                {
                    t1 = qhr[j - 1];
                    t2 = qhi[j - 1];
                    hr[j] = tr * t1 - ti * t2 + qpr[j];
                    hi[j] = tr * t2 + ti * t1 + qpi[j];
                }

                hr[0] = qpr[0];
                hi[0] = qpi[0];
            }
            else
            {
                /* if h(s) is zero replace h with qh. */
                for (j = 1; j < n; j++)
                {
                    hr[j] = qhr[j - 1];
                    hi[j] = qhi[j - 1];
                }

                hr[0] = 0.0;
                hi[0] = 0.0;
            }
        }

        /*--------------------- Independent Complex Polynomial Utilities ----------*/

        void polyev(int n,
                double s_r, double s_i,
                double[] p_r, double[] p_i,
                double[] q_r, double[] q_i,
                out double v_r, out double v_i)
        {
            /* 
             * evaluates a polynomial  p  at  s	 by the horner recurrence
             * placing the partial sums in q and the computed value in v_.
             * 
             */
            int i;
            double t;
            q_r[0] = p_r[0];
            q_i[0] = p_i[0];
            v_r = q_r[0];
            v_i = q_i[0];
            for (i = 1; i < n; i++)
            {
                t = v_r * s_r - v_i * s_i + p_r[i];
                q_i[i] = v_i = v_r * s_i + v_i * s_r + p_i[i];
                q_r[i] = v_r = t;
            }
        }

        double errev(int n, double[] qr, double[] qi, double ms, double mp, double a_re, double m_re)
        {
            /*	
             *	bounds the error in evaluating the polynomial by the horner
             *	recurrence.
             *
             *	qr,qi	 - the partial sum vectors
             *	ms	 - modulus of the point
             *	mp	 - modulus of polynomial value
             *  a_re,m_re - error bounds on complex addition and multiplication
             *  
             */

            double e;
            int i;
            e = RVaria.hypot(qr[0], qi[0]) * m_re / (a_re + m_re);
            for (i = 0; i < n; i++)
            {
                e = e * ms + RVaria.hypot(qr[i], qi[i]);
            }

            return e * (a_re + m_re) - mp * m_re;
        }

        double cpoly_cauchy(int n, double[] pot, double[] q)
        {
            /* 
             * Computes a lower bound on the moduli of the zeros of a polynomial
             * pot[1:nn] is the modulus of the coefficients.
             * 
             */
            double f, x, delf, dx, xm;
            int i, n1 = n - 1;
            pot[n1] = -pot[n1];

            /* compute upper estimate of bound. */
            x = Math.Exp((Math.Log(-pot[n1]) - Math.Log(pot[0])) / (double)n1);

            /* if newton step at the origin is better, use it. */
            if (pot[n1 - 1] != 0.0)
            {
                xm = -pot[n1] / pot[n1 - 1];
                if (xm < x)
                {
                    x = xm;
                }
            }

            /* chop the interval (0,x) unitl f le 0.0 */
            for (;;)
            {
                xm = x * 0.1;
                f = pot[0];
                for (i = 1; i < n; i++)
                {
                    f = f * xm + pot[i];
                }

                if (f <= 0.0)
                {
                    break;
                }

                x = xm;
            }

            dx = x;

            /* do Newton iteration until x converges to two decimal places. */
            while (Math.Abs(dx / x) > 0.005)
            {
                q[0] = pot[0];
                for (i = 1; i < n; i++)
                {
                    q[i] = q[i - 1] * x + pot[i];
                }

                f = q[n1];
                delf = q[0];
                for (i = 1; i < n1; i++)
                {
                    delf = delf * x + q[i];
                }

                dx = f / delf;
                x -= dx;
            }

            return x;
        }

        double cpoly_scale(int n, double[] pot,
                   double eps, double BIG, double small, double base_)
        {
            /* 
             * Returns a scale factor to multiply the coefficients of the polynomial.
             * The scaling is done to avoid overflow and to avoid
             *	undetected underflow interfering with the convergence criterion.
             * The factor is a power of the base_.

             * pot[1:n] : modulus of coefficients of p
             * eps,BIG,
             * small,base_ - constants describing the floating point arithmetic.
             * 
             */

            int i, ell;
            double x, high, sc, lo, min_, max_;

            /* find largest and smallest moduli of coefficients. */
            high = Math.Sqrt(BIG);
            lo = small / eps;
            max_ = 0.0;
            min_ = BIG;
            for (i = 0; i < n; i++)
            {
                x = pot[i];
                if (x > max_)
                {
                    max_ = x;
                }

                if ((x != 0.0) && (x < min_))
                {
                    min_ = x;
                }
            }

            /* scale only if there are very large or very small components. */

            if (min_ < lo || max_ > high)
            {
                x = lo / min_;
                if (x <= 1.0)
                {
                    sc = 1.0 / (Math.Sqrt(max_) * Math.Sqrt(min_));
                }
                else
                {
                    sc = x;
                    if (BIG / sc > max_)
                    {
                        sc = 1.0;
                    }
                }

                ell = (int)(Math.Log(sc) / Math.Log(base_) + 0.5);
                return DiGamma.R_pow_di(base_, ell);
            }
            else
            {
                return 1.0;
            }
        }



        void cdivid(double ar, double ai, double br, double bi, out double cr, out double ci)
        {

            /* 
             * complex division c = a/b, i.e., (cr +i*ci) = (ar +i*ai) / (br +i*bi), avoiding overflow. 
             * 
             */
            double d, r;
            if (br == 0.0 && bi == 0.0)
            {
                /* division by zero, c = infinity. */
                cr = ci = double.PositiveInfinity;
            }
            else if (Math.Abs(br) >= Math.Abs(bi))
            {
                r = bi / br;
                d = br + r * bi;
                cr = (ar + ai * r) / d;
                ci = (ai - ar * r) / d;
            }
            else
            {
                r = br / bi;
                d = bi + r * br;
                cr = (ar * r + ai) / d;
                ci = (ai * r - ar) / d;
            }
        }


        /* 
         * static double cpoly_cmod(double *r, double *i)
         * --> replaced by hypot() everywhere
         * 
        */
    }
}
