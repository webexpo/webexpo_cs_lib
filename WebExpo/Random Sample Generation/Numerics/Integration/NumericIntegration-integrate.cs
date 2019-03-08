using System;
using Zygotine.Statistics.Util;

namespace Zygotine.Numerics.Integration
{
    public partial class NumericIntegration
    {
        /* R file: integrate.c
 */
        private static void rdqagse(Parameters parameters, Workspace ws, Results results)
        {
            /* Local variables */
            bool noext, extrap;
            bool goto140 = false;
            int k, ksgn, nres;
            int ierro;
            int ktmin, nrmax;
            int iroff1, iroff2, iroff3;
            int id;
            int numrl2;
            int jupbnd;
            int maxerr;
            int last;
            int ier;
            int neval;
            int rdqk21_count = 0;
            double[] res3la = new double[4];
            double[] rlist2 = new double[53];
            double
                result = 0.0,
                abserr = 0.0,
                abseps,
                area,
                area1,
                area2,
                area12,
                dres,
                epmach;
            double
                a1,
                a2,
                b1,
                b2,
                defabs,
                defab1,
                defab2,
                oflow,
                uflow,
                resabs,
                reseps;
            double
                error1,
                error2,
                erro12,
                errbnd,
                erlast,
                errmax,
                errsum;

            double
                correc = 0.0,
                erlarg = 0.0,
                ertest = 0.0,
                small = 0.0;

            //Let's go
            epmach = RVaria.DBL_EPSILON;
            ier = 0;
            neval = 0;
            last = 0;

            setArray(ws.alist);
            setArray(ws.blist);
            setArray(ws.rlist);
            setArray(ws.elist);
            setArray(ws.iord);


            ws.alist[1] = parameters.a;
            ws.blist[1] = parameters.b;
            ws.rlist[1] = 0.0;
            ws.elist[1] = 0.0;

            //TODO déplacer(?) dans compute
            if ((parameters.epsabs <= 0.0) && (parameters.epsrel < RVaria.fmax2(epmach * 50.0, 5e-29)))
            {
                ier = 6;
                return;
            }
            /*           first approximation to the integral */
            /*           ----------------------------------- */

            uflow = RVaria.DBL_MIN;
            oflow = RVaria.DBL_MAX;
            ierro = 0;

            rdqk21(parameters.f, parameters.a, parameters.b, out result, out abserr, out defabs, out resabs, ref rdqk21_count);

            /*           test on accuracy. */

            dres = Math.Abs(result);
            errbnd = RVaria.fmax2(parameters.epsabs, parameters.epsrel * dres);
            last = 1;
            ws.rlist[1] = result;
            ws.elist[1] = abserr;
            ws.iord[1] = 1;
            if ((abserr <= (epmach * 100.0 * defabs)) && (abserr > errbnd))
                ier = 2;
            if (ws.limit == 1)
                ier = 1;
            if ((ier != 0) || (abserr <= errbnd && abserr != resabs)
            || (abserr == 0.0))
            {
                goto140 = true;
                goto L139;
            }

            /*           initialization */
            /*           -------------- */

            rlist2[0] = result;
            errmax = abserr;
            maxerr = 1;
            area = result;
            errsum = abserr;
            abserr = oflow;
            nrmax = 1;
            nres = 0;
            numrl2 = 2;
            ktmin = 0;
            extrap = false;
            noext = false;
            iroff1 = 0;
            iroff2 = 0;
            iroff3 = 0;
            ksgn = -1;
            if (dres >= (1.0 - epmach * 50.0) * defabs)
            {
                ksgn = 1;
            }

            /*------------------------*/


            for (last = 2; last <= parameters.limit; ++(last))
            {

                /*           bisect the subinterval with the nrmax-th largest error estimate. */

                a1 = ws.alist[maxerr];
                b1 = (ws.alist[maxerr] + ws.blist[maxerr]) * .5;
                a2 = b1;
                b2 = ws.blist[maxerr];
                erlast = errmax;
                rdqk21(parameters.f, a1, b1, out area1, out error1, out resabs, out defab1, ref rdqk21_count);
                rdqk21(parameters.f, a2, b2, out area2, out error2, out resabs, out defab2, ref rdqk21_count);

                /*           improve previous approximations to integral
                         and error and test for accuracy. */

                area12 = area1 + area2;
                erro12 = error1 + error2;
                errsum = errsum + erro12 - errmax;
                area = area + area12 - ws.rlist[maxerr];
                if (!(defab1 == error1 || defab2 == error2))
                {

                    if (Math.Abs(ws.rlist[maxerr] - area12) <= Math.Abs(area12) * 1e-5 &&
                    erro12 >= errmax * .99)
                    {
                        if (extrap)
                            ++iroff2;
                        else /* if(! extrap) */
                            ++iroff1;
                    }
                    if (last > 10 && erro12 > errmax)
                        ++iroff3;
                }
                ws.rlist[maxerr] = area1;
                ws.rlist[last] = area2;
                errbnd = RVaria.fmax2(parameters.epsabs, parameters.epsrel * Math.Abs(area));

                /*           test for roundoff error and eventually set error flag. */

                if (iroff1 + iroff2 >= 10 || iroff3 >= 20)
                    ier = 2;
                if (iroff2 >= 5)
                    ierro = 3;

                /* set error flag in the case that the number of subintervals equals limit. */
                if (last == ws.limit)
                    ier = 1;

                /*           set error flag in the case of bad integrand behaviour
                         at a point of the integration range. */

                if (RVaria.fmax2(Math.Abs(a1), Math.Abs(b2)) <=
                    (epmach * 100.0 + 1.0) * (Math.Abs(a2) + uflow * 1e3))
                {
                    ier = 4;
                }

                /*           append the newly-created intervals to the list. */

                if (error2 > error1)
                {
                    ws.alist[maxerr] = a2;
                    ws.alist[last] = a1;
                    ws.blist[last] = b1;
                    ws.rlist[maxerr] = area2;
                    ws.rlist[last] = area1;
                    ws.elist[maxerr] = error2;
                    ws.elist[last] = error1;
                }
                else
                {
                    ws.alist[last] = a2;
                    ws.blist[maxerr] = b1;
                    ws.blist[last] = b2;
                    ws.elist[maxerr] = error1;
                    ws.elist[last] = error2;
                }

                /*           call subroutine dqpsrt to maintain the descending ordering
                         in the list of error estimates and select the subinterval
                         with nrmax-th largest error estimate (to be bisected next). */

                /*L30:*/

                rdqpsrt(ws.limit, last, ref maxerr, out errmax, ws.elist, ws.iord, ref nrmax);
                if (errsum <= errbnd)
                    goto L115;/* ***jump out of do-loop */
                if (ier != 0) break;
                if (last == 2)
                { /* L80: */
                    small = Math.Abs(parameters.b - parameters.a) * .375;
                    erlarg = errsum;
                    ertest = errbnd;
                    rlist2[1] = area; continue;
                }

                if (noext) continue;
                erlarg -= erlast;
                if (Math.Abs(b1 - a1) > small)
                {
                    erlarg += erro12;
                }
                if (!extrap)
                {

                    /*          test whether the interval to be bisected next is the
                            smallest interval. */

                    if (Math.Abs(ws.blist[maxerr] - ws.alist[maxerr]) > small)
                    {
                        continue;
                    }
                    extrap = true;
                    nrmax = 2;
                }

                if (ierro != 3 && erlarg > ertest)
                {

                    /*           the smallest interval has the largest error.
                             before bisecting decrease the sum of the errors over the
                             larger intervals (erlarg) and perform extrapolation. */

                    id = nrmax;
                    jupbnd = last;
                    if (last > parameters.limit / 2 + 2)
                    {
                        jupbnd = parameters.limit + 3 - last;
                    }
                    for (k = id; k <= jupbnd; ++k)
                    {
                        maxerr = ws.iord[nrmax];
                        errmax = ws.elist[maxerr];
                        if (Math.Abs(ws.blist[maxerr] - ws.alist[maxerr]) > small)
                        {
                            goto L90;
                        }
                        ++nrmax;
                        /* L50: */
                    }
                }
                /*           perform extrapolation.  L60: */

                ++numrl2;
                rlist2[numrl2 - 1] = area;

                rdqelg(ref numrl2, rlist2, out reseps, out abseps, res3la, ref nres);


                ++ktmin;
                if (ktmin > 5 && abserr < errsum * .001)
                {
                    ier = 5;
                }
                if (abseps < abserr)
                {
                    ktmin = 0;
                    abserr = abseps;
                    result = reseps;
                    correc = erlarg;
                    ertest = RVaria.fmax2(parameters.epsabs, parameters.epsrel * Math.Abs(reseps));
                    if (abserr <= ertest)
                    {
                        break;
                    }
                }

                /*           prepare bisection of the smallest interval.  L70: */

                if (numrl2 == 1)
                {
                    noext = true;
                }
                if (ier == 5)
                {
                    break;
                }
                maxerr = ws.iord[1];
                errmax = ws.elist[maxerr];
                nrmax = 1;
                extrap = false;
                small *= 0.5;
                erlarg = errsum;
            L90:
                ;
            }

            /* L100:	set final result and error estimate. */
            /*		------------------------------------ */

            if (abserr == oflow) goto L115;
            if (ier + ierro == 0) goto L110;
            if (ierro == 3)
                abserr += correc;
            if (ier == 0)
                ier = 3;
            if (result == 0.0 || area == 0.0)
            {
                if (abserr > errsum) goto L115;
                if (area == 0.0) goto L130;
            }
            else
            { /* L105:*/
                if (abserr / Math.Abs(result) > errsum / Math.Abs(area))
                    goto L115;
            }

        L110:/*		test on divergence. */
            if (ksgn == -1 && RVaria.fmax2(Math.Abs(result), Math.Abs(area)) <= defabs * .01)
            {
                goto L130;
            }
            if (.01 > result / area || result / area > 100.0 || errsum > Math.Abs(area))
            {
                ier = 5;
            }
            goto L130;

        L115:/*		compute global integral sum. */
            result = 0.0;
            for (k = 1; k <= last; ++k)
                result += ws.rlist[k];
            abserr = errsum;
        L130:
        L139:
            if ((ier > 2) || goto140)
            {
            /*L140:*/
                neval = last * 42 - 21;
            }
            results.AbsErr = abserr;
            results.IEr = ier;
            results.Last = last;
            results.NEval = neval;
            results.Value = result;
            return;

            // ---------------------------*/


        }

        /* R file: integrate.c
         */
        private static void rdqk21(Func<double[], double[]> f, double a, double b, out double result, out double abserr, out double resabs, out double resasc, ref int rdqk21_count)
        {
            //static void  rdqk21(integr_fn f, void *ex, double *a, double *b, double *result,
            // double *abserr, double *resabs, double *resasc)
            double[]
                fv1 = new double[10],
                fv2 = new double[10],
                vec = new double[21];
            double absc, resg, resk, fsum, fval1, fval2;
            double hlgth, centr, reskh, uflow;
            double fc, epmach, dhlgth;
            int j, jtw, jtwm1;

            rdqk21_count++;
            epmach = RVaria.DBL_EPSILON;
            uflow = RVaria.DBL_MIN;

            centr = (a + b) * 0.5;
            hlgth = (b - a) * 0.5;
            dhlgth = Math.Abs(hlgth);

            /*           compute the 21-point kronrod approximation to
	     the integral, and estimate the absolute error. */

            resg = 0.0;
            vec[0] = centr;
            for (j = 1; j <= 5; ++j)
            {
                jtw = j << 1;
                absc = hlgth * xgk[jtw - 1];
                vec[(j << 1) - 1] = centr - absc;
            /*L5:*/
                vec[j * 2] = centr + absc;
            }

            for (j = 1; j <= 5; ++j)
            {
                jtwm1 = (j << 1) - 1;
                absc = hlgth * xgk[jtwm1 - 1];
                vec[(j << 1) + 9] = centr - absc;
                vec[(j << 1) + 10] = centr + absc;
            }

            double[] vec_y = f(vec);//, 21);

            fc = vec_y[0];
            resk = wgk[10] * fc;
            resabs = Math.Abs(resk);

            for (j = 1; j <= 5; ++j)
            {
                jtw = j << 1;
                absc = hlgth * xgk[jtw - 1];
                fval1 = vec_y[(j << 1) - 1];
                fval2 = vec_y[j * 2];
                fv1[jtw - 1] = fval1;
                fv2[jtw - 1] = fval2;
                fsum = fval1 + fval2;
                resg += wg[j - 1] * fsum;
                resk += wgk[jtw - 1] * fsum;
                resabs += wgk[jtw - 1] * (Math.Abs(fval1) + Math.Abs(fval2));
                //L10: 
            }


            for (j = 1; j <= 5; ++j)
            {
                jtwm1 = (j << 1) - 1;
                absc = hlgth * xgk[jtwm1 - 1];
                fval1 = vec_y[(j << 1) + 9];
                fval2 = vec_y[(j << 1) + 10];
                fv1[jtwm1 - 1] = fval1;
                fv2[jtwm1 - 1] = fval2;
                fsum = fval1 + fval2;
                resk += wgk[jtwm1 - 1] * fsum;
                resabs += wgk[jtwm1 - 1] * (Math.Abs(fval1) + Math.Abs(fval2));
                // L15:
            }


            reskh = resk * .5;
            resasc = wgk[10] * Math.Abs(fc - reskh);
            for (j = 1; j <= 10; ++j)
            {
                resasc += wgk[j - 1] * (Math.Abs(fv1[j - 1] - reskh) +
                             Math.Abs(fv2[j - 1] - reskh));
                /* L20: */
            }
            //vec[0] = 0;

            result = resk * hlgth;
            resabs *= dhlgth;
            resasc *= dhlgth;
            abserr = Math.Abs((resk - resg) * hlgth);
            if ((resasc != 0.0) && (abserr != 0.0))
            {
                abserr = resasc * RVaria.fmin2(1.0, Math.Pow(abserr * 200.0 / resasc, 1.5));
            }

            if (resabs > uflow / (epmach * 50.0))
            {
                abserr = RVaria.fmax2(epmach * 50.0 * resabs, abserr);
            }


            return;
            // throw new System.NotImplementedException();
        }

        /* R file: integrate.c
         */
        #region constants defined in rdqk21
        private static readonly double[] wg = new double[5]{
            .066671344308688137593568809893332,
            .149451349150580593145776339657697,
            .219086362515982043995534934228163,
            .269266719309996355091226921569469,
            .295524224714752870173892994651338
        };

        private static readonly double[] xgk = new double[11]{
            .995657163025808080735527280689003,
            .973906528517171720077964012084452,
            .930157491355708226001207180059508,
            .865063366688984510732096688423493,
            .780817726586416897063717578345042,
            .679409568299024406234327365114874,
            .562757134668604683339000099272694,
            .433395394129247190799265943165784,
            .294392862701460198131126603103866,
            .14887433898163121088482600112972,
            0.0
        };

        private static readonly double[] wgk = new double[11]{
            .011694638867371874278064396062192,
            .03255816230796472747881897245939,
            .05475589657435199603138130024458,
            .07503967481091995276704314091619,
            .093125454583697605535065465083366,
            .109387158802297641899210590325805,
            .123491976262065851077958109831074,
            .134709217311473325928054001771707,
            .142775938577060080797094273138717,
            .147739104901338491374841515972068,
            .149445554002916905664936468389821
        };
        #endregion

        /* R file: integrate.c
         */
        private static void rdqelg(ref int n, double[] epstab, out double result, out double abserr, double[] res3la, ref int nres)
        {
            /* Local variables */
            int i__, indx, ib, ib2, ie, k1, k2, k3, num, newelm, limexp;
            double delta1, delta2, delta3, e0, e1, e1abs, e2, e3, epmach, epsinf;
            double oflow, ss, res;
            double errA, err1, err2, err3, tol1, tol2, tol3;

            #region prologue
            /*  ***begin prologue  dqelg
		        ***refer to  dqagie,dqagoe,dqagpe,dqagse
		        ***revision date  830518   (yymmdd)
		        ***keywords  epsilon algorithm, convergence acceleration,
		        extrapolation
		        ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
		        de doncker,elise,appl. math & progr. div. - k.u.leuven
		        ***purpose  the routine determines the limit of a given sequence of
		        approximations, by means of the epsilon algorithm of
		        p.wynn. an estimate of the absolute error is also given.
		        the condensed epsilon table is computed. only those
		        elements needed for the computation of the next diagonal
		        are preserved.
		        ***description

		        epsilon algorithm
		        standard fortran subroutine
		        double precision version

		        parameters
		        n      - int
		        epstab(n) contains the new element in the
		        first column of the epsilon table.

		        epstab - double precision
		        vector of dimension 52 containing the elements
		        of the two lower diagonals of the triangular
		        epsilon table. the elements are numbered
		        starting at the right-hand corner of the
		        triangle.

		        result - double precision
		        resulting approximation to the integral

		        abserr - double precision
		        estimate of the absolute error computed from
		        result and the 3 previous results

		        res3la - double precision
		        vector of dimension 3 containing the last 3
		        results

		        nres   - int
		        number of calls to the routine
		        (should be zero at first call)

		        ***end prologue  dqelg


		        list of major variables
		        -----------------------

		        e0     - the 4 elements on which the computation of a new
		        e1       element in the epsilon table is based
		        e2
		        e3                 e0
		        e3    e1    new
		        e2

		        newelm - number of elements to be computed in the new diagonal
		        errA   - errA = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
		        result - the element in the new diagonal with least value of errA

		        machine dependent constants
		        ---------------------------

		        epmach is the largest relative spacing.
		        oflow is the largest positive magnitude.
		        limexp is the maximum number of elements the epsilon
		        table can contain. if this number is reached, the upper
		        diagonal of the epsilon table is deleted. */
            #endregion
            // -res3la;
            //--epstab;
            /* Function Body */
            epmach = RVaria.DBL_EPSILON;
            oflow = RVaria.DBL_MAX;
            ++nres;
            abserr = oflow;
            result = epstab[n - 1]; // modifié 
            if (n < 3)
            {
                goto L100;
            }
            limexp = 50;
            epstab[n + 2 - 1] = epstab[n - 1];// modifiés 
            newelm = (n - 1) / 2;
            epstab[n - 1] = oflow;// modifié 
            num = n;
            k1 = n;
            for (i__ = 1; i__ <= newelm; ++i__)
            {
                k2 = k1 - 1;
                k3 = k1 - 2;
                res = epstab[k1 + 2 - 1];// modifié 
                e0 = epstab[k3 - 1]; // modifié 
                e1 = epstab[k2 - 1];// modifié 
                e2 = res;
                e1abs = Math.Abs(e1);
                delta2 = e2 - e1;
                err2 = Math.Abs(delta2);
                tol2 = RVaria.fmax2(Math.Abs(e2), e1abs) * epmach;
                delta3 = e1 - e0;
                err3 = Math.Abs(delta3);
                tol3 = RVaria.fmax2(e1abs, Math.Abs(e0)) * epmach;
                if (err2 <= tol2 && err3 <= tol3)
                {
                    /*           if e0, e1 and e2 are equal to within machine
                    accuracy, convergence is assumed. */
                    result = res;/*		result = e2 */
                    abserr = err2 + err3;/*	abserr = fabs(e1-e0)+fabs(e2-e1) */

                    goto L100;	/* ***jump out of do-loop */
                }

                e3 = epstab[k1 - 1];// modifié 
                epstab[k1 - 1] = e1;// modifié 
                delta1 = e1 - e3;
                err1 = Math.Abs(delta1);
                tol1 = RVaria.fmax2(e1abs, Math.Abs(e3)) * epmach;

                /*           if two elements are very close to each other, omit
                a part of the table by adjusting the value of n */

                if (err1 > tol1 && err2 > tol2 && err3 > tol3)
                {
                    ss = 1.0 / delta1 + 1.0 / delta2 - 1.0 / delta3;
                    epsinf = Math.Abs(ss * e1);

                    /* 
                     * test to detect irregular behaviour in the table, and
                     * eventually omit a part of the table adjusting the value of n. 
                     */

                    if (epsinf > 1e-4)
                    {
                        goto L30;
                    }
                }

                n = i__ + i__ - 1;
                goto L50;
                /* ***jump out of do-loop */


                L30:/* compute a new element and eventually adjust the value of result. */

                res = e1 + 1.0 / ss;
                epstab[k1 - 1] = res; // modifié 
                k1 += -2;
                errA = err2 + Math.Abs(res - e2) + err3;
                if (errA <= abserr)
                {
                    abserr = errA;
                    result = res;
                }
            }

            /*           shift the table. */

            L50:
            if (n == limexp)
            {
                n = (limexp / 2 << 1) - 1;
            }

            if (num / 2 << 1 == num) ib = 2; else ib = 1;
            ie = newelm + 1;
            for (i__ = 1; i__ <= ie; ++i__)
            {
                ib2 = ib + 2;
                epstab[ib - 1] = epstab[ib2 - 1];// modifiés
                ib = ib2;
            }
            if (num != n)
            {
                indx = num - n + 1;
                for (i__ = 1; i__ <= n; ++i__)
                {
                    epstab[i__ - 1] = epstab[indx - 1]; // modifiés
                    ++indx;
                }
            }
            /*L80:*/
            if (nres >= 4)
            {
                /* L90: */
                abserr = Math.Abs(result - res3la[3 - 1]) +
                    Math.Abs(result - res3la[2 - 1]) +
                    Math.Abs(result - res3la[1 - 1]); //modifiés  3x
                res3la[1 - 1] = res3la[2 - 1];
                res3la[2 - 1] = res3la[3 - 1];
                res3la[3 - 1] = result;//modifiés 3x
            }
            else
            {
                res3la[nres - 1] = result;
                abserr = oflow;
            }

            L100:/* compute error estimate */
            abserr = RVaria.fmax2(abserr, epmach * 5.0 * Math.Abs(result));
            return;
        }
        //(out numrl2, rlist2, out reseps, out abseps, res3la, out nres);
 
        /* R file: integrate.c
         */
        private static void rdqpsrt(int limit, int last, ref int maxerr, out double ermax, double[] elist, int[] iord, ref int nrmax)
        {
            /* Local variables */
            int i, j, k, ido, jbnd, isucc, jupbn;
            double errmin, errmax;
            #region prologue
            /* ***begin prologue  dqpsrt
             ***refer to  dqage,dqagie,dqagpe,dqawse
             ***routines called  (none)
             ***revision date  810101   (yymmdd)
             ***keywords  sequential sorting
             ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
                       de doncker,elise,appl. math. & progr. div. - k.u.leuven
             ***purpose  this routine maintains the descending ordering in the
                        list of the local error estimated resulting from the
                        interval subdivision process. at each call two error
                        estimates are inserted using the sequential search
                        method, top-down for the largest error estimate and
                        bottom-up for the smallest error estimate.
             ***description

                       ordering routine
                       standard fortran subroutine
                       double precision version

                       parameters (meaning at output)
                          limit  - int
                                   maximum number of error estimates the list
                                   can contain

                          last   - int
                                   number of error estimates currently in the list

                          maxerr - int
                                   maxerr points to the nrmax-th largest error
                                   estimate currently in the list

                          ermax  - double precision
                                   nrmax-th largest error estimate
                                   ermax = elist(maxerr)

                          elist  - double precision
                                   vector of dimension last containing
                                   the error estimates

                          iord   - int
                                   vector of dimension last, the first k elements
                                   of which contain pointers to the error
                                   estimates, such that
                                   elist(iord(1)),...,  elist(iord(k))
                                   form a decreasing sequence, with
                                   k = last if last <= (limit/2+2), and
                                   k = limit+1-last otherwise

                          nrmax  - int
                                   maxerr = iord(nrmax)

            ***end prologue  dqpsrt
            */
            #endregion

            /* Parameter adjustments */
            // mis en commentaires par zygotine
            //--iord;
            //--elist;

            /* Function Body */

            /*           check whether the list contains more than
                     two error estimates. */
            if (last <= 2)
            {
                iord[1] = 1;
                iord[2] = 2;
                goto Last;
            }
            /*           this part of the routine is only executed if, due to a
                     difficult integrand, subdivision increased the error
                     estimate. in the normal case the insert procedure should
                     start after the nrmax-th largest error estimate. */

            errmax = elist[maxerr];
            if (nrmax > 1)
            {
                ido = nrmax - 1;
                for (i = 1; i <= ido; ++i)
                {
                    isucc = iord[nrmax - 1];
                    if (errmax <= elist[isucc])
                        break; /* out of for-loop */
                    iord[nrmax] = isucc;
                    --nrmax;
                    /* L20: */
                }
            }

            /*L30:       compute the number of elements in the list to be maintained
                     in descending order. this number depends on the number of
                     subdivisions still allowed. */
            if (last > limit / 2 + 2)
                jupbn = limit + 3 - last;
            else
                jupbn = last;

            errmin = elist[last];

            /*           insert errmax by traversing the list top-down,
                     starting comparison from the element elist(iord(nrmax+1)). */

            jbnd = jupbn - 1;
            for (i = nrmax + 1; i <= jbnd; ++i)
            {
                isucc = iord[i];
                if (errmax >= elist[isucc])
                {/* ***jump out of do-loop */
                    /* L60: insert errmin by traversing the list bottom-up. */
                    iord[i - 1] = maxerr;
                    for (j = i, k = jbnd; j <= jbnd; j++, k--)
                    {
                        isucc = iord[k];
                        if (errmin < elist[isucc])
                        {
                            /* goto L80; ***jump out of do-loop */
                            iord[k + 1] = last;
                            goto Last;
                        }
                        iord[k + 1] = isucc;
                    }
                    iord[i] = last;
                    goto Last;
                }
                iord[i - 1] = isucc;
            }

            iord[jbnd] = maxerr;
            iord[jupbn] = last;

        Last:/* set maxerr and ermax. */

            maxerr = iord[nrmax];
            ermax = elist[maxerr];
            return;
        } /* rdqpsrt_ */
    }
}
