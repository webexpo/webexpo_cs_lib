using System;
using Zygotine.Statistics.Util;

namespace Zygotine.Numerics.Integration
{

    /*
     * R file: integrate.c
     * 
     * La classe NumericIntegration implante en C# la fonction Rdqagi définie dans le fichier integrate.c.
     * 
     */

    public partial class NumericIntegration
    {

        public NumericIntegration(Func<double[], double[]> f, double lowerBound, double upperBound)
        {
            this.a = lowerBound;
            this.b = upperBound;
            this.fn = f;
        }

        public static NumericIntegration.Results Integrate(
                Func<double,
                double> f,
                double lowerBound,
                double upperBound,
                double epsabs = Accuracy,
                double epsrel = Accuracy,
                int limit = 100)
        {
            Func<double[],
                double[]> f2 = GetF4Integration(f);
            return Integrate(f2, lowerBound, upperBound, epsabs, epsrel, limit);
        }

        private static Func<double[], double[]> GetF4Integration(Func<double, double> f)
        {
            Func<double[], double[]> _dummy = delegate (double[] x)
            {
                double[] y = new double[x.Length];
                for (int i = 0; i < x.Length; i++)
                {
                    y[i] = f(x[i]);
                }
                return y;
            };
            return _dummy;
        }


        public static NumericIntegration.Results Integrate(
                Func<double[],
                double[]> f,
                double lowerBound,
                double upperBound,
                double epsabs = Accuracy,
                double epsrel = Accuracy,
                int limit = 100)
        {
            NumericIntegration ni = new NumericIntegration(f, lowerBound, upperBound);
            return ni.Compute(epsabs, epsrel, limit);
        }

        // fields
        private const double Accuracy = 0.0001220703125; // 2.220446049250313e-16 ^ .25 == .Machine$double.eps ^ .25
        private Func<double[], double[]> fn;
        private double a, b; // bounds

        private Results Compute(
            double epsabs = Accuracy,
            double epsrel = Accuracy,
            int limit = 100
            )
        {
            /*
                f : double, function subprogram defining the integrand
                epsabs - double, absolute accuracy requested, > 0.0
                epsrel - double precision, relative accuracy requested


                on return
                result - double, approximation to the integral
                abserr - double, estimate of the modulus of the absolute error, which should equal or exceed abs(i-result)
                neval  - int, number of integrand evaluations
                ier     int, 
                        ier = 0 normal and reliable termination of the routine...
                        ier = 1 maximum number of subdivisions allowed has been achieved. 
                        ier = 2 the occurrence of roundoff error is detected ...
                        ier = 3 extremely bad integrand behaviour occurs at some points of the integration interval.
                        ier = 4 the algorithm does not converge...
                        ier = 5 the integral is probably divergent, or slowly convergent. ...
                        ier = 6 the input is invalid
                limit - int, limit determines the maximum number of subintervals
                lenw  - int, >= limit * 4
                last  - int, on return, last equals the number of subintervals ...

                work arrays
                    iwork - int
                            vector of dimension at least limit, the first k
                            elements of which contain pointers
                            to the error estimates over the subintervals
                            such that work(limit*3+iwork(1)),... ,
                            work(limit*3+iwork(k)) form a decreasing
                            sequence, with k = last if last <= (limit/2+2),
                            and k = limit+1-last otherwise

                    work  - double precision
                            vector of dimension at least lenw
                            on return
                            work(1), ..., work(last) contain the left
                            end-points of the subintervals in the
                            partition of (a,b),
                            work(limit+1), ..., work(limit+last) contain
                            the right end-points,
                            work(limit*2+1), ..., work(limit*2+last) contain
                            the integral approximations over the subintervals,
                            work(limit*3+1), ..., work(limit*3+last)
                            contain the error estimates.
            */

            double a = this.a;
            double b = this.b;

            Results results = new Results();

            results.Value = 0.0;
            results.NEval = 0;
            results.IEr = 6;
            results.Last = 0;
            results.AbsErr = 0.0;

            if (limit < 1)
            {
                return results;
            }

            Workspace workspace = new Workspace(limit);

            //int
            //    l1 = limit,
            //    l2 = limit + l1,
            //    l3 = limit + l2;
            Parameters parameters = new Parameters(fn, a, b, epsabs, epsrel, limit);

            rdqagse(parameters, workspace, results);

            return results;
        }


        internal class Workspace
        {
            public int limit { get; internal set; }
            public double[] alist { get; internal set; }
            public double[] blist { get; internal set; }
            public double[] rlist { get; internal set; }
            public double[] elist { get; internal set; }
            public int[] iord { get; internal set; }

            internal Workspace(int limit)
            {
                this.limit = limit;
                this.alist = new double[limit + 1];
                this.blist = new double[limit + 1];
                this.rlist = new double[limit + 1];
                this.elist = new double[limit + 1];
                this.iord = new int[limit + 1];
            }

        }

        public class Results
        {

            internal Results() { }

            public double Value { get; internal set; }
            public int NEval { get; internal set; }
            public int IEr { get; internal set; }
            public int Last { get; internal set; }
            public double AbsErr { get; internal set; }
            public override string ToString()
            {
                return string.Format("result={0},  absoluteError={1}, errorNo={2}, nSubdivisions={3}", Value, AbsErr, IEr, Last);
            }
        }

        internal class Parameters
        {
            internal Func<double[], double[]> f;
            internal double a;
            internal double b;
            internal double epsabs;
            internal double epsrel;
            internal int limit;

            internal Parameters(Func<double[], double[]> f, double a, double b, double epsabs, double epsrel, int limit)
            {
                this.f = f;
                this.a = a;
                this.b = b;
                this.epsabs = epsabs;
                this.epsrel = epsrel;
                this.limit = limit;
            }
        }

        private const double _stopd = 8888d;
        private const int _stopi = 8888;

        private static void setArray(double[] arr)
        {
            for (int i = 0; i < arr.Length; i++)
                arr[i] = _stopd;
        }

        private static void setArray(int[] arr)
        {
            for (int i = 0; i < arr.Length; i++)
                arr[i] = _stopi;
        }


        public static Func<double[], double[]> GetFunc2(int n, double b, double m, double s, double k)
        {
            Func<double[], double[]> _dummy = delegate (double[] x)
           {
               double[] y = new double[x.Length];

               Func<double, double> log_f = delegate (double x2)
               {
                   return -(n + 1) * Math.Log(x2) - b / (x2 * x2) - Math.Pow((Math.Log(x2) - m), 2) / (2 * s * s);
               };

               for (int i = 0; i < x.Length; i++)
               {

                   y[i] = Math.Exp(log_f(x[i]) + k);
               }
               return y;
           };
            return _dummy;

        }

        public static NumericIntegration.Results TestMe(Func<double[], double[]> fn, double a, double b)
        {
            NumericIntegration integration = new NumericIntegration(fn, a, b);
            Results res = integration.Compute(Math.Pow(RVaria.DBL_EPSILON, 0.25), Math.Pow(RVaria.DBL_EPSILON, 0.25));
            return res;
        }


    }

}

