namespace Zygotine.WebExpo
{
    using System;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    public static class RNorm4CensoredMeasures
    {
        //nouveau 24 avril
        internal static double[] RNormCensored(double[] mu, double[] sd, double[] lower = null, double[] upper = null, bool negativeDisallowed = false)
        {
            //# Beginning of calculations
            //# new_0.4
            //# (note that negative.values.disallowed was added to fct arguments)
            if (mu.Length == 1)
            {
                return RNormCensored(mu[0], sd[0], lower, upper, negativeDisallowed);
            }

            double[] upperBound = upper == null ? new double[0] : upper;
            double[] lowerBound = lower == null ? new double[0] : lower;
            int len = mu.Length;
            if (negativeDisallowed)
            {
                lowerBound = new double[Math.Max(1, upperBound.Length)]; // lowerbound ne contient que des zéros
            }

            bool leftCensored = lowerBound.Length > 0;
            double[] z = new double[len];
            if (leftCensored)
            {
                bool rightCensored = upperBound.Length > 0;
                if (rightCensored)
                {
                    for (int i = 0; i < len; i++)
                    {
                        //# it is interval-censored
                        z[i] = RNormIntervalCensored(mu[i], sd[i], lowerBound[i], upperBound[i]);
                    }
                }
                else
                {
                    for (int i = 0; i < len; i++)
                    {
                        //# is is left-censored
                        z[i] = RNormLeftCensored(mu[i], sd[i], lowerBound[i]);
                    }
                }
            }
            else
            {
                for (int i = 0; i < len; i++)
                {
                    //# then it is right-censored
                    z[i] = RNormRightCensored(mu[i], sd[i], upperBound[i]);
                }
            }
            return z;
        } //# end of rnorm.censored

        //nouveau 24 avril
        internal static double RNormCensored(double mu, double sd, double lower = double.NaN, double upper = double.NaN, bool negativeDisallowed = false)
        {
            return RNormCensored(
                mu: Tools.Combine(mu),
                sd: Tools.Combine(sd),
                lower: double.IsNaN(lower) ? null : Tools.Combine(lower),
                upper: double.IsNaN(upper) ? null : Tools.Combine(upper),
                negativeDisallowed: negativeDisallowed)[0];
        }

        internal static double[] RNormCensored(double mu, double sd, double[] lower = null, double[] upper = null, bool negativeDisallowed = false)
        {
            //# Beginning of calculations
            //# new_0.4
            //# (note that negative.values.disallowed was added to fct arguments)

            double[] upperBound = upper == null ? new double[0] : upper;
            double[] lowerBound = lower == null ? new double[0] : lower;
            if (negativeDisallowed)
            {
                int len = Math.Max(1, upperBound.Length);
                lowerBound = new double[len];
            }

            bool leftCensored = lowerBound.Length > 0;
            double[] z = null;
            if (leftCensored)
            {
                bool rightCensored = upperBound.Length > 0;
                if (rightCensored)
                {
                    //# it is interval-censored
                    z = RNormIntervalCensored(mu, sd, lowerBound, upperBound);
                }
                else
                {
                    //# is is left-censored
                    z = RNormLeftCensored(mu, sd, lowerBound);
                }
            }
            else
            {
                //# then it is right-censored
                z = RNormRightCensored(mu, sd, upper);
            }

            return z;
        } //# end of rnorm.censored

        internal static double RUnifLogP1(double[] logPLim, bool lowerTail = true, double u = double.NaN)
        {
            // On suppose logPLim est de longueur 2, et que logPLim[1] > logPLim[0]
            if (double.IsNaN(u))
            {
                u = UniformDistribution.RUnif();
            }
            //# To sample 'size' values uniformly in c(exp(logp.lim[1]), exp(logp.lim[2]))

            double w = logPLim[1] - logPLim[0];
            if (!lowerTail)
            {
                w = Math.Abs(w);
            }

            double logP = Math.Log(u) + ExponentialDistribution.PExp(w, logP: true);
            double x = ExponentialDistribution.QExp(logP, logP: true);
            x = logPLim.Max() - x;
            return x;
        }

        //nouveau 24 avril
        private static double RNormIntervalCensored(double mu, double sd, double lower, double upper)
        {
            //# mu, sd: same length (1, or same as length(lower))
            //# lower and upper: same length, but length can be > 1

            double logPLower = NormalDistribution.PNorm(lower, mu: mu, sigma: sd, log_p: true);
            double logPUpper = NormalDistribution.PNorm(upper, mu: mu, sigma: sd, log_p: true);
            double logP = RUnifLogP1(Tools.Combine(logPLower, logPUpper));
            double y = NormalDistribution.QNorm(logP, mu: mu, sigma: sd, log_p: true);
            return y;
        }//# end of rnorm.interval.censored

        private static double[] RNormIntervalCensored(double mu, double sd, double[] lower, double[] upper)
        {
            //# mu, sd: same length (1, or same as length(lower))
            //# lower and upper: same length, but length can be > 1

            double[] logPLower = NormalDistribution.PNorm(lower, mu: mu, sigma: sd, log_p: true);
            double[] logPUpper = NormalDistribution.PNorm(upper, mu: mu, sigma: sd, log_p: true);
            double[] logP = new double[lower.Length];
            for (int i = 0; i < lower.Length; i++)
            {
                logP[i] = RUnifLogP1(Tools.Combine(logPLower[i], logPUpper[i]));
            }

            double[] y = NormalDistribution.QNorm(logP, mu: mu, sigma: sd, log_p: true);
            return y;
        }//# end of rnorm.interval.censored
        
        //nouveau le 24 avril
        private static double RNormRightCensored(double mu, double sd, double upper)
        {
            //# This function was incorrectly named rnorm.left.censored in previous versions
            //# mu, sd: same length (1, or same as length(upper))

            double logP =  Math.Log(UniformDistribution.RUnif()) + NormalDistribution.PNorm(upper, mu: mu, sigma: sd, log_p: true);
            double y = NormalDistribution.QNorm(logP, mu: mu, sigma: sd, log_p: true);
            return y;
        } //# end of rnorm.right.censored

        private static double[] RNormRightCensored(double mu, double sd, double[] upper)
        {
            //# This function was incorrectly named rnorm.left.censored in previous versions
            //# mu, sd: same length (1, or same as length(upper))

            double[] logP = UniformDistribution.RUnif(upper.Length).Log().Add(NormalDistribution.PNorm(upper, mu: mu, sigma: sd, log_p: true));
            double[] y = NormalDistribution.QNorm(logP, mu: mu, sigma: sd, log_p: true);
            return y;
        } //# end of rnorm.right.censored

        //nouveau 24 avril
        private static double RNormLeftCensored(double mu, double sd, double lower)
        {
            double y = -RNormRightCensored(-mu, sd, -lower);
            return y;
        } //# end of rnorm.left.censored

        private static double[] RNormLeftCensored(double mu, double sd, double[] lower)
        {
            double[] y = RNormRightCensored(-mu, sd, lower.Multiply(-1)).Multiply(-1);
            return y;
        } //# end of rnorm.left.censored
    }
}
