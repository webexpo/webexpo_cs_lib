namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;


    public static class WebExpoFunctions3
    {
        internal static readonly double[] LPQAC2 = Tools.Combine(-0.6931471805599452862268, 0.7978845608028654057264, -0.3183098861837906912164);
            //Tools.Combine(-0.693147180559945, 0.797884560802865, -0.318309886183791); // logPhi.quadratic.approx.coeff
        internal static readonly double[] LPQAC3 = Tools.Combine(-0.6931471805599452862268, 0.7978845608028654057264, -0.3183098861837906912164, 0.03633560235749837274977);
            //Tools.Combine(-0.693147180559945, 0.797884560802865, -0.318309886183791, 0.0363356023574984); // logPhi.quadratic.approx.coeff

        internal static InitialValues DefaultInits(DataSummary data, bool logNormalDistrn, double[] muRange, double[] sigmaRange, bool includeCensoredData)
        {

            //# modif_0.13
            //# added argument include.censored.data to this fct
            double mu;
            double sigma;
            double sigmaLower = sigmaRange.Min();
            double sigmaUpper = sigmaRange.Max();
            List<double> normalizedY = new List<double>();
            if (data.YLength > 0)
            {
                normalizedY.AddRange(logNormalDistrn ? data.LogY : data.Y);
            }

            // new_0.12
            if (includeCensoredData)
            {
                if (data.GTLength > 0)
                {
                    normalizedY.AddRange(logNormalDistrn ? data.LogGT : data.GT);
                }

                if (data.LTLength > 0)
                {
                    normalizedY.AddRange(logNormalDistrn ? data.LogLT : data.LT);
                }

                if (data.IntervalLength > 0)
                {
                    double[] midPoints = new double[data.IntervalLength];
                    double[] lt = logNormalDistrn ? data.LogIntervalLT : data.IntervalLT;
                    double[] gt = logNormalDistrn ? data.LogIntervalGT : data.IntervalGT;
                    for (int i = 0; i < data.IntervalLength; i++)
                    {
                        midPoints[i] = gt[i] + (lt[i] - gt[i]) / 2.0;
                    }
                    normalizedY.AddRange(midPoints);
                }
            }

            if (data.YLength > 0)
            {
                mu = normalizedY.Mean();
                double lower = muRange.Min();

                if (mu < lower)
                {
                    mu = lower;
                }
                else
                {
                    double muUpper = muRange.Max();
                    if (mu > muUpper)
                    {
                        mu = muUpper;
                    }
                }

                if (normalizedY.Count > 1)
                {
                    sigma = Math.Sqrt(normalizedY.Variance());
                }
                else if (!includeCensoredData)
                {
                    InitialValues tmp = DefaultInits(data, logNormalDistrn, muRange, sigmaRange, includeCensoredData: true);
                    sigma = tmp.SigmaWithin;
                }
                else
                {
                    //# new_0.13
                    sigma = 0; //# will be corrected below
                }

                if (sigma < sigmaLower)
                {
                    sigma = sigmaLower;
                }

                // sigma >= 0 && sigma >= sigmaLower
                if (sigma == 0.0)
                {
                    // sigmaLower <= 0 !!
                    sigma = sigmaUpper / 10.0;
                }
            }
            else
            {
                sigma = sigmaLower > 0 ? sigmaLower : sigmaUpper / 10.0;
                mu = Tools.Mean(muRange);
            }

            return new InitialValues(mu, sigma);
        } // Default.inits


        internal static double[] Diff(this double[] x)
        {
            if (x == null)
            {
                return null;
            }

            if (x.Length <= 1)
            {
                return new double[0];
            }

            double[] rep = new double[x.Length - 1];
            for (int i = 0; i < x.Length - 1; i++)
            {
                rep[i] = x[i + 1] - x[i];
            }

            return rep;
        }

       internal static LPhiObject LPhi(double z)
        {
            //# z: can be scalar or vector
            return LPhi(Tools.Combine(z));
        }

        internal static LPhiObject LPhi(double[] z)
        {
            //# z: can be scalar or vector
            LPhiObject rep = new LPhiObject(z.Length);
            /*
                lphi <- function(z)
                {
                  # z: can be scalar or vector
  
                  log.phi <- dnorm(z, log=T)      # log(phi(z))
                  log.Phi <- pnorm(z, log.p=T)    # log(Phi(z))
  
                  r  <- exp(log.phi-log.Phi)      # phi(z)/Phi(z)
                  r2 <- exp(2*(log.phi-log.Phi))  # phi^2(z) / Phi^2 (z)
  
                  list(r=r, r2=r2) # r and r2 are of same length as z
                } # end of lphi
            */

            double[] tmp = NormalDistribution.DNorm(z, give_log: true).Substract(NormalDistribution.PNorm(z, log_p: true));
            rep.R = tmp.Exp();  // phi(z)/Phi(z)
            rep.R2 = tmp.Multiply(2).Exp();  // phi^2(z) / Phi^2 (z)
            return rep;
        }


        internal class SAEList
        {
            internal bool[] Ok = Tools.Rep(false, 2);
            internal double[] X = Tools.Rep(Tools.NA, 2);
            internal double[] A = Tools.Rep(Tools.NA, 2);
            internal double[] F = Tools.Rep(Tools.NA, 2);
        }
        //smoothed.area.estimate
        internal static double SmoothedAreaEstimate(Func<double, SGNFnAParam, double> area, double x, SGNFnAParam A, double mathLowerLimit, double[] refRange, bool inestimableLowerLimit, bool remoteRight = false, int hpMult = 1, int maxCount = 10000)
        {
            double expectedSlope = Tools.NA; // non défini
            double fX = A.F(x, A);
            double[] a = Tools.Rep(Tools.NA, 3);
            // # values around x
            double d = 1e-4;
            double tmp = (refRange[1] - refRange[0]) / 1000.0;
            if (tmp < d)
            {
                d = tmp;
            }

            double[] xStar = Tools.Combine(x - d, x, x + d);

            if (xStar[0] < mathLowerLimit)
            {
                //xStar[0] = inestimableLowerLimit ? ((x + mathLowerLimit) / 2.0) : mathLowerLimit;
                if (inestimableLowerLimit)
                {
                    xStar[0] = (x + mathLowerLimit) / 2.0;
                }
                else
                {
                    xStar[0] = mathLowerLimit;
                }
            }

            for (int i = 0; i < 3; i++)
            {
                a[i] = area(xStar[i], A);
            }

            double aX = a[1];
            double[] observedSlopes = a.Diff().Divide(xStar.Diff()).Multiply(hpMult);
            expectedSlope = fX;
            double[] observedSlopeAngles = Tools.Combine(Math.Atan(observedSlopes[0]), Math.Atan(observedSlopes[1]));
            double expectedSlopeAngle = Math.Atan(expectedSlope);
            double angleDiff = Math.Abs(expectedSlopeAngle - observedSlopeAngles[0]);
            tmp = Math.Abs(expectedSlopeAngle - observedSlopeAngles[1]);
            if (tmp > angleDiff)
            {
                angleDiff = tmp;
            }

            bool estimateOk = (angleDiff < 0.035) && observedSlopes.All(z => z > 0);
            //# 0.035 radians is approximately 2 degrees
            bool cntn = !estimateOk;
            SAEList found = new SAEList();
            if (cntn)
            {
                /*
                 * # We look for a point on the left side of x and one on its right side
                 * # where the integral [area] seems to be correctly estimated
                 * # (with positive slopes of value close to expected value [density])
                 * # If we are estimating the area on the remote right end of the distrn (remote.right=T)
                 * # then we will be happy if an estimate to the left of x is obtained 
                 * 
                 */
                int side = -1;
                //# Store starting points
                double[] tmpX = xStar.Copy();//?
                double[] tmpA = a.Copy();//?
                while (cntn)
                {
                    side = side + 1;    // # side = 0 -> Left  of x
                                        // # side = 1 -> Right of x
                    int direction = side == 0 ? -1 : 1;
                    if (side == 1)
                    {
                        xStar = tmpX.Copy(); //?
                        a = tmpA.Copy();
                    }

                    int j = side == 0 ? 0 : 2;
                    double z = xStar[j];
                    int count = 0;
                    bool onBorder = false;
                    double aZ = Tools.NA;
                    while (cntn)
                    {
                        z = z + direction * d;
                        //# make sure that we don't step to the left of mathematical lower limit
                        if (z <= mathLowerLimit)
                        {
                            if (inestimableLowerLimit)
                            {
                                z = (mathLowerLimit + z + d) / 2.0;
                            }
                            else
                            {
                                z = mathLowerLimit;
                                onBorder = true;
                            }
                        }

                        if ((side == 0) && z == mathLowerLimit)
                        {
                            //# This is still possible even with the above 'protection' against
                            //# this scenario, due to numerical imprecision
                            aZ = 0;
                            cntn = false;
                            estimateOk = true;
                        }
                        else
                        {
                            aZ = area(z, A);
                            if (side == 0)
                            {
                                xStar = Tools.Combine(z, xStar[0], xStar[1]); //[-3]
                                a = Tools.Combine(aZ, a[0], a[1]);
                            }
                            else
                            {
                                xStar = Tools.Combine(xStar[1], xStar[1], z); //[-1]
                                a = Tools.Combine(a[1], a[2], aZ);
                            }

                            //a.Diff().Divide(xStar.Diff()).Multiply(hpMult);
                            observedSlopes = a.Diff().Divide(xStar.Diff()).Multiply(hpMult);
                            expectedSlope = A.F(z, A);
                            observedSlopeAngles = Tools.Combine(Math.Atan(observedSlopes[0]), Math.Atan(observedSlopes[1]));
                            expectedSlopeAngle = Math.Atan(expectedSlope);
                            angleDiff = Math.Abs(expectedSlopeAngle - observedSlopeAngles[0]);
                            tmp = Math.Abs(expectedSlopeAngle - observedSlopeAngles[1]);
                            if (tmp > angleDiff)
                            {
                                angleDiff = tmp;
                            }
                            estimateOk = ((angleDiff < 0.035) && observedSlopes.All(w => w > 0)) || (angleDiff < 1E-6);
                            //# 0.035 radians is approximately 2 degrees
                            cntn = !estimateOk && (!onBorder);
                        }
                        if (cntn)
                        {
                            if (remoteRight && (side == 0))
                            {
                                cntn = (expectedSlope / fX) < 10000;
                            }
                            else
                            {
                                count = count + 1;
                                cntn = count <= maxCount;
                            }
                        }
                    }

                    if ((side == 0) && remoteRight)
                    {
                        cntn = !estimateOk;
                        if (estimateOk)
                        {
                            aX = aZ;
                        }
                    }
                    else
                    {
                        cntn = side == 0;
                    }

                    //# Store results
                    found.Ok[side] = estimateOk;
                    found.X[side] = z;
                    found.A[side] = aZ;
                    found.F[side] = expectedSlope;
                } //# end of while-continue

                if (!remoteRight || !found.Ok[0])
                {
                    if (found.Ok[0] && found.Ok[1])
                    {
                        double fDiff = found.F.Diff()[0];
                        bool extrapolate = fDiff == 0;
                        if (!extrapolate)
                        {
                            //- (diff(found$a) - diff(found$f*found$x)) / f.diff
                            double t1 = (found.A[1] - found.A[0]); //diff(found$a)
                            double t2 = (found.F[1] * found.X[1]) - (found.F[0] * found.X[0]); //diff(found$f*found$x)
                            double xIntersect = -(t1 - t2) / fDiff;
                            extrapolate = (xIntersect < found.X[0]) || (xIntersect > found.X[1]);
                            if (!extrapolate)
                            {
                                int j = x <= xIntersect ? 0 : 1;
                                aX = found.A[j] + found.F[j] * (x - found.X[j]);
                            }
                        }

                        if (extrapolate)
                        {
                            aX = found.A[0] + Diff(found.A)[0] / Diff(found.X)[0] * (x - found.X[0]);
                        }
                    }
                    else if (found.Ok[0] || found.Ok[1])
                    {
                        int j = found.Ok[0] ? 0 : 1;
                        aX = found.A[j] + found.F[j] * (x - found.X[j]);
                    }
                    else
                    {
                        //TODO Exception
                        throw new WEException("Could not get a smoothed cumulative density estimate.");
                    }
                }
            }
            return aX;
        } // end of smoothed.area.estimate                    

        //random.pow < -function(a, range)
        public static double RandomPow(int a, double[] range)
        {
            //# sample a value from f(x) = 1/x^a on the range specified

            double z;
            double U = UniformDistribution.RUnif();
            if (a == 1)
            {
                if (double.IsInfinity(range[1]))
                {
                    range[1] = 1E8;
                }

                z = Diff(Tools.Combine(U - 1, U).Multiply(range.Log())).Exp()[0];
            }
            else
            {
                double[] fCum = Tools.Combine(Math.Pow(range[0], (1 - a)), Math.Pow(range[1], (1 - a)));
                double tmp = fCum.Multiply(Tools.Combine(1.0 - U, U)).Sum();
                z = Math.Pow(tmp, (1.0 / (1.0 - a)));
            }

            return z;
        } //# end of random.pow

        public static double SqrtInvertedGammaGen(int n, double beta, double[] xRange, double betaMin = 1E-8)
        {
            return SqrtInvertedGammaGen(n, beta, xRange, null, betaMin);
        }

        //sqrt.invertedGamma.gen <- function(n, beta, xrange, o=list(), beta.min=1e-8)
        internal static double SqrtInvertedGammaGen(int n, double beta, double[] xRange, GenObject o, double betaMin = 1E-8)
        {
            //# Sample a value from the distrn 
            //#
            //# f(x) =     1
            //#         ------- . exp(-beta/x^2)
            //# x^n

            double x = Tools.ND;
            if (n <= 1)
            {
                SGNFnAParam A = o.A.Clone();
                A.B = beta;
                Icdf icdf = new Icdf(o, A, xRange);
                x = icdf.Bidon(start: Math.Sqrt(2.0 * beta), inestLowerLim: xRange[0] == 0);

            }
            else if (beta < betaMin)
            {
                if (xRange[0] == 0)
                {
                    xRange[0] = 0.0001;
                }
                x = RandomPow(n, xRange);
            }
            else
            {
                double alpha = (n - 1) / 2.0;
                double[] invX2Range = xRange.Sqr().Reverse().ToArray().Reciprocal();
                double invX2 = RGammaTruncated(alpha, beta, invX2Range);
                x = 1 / Math.Sqrt(invX2);
            }

            return x;
        } //# end of sqrt.invertedGamma.gen

        //rgamma.truncated <- function(alpha, beta, range)
        public static double RGammaTruncated(double alpha, double beta, double[] range)
        {
            //# We first look on which side of the distribution (comparing to its mode)
            //# the range lies, in order to decide on which side we will compute the log.prob
            //# (choosing the correct side gives more accuracy in the sampling in the eventuality of a remote range)
            bool lowerTail = false;
            if (alpha < 1)
            {
                lowerTail = false;
            }
            else
            {
                double mode = (alpha - 1) / beta;
                double minRange = range.Min();
                double maxRange = range.Max();

                if (minRange > mode)
                {
                    lowerTail = false;
                }
                else if (maxRange < mode)
                {
                    lowerTail = true;
                }
                else
                {
                    double logPLeft = GammaDistribution.PGammaFromRateParameter(x: minRange, alph: alpha, rate: beta, logP: true, lowerTail: true);
                    double logPRight = GammaDistribution.PGammaFromRateParameter(x: maxRange, alph: alpha, rate: beta, logP: true, lowerTail: false);
                    lowerTail = logPLeft < logPRight;
                }
            }

            double[] logPLim = GammaDistribution.PGammaFromRateParameter(x: range, alph: alpha, rate: beta, lowerTail: lowerTail, logP: true);
            double logP = RNorm4CensoredMeasures.RUnifLogP1(logPLim, lowerTail: lowerTail);
            return GammaDistribution.QGammaFromRate(logP, shape: alpha, rate: beta, logP: true, lowerTail: lowerTail);

        }


    }
}

