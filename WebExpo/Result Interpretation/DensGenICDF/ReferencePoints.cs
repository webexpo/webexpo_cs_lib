namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Zygotine.Numerics;
    using Zygotine.Util;

    internal class ReferencePoints
    {
        internal SGNFnAParam A { get; private set; }
        internal GenObject O { get; private set; }
        internal double[] Start { get; private set; } = new double[0];
        internal double[] Range { get; private set; }
        internal bool InestimableLowerLimit { get; private set; }
        internal double FRatioRemote { get; private set; }
        internal double Epsilon { get; private set; }
        internal int MaxNIter { get; private set; }

        internal ReferencePoints(
                GenObject o,
                SGNFnAParam a,
                double[] start,
                double[] range,
                bool inestimableLowerLimit = false,
                double fRatioRemote = 1e-8,
                double epsilon = 1e-6,
                int maxNIter = 200)
        {
            this.A = a;
            this.O = o;
            this.Start = start;
            this.Range = range;
            this.InestimableLowerLimit = inestimableLowerLimit;
            this.FRatioRemote = fRatioRemote;
            this.Epsilon = epsilon;
            this.MaxNIter = 200;
        }


        private class MaxFastTrackObject
        {
            internal double X { get; set; } = Tools.NA;
            internal double H { get; set; } = Tools.NA;
            internal bool Converged { get; set; } = false;

            private MaxFastTrackObject()
            {
            }

            internal MaxFastTrackObject(double x, double h, bool converged)
            {
                X = x;
                H = h;
                Converged = converged;
            }
        }

        private double HigherLocalMax(GenObject o, SGNFnAParam A, double dipX, double[] range, double epsilon)
        {
            //# Find a starting point on each side of dip.x
            double[] start = o.Start(A);
            double tmp = start.Substract(dipX).Abs().Max();
            double startLeft = start.Min();
            if (startLeft > dipX)
            {
                startLeft = dipX - tmp;
            }

            double startRight = start.Max();
            if (startRight < dipX)
            {
                startRight = dipX + tmp;
            }

            //# Look for starting points with negative values for h''
            //# i) on left side
            bool cntn = o.LogFSecond(startLeft, A) > 0;
            while (cntn)
            {
                startLeft = startLeft - tmp;
                if (startLeft < range[0])
                {
                    startLeft = (range[0] + startLeft + tmp) / 2;
                }
                cntn = o.LogFSecond(startLeft, A) > 0;
            }

            //# ii) on right side
            cntn = o.LogFSecond(startRight, A) > 0;
            while (cntn)
            {
                startRight = startRight + tmp;
                cntn = o.LogFSecond(startRight, A) > 0;
            }

            //# Run pure Newton-Raphson to find local max on each side
            //# i) left side
            double x = startLeft;
            cntn = true;
            double hp = o.LogFPrime(x, A);
            while (cntn)
            {
                double change = hp / o.LogFSecond(x, A);
                x = x - change;
                hp = o.LogFPrime(x, A);
                cntn = Math.Abs(hp) > epsilon;
            }

            double localModeLeft = x;
            //# ii) right-side
            x = startRight;
            cntn = true;
            hp = o.LogFPrime(x, A);
            while (cntn)
            {
                double change = hp / o.LogFSecond(x, A);
                x = x - change;
                hp = o.LogFPrime(x, A);
                cntn = Math.Abs(hp) > epsilon;
            }
            double localModeRight = x;
            double hLeft = o.LogF(localModeLeft, A);
            double hRight = o.LogF(localModeRight, A);
            x = hLeft > hRight ? localModeLeft : localModeRight;
            return x;
        } //# end of higher.local.max

        //# new_0.11
        private double[] SplitStartingValuesLeftRight(GenObject o, double target, SGNFnAParam A, double[] range, double[] startValues)
        {
            double[] start = new double[2] { Tools.NA, Tools.NA };
            double[] lsStart = startValues.Where(number => (number < A.Mode) && (number > range[0])).ToArray();
            double[] rsStart = startValues.Where(number => (number > A.Mode) && (number < range[1])).ToArray();

            if (lsStart.Length == 1)
            {
                start[0] = lsStart[0];
            }
            else if (lsStart.Length > 1)
            {
                double[] f = new double[lsStart.Length];
                for (int i = 0; i < lsStart.Length; i++)
                {
                    f[i] = o.LogF(lsStart[i], A);
                }

                double[] d = f.Substract(target).Abs();
                int w = d.WhichMin();
                start[0] = lsStart[w];
            }

            if (rsStart.Length == 1)
            {
                start[1] = rsStart[0];
            }
            else if (rsStart.Length > 1)
            {
                double[] f = new double[rsStart.Length];
                for (int i = 0; i < rsStart.Length; i++)
                {
                    f[i] = o.LogF(rsStart[i], A);
                }

                double[] d = f.Substract(target).Abs();
                int w = d.WhichMin();
                start[1] = rsStart[w];
            }

            return start;
        } //# end of split.starting.values.leftRight  

        //# new_0.11
        private MaxFastTrackObject MaxFastTrack(GenObject o, SGNFnAParam A, double lowerLimit, double rightsideX, double rightsideHp, double epsilon = 1e-4, double epsilonX = 1e-20)
        {
            double x = rightsideX;
            double h = o.LogF(x, A);
            double hp = rightsideHp;
            //b = list(x = c(NA, x), h = c(NA, h), hp = c(NA, hp))
            double[] bX = new double[2] { x, Tools.NA };
            double[] bH = new double[2] { h, Tools.NA };
            double[] bHp = new double[2] { hp, Tools.NA };

            x = (x + lowerLimit) / 2;
            bool cntn = true;
            //# Find a point with positive slope (that is, to the left of the mode) 
            while (cntn)
            {
                hp = o.LogFPrime(x, A);
                if (hp > 0)
                {
                    cntn = double.IsInfinity(hp);

                    if (cntn)
                    {
                        double nextX = (3 * x - lowerLimit) / 2.0;
                        lowerLimit = x;
                        cntn = Math.Abs(nextX - x) > epsilonX;
                        if (cntn)
                        {
                            x = nextX;
                        }
                    }
                }
                else
                {
                    x = (x + lowerLimit) / 2.0;
                }
            }

            hp = o.LogFPrime(x, A);
            bool converged = Tools.IsFinite(hp) & (hp > 0);
            cntn = converged;

            while (cntn)
            {
                int side = hp > 0 ? 0 : 1;
                bX[side] = x;
                bH[side] = o.LogF(x, A);
                bHp[side] = hp;

                double bDiff = bHp.Diff()[0];

                if (bDiff == 0)
                {
                    x = bX.Mean();
                }
                else
                {
                    x = (bHp.Multiply(bX).Diff()[0] - bH.Diff()[0]) / bDiff;
                    if ((x <= bX[0]) || (x >= bX[1]))
                    {
                        x = bX.Mean();
                    }

                    hp = o.LogFPrime(x, A);
                    cntn = Math.Abs(hp) > epsilon;
                }
            }
            h = o.LogF(x, A);

            return new MaxFastTrackObject(x, h, converged);

        } //# end of max.fastTrack

        internal LModeExt GetPoints()
        {
            //# Arguments:
            //# ----------------------------------------------
            //# o: same nature as o in dens.gen.icdf arguments
            //# start: EITHER a vector of length 2 [two potential starting points]
            //#        or a scalar
            //# f.ratio.remote: scalar
            //# epsilon: scalar
            //# range: vector of length 2
            //# --------------------------------------------------------------------------
            //# If range is finite and log.f not inestimable at lower end, 
            //# then we first look for values at both ends before searching for a mode

            //mode = list(x = numeric(0), h = numeric(0), found = false) # modif_0.11
            LMode mode = new LMode();
            if ((this.Range[0].IsFinite() && this.Range[1].IsFinite() && !this.InestimableLowerLimit))
            {
                double[] logFPrime = new double[] { this.O.LogFPrime(this.Range[0], A), this.O.LogFPrime(this.Range[1], A) };
                //# If slope (f') is the same at both ends, then there is no mode between the 2 endpoints (assuming a unimodal distrn)
                //# and we then assume that the highest end value for log.f is hence the local mode
                double[] logFPrimeSign = new double[] { Math.Sign(logFPrime[0]), Math.Sign(logFPrime[1]) };
                if (logFPrime.Prod() > 0)
                {
                    double[] logF = new double[] { this.O.LogF(this.Range[0], A), this.O.LogF(this.Range[1], A) };
                    int w = (logF[0] >= logF[1]) ? 0 : 1;
                    mode = new LMode(w == 0 ? this.Range[0] : this.Range[1], logF[w], true);
                } 
            } // if ((this.Range[0]

            SecuredNRSearch nr = null;
            double startingPoint;
            if (!mode.Found)
            {
                //# Choose a starting point among suggested starting points, if more than one point was suggested
                if (this.Start.Length > 1)
                {
                    this.Start = Tools.Combine(this.Start, this.Start.Mean());
                    //# Add the middle point to the two potential starting points
                    double[] logF = new double[this.Start.Length];
                    for(int i=0; i<this.Start.Length; i++)
                    {
                        logF[i] = this.O.LogF(this.Start[i], A);
                    }
                    //double[] logF = new double[] { this.O.LogF(this.Start[0], A), this.O.LogF(this.Start[1], A), this.O.LogF(this.Start.Mean(), A) };
                    //TODO peut planter
                    startingPoint = this.Start[logF.WhichMax()];
                } //if (this.Start.Length > 1)
                else
                {
                    startingPoint = this.Start[0];
                }
                //# Find mode by Newton-Raphson
                SecuredNRA oMode = new SecuredNRA(this.O.LogFPrime, this.O.LogFSecond, this.Range);
                nr = new SecuredNRSearch(oMode, A, startingPoint, this.Epsilon, this.Range); 
                // 4 derniers paramètres prennent les valeurs par défaut: 
                // maxPoint = null, target = 0, inestimableLowerLimit = false,  expectedHpSign = -1
                if (nr.Converged)
                {
                    double hs = this.O.LogFSecond(nr.X, A);
                    if (hs > 0)
                    {
                        //# we have found a mimimum point (between two local modes):
                        //# we redo the search on both sides of this dip
                        mode.X = HigherLocalMax(this.O, A, nr.X, this.Range, this.Epsilon);
                    } //if (hs > 0)
                    else
                    {
                        mode.X = nr.X;
                    } //if (hs > 0)

                    mode.Found = true;
                } //if (nr.Converged)
                else
                {
                    double[] xA = nr.Bounds.X;
                    double[] hA = nr.Bounds.H;
                    bool hOnBothSides = hA.Aggregate(1, (x, y) => x * Math.Sign(y)) < 0;
                    double diffX = xA.Diff()[0];
                    if ((diffX < this.Epsilon) && hOnBothSides)
                    {
                        mode.X = xA.Mean();
                        mode.Found = true;
                    } //if ((diffX < this.Epsilon) && hOnBothSides)
                    else
                    {
                        bool converged = false;
                        if ((nr.Bounds.X.Diff()[0] < this.Epsilon) && nr.Bounds.H.All(a => a < 0))
                        {
                            if (this.Range[0].IsFinite())
                            {
                                //# We try a little bit further to find the actual mode, 
                                //# since it looks like we are almost there...
                                converged = true;
                                MaxFastTrackObject tmp = MaxFastTrack(this.O, A, this.Range[1], nr.Bounds.X[0], nr.Bounds.H[0]); //# yes: 4th argument is bounds$h (and not hp!)

                                if (tmp.Converged)
                                {
                                    mode.X = tmp.X;
                                    mode.Found = true;
                                } //if (tmp.Converged)
                                else
                                {
                                    mode.X = this.Range[1];
                                    mode.Found = !this.InestimableLowerLimit;
                                } //if (tmp.Converged)
                            } //if (this.Range[0].IsFinite())
                            else
                            {
                                double x = nr.Bounds.X[0];
                                double h = nr.Bounds.H[0];
                                double tmp = Math.Abs(h / this.O.LogFSecond(x, this.A));
                                if (tmp < 1e-12)
                                {
                                    mode = new LMode(x, h, true);
                                    converged = true;
                                } //if (tmp < 1e-12)
                            } //if (this.Range[0].IsFinite())
                        }

                        if (!converged)
                        {
                            throw new WEException("Algorithm did not converge.");
                        } //if (!converged)

                    } //if ((diffX < this.Epsilon) && hOnBothSides)

                } //if (nr.Converged)

                if (mode.Found)
                {
                    mode.H = this.O.LogF(mode.X, this.A);
                }
            } //if (!mode.Found)

            double logFRatioRemote = Math.Log(this.FRatioRemote);// # new_0.11
            double target = Tools.ND; // R n'a pas vraiment de notion de scope
            if (mode.Found)
            {
                target = mode.H + logFRatioRemote;
            }

            this.A.Mode = mode.X; //  # some 'remote.*' functions need that information
                                  //# modif_0.11
                                  //# Find remote points on both sides of mode
            double remoteLeft = 0, remoteRight = 0;
            bool[] modeOnBorder = new bool[2];
            double fs = Tools.ND;
            if (mode.Found)
            {
                //TODO pourquoi la définition et l'affectation de target ne serait-elle pas ici?
                double[] x = new double[] { Tools.NA, Tools.NA };
                double[] f = new double[] { Tools.NA, Tools.NA };
                double[] fp = new double[] { Tools.NA, Tools.NA };
                fs = this.O.LogFSecond(mode.X, A);// # scalar
                double d = 10 / Math.Abs(fs);
                //# new_0.11
                double[] remoteStart = this.O.LogFInvRemote(target, this.A).ToArray();//  # vector
                                                                                      // TODO start est un paramètre de la fonction en R, vérifier l'impact de modifier la référence d'objet ici.
                this.Start = SplitStartingValuesLeftRight(this.O, target, A, this.Range, remoteStart);// # vector of length 2

                modeOnBorder = new bool[2] { this.Range[0] == mode.X, this.Range[1] == mode.X }; //# new_0.11
                List<int> jSeq = new List<int>();
                for (int i = 0; i < 2; i++)
                {
                    if (!modeOnBorder[i] && this.Start[i].IsNA())
                    {
                        jSeq.Add(i);
                    }
                }

                foreach (int j in jSeq)
                {
                    //# modif_0.11 (first block that appeared here was moved above)
                    // condition inutile
                    if (Tools.IsNA(this.Start[j]))
                    {
                        int dir = (j == 0) ? -1 : 1;
                        double tmp = mode.X + dir * d;
                        bool accepted = (tmp > this.Range[0]) && (tmp < this.Range[1]);
                        if (!accepted)
                        {
                            double myD = d;
                            while (!accepted)
                            {
                                myD = myD / 2;
                                tmp = mode.X + dir * myD;
                                accepted = (tmp > this.Range[0]) && (tmp < this.Range[1]);
                            }
                        }

                        fp[0] = this.O.LogFPrime(tmp, A);
                        accepted = fp[0].IsFinite();
                        while (!accepted)
                        {
                            tmp = (tmp + mode.X) / 2;
                            fp[0] = this.O.LogFPrime(tmp, A);
                            accepted = fp[0].IsFinite();
                        }

                        double a = fp[0] / d; // # estimated rate of slope acceleration
                        double d2 = Math.Sqrt(2 * Math.Abs(Math.Log(this.FRatioRemote) / a));// # estimated distance we need to get from mode, 
                                                                                             //# at the rate above, to reach the f.ratio.remote distance
                        double x2 = mode.X + dir * d2;
                        accepted = (x2 > this.Range[0]) && (x2 < this.Range[1]);
                        while (!accepted)
                        {
                            d2 = d2 / 2;
                            x2 = mode.X + dir * d2;
                            accepted = x2 > this.Range[0] & x2 < this.Range[1];
                        }

                        double fp2 = this.O.LogFPrime(x2, A);
                        accepted = fp2.IsFinite();

                        while (!accepted)
                        {
                            d2 = d2 / 2;
                            x2 = mode.X + dir * d2;
                            fp2 = this.O.LogFPrime(x2, A);
                            accepted = Tools.IsFinite(fp2);
                        }

                        x = new double[2] { tmp, x2 };
                        f = new double[2] { Tools.NA, Tools.NA };
                        for (int i = 0; i < 2; i++)
                        {
                            f[i] = this.O.LogF(x[i], A);
                        }
                        fp[1] = fp2;

                        //# Approximating the shape of log.f.prime by a straight line, 
                        //# we can estimate the point where the distance
                        //# prescribed by f.ratio.remote is reached 

                        double[,] xSolved = Matrix2x2.Solve2x2(new double[] { x[0], x[1], 1, 1 }, byCol: true); // # 2 x 2 matrix
                        double[] theta = Matrix2x2.Product(xSolved, fp);
                        double qA = theta[0];
                        double qB = theta[1];
                        // ifelse(xor(j==1, x[1]>x[2]), 2, 1)
                        int innerSide = (j == 0) ^ (x[0] > x[1]) ? 1 : 0;
                        double innerX = x[innerSide];
                        double qC = -qA * Math.Pow(innerX, 2) - qB * innerX - (target - f[innerSide]);

                        double l;
                        double u;
                        if (j == 0)
                        {
                            l = this.Range[0];
                            u = mode.X;
                        }
                        else
                        {
                            l = mode.X;
                            u = this.Range[1];
                        }

                        List<double> roots = QuadraticEquation.GetRealRoots(new double[] { qC, qB, qA }, l: l, u: u, ROrder: true);// # modif_0.11
                        switch (roots.Count)
                        {
                            case 0:
                                this.Start[j] = x[1 - innerSide];
                                break;
                            case 1:
                                this.Start[j] = roots[0];
                                break;
                            default:
                                if (j == 1)
                                {
                                    this.Start[j] = roots.Max();
                                }
                                else
                                {
                                    this.Start[j] = roots.Min();
                                }
                                break;
                        }
                    }
                }
                //# new_0.11
                //# Suggest another starting point (on both sides)
                //# based on the 2nd degree polynomial approximation (Taylor series devpmt) of log.f

                f = new double[0];
                if (!modeOnBorder.Any(pX => pX))
                {
                    double delta = Math.Sqrt(2 * logFRatioRemote / fs);
                    double[] altStart;

                    //# Left 
                    double xScalar = mode.X - delta;
                    if (xScalar > this.Range[0])
                    {
                        altStart = new double[] { this.Start[0], xScalar };
                        f = new double[] { Math.Abs(this.O.LogF(altStart[0], A) - target), Math.Abs(this.O.LogF(altStart[1], A) - target) };
                        IEnumerable<double> lstD = f.Substract(target).Abs();
                        int w = (lstD.First() < lstD.Last()) ? 0 : 1;
                        this.Start[0] = altStart[w];
                    }

                    //# Right side
                    xScalar = mode.X + delta;
                    if (xScalar < this.Range[1])
                    {
                        altStart = new double[] { this.Start[1], xScalar };
                        f = new double[] { Math.Abs(this.O.LogF(altStart[0], A) - target), Math.Abs(this.O.LogF(altStart[1], A) - target) };
                        IEnumerable<double> lstD = f.Substract(target).Abs();
                        int w = (lstD.First() < lstD.Last()) ? 0 : 1;
                        this.Start[1] = altStart[w];
                    }
                }
                //oh <- list(h=o$log.f, hp=o$log.f.prime, hs=o$log.f.second)
                SecuredNRA oh = new SecuredNRA(this.O.LogF, this.O.LogFPrime, this.O.LogFSecond);// # modif_0.11
                if (Tools.IsNA(this.Start[0]))
                {
                    remoteLeft = mode.X;
                }
                else
                {
                    double[] range = new double[] { this.Range[0], mode.X };
                    //LMaxPoint lm = new LMaxPoint(mode);
                    nr = new SecuredNRSearch(oh, A, this.Start[0], this.Epsilon, range, 
                    maxPoint: mode, target: target, inestimableLowerLimit: this.InestimableLowerLimit);
                    // le dernier paramètre a une valeur par défaut: expectedHpSign = -1

                    if (nr.Converged)
                    {
                        remoteLeft = nr.X;
                    }
                    else if (this.Range[0].IsFinite())
                    {
                        remoteLeft = this.Range[0];
                    }
                    else
                    {
                        //TODO Exception
                        throw new WEException("Fin temporaire!");
                    }
                }

                if (Tools.IsNA(this.Start[1]))
                {
                    remoteRight = mode.X;
                } // if
                else
                {
                    this.Range = new double[] { mode.X, this.Range[1] };
                    nr = new SecuredNRSearch(oh, A, this.Start[1], this.Epsilon, this.Range,
                    maxPoint: mode, target: target);
                    // les 2 derniers paramètres prennent des valeurs par défaut: 
                    // inestimableLowerLimit = false, expectedHpSign = -1
                    if (nr.Converged)
                    {
                        remoteRight = nr.X;
                    }
                    else
                    {
                        throw new WEException("Algorithm did not converge.");// # should not happen
                    }
                }
            }
            else
            {
                //# new_0.11
                //# mode was not found
                remoteLeft = Tools.NA;
                remoteRight = Tools.NA;
            }

            //list(mode = mode, remote = c(remote.left, remote.right))

            LModeExt lmext = new LModeExt(mode, remoteLeft, remoteRight);
            return lmext;
        } //# end of ref.points
    }
} //# end of ref.points
