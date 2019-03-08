namespace Zygotine.WebExpo
{
    using System;
    using System.Linq;
    using System.Collections.Generic;
    using Zygotine.Numerics;
    using Zygotine.Util;

    internal class SecuredNRSearch
    {
        //Semblable à ParamB01, X ici est un tableau ...
        internal double X { get; private set; }
        internal bool Converged { get; private set; }
        internal Bounds Bounds { get; private set; }

        private class L1
        {
            public double X { get; internal set; } = Tools.NA;
            public double H { get; internal set; } = Tools.NA;
            public double Hp { get; internal set; } = Tools.NA;
            public L1()
            {
            }

            public L1(double x, double h, double hp)
            {
                this.X = x;
                this.H = h;
                this.Hp = hp;
            }

        }

        internal class L2
        {
            public double X { get; internal set; }
            public double H { get; internal set; }
            public double HMax { get; internal set; }
            public bool Converged { get; internal set; }

            public L2(double x, double hMax)
            {
                this.X = x;
                this.HMax = hMax;
            }
        }

        private class LWhere
        {
            public bool Within { get; internal set; }
            public bool Below { get; internal set; }
            public bool Above { get; internal set; }

            public LWhere(bool within, bool below, bool above)
            {
                this.Within = within;
                this.Below = below;
                this.Above = above;
            }
        }

        private class LLimit
        {
            public double[] X = null;
            internal double[] H { get; set; } = new double[] { Tools.ND, Tools.ND };
            internal double[] Hp { get; set; } = new double[] { Tools.ND, Tools.ND };
            internal bool[] Usable { get; set; } = new bool[2] { false, false };
            internal bool[] Chckd { get; set; } = new bool[2] { false, false };
        }

        //leftbump.characteristics
        private L2 LeftBumpCharacteristics(
                ISecuredNR o,
                double start,
                double hpStart,
                double target,
                SGNFnAParam A,
                double lowerSideLimit,
                double h2Modeless,
                L1 lsStart,
                double epsilon = 1e-8)
        {
            double[] lim = new double[] { lowerSideLimit, start };
            L2 repOut = new L2(0.0, 0.0);

            //# Find maximum point

            //double[] lim = new double[] { lowerSideLimit, start };
            double h;

            if (!Tools.IsNA(lsStart.X))
            {
                lim[0] = lsStart.X;
            }

            double x = start;
            double hp = hpStart;
            bool cntn = true;
            bool converged = false;
            while (cntn)
            {
                if (hp < 0)
                {
                    lim[1] = x;
                }
                else
                {
                    lim[0] = x;
                }

                double hs = -Math.Abs(o.HSecond(x, A));
                double change = hp / hs;
                x = x - change;

                if (x > lim[1])
                {
                    x = (x + change + lim[1]) / 2;
                }
                else if (x < lim[0])
                {
                    x = (x + change + lim[0]) / 2;
                }

                hp = o.HPrime(x, A);
                converged = Math.Abs(hp) < epsilon;
                cntn = !converged;

                if (cntn)
                {
                    h = o.H(x, A);
                    cntn = h < h2Modeless;
                }
            }

            if (converged)
            {
                repOut.HMax = o.H(x, A);
                cntn = repOut.HMax > target;
            }

            if (cntn)
            {
                //# Find bump-left-side solution
                lim = new double[] { lowerSideLimit, x };
                if (!Tools.IsNA(lsStart.X))
                {
                    x = lsStart.X;
                    h = lsStart.H;
                }
                else
                {
                    double hs = o.HSecond(x, A);
                    double change = (hp - 1) / hs;
                    x = x - change;
                    if (x < lim[0])
                    {
                        x = (x + change + lim[0]) / 2;
                    }

                    h = o.H(x, A);
                }

                while (cntn)
                {
                    hp = o.HPrime(x, A);
                    double change = (h - target) / hp;
                    x = x - change;
                    if (x > lim[1])
                    {
                        x = (x + change + lim[1]) / 2;
                    }
                    else if (x < lim[0])
                    {
                        x = (x + change + lim[0]) / 2;
                    }

                    h = o.H(x, A);
                    cntn = Math.Abs(h - target) > epsilon;
                }

                repOut.X = x;
                repOut.H = h;
            }

            repOut.Converged = converged;
            return repOut;
        } //# end of leftbump.characteristics


        //# new_0.11
        private ParamB01 LeftSideFastScan(Bounds bounds, ISecuredNR o, SGNFnAParam A, double target, double lowerSideLimit, double h2Modeless, L1 lsStart = null, double epsilon = 1e-8, int lowerSide = 0) 
        {
            //# in input bounds, the slopes (hp values) must be negative & positive, respectively
            //# (thus pointing towards a 'dip', or local minimum)

            if (lsStart == null)
            {
                lsStart = new L1();
            }
            int higherSide = 1 - lowerSide; // correction juillet 2017
            L1 dip = DipCharacteristics(bounds.Copy(), o, A, epsilon);
            L2 leftBump = LeftBumpCharacteristics(o, bounds.X[lowerSide], bounds.Hp[lowerSide], target, A, lowerSideLimit, h2Modeless, lsStart);

            double x, h, hp;
            bool cntn = false;

            if (leftBump.Converged)
            {
                bool foundSoln = leftBump.HMax > target;
                cntn = !foundSoln;
                //TODO SUIVRE LE NA
                if (foundSoln)
                {
                    x = leftBump.X;
                    h = leftBump.H;
                    //TODO SUIVRE LE NA
                    hp = Tools.NA;
                }
                else
                {
                    x = dip.X;
                    h = dip.H;
                    hp = (bounds.X[higherSide] - h) / (bounds.X[higherSide] - x);
                }
            }
            else
            {
                x = lowerSideLimit;
                h = target;
                hp = Tools.NA;
                cntn = false;
            }

            return new ParamB01(x, h, hp, cntn);
        } //# end of leftSide.fastScan

        private L1 DipCharacteristics(Bounds b, ISecuredNR o, SGNFnAParam A, double epsilon = 1e-8)
        {
            bool cntn = true;
            double x = 0, h = 0, hp = 0;
            while (cntn)
            {
                //# intersection between the tangents at both ends of bounds (b)
                double bDiff = b.Hp.Diff()[0];
                if (bDiff == 0)
                {
                    x = (b.X.Mean());
                }
                else
                {
                    x = (b.Hp.Multiply(b.X).Diff()[0] - b.H.Diff()[0]) / bDiff;
                    if (x <= b.X[0] | x >= b.X[1])
                    {
                        x = b.X.Mean();
                    }
                }
                h = o.H(x, A);
                hp = o.HPrime(x, A);
                int side = hp < 0 ? 0 : 1;
                b.X[side] = x;
                b.H[side] = h;
                b.Hp[side] = hp;
                cntn = b.X.Diff()[0] > epsilon;
            }
            return new L1(x, h, hp);
            //TODO Voir si l'on peut partager a l'intérieur de SecuredNRSearch, en utilisant des champs, x, h et hp et possiblement d'autres
            //list(x = x, h = h, hp = hp)
        } //# end of dip.characteristics

        private ParamB01 NewBound(double start, Bounds bounds, bool monotonic, ISecuredNR o, SGNFnAParam A, double target, double[] range, bool above)
        {
            double x = start;
            double h = o.H(x, A);
            double hp = o.HPrime(x, A);
            bool accepted = h.IsFinite() && hp.IsFinite();

            if (!accepted)
            {
                //# Try a cubic extrapolation if we have two points in bounds
                int nBounds = bounds.X.Where(z => !Tools.IsNA(z)).Count();
                if (nBounds == 2)
                {
                    double tmp = CubicExtrapolation.Extrapolate(target, bounds, monotonic, range);

                    if (!Tools.IsNA(tmp))
                    {
                        x = tmp;
                        h = o.H(x, A);
                        hp = o.HPrime(x, A);
                        accepted = h.IsFinite() && hp.IsFinite();
                    }
                }
            }

            if (!accepted)
            {
                double bound = 1;
                if (above)
                {
                    bound = bounds.X.Where(a => !Tools.IsNA(a)).Max(); // peut-être double.NegativeInfinity comme en R si la liste est vide.
                }
                else
                {
                    bound = bounds.X.Where(a => !Tools.IsNA(a)).Min(); // peut-être double.PositiveInfinity comme en R si la liste est vide.
                }

                while (!accepted)
                {
                    x = (x + bound) / 2;
                    h = o.H(x, A);
                    hp = o.HPrime(x, A);
                    accepted = h.IsFinite() && hp.IsFinite();
                }
            }


            return new ParamB01(x, h, hp);
        } //# end of new.bound

        private bool UnreachedTarget(double target, double h, double higherSide, bool monotonic, bool modeSearch, bool leftSide = false)
        {
            bool unreached = false;
            if (!monotonic || modeSearch)
            {
                unreached = false;
            }
            else
            {
                bool posSlope = higherSide == 1; // probablement 1
                unreached = posSlope ^ (h > target);
                if (leftSide)
                {
                    unreached = !unreached;
                }
            }
            return unreached;
        } // end of unreachedTarget

        private LWhere WithinBounds(double x, Bounds bounds)
        {
            IEnumerable<double> bXLst = bounds.X.Where(a => !Tools.IsNA(a));
            int nBounds = bXLst.Count();
            bool above, below;
            if (nBounds == 2)
            {
                above = x > bounds.X[1];
                below = x < bounds.X[0];
            }
            else
            {
                //TODO vérifier les implication si bXLst soit vide
                // on semble croire que si nBounds n'est pas deux il est 1 et non pas 0.
                /*
                 * 
                w =  which(!is.na(bounds.x))
                bound =  bounds.x[w]
                *
                * Dans ce qui précède bound correspond à bXLst
                */
                above = x > bXLst.First();
                below = !above;
            }

            bool within = !above & !below;
            return new LWhere(within, below, above);

        } //# end of within.bounds

        internal static int MaxNIter { get; set;} = 500;
        internal SecuredNRSearch(
            ISecuredNR o,
            SGNFnAParam A,
            double start,
            double epsilon,
            double[] range,
            LMode maxPoint = null,
            double target = 0,
            bool inestimableLowerLimit = false,
            // int maxNIter = 500,
            int expectedHpSign = -1)

        //max.point=list(x= numeric(0), h= numeric(0)),
        {
            //# Throughout this function, 'bounds' is a list with the following elements/dimensions:
            //# x, h, hp, range:   numeric (vectors of length 2)
            //# higherSide:       numeric (1)
            //# IncludeSoln:      logical (1)

            //TODO l'accès à h2Modeless est difficile à comprendre en R qui permet l'accès à des variables non initialisé ... quitte à planter!
            //pour ce faire on introduit ici h2Modeless ...

            //-->secured.NR.search ,start=9.24215468292111 ,target=0 ,range=( 0, Inf ) ,inestimable.lower.limit=FALSE ,max.point$x= ,max.point$h= ,expected.hpsign=-1

            double h2Modeless = Tools.NA;

            if (maxPoint == null)
            {
                maxPoint = new LMode();
            }
            bool modeSearch = Tools.IsNA(maxPoint.X);
            bool monotonic = !modeSearch;
            LLimit limit = new LLimit();
            limit.X = range.Copy(); // les autres champs prennent les valeurs par défaut désirée
            limit.Chckd[0] = inestimableLowerLimit; //# new_0.11
            Bounds bounds = new Bounds(
                Tools.Rep(Tools.NA, 2), //x
                Tools.Rep(Tools.NA, 2), //h 
                Tools.Rep(Tools.NA, 2), //hp
                null, // range
                -1, // higherSide
                false); // includeSoln)

            bool found2Bounds;
            bool correct4HpSign;
            int higherSide;
            if (!Tools.IsNA(maxPoint.X))
            {
                higherSide = (start < maxPoint.X) ? 1 : 0; // higherSide est un indice
                bounds.X[higherSide] = maxPoint.X;
                bounds.H[higherSide] = maxPoint.H;
                bounds.Hp[higherSide] = 0.0;
                found2Bounds = true;
                expectedHpSign = higherSide == 0 ? -1 : 1;
                correct4HpSign = start > maxPoint.X;
                h2Modeless = (maxPoint.H * 2.0) - target;
            }
            else
            {
                found2Bounds = false;
                higherSide = 0;
                correct4HpSign = true;
            }

            int lowerSide = 1 - higherSide;
            double x = start;
            //# make sure that x is not out of range
            if (x < limit.X[0])
            {
                x = (maxPoint.X + limit.X[0]) / 2;
            }


            bounds.Range = range.Copy();
            bounds.HigherSide = higherSide;
            //# Register initial x into bounds
            double h = o.H(x, A);
            double hp = o.HPrime(x, A);
            int side;
            if (!Tools.IsNA(maxPoint.X))
            {
                bounds.IncludeSoln = h < target;
                side = lowerSide;
            }
            else
            {
                side = 0;
            }

            bounds.X[side] = x;
            bounds.H[side] = h;
            bounds.Hp[side] = hp;
            Visits visitedX = new Visits(SecuredNRSearch.MaxNIter);
            visitedX.Add(x);

            //# Do a first step
            bool cntn;
            if (correct4HpSign)
            {
                hp = expectedHpSign * Math.Abs(hp);
            }

            double change = (h - target) / hp;
            x = x - change;
            //# cntn until convergence
            bool unreachableMode = false;
            cntn = true;
            int count = 0;
            bool converged = false;
            int debCount = 0;
            while (cntn)
            {
                debCount++;
                if (x.IsNaN())
                {
                    if (1 == 1) { }
                }
                bool computedH = false;
                bool accepted = false;
                count = count + 1;
                LWhere wb; // = WithinBounds(x, bounds);
                while (!accepted)
                {
                     wb = WithinBounds(x, bounds);
                    if (bounds.IncludeSoln)
                    {
                        if (wb.Within)
                        {
                            accepted = true;
                        }
                        else
                        {
                            x = CubicExtrapolation.Extrapolate(target, bounds, monotonic, range);
                            accepted = !Tools.IsNA(x);
                        }
                    }
                    else if (wb.Above)
                    {
                        //# bounds do not include solution and we are looking above bounds
                        if (x > limit.X[1])
                        {
                            if (limit.Chckd[1])
                            {
                                accepted = true;
                                h = limit.H[1];

                                if (UnreachedTarget(target, h, higherSide, monotonic, modeSearch))
                                {
                                    //# force end of while-loop, as target is not reached even at this end
                                    x = limit.X[1];
                                    h = target;
                                    computedH = true;
                                }
                                else if (limit.Usable[1])
                                {
                                    x = bounds.Range[1];
                                    hp = limit.Hp[1];
                                    computedH = true;
                                }
                                else
                                {
                                    x = (x + change + bounds.Range[1]) / 2;
                                }
                            }
                            else
                            {
                                limit.Chckd[1] = true;
                                h = o.H(limit.X[1], A);
                                if (UnreachedTarget(target, h, higherSide, monotonic, modeSearch))
                                {
                                    //# force end of while-loop, as target is not reached even at this end
                                    x = limit.X[1];
                                    h = target;
                                    accepted = true;
                                    computedH = true;
                                }
                                else
                                {
                                    limit.H[1] = h;
                                    if (h.IsFinite())
                                    {
                                        hp = o.HPrime(limit.X[1], A);
                                        limit.Hp[1] = hp;
                                        limit.Usable[1] = hp.IsFinite();
                                        accepted = limit.Usable[1];
                                    }
                                    else
                                    {
                                        accepted = false;
                                    }

                                    if (accepted)
                                    {
                                        x = bounds.Range[1];
                                        computedH = true;
                                        unreachableMode = modeSearch && h > 0;//# new_0.11
                                    }
                                    else
                                    {
                                        x = (x + change + bounds.Range[1]) / 2;
                                    }
                                }
                            }
                        }
                        else
                        {
                            ParamB01 tmp = NewBound(x, bounds, monotonic, o, A, target, range, wb.Above);
                            x = tmp.X;
                            h = tmp.H;
                            hp = tmp.Hp;
                            computedH = true;
                            accepted = true;
                        }
                    }
                    else
                    {
                        //# bounds do not include solution and we are looking below bounds
                        if (x <= limit.X[0])
                        {
                            //# we are looking out of the variable domain (bad!)
                            if (!limit.Chckd[0])
                            {
                                limit.Chckd[0] = true;
                                h = o.H(limit.X[0], A);
                                limit.H[0] = h;
                                if (h.IsFinite())
                                {
                                    if (UnreachedTarget(target, h, higherSide, monotonic, modeSearch, leftSide: true))
                                    {
                                        //# force end of while-loop, as target is not reached even at this end
                                        x = limit.X[0];
                                        h = target;
                                        accepted = true;
                                    }
                                    else
                                    {
                                        hp = o.HPrime(limit.X[0], A);
                                        limit.Hp[0] = hp;
                                        limit.Usable[0] = hp.IsFinite();
                                        accepted = limit.Usable[0];
                                    }

                                    computedH = accepted;
                                    if (accepted)
                                    {
                                        x = limit.X[0];
                                    }
                                }
                                else
                                {
                                    limit.Usable[0] = false;
                                    x = (bounds.X.Where(a => !Tools.IsNA(a)).Min() + limit.X[0]) / 2;
                                }
                            }
                            else if (limit.Usable[0])
                            {
                                x = limit.X[0];
                                h = limit.H[0];
                                hp = limit.Hp[0];
                                computedH = true;
                                accepted = true;
                            }
                            else
                            {
                                //# not limit$usable[1]
                                x = (bounds.X.Min() + limit.X[0]) / 2;
                            }
                        }
                        else
                        {
                            ParamB01 tmp = NewBound(x, bounds, monotonic, o, A, target, range, wb.Above);
                            x = tmp.X;
                            h = tmp.H;
                            hp = tmp.Hp;
                            computedH = true;
                            accepted = true;
                            //    # new_0.11
                            unreachableMode = modeSearch && (((x == limit.X[1]) && (h > 0)) || ((x == limit.X[0]) && (h < 0)));
                        }
                    }
                } // while(!accepted)

                if (!computedH)
                {
                    h = o.H(x, A);
                    hp = o.HPrime(x, A);
                }

                converged = (Math.Abs(h - target) < epsilon) || unreachableMode;// # modif_0.11
                cntn = !converged;
                if (cntn)
                {
                    //# register results in bounds
                    bool hLtTarget = h < target;
                    if (bounds.IncludeSoln)
                    {
                        side = hLtTarget ? lowerSide : higherSide;
                    }
                    else
                    {
                        bounds.IncludeSoln = hLtTarget ^ (bounds.H[0] < target);

                        if (found2Bounds)
                        {
                            bool changeOppositeSideBounds = true;
                            if (modeSearch)
                            {
                                side = hLtTarget ? lowerSide : higherSide;
                                changeOppositeSideBounds = true;
                            }
                            else if (bounds.IncludeSoln)
                            {
                                //# bounds newly include solution (happening for the first time with current x)
                                wb = WithinBounds(x, bounds);
                                if (wb.Within)
                                {
                                    if (higherSide == 1)
                                    {
                                        side = lowerSide;
                                        changeOppositeSideBounds = false;
                                        if (hp < 0)
                                        {
                                            bounds.X[lowerSide] = x;
                                            bounds.H[lowerSide] = h;
                                            bounds.Hp[lowerSide] = hp;
                                            ParamB01 tmp = LeftSideFastScan(bounds, o, A, target, range[0], h2Modeless);
                                            x = tmp.X;
                                            h = tmp.H;
                                            hp = tmp.Hp;
                                            cntn = tmp.Cntn;
                                        }
                                    }
                                    else
                                    {
                                        throw new WEException("Scénario imprévu.");
                                    }
                                }
                                else if (wb.Below)
                                {
                                    side = 0;
                                    changeOppositeSideBounds = true;
                                }
                                else
                                {
                                    //# wb.above
                                    side = 1;
                                    changeOppositeSideBounds = true;
                                }
                            }
                            else
                            {
                                //# bounds still do not include solution
                                wb = WithinBounds(x, bounds);
                                if (wb.Within)
                                {
                                    changeOppositeSideBounds = false;
                                    if (h > bounds.H.Max())
                                    {
                                        side = higherSide;
                                    }
                                    else if (h < bounds.H.Min())
                                    {
                                        side = hp < 0 ? lowerSide : higherSide;
                                    }
                                    else
                                    {
                                        side = higherSide;
                                    }
                                }
                                else if (wb.Below)
                                {
                                    //Modifié à 1 
                                    if (higherSide == 1)
                                    {
                                        if (hp < 0)
                                        {
                                            bounds.X[lowerSide] = x;
                                            bounds.H[lowerSide] = h;
                                            bounds.Hp[lowerSide] = hp;
                                            ParamB01 tmp = LeftSideFastScan(bounds, o, A, target, range[lowerSide], h2Modeless);
                                            x = tmp.X;
                                            h = tmp.H;
                                            hp = tmp.Hp;
                                            cntn = tmp.Cntn;
                                            changeOppositeSideBounds = true;
                                        }
                                        else if (h < bounds.H.Min())
                                        {
                                            side = lowerSide;
                                            changeOppositeSideBounds = true;
                                        }
                                        else
                                        {
                                            L1 leftSidePotentialStartPoint = new L1(x, h, hp);
                                            double[] lim = new double[] { x, bounds.X[lowerSide] };
                                            x = lim.Mean();
                                            hp = o.HPrime(x, A);
                                            while (hp > 0)
                                            {
                                                h = o.H(x, A);
                                                side = h > bounds.H[lowerSide] ? 0 : 1;
                                                lim[side] = x;
                                                x = lim.Mean();
                                                hp = o.HPrime(x, A);
                                            }
                                            h = o.H(x, A);

                                            bounds.X[higherSide] = bounds.X[lowerSide];
                                            bounds.H[higherSide] = bounds.H[lowerSide];
                                            bounds.Hp[higherSide] = bounds.Hp[lowerSide];
                                            bounds.X[lowerSide] = x;
                                            bounds.H[lowerSide] = h;
                                            bounds.Hp[lowerSide] = hp;
                                            ParamB01 tmp = LeftSideFastScan(bounds, o, A, target, range[lowerSide], h2Modeless, leftSidePotentialStartPoint);
                                            x = tmp.X;
                                            h = tmp.H;
                                            hp = tmp.Hp;
                                            cntn = tmp.Cntn;
                                            side = higherSide;
                                            changeOppositeSideBounds = false;
                                        }
                                    }
                                    else
                                    {
                                        throw new WEException("Scénario imprévu.");
                                    }
                                }
                                else
                                {
                                    //# wb.above
                                    side = lowerSide;
                                    changeOppositeSideBounds = true;
                                }
                            }

                            if (changeOppositeSideBounds)
                            {
                                int oppositeSide = 1 - side;
                                bounds.X[oppositeSide] = bounds.X[side];
                                bounds.H[oppositeSide] = bounds.H[side];
                                bounds.Hp[oppositeSide] = bounds.Hp[side];
                            }
                        }//if(found2bounds)
                        else
                        {
                            //# !found2bounds
                            found2Bounds = true;
                            side = h < bounds.H[0] ? lowerSide : higherSide;
                            if (side == 0)
                            {
                                bounds.X[1] = bounds.X[0];
                                bounds.H[1] = bounds.H[0];
                                bounds.Hp[1] = bounds.Hp[0];
                            }
                        }
                    }

                    bounds.X[side] = x;
                    bounds.H[side] = h;
                    bounds.Hp[side] = hp;
                    visitedX.Add(x);
                    converged = (Math.Abs(h - target) < epsilon) || unreachableMode;// # modif_0.11
                    if (!converged)
                    {
                        cntn = count < SecuredNRSearch.MaxNIter;// # new_0.11: changed <= for <

                        if (!cntn)
                        {
                            //# Before we give up, we check one last thing:
                            //# if the last few steps were all leaning in the same direction,
                            //# then convergence may only be a question of time & patience!
                            //# Give it (yet) another chance!

                            //TODO BIEN VÉRIFIER LA LOGIQUE
                            Visits.DirectionChange[] directionChanges = visitedX.GetDirectionChanges();
                            if (directionChanges.Length > 0)
                            {
                                Visits.DirectionChange lastChangeDirection = directionChanges.Last(); //
                                IEnumerable<Visits.DirectionChange> inOppositeDirection = directionChanges.Where(a => a.Direction == -lastChangeDirection.Direction);
                                if (inOppositeDirection.Count() == 0)
                                {
                                    cntn = true;
                                }
                                else
                                {
                                    Visits.DirectionChange lastInOppositeDirection = inOppositeDirection.Last();
                                    if ((count - lastInOppositeDirection.Index) > 20)
                                    {
                                        cntn = true;
                                    }
                                }
                            }
                            else
                            {
                                //TODO ne devrait pas se produire.
                            }

                            if (cntn)
                            {
                                count = 0;
                                visitedX = new Visits(SecuredNRSearch.MaxNIter);
                            }
                        } //if(!cntn)
                    } // if (!converged)

                    if (cntn)
                    {
                        found2Bounds = true;
                        double absHp = Math.Abs(hp);
                        if (correct4HpSign)
                        {
                            hp = expectedHpSign * absHp;
                        }
                        change = (h - target) / hp;
                        double pctChange = Math.Abs(change) / bounds.X.Diff()[0];

                        if ((pctChange < 1e-4) && absHp > 1000)
                        {
                            double[] p;
                            if (change < 0)
                            {
                                p = new double[] { 0.9, 0.1 };
                            }
                            else
                            {
                                p = new double[] { 0.1, 0.9 };
                            }
                            x = p.Multiply(bounds.X).Sum();
                        }
                        else
                        {
                            x = x - change;
                        }
                    } //if (cntn)

                    //#cat('count=', count, '\n')
                    //#cat('x = ', x, '\n')
                    //#cat('h = ', h, '\n')
                    //#cat('target = ', target, '\n')
                    //#cat('abs(h-target)', abs(h-target), '\n')
                    //#cat('bounds.x = ', bounds.x, '\n')
                    //#cat('diff(bounds.x) = ', diff(bounds.x), '\n')
                    //#if (diff(bounds.x) < 0) cat('***************\n')
                } // if (cntn)
            } //while (cntn)

            this.X = x;
            this.Converged = converged;
            this.Bounds = bounds;
            //list(x = x, converged = converged, bounds = bounds);
        } //# end of secured.NR.search

        public override string ToString()
        {
            return string.Format("X ={0}, Converged ={1} Bounds: H ={2}, Hp ={3}, IncludeSoln ={4}, Range ={5}, HigherSide ={6}",
                    this.X,
                    this.Converged,
                    Tools.Show(this.Bounds.H),
                    Tools.Show(this.Bounds.Hp),
                    this.Bounds.IncludeSoln,
                    this.Bounds.Range,
                    this.Bounds.HigherSide);
        }
    }
}



