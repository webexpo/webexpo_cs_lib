namespace Zygotine.WebExpo
{
    using System;
    using Zygotine.Util;
    using Zygotine.Statistics.Distribution;

    public class Icdf
    {
        internal SGNFnAParam A { get; private set; }
        // TODO s'assurer que la lazy evaluation ne vient pas fucker le chien dans tout cela.
        public bool InestimableLowerLimit { get; private set; } = false;
        public int NRMaxIter { get; private set; } = 200;
        internal GenObject Ob01 { get; private set; }
        public double Precision { get; private set; } = 1E-6;
        public double[] Range { get; private set; }
        public double[] Start { get; private set; }
        public int UCount = 0;
        // TODO localiser l'utilisation de U, et si utiliser l'initialiser à runif()
        // On assume qu'une valeur réelle de u n'est jamais passée en paramètre.
        // À vérifier.
        private double theRealU = Tools.ND; // not defined
        private double U
        {
            get
            {
                if (theRealU.IsND())
                {
                    theRealU = UniformDistribution.RUnif();
                }

                return theRealU;
            }
        }


        internal Icdf(GenObject o, SGNFnAParam a, double[] range, double precision = 1E-6, int NRMaxIteration = 200)
        {
            this.Ob01 = o;
            this.A = a == null ? o.A : a;
            this.Range = (range == null) ? new double[0] : range;
            this.Start = Tools.Combine(Range.Mean()); // donne NaN si range est de longueur 0.
            this.Precision = precision;
            this.NRMaxIter = NRMaxIteration;
        }

        public double Bidon(double start)
        {
            this.Start = Tools.Combine(start);
            return Bidon();
        }
        public double Bidon(double[] start)
        {
            this.Start = start;
            return Bidon();
        }

        public double Bidon(bool inestLowerLim)
        {
            this.InestimableLowerLimit = inestLowerLim;
            return Bidon();
        }

        public double Bidon(double start, bool inestLowerLim)
        {
            this.InestimableLowerLimit = inestLowerLim;
            return Bidon(start);
        }

        public double Bidon(double[] start, bool inestLowerLim)
        {
            this.InestimableLowerLimit = inestLowerLim;
            return Bidon(start);
        }

        public double Bidon()
        {

            ReferencePoints rp = new ReferencePoints(this.Ob01, this.A, this.Start.Copy(), this.Range, this.InestimableLowerLimit, epsilon: this.Precision);
            LModeExt refPtsMode = rp.GetPoints();
            bool cntn = false;
            double uMode;
            double x;
            double target;
            double areaMode;
            bool distrnLeftSide = false;
            double mathLowerLimit = this.Ob01.MathLowerLimit;
            bool converged = false;
            double areaX;
            double precision = this.Precision;
            Func<double, SGNFnAParam, double> area;
            if (refPtsMode.Found)
            {
                // A$f < -o$f
                this.A.F = this.Ob01.F;
                this.A.LowerLimit = refPtsMode.Remote[0];
                this.A.M = refPtsMode.H; // # (kind of a) standardizing constant
                area = this.A.AreaFnFromALLimToX;
                double areaTot = WebExpoFunctions3.SmoothedAreaEstimate(
                    area,
                    refPtsMode.Remote[1],
                    this.A,
                    mathLowerLimit,
                    refPtsMode.Remote,
                    this.InestimableLowerLimit,
                    remoteRight: true);

                target = this.U * areaTot;
                precision = this.Precision * areaTot;
                x = refPtsMode.X;
                areaX = 0;
                areaMode = WebExpoFunctions3.SmoothedAreaEstimate(area, x, this.A, mathLowerLimit, refPtsMode.Remote, this.InestimableLowerLimit);
                uMode = areaMode / areaTot;
                converged = Math.Abs(areaMode - target) < precision;
                cntn = !converged;
            }
            else
            {
                throw new WEException("Mode not found");
            }

            if (cntn)
            {
                distrnLeftSide = this.U < uMode;
                if (distrnLeftSide)
                {
                    area = this.A.AreaFn_FromXToAULim;
                    this.A.UpperLimit = refPtsMode.X;
                    this.Range = Tools.Combine(refPtsMode.Remote[0], refPtsMode.X);
                }
                else
                {
                    this.A.LowerLimit = refPtsMode.X;
                    this.Range = Tools.Combine(refPtsMode.X, refPtsMode.Remote[1]);
                }
            }

            //converged = true;
            if (cntn && this.Ob01.PotentiallyBimodal)
            {
                // Slower but safer algorithm for posterior distrns that are potentially bimodal
                double fMode = A.F(refPtsMode.X, A);
                double start = refPtsMode.X + (target - areaMode) / fMode;
                cntn = false;
                converged = true; // voir l'impact !!!

                double highestPoint;
                double hMax;
                if (distrnLeftSide)
                {
                    highestPoint = refPtsMode.Remote[0];
                    hMax = areaMode;
                    target = (uMode - this.U) / uMode * hMax;
                }
                else
                {
                    highestPoint = refPtsMode.Remote[1];
                    hMax = WebExpoFunctions3.SmoothedAreaEstimate(area, refPtsMode.Remote[1], this.A, mathLowerLimit, this.Range, this.InestimableLowerLimit, remoteRight: true);
                    target = (this.U - uMode) / (1 - uMode) * hMax;
                }

                SecuredNRA oCum = new SecuredNRA(area, this.Ob01.F);
                SecuredNRSearch nr = new SecuredNRSearch(oCum, this.A, start, precision, this.Range,
                    maxPoint: new LMode(highestPoint, hMax), target: target, inestimableLowerLimit: this.InestimableLowerLimit);
                // le dernier paramètre a une valeur par défaut: expectedHpSign = -1
                if (nr.Converged)
                {
                    x = nr.X;
                }
                else
                {
                    //# it seems we were caught in an infinite loop, had we not limited the number of iterations;
                    //# see if we can rule that problem out
                    x = nr.Bounds.X[0];
                    double h = nr.Bounds.H[0];
                    x = ICDFFunctions.PingpongTieBreaker(area, this.A, x, h, target, mathLowerLimit, this.Range, this.InestimableLowerLimit, precision, distrnLeftSide);
                }
            }

            int direction = 0;
            int count = 0;
            if (cntn)
            {

                direction = distrnLeftSide ? -1 : 1;
                if (distrnLeftSide)
                {
                    target = (uMode - this.U) / uMode * areaMode;
                }
                else
                {
                    double areaRemoteRight = WebExpoFunctions3.SmoothedAreaEstimate(
                        area,
                        refPtsMode.Remote[1],
                        this.A,
                        mathLowerLimit,
                        this.Range,
                        this.InestimableLowerLimit,
                        remoteRight: true);
                    target = (this.U - uMode) / (1 - uMode) * areaRemoteRight;
                }
            }

            while (cntn)
            {
                double hp = this.Ob01.F(x, A);
                double change = (target - areaX) / hp * direction;
                x = x + change;
                areaX = area(x, A);
                count = count + 1;
                converged = Math.Abs(areaX - target) < precision;
                cntn = !converged && count <= this.NRMaxIter;
            }

            if (!converged)
            {
                if (x > refPtsMode.Remote[1] || x < refPtsMode.Remote[0])
                {
                    
                    throw new WEException("Stepped out of bounds -- due to numerical imprecision?");
                }
                else
                {
                    //# it seems we were caught in an infinite loop, had we not limited the number of iterations;
                    //# see if we can rule that problem out
                    x = ICDFFunctions.PingpongTieBreaker(area, this.A, x, areaX, target, mathLowerLimit, this.Range, this.InestimableLowerLimit, precision, distrnLeftSide);
                }
            }
            return x;
        }
    }
}


