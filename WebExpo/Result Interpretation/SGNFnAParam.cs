namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using Zygotine.Util;

    internal class SGNFnAParam
    {
        // A parameter for sigma.gen.object
        public double B { get; set; } = Tools.ND;
        public double LowerLimit { get; internal set; } = Tools.ND;
        public double UpperLimit { get; internal set; } = Tools.ND;
        internal double LS2 { get; set; } = Tools.ND;
        internal double M { get; set; } = Tools.ND;
        internal double Mode { get; set; } = Tools.ND;
        internal int N { get; set; } = 0;
        internal double S { get; set; } = Tools.ND;
        internal double[] Range { get; set; } = new double[0];
        internal double[] Truevalues { get; set; } = new double[0];
        public Func<double, SGNFnAParam, double> F { get; internal set; }
        public double LM { get; internal set; } = Tools.ND;
        public double Mu { get; internal set; } = Tools.ND;
        public double[] Muk { get; internal set; } = new double[0];
        public double S2 { get; internal set; } = Tools.ND;
        public double Y { get; internal set; } = Tools.ND;
        public double Ksi { get; internal set; } = Tools.ND;
        public double Ksi2 { get; internal set; } = Tools.ND;
        public double Sigma2 { get; internal set; } = Tools.ND;
        public double LogY { get; internal set; } = Tools.ND;
        public double CV2 { get; internal set; } = Tools.ND;
        public double Y2 { get; internal set; } = Tools.ND;
        public double MuMean { get; internal set; } = Tools.ND;

        /**
         * Property F (Func<double, SGNFnAParam, double>) must be defined
         */
        public double AreaFnFromALLimToX(double x, SGNFnAParam A)
        {
            Func<double, double> preArea = delegate (double xParm)
            {
                return this.F(xParm, A);
            };
            Numerics.Integration.NumericIntegration.Results res;
            res = Numerics.Integration.NumericIntegration.Integrate(preArea, A.LowerLimit, x);
            return res.Value;
        }

        public double AreaFn_FromXToAULim(double x, SGNFnAParam A)
        {
            Func<double, double> preArea = delegate (double xParm)
             {
                 return this.F(xParm, A);
             };
            Numerics.Integration.NumericIntegration.Results res;
            res = Numerics.Integration.NumericIntegration.Integrate(preArea, x, A.UpperLimit);
            return res.Value;
        }

        public Func<double, double> AreaFromXToUpperLimitFn(SGNFnAParam A)
        {
            Func<double, double> dummy = null;
            if (this.F == null)
            {
                dummy = delegate (double x)
                {
                    return Tools.NA;
                };
            }
            else
            {
                Func<double, double> area = delegate (double xParm)
                {
                    return F(xParm, A);
                };
                dummy = delegate (double x)
                {
                    Numerics.Integration.NumericIntegration.Results res;
                    res = Numerics.Integration.NumericIntegration.Integrate(area, x, A.UpperLimit);
                    return res.Value;
                };
            }
            return dummy;
        }

        internal SGNFnAParam()
        {
        }

        internal SGNFnAParam Clone()
        {
            SGNFnAParam a = (SGNFnAParam)this.MemberwiseClone();

            if (this.Range == null)
            {
                a.Range = null;
            }
            else if (this.Range.Length == 0)
            {
                a.Range = new double[0];
            }
            else
            {
                a.Range = this.Range.Copy();
            }

            if (this.Truevalues == null)
            {
                a.Truevalues = null;
            }
            else if (this.Truevalues.Length == 0)
            {
                a.Truevalues = new double[0];
            }
            else
            {
                a.Truevalues = this.Truevalues.Copy();
            }

            return a;
        }

        public override string ToString()
        {
            return this.Show();
        }

        public string Show()
        {
            List<string> ls = new List<string>();
            if (!double.IsNaN(this.B))
            {
                ls.Add(string.Format("B={0:R}", this.B));
            }

            if (!double.IsNaN(this.CV2))
            {
                ls.Add(string.Format("CV2={0:R}", this.CV2));
            }

            if (!double.IsNaN(this.Ksi))
            {
                ls.Add(string.Format("Ksi={0:R}", this.Ksi));
            }

            if (!double.IsNaN(this.Ksi2))
            {
                ls.Add(string.Format("Ksi2={0:R}", this.Ksi2));
            }

            if (!double.IsNaN(this.LM))
            {
                ls.Add(string.Format("LM={0:R}", this.LM));
            }

            if (!double.IsNaN(this.LogY))
            {
                ls.Add(string.Format("LogY={0:R}", this.LogY));
            }


            if (!double.IsNaN(this.LowerLimit))
            {
                ls.Add(string.Format("LowerLimit={0:R}", this.LowerLimit));
            }

            if (!double.IsNaN(this.LS2))
            {
                ls.Add(string.Format("LS2={0:R}", this.LS2));
            }

            if (!double.IsNaN(this.M))
            {
                ls.Add(string.Format("M={0:R}", this.M));
            }

            if (!double.IsNaN(this.Mode))
            {
                ls.Add(string.Format("Mode={0:R}", this.Mode));
            }

            if (!double.IsNaN(this.Mu))
            {
                ls.Add(string.Format("Mu={0:R}", this.Mu));
            }

            if (this.Muk.Length > 0)
            {
                ls.Add(string.Format("Muk={0}", this.Muk.ShowR()));
            }

            if (this.N > 0)
            {
                ls.Add(string.Format("N={0}", this.N));
            }

            if (this.Range.Length > 0)
            {
                ls.Add(string.Format("Range={0}", this.Range.ShowR()));
            }

            if (!double.IsNaN(this.S))
            {
                ls.Add(string.Format("S={0:R}", this.S));
            }

            if (!double.IsNaN(this.S2))
            {
                ls.Add(string.Format("S2={0:R}", this.S2));
            }

            if (this.Truevalues.Length > 0)
            {
                ls.Add(string.Format("Truevalues={0}", this.Truevalues.ShowR()));
            }

            if (!double.IsNaN(this.UpperLimit))
            {
                ls.Add(string.Format("UpperLimit={0:R}", this.UpperLimit));
            }

            if (!double.IsNaN(this.Y))
            {
                ls.Add(string.Format("Y={0:R}", this.Y));
            }

            if (!double.IsNaN(this.Y2))
            {
                ls.Add(string.Format("Y2={0:R}", this.Y2));
            }

            return string.Join(", ", ls);
        }
    }
}
