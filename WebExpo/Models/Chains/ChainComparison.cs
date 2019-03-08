namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using Zygotine.Statistics;

    public class ChainComparison
    {
        private static readonly double[] defaultQuantileDefinition = new double[] { 0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99 };
        private static readonly double[] defaultTolerance = new double[] { 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.02, 0.03, 0.04 };
        public static double[] GetDefaultQuantile()
        {
            double[] rep = new double[ChainComparison.defaultQuantileDefinition.Length];
            Array.Copy(ChainComparison.defaultQuantileDefinition, rep, ChainComparison.defaultQuantileDefinition.Length);
            return rep;
        }

        private double[] quantileDefinition;
        private double[] tolerance;
        private double[] referenceChain;
        private double[] referenceQuantiles;

        public Quantile Quantile { get; private set; }

        public ChainComparison()
        {
            this.quantileDefinition = ChainComparison.defaultQuantileDefinition;
            this.tolerance = ChainComparison.defaultTolerance;
            this.Quantile = new Quantile(this.quantileDefinition);
        }

        private ChainComparison(double[] probs, double[] criteria)
        {
            if ((probs == null) || (criteria == null))
            {
                throw new ArgumentNullException("probs or criteria cannot be null array");
            }

            if ((probs.Length == 0) || (criteria.Length == 0))
            {
                throw new ArgumentException("probs or/and criteria array is empty.");
            }

            if (probs.Length != criteria.Length)
            {
                throw new ArgumentException("probs and criteria must be of the same length.");
            }

            for (int i = 0; i < probs.Length; i++)
            {
                if (probs[i] <= 0 || probs[i] >= 1)
                {
                    throw new ArgumentException("each probability must be within the interval (0, 1).");
                }

                if (criteria[i] <= 0)
                {
                    throw new ArgumentException("each criteria must be greater than 0.");
                }
            }

            this.quantileDefinition = new double[probs.Length];
            Array.Copy(probs, this.quantileDefinition, probs.Length);

            this.tolerance = new double[probs.Length];
            Array.Copy(criteria, this.tolerance, probs.Length);

            this.Quantile = new Quantile(this.quantileDefinition);
        }

        public ChainComparison(double[] probs, double[] criteria, double[] reference) : this(probs,criteria)
        {
            this.referenceChain = reference;
            this.referenceQuantiles = this.Quantile.Compute(this.referenceChain);
        }

        public ChainComparison.Result[] Compare( double[] chainToCompare)
        {
            ChainComparison.Result[] rep = new ChainComparison.Result[this.tolerance.Length];
            double[] qForChainToCompare = this.Quantile.Compute(chainToCompare);
            for (int i = 0; i < this.referenceQuantiles.Length ; i++)
            {
                rep[i] = new ChainComparison.Result(this.quantileDefinition[i], this.tolerance[i], Math.Abs((this.referenceQuantiles[i] - qForChainToCompare[i]) / this.referenceQuantiles[i]));
            }
            return rep;
        }

        public string GetHeaders()
        {
            System.Text.StringBuilder sb = new System.Text.StringBuilder();
            for (int i = 0; i < this.quantileDefinition.Length; i++)
            {
                sb.AppendFormat(System.Globalization.CultureInfo.InvariantCulture, "\t Q{0} p={1}/Δ={2:P0}", i, this.quantileDefinition[i], this.tolerance[i]);
            }

            return sb.ToString();
        }

        public static ChainComparison[] GetComparator(double[] probs, double[] criteria, string pathToReferenceFile)
        {
            System.IO.StreamReader srReference = new System.IO.StreamReader(pathToReferenceFile);

            List<List<double>> referenceChains = new List<List<double>>();
            if (srReference.EndOfStream)
            {
                srReference.Close();
                throw new IOException(string.Format("\"{0}\" is empty!", pathToReferenceFile));
            }

            string lineA = srReference.ReadLine();
            string[] tmpHeader = lineA.Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);

            int columnCount = tmpHeader.Length;
            referenceChains = new List<List<double>>();
            for (int columnIndex = 0; columnIndex < columnCount; columnIndex++)
            {
                referenceChains.Add(new List<double>());
            }

            if (srReference.EndOfStream)
            {
                srReference.Close();
                throw new IOException(string.Format("\"{0}\" contains only headers!", pathToReferenceFile));
            }

            while (!srReference.EndOfStream)
            {
                lineA = srReference.ReadLine();
                string[] tmpA = lineA.Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                for (int columnIndex = 0; columnIndex < columnCount; columnIndex++)
                {
                    referenceChains[columnIndex].Add(double.Parse(tmpA[columnIndex]));
                }
            }

            ChainComparison[] comp = new ChainComparison[columnCount];
            for(int columnIndex = 0; columnIndex < columnCount; columnIndex++)
            {
                comp[columnIndex] = new ChainComparison(probs, criteria, referenceChains[columnIndex].ToArray());
            }

            return comp;
        }

        public class Result
        {

            public double Probability { get; private set; }
            public double Value { get; private set; }
            public double Tolerance { get; private set; }
            public bool Valid { get; private set; } = false;

            internal Result(double qdef, double tolerance, double value)
            {

                this.Probability = qdef;
                this.Tolerance = tolerance;
                this.Value = value;
                this.Valid = this.Value <= this.Tolerance;
            }
        }
    }
}
