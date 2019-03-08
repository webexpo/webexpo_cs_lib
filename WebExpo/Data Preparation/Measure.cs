namespace Zygotine.WebExpo
{
    using System.Text.RegularExpressions;
 
    public abstract class Measure
    {
        public double A { get; internal set; } = double.NaN;
        public double B { get; internal set; } = double.NaN;
        //public double GenVal { get; internal set; } = double.NaN;

        public MeasureType Censoring { get; protected set; } = MeasureType.Uncensored;
        public Worker Worker { get; internal set; } = null;
        //public double Uncertainty { get; internal set; } = 0.0;
        public int Ordinal { get; internal set; }
        public int WorkerDataOrdinal { get; internal set; } // pour assurer le lien avec R dont les données sont 'ordonnées' comme suit:
                                                            // Uncensored, GT, LT  et Interval. Voir WorkerData

        internal void SetWorker(string worker)
        {
            if (worker != null)
            {
                worker = Regex.Replace(worker, @"[\s]", "");
                if (worker != string.Empty)
                    this.Worker = new Worker(worker);
            }
        }

        internal abstract Measure Duplicate();
 
        public enum MeasureType { Uncensored = 0, GreaterThan, LessThan, Interval };

    }





}
