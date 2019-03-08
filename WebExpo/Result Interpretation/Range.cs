namespace Zygotine.WebExpo
{
    public class Range
    {
        public double Minimum { get; set; } = double.NaN;
        public double Maximum { get; set; } = double.NaN;
        public double InitialValue { get; internal set; } = double.NaN;

        protected Range()
        {
        }

        public Range(double min, double max)
        {
            this.Minimum = min;
            this.Maximum = max;
        }

        public Range(double min, double max, double initValue) : this(min,max)
        {
            this.InitialValue = initValue;
        }

        public static MeasurementErrorRange GetCVRange(double min, double max)
        {
            return new MeasurementErrorCVRange(min, max);
        }
    }

    public abstract class MeasurementErrorRange : Range
    {
        public bool Exists { get; internal set; } = false;
        public bool ThroughSD { get; internal set; } = false;
        public bool ThroughCV { get; internal set; } = false;
        public bool Known { get; internal set; } = false;

        internal MeasurementErrorRange(double min, double max) : base(min, max)
        {
            this.Known = (max - min) == 0;
            this.InitialValue = (max + min) / 2;
        }
    }

    public class MeasurementErrorCVRange : MeasurementErrorRange
    {
        internal MeasurementErrorCVRange(double min, double max) : base(min, max)
        {
            this.Exists = true;
            this.ThroughCV = true;
        }
    }

    public class CoefficientOfVariationRange : Range
    {
        public CoefficientOfVariationRange(double min, double max) : base(min, max)
        {
        }
    }
}
