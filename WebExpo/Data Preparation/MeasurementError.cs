
namespace Zygotine.WebExpo
{
    public class MeasurementError : Range
    {
        public bool Any { get; internal set; } = false;
        public bool ThroughSD { get; internal set; } = false;
        public bool ThroughCV { get; internal set; } = false;
        public bool Known { get; internal set; } = false;
        public double Parm { get; internal set; }
        public double[] GetRange()
        {
            return new double[] { this.Minimum, this.Maximum };
        }

        public static readonly MeasurementError None = new MeasurementError();

        private MeasurementError() : base()
        {
        }

        public override string ToString()
        {
            if (this == MeasurementError.None)
            {
                return "";
            }
            return (this.ThroughCV ? "CV(" : "SD(") + this.Minimum.ToString() + (this.Known ? ")" : ", " + this.Maximum.ToString() + ")");
        }

        private static MeasurementError NewError(double min, double max)
        {
            MeasurementError me = new MeasurementError();
            me.Any = true;
            me.Minimum = min;
            me.Maximum = max;
            me.Known = min == max;
            me.InitialValue = (me.Minimum + me.Maximum) / 2.0;
            return me;
        }

        internal static MeasurementError NewCVError(double min, double max)
        {
            MeasurementError me = NewError(min, max);
            me.ThroughCV = true;
            return me;
        }

        internal static MeasurementError NewSDError(double min, double max)
        {
            MeasurementError me = NewError(min, max);
            me.ThroughSD = true;
            return me;
        }
    }
}
