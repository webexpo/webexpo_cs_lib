
namespace Zygotine.WebExpo
{
    using System.Globalization;
    using Zygotine.Util;

    internal class IntervalMeasure : UncensoredMeasure
    {
        internal IntervalMeasure(double a, double b, string workerTag = null) : base(a, workerTag)
        {
            if (!Tools.IsFinite(b))
                throw new WEException("Parameter b must be a finite number.");

            B = b;
            Censoring = MeasureType.Interval;
        }

        internal IntervalMeasure(double a, double b, Worker worker) : base(a, worker)
        {
            if (!Tools.IsFinite(b))
                throw new WEException("Parameter b must be a finite number.");

            this.B = b;
            this.Censoring = MeasureType.Interval;
        }

        internal override Measure Duplicate()
        {
            return new IntervalMeasure(A, B, Worker);
        }

        public override string ToString()
        {
            return ToString(MeasureList.FieldSeparator.Semicolon);
        }

        public override string ToString(MeasureList.FieldSeparator fs)
        {
            string fieldSeparator = MeasureList.GetFieldSeparator(fs);
            return string.Format(CultureInfo.InvariantCulture, "[{0},{1}]{3}{2}",
                A,
                B,
                Worker == null ? "" : Worker.Tag,
                fieldSeparator);
        }
    } // IntervalMeasure
}
