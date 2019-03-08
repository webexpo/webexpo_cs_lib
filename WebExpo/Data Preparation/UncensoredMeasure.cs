using System;
using System.Globalization;
using Zygotine.Util;
namespace Zygotine.WebExpo
{
    internal class UncensoredMeasure : Measure
    {

        internal UncensoredMeasure(double a, string workerTag = null)
        {
            if (!Tools.IsFinite(a)) 
                throw new WEException("Parameter a must be a finite number."); 
            A = a;
            SetWorker(workerTag);
        }
        internal UncensoredMeasure(double a, Worker worker)
        {
            if (!Tools.IsFinite(a))
                throw new WEException("Parameter a must be a finite number.");
            A = a;
            this.Worker = worker;
        }

        internal override Measure Duplicate()
        {

            UncensoredMeasure m = new UncensoredMeasure(A, Worker);
            m.Ordinal = this.Ordinal;
            return m;
        }

         
        public override string ToString()
        {
            return ToString(MeasureList.FieldSeparator.Semicolon);
        }
        public virtual string ToString(MeasureList.FieldSeparator fs)
        {
            string fieldSeparator = MeasureList.GetFieldSeparator(fs);
            return string.Format(CultureInfo.InvariantCulture, "{0}{2}{1}",
                A,
                Worker == null ? "" : Worker.Tag,
                fieldSeparator);
        }
    }
}
