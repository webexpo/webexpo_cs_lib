using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Zygotine.WebExpo
{
    internal class LTMeasure : UncensoredMeasure
    {
        internal LTMeasure(double a, string workerTag = null) : base(a,  workerTag)
        {
            this.Censoring = MeasureType.LessThan;
        }

        internal LTMeasure(double a, Worker worker) : base(a, worker)
        {
            this.Censoring = MeasureType.LessThan;
        }

        internal override Measure Duplicate()
        {
            return new LTMeasure(A, Worker);
        }

       public override string ToString()
        {
            return ToString(MeasureList.FieldSeparator.Semicolon);
        }

        public override string ToString(MeasureList.FieldSeparator fs)
        {
            return "<" + base.ToString(fs);
        }
    }

}
