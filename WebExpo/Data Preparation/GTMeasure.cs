namespace Zygotine.WebExpo
{
    internal class GTMeasure : UncensoredMeasure
    {
        internal GTMeasure(double a, string workerTag = null) : base(a, workerTag)
        {
            this.Censoring = MeasureType.GreaterThan;
        }

        internal GTMeasure(double a, Worker worker) : base(a, worker)
        {
            this.Censoring = MeasureType.GreaterThan;
        }

        internal override Measure Duplicate()
        {

            return new GTMeasure(A, Worker);
        }

         public override string ToString()
        {
            return ToString(MeasureList.FieldSeparator.Semicolon);
        }

        public override string ToString(MeasureList.FieldSeparator fs)
        {
            return ">" + base.ToString(fs);
        }
    } //class GTMeasure


}
