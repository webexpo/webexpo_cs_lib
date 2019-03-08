
namespace Zygotine.WebExpo
{
    public class Chain
    {
        public int Length { get { return Values.Length; } }
        public double[] Values { get; internal set; }
        public string Label { get; internal set; }
        public bool IsSample { get; internal set; }
        public string Header { get; internal set; }

        internal Chain(string label, string header, int length, bool isSample)
        {
            Values = new double[length];
            Label = label;
            Header = header;
            IsSample = isSample;
        }
    }
}
