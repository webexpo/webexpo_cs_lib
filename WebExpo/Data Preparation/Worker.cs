namespace Zygotine.WebExpo
{
    using System;

    public class Worker : IComparable<Worker>
    {
        public string Tag { get; private set; } = null;
        public int Ordinal { get; internal set; } = -1;

        internal Worker(string tag)
        {
            Tag = tag;
        }

        public int CompareTo(Worker other)
        {
            return this.Tag.CompareTo(other.Tag);
        }
    }
}
