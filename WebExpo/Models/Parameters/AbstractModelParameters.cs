namespace Zygotine.WebExpo
{
    public abstract class AbstractModelParameters
    {
        public bool LogNormalDstrn { get; set; }
        protected AbstractModelParameters(bool logNormalDstrn = true)
        {
            this.LogNormalDstrn = logNormalDstrn;
        }
    }
}
