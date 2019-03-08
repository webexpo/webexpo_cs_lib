
namespace Zygotine.WebExpo
{
    /// <summary>
    /// The total samples count for each posterior will be : 
    /// <para><code> NBurnin + (NIter * NThin)</code></para>
    /// <para>If MonitorBurnin is true, then for each posterior, NBurnin samples will be collected.</para>
    /// <para>Of the remaining (NIter * NThin) samples, only NIter (i.e. 1 out of NThin) will be retained.</para>
    /// </summary>
    public class McmcParameters
    {
        public int NIter { get; set; } = 15000;
        public bool MonitorBurnin { get; set; } = false;
        public int NBurnin { get; set; } = 500;
        public int NThin { get; set; } = 1;
        internal Messages Error { get; private set; }
        public string ClassName { get; private set; }

        // These members are not technically MCMC parameters but I keep them here for convenince
        public string CsvSeparator { get; set; } = ";";
        public int NumHistoClasses { get; set; } = -1;

        public McmcParameters(int nIter = 15000, int nBurnin = 500, int nThin = 1, bool monitorBurnin = false)
        {
            this.ClassName = this.GetType().Name;
            Error = new Messages();
            this.NIter = nIter;
            this.NBurnin = nBurnin;
            this.NThin = nThin;
            this.MonitorBurnin = monitorBurnin;

            if (this.NIter <= 0)
            {
                this.Error.AddError("Parameter NIter must be > 0.", this.ClassName);
            }

            if (this.NBurnin < 0)
            {
                this.Error.AddWarning(string.Format("Parameter nBurnin (={0}) was set to 0.", this.NBurnin), this.ClassName);
                this.NBurnin = 0;
            }

            if (this.NThin <= 0)
            {
                this.Error.AddWarning(string.Format("Parameter nThin (={0}) was set to 1.", this.NThin), this.ClassName);
                this.NThin = 1;
            }
        }

        public override string ToString()
        {
            return string.Format("sampleSize={0}, burninSize={1}, burninMonitoring={2}, thinning={3}", NIter, NBurnin, MonitorBurnin, NThin);
        }

        public McmcParameters Clone()
        {
            return new McmcParameters(nIter: this.NIter, nBurnin: this.NBurnin, nThin: this.NThin, monitorBurnin: this.MonitorBurnin);
        }
    }
}
