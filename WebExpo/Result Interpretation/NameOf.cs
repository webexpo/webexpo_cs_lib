namespace Zygotine
{
    public static class NameOf
    {
        public static string GetCallerName([System.Runtime.CompilerServices.CallerMemberName] string caller = null)
        {
            return caller;
        }
    }
}
