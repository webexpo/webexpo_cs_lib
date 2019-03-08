namespace Zygotine.WebExpo
{
    public static class GlobalParameters
    {
        /*
         * Cette méthode permet d'ajuster 
         * le paramètre MaxNIter de la classe SecuredNRSearch.
         * A des fins de test
         */
        public static int SetNewtonRaphsonMaxIter( int maxIter)
        {
            if( maxIter <= 200)
            {
                maxIter = 200;
            }

            SecuredNRSearch.MaxNIter = maxIter;
            return SecuredNRSearch.MaxNIter;
        }

        public static int GetNewtonRaphsonMaxIter(int maxIter)
        {
            return SecuredNRSearch.MaxNIter;
        }
    }
}
