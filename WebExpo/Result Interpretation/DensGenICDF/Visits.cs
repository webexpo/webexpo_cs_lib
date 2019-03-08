
namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Zygotine.Util;

    internal class Visits
    {
        internal class DirectionChange
        {
            public int Index { get; private set; }
            public int Direction { get; private set; } 
            internal DirectionChange(int dir, int index)
            {
                this.Index = index;
                this.Direction = dir;
            }
        }

        
        List<double> list; // la liste des points visités selon l'ordre des visites
        Dictionary<double, List<int>> dict = new Dictionary<double, List<int>>();   // pour chaque x visité, la liste des index de list.
                                                                                    // normalement un point ne devrait être visité qu'une fois ... s'il est fini.
        internal int Count { get { return this.list.Count; } } // Le nombre de visites
        internal int IndexOfLastDuplicateFound { get; private set; } = -1; // L'index du dernier doublon visité 
        internal int IndexOfLastFiniteDuplicateFound { get; private set; } = -1; // L'index du dernier doublon 'fini' visité 
        internal bool caughtInLoop = false;
        internal bool CaughtInLoop { get { return caughtInLoop; } }

        //
        internal DirectionChange[] GetDirectionChanges()
        {
            int lc = list.Count - 1;

            DirectionChange[] dirChanges = new DirectionChange[lc];
            for (int i = 0; i < lc; i++)
            {
                dirChanges[i] = new DirectionChange( Math.Sign(list[i + 1] - list[i]), i);
            }

            return dirChanges;
        }

        internal double[] GetFromTail(double value)
        {
            if (!dict.ContainsKey(value))
            {
                return new double[0];
            }
            else
            {
                int i = dict[value][0]; // l'indice de la première visite faite à value
                return list.Skip(i).ToArray();
            }
        }

        internal Visits(int size)
        {
            this.list = new List<double>(size);
        }

        internal void Add(double x)
        {
            bool test = false;
            if (dict.ContainsKey(x))
            {
                this.IndexOfLastDuplicateFound = dict[x][0];
                if (x.IsFinite())
                {
                    this.IndexOfLastFiniteDuplicateFound = dict[x][0];
                    test = true;
                }
            }
            else
            {
                this.IndexOfLastDuplicateFound = -1;
                this.dict.Add(x, new List<int>());
            }

            this.list.Add(x);
            this.dict[x].Add(Count);
            caughtInLoop = caughtInLoop || test;
        }
    }

}
