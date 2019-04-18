using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.Diagnostics;
using System.Numerics;
using System.IO;
using System.Collections.Concurrent;

namespace ConsoleApplication1
{
    class Program
    {

        // Shenks linear
        public static BigInteger Shenks_linear(BigInteger a, BigInteger b, BigInteger p)
        {
            BigInteger m = (int)Math.Pow(Math.E, BigInteger.Log(p) / 2) + 1; //sqrt            
            Dictionary<BigInteger, BigInteger> vals = new Dictionary<BigInteger, BigInteger>();

            for (BigInteger i = 1; i <= m; i++)
            {
                BigInteger cur_a = BigInteger.ModPow(a, i + m - 1, p);
                if (!vals.ContainsKey(cur_a))
                    vals[cur_a] = i;
            }

            for (BigInteger i = 0, cur = b; i <= m; ++i)
            {
                if (vals.ContainsKey(cur))
                {
                    BigInteger ans = vals[cur] * m - i;
                    if (ans < p)
                        return ans;
                }
                cur = (cur * a) % p;
            }
            return -1;
        }
        static void fill_table(ref ConcurrentDictionary<BigInteger, BigInteger> vals, BigInteger a, BigInteger m, BigInteger p, BigInteger start, BigInteger end)
        {
            for (BigInteger i = start; i <= end; i++)
            {
                BigInteger cur_a = BigInteger.ModPow(a, i + m - 1, p);
                if (!vals.ContainsKey(cur_a))
                    vals[cur_a] = i;
            }
        }
        public static BigInteger Shenks_parallel(BigInteger a, BigInteger b, BigInteger p)
        {
            const int n = 3; //число потоков
            BigInteger m = (int)Math.Pow(Math.E, BigInteger.Log(p) / 2) + 1; //sqrt

            ConcurrentDictionary <BigInteger, BigInteger> vals = new ConcurrentDictionary <BigInteger, BigInteger>();
            
            Thread[] th = new Thread[n];

            Parallel.For(0, n, j => {
                BigInteger start = 1 + j * (m / n);
                BigInteger end = (j + 1) * (m / n);
                fill_table(ref vals, a, m, p, start, end);
            });
            
            for (BigInteger i = 0, cur = b; i <= m; i++)
            {
                if (vals.ContainsKey(cur))
                {
                    BigInteger ans = vals[cur] * m - i;
                    if (ans < p)
                        return ans;
                }
                cur = (cur * a) % p;
            }
            return -1;
        }


        public static void read_file(StreamReader file, string model)
        {
            double time_with_log = 0; // для случая, когда решение существует
            double time_without_log = 0; // для случая, когда решения не существует
            int count_with_log = 0;
            int count_without_log = 0;
            DateTime date1 = new DateTime();
            DateTime date2 = new DateTime();

            string q =  file.ReadLine();
            while (q != null)
            {
                string[] temp = q.Split(' ');
                BigInteger a = BigInteger.Parse(temp[0]);
                BigInteger b = BigInteger.Parse(temp[1]);
                BigInteger p = BigInteger.Parse(temp[2]);
                BigInteger x = 0;
                if (model=="linear")
                {
                    date1 = DateTime.Now;
                    x = Shenks_linear(a, b, p);
                    date2 = DateTime.Now;
                }
                else
                {
                    date1 = DateTime.Now;
                    x = Shenks_parallel(a, b, p);
                    date2 = DateTime.Now;
                }

                if (x == -1)
                {
                    time_without_log = time_without_log + (date2 - date1).TotalMilliseconds;
                    count_without_log++;
                }
                else
                {
                    time_with_log = time_with_log + (date2 - date1).TotalMilliseconds;
                    count_with_log++;
                }
                q = file.ReadLine();
            }
            file.Close();
            Console.WriteLine("model = "+ model + "; time with log: " + time_with_log / count_with_log);
            Console.WriteLine("model = " + model + "; time without log: " + time_without_log / count_without_log);

        }

 
        static void Main(string[] args)
        {
            
            StreamReader file = new StreamReader(@"../../../../data/example.txt");
            read_file(file, "linear");
            file = new StreamReader(@"../../../../data/example.txt");
            read_file(file, "parallel");

            Console.ReadKey();
        }
    }
}

  

