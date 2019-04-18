using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using System.Threading;
using System.Diagnostics;
using System.IO;



namespace Pollard
{
    class Program
    {
        static void ExtendedGCD(BigInteger a, BigInteger b, ref BigInteger d, ref BigInteger x, ref BigInteger y)
        {
            if (b != 0)
            {
                //return a , 1, 0
                ExtendedGCD(b, a % b, ref d, ref x, ref y);
                BigInteger c = x;
                x = y;
                y = c - (a / b) * y;
            }
            else
            {
                d = a; x = 1; y = 0;
            }
        }

        static BigInteger EulerPhi(BigInteger input)
        {
            BigInteger res = 1;
            for (int i = 2; i * i <= input; i++)
            {
                BigInteger p = 1;
                while (input % i == 0)
                {
                    input = input / i;      // если не взаимно просты, делим
                    p = p * i;			// произведение делителей i втч и кратных
                }
                p = p / i;
                if (p != 0)
                    // если мы хоть раз делили на текущее i
                    // то общее произведение делителей
                    // умножаем на (i - 1)*i^(число раз - 1)
                    res = res * (p * (i - 1));
            }

            BigInteger n = input - 1;
            if (n == 0)
                return res;
            else
                // умножаем на (input - 1)*input^(число раз - 1)
                // но число раз = 1
                return n * res;
        }
        static void new_state(ref List<Tuple<BigInteger, BigInteger, BigInteger>> xes,
            int index, int w, BigInteger a, BigInteger g, BigInteger p, BigInteger n)
        {
            for (int i = 0; i < w; i++)
                for (int j = 0; j < index; j++)
                {
                    if (xes.Last().Item1 < p / 3)
                        xes.Add(Tuple.Create(a * xes.Last().Item1 % p, xes.Last().Item2, (xes.Last().Item3 + 1) % n));

                    if ((xes.Last().Item1 >= p / 3) && (xes.Last().Item1 < 2 * p / 3))
                        xes.Add(Tuple.Create(xes.Last().Item1 * xes.Last().Item1 % p, 2 * xes.Last().Item2 % n, 2 * xes.Last().Item3 % n));

                    if (xes.Last().Item1 >= 2 * p / 3)
                        xes.Add(Tuple.Create(g * xes.Last().Item1 % p, (xes.Last().Item2 + 1) % n, xes.Last().Item3));
                }

        }
        static BigInteger Pollard_parallel(BigInteger g, BigInteger a, BigInteger p)
        {
            List<Tuple<BigInteger, BigInteger, BigInteger>> xes1 = new List<Tuple<BigInteger, BigInteger, BigInteger>>();
            List<Tuple<BigInteger, BigInteger, BigInteger>> xes2 = new List<Tuple<BigInteger, BigInteger, BigInteger>>();
            BigInteger x0 = 1, a0 = 0, b0 = 0;
            xes1.Add(Tuple.Create(x0, a0, b0));
            xes2.Add(Tuple.Create(x0, a0, b0));
            xes2.Add(Tuple.Create(x0, a0, b0));
            const int N = 2, W = 4;
            var tasks = new List<Thread>();
            //var tasks = new Task[W];
            BigInteger n = EulerPhi(p);
            bool stop = false;
            int global_step = 1, nom = 0;

            if (a == g)
                return 1;
            do
            {
                tasks.Add(new Thread(delegate () { new_state(ref xes2, 2, W, a, g, p, n); }));
                tasks.Last().Start();
                tasks.Add(new Thread(delegate () { new_state(ref xes1, 1, W, a, g, p, n); }));
                tasks.Last().Start();
                foreach (Thread task in tasks) 
                {
                    task.Join();
                }
                // tasks.Clear();
                Parallel.For(0, W, i =>
                {
                    if (xes1[global_step + i].Item1 == xes2[2 * global_step + 2 * i].Item1) { nom = global_step + i; stop = true; }
                });

                global_step = global_step + W;
            } while (stop == false);

            BigInteger a1 = xes1[nom].Item2, a2 = xes2[2 * nom].Item2, b1 = xes1[nom].Item3, b2 = xes2[2 * nom].Item3;
            BigInteger x1 = xes1[nom].Item1, x2 = xes2[2 * nom].Item1;
            BigInteger u = (a1 - a2) % n;
            BigInteger v = (b2 - b1) % n;

            //if (v % n == 0)
            //    return -1;

            BigInteger d = 0, nu = 0, mu = 0;
            ExtendedGCD(v, n, ref d, ref nu, ref mu);
            //if (d < 0) d = d + p;
            BigInteger x = 1;
            for (int i = 0; i < (d + 1); i++)
            {
                BigInteger w = i;
                x = ((u * nu + w * n) / d) % n;
                if (x < 0) x = x + n;
                if ((BigInteger.ModPow(g, x, p) - a) % p == 0)
                    return x;
            }
            return -1;
        }
        static BigInteger Pollard_linear(BigInteger g, BigInteger a, BigInteger p)
        {
            List<Tuple<BigInteger, BigInteger, BigInteger>> xes1 = new List<Tuple<BigInteger, BigInteger, BigInteger>>();
            List<Tuple<BigInteger, BigInteger, BigInteger>> xes2 = new List<Tuple<BigInteger, BigInteger, BigInteger>>();
            BigInteger x0 = 1, a0 = 0, b0 = 0;
            xes1.Add(Tuple.Create(x0, a0, b0));
            xes2.Add(Tuple.Create(x0, a0, b0));
            xes2.Add(Tuple.Create(x0, a0, b0));
            const int N = 2, W = 4;
            BigInteger n = EulerPhi(p);
            bool stop = false;
            int global_step = 1, nom = 0;

            if (a == g)
                return 1;
            do
            {
                new_state(ref xes1, 2, W, a, g, p, n);
                new_state(ref xes2, 2, W, a, g, p, n);
                new_state(ref xes2, 2, W, a, g, p, n);
                for (int i = 0; i < W; i++)
                {
                    if (xes1[global_step + i].Item1 == xes2[2 * global_step + 2 * i].Item1) { nom = global_step + i; stop = true; }
                }
                global_step = global_step + W;
            } while (stop == false);

            BigInteger a1 = xes1[nom].Item2, a2 = xes2[2 * nom].Item2, b1 = xes1[nom].Item3, b2 = xes2[2 * nom].Item3;
            BigInteger x1 = xes1[nom].Item1, x2 = xes2[2 * nom].Item1;
            BigInteger u = (a1 - a2) % n;
            BigInteger v = (b2 - b1) % n;

            //if (v % n == 0)
             //   return -1;

            BigInteger d = 0, nu = 0, mu = 0;
            ExtendedGCD(v, n, ref d, ref nu, ref mu);
            //if (d < 0) d = d + p;
            BigInteger x = 1;
            for (int i = 0; i < (d + 1); i++)
            {
                BigInteger w = i;
                x = ((u * nu + w * n) / d) % n;
                if (x < 0) x = x + n;
                if ((BigInteger.ModPow(g, x, p) - a) % p == 0)
                    return x;
            }
            return -1;
        }

        static BigInteger premutive_log(BigInteger g, BigInteger a, BigInteger b)
        {
            for (BigInteger x = 0; x != b; x++)
                if ((BigInteger.ModPow(g, x, b) - a) % b == 0)
                    return x;
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

            string q = file.ReadLine();
            while (q != null)
            {
                string[] temp = q.Split(' ');
                BigInteger a = BigInteger.Parse(temp[0]);
                BigInteger b = BigInteger.Parse(temp[1]);
                BigInteger p = BigInteger.Parse(temp[2]);
                BigInteger x = 0;
                if (model == "linear")
                {
                    date1 = DateTime.Now;
                    x = Pollard_linear(a, b, p);
                    date2 = DateTime.Now;
                }
                else
                {
                    date1 = DateTime.Now;
                    x = Pollard_parallel(a, b, p);
                    date2 = DateTime.Now;
                }
                // Console.WriteLine("p =" + p + " time = " + (date2 - date1).TotalMilliseconds + " x = " + x);

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
            Console.WriteLine("model = " + model + "; time with log: " + time_with_log / count_with_log);
            Console.WriteLine("model = " + model + "; time without log: " + time_without_log / count_without_log);
           
        }
        static void Main(string[] args)
        {

            StreamReader file = new StreamReader(@"../../../../data/example.txt");
            read_file(file, "linear");
            
            file = new StreamReader(@"../../../../data/example.txt");
            read_file(file, "parallel");
            file = new StreamReader(@"../../../../data/example.txt");
            read_file(file, "linear");

            BigInteger g = 4, a = 5, m = 11;
            Console.WriteLine(Pollard_linear(g, a, m));
            Console.WriteLine(Pollard_parallel(g, a, m));
            Console.WriteLine(premutive_log(g, a, m));
            Console.ReadKey();
        }
    }
}
