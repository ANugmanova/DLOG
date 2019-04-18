using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using System.Threading;
using System.Diagnostics;
using System.IO;

namespace ConsoleApplication1
{
    class Program
    {
         //"""Extended GCD:
    //Returns (gcd, x, y) where gcd is the greatest common divisor of a and b
    //with the sign of b if b is nonzero, and with the sign of a if b is 0.
    //The numbers x,y are such that gcd = ax+by."""
        static void xgcd(BigInteger a, BigInteger b, out BigInteger gcd, out BigInteger prevx, out BigInteger prevy)
        {
            prevx = 1; prevy = 0;
            BigInteger x = 0, y = 1, c;
            
            while (b != 0)
            {
                BigInteger q = a / b;
                BigInteger r = a % b;
                c = x; x = prevx - q * x; prevx = c;
                c = y; y = prevy - q * y; prevy = c;
                a = b; b =  r;
            }
            gcd = a;
        }

        //Returns x so pow(a, x, p) is b mod p, or None if no solution. // brute force version
        static BigInteger discreteLogModP(BigInteger a, BigInteger b, BigInteger p)
        {
            BigInteger a_x = 1;
            b = b % p;
            for (BigInteger x = 0; x < p - 1; x++)
            {
                if (a_x == b)
                    return x;
                a_x = (a_x * a) % p;
            }
            return 0;
        }
        
        //returns a list of prime factors of n, knowing min possible >= startFrom.
        static List<BigInteger> factor(BigInteger n, int startFrom = 2) //=((((
        {
            List<BigInteger> factors = new List<BigInteger>();
            if (n <= 1) return factors;
            BigInteger d = startFrom;
            while (n >= d*d)
            {
                if (n % d == 0)
                {
                    factors.Add(d);
                    n = n/d;
                }
                else
                    d += 1 + d % 2;  // 2 -> 3, odd -> odd + 2
            }            
            factors.Add(n);
            return factors;
        }

        //Given a sequence, return a list of (item, consecutive_repetitions)      
        static List<Tuple<BigInteger, BigInteger>> countConsecutiveSame(List<BigInteger> seq)
        {
            BigInteger n = 0, current = -1;
            List<Tuple<BigInteger, BigInteger>> pairs = new List<Tuple<BigInteger, BigInteger>>();
            if (seq == null) return pairs;
            foreach (BigInteger e in seq)
            {
                if (e == current)
                    n += 1;
                else
                {
                    if (n > 0)
                        pairs.Add(Tuple.Create(current, n));
                    n = 1;
                    current = e;
                }
            }
            pairs.Add(Tuple.Create(current, n));            
            return pairs;
        }

        static List<Tuple<BigInteger, BigInteger>> factorMultiplicity(BigInteger n)
        {
            return countConsecutiveSame(factor(n));
        }
    
        //Return the solution to the Chinese Remainder Theorem, (x, M)
        //Pairs contains tuples (a, m) with all m's positive and coprime.
        //Return the smallest nonnegative integer x so 
        //x mod m  = a mod m for each (a, m) in pairs.
        //M is the product ofthe m's.
        static Tuple<BigInteger, BigInteger> ChineseRemainder(List<Tuple<BigInteger, BigInteger>> pairs)
        {
            BigInteger a = pairs[0].Item1, m = pairs[0].Item2;
            pairs.Remove(pairs[0]);
            foreach (Tuple<BigInteger, BigInteger> temp in pairs)
            {
                BigInteger b = temp.Item1, p = temp.Item2;
                BigInteger gcd_x, gcd_y, gcd;
                xgcd(m, p, out gcd, out gcd_x, out gcd_y);
                BigInteger k=((b-a)*gcd_x) % p; //moduli coprime so inverse exists for m mod p
                a=(a+m*k) % (m*p);// joining a mod m and b mod p gives a mod(mp)
                m = m * p; // mod mp
            }
            return Tuple.Create(a, m);
        }

        // return (x, q**r) with (p-1)/q**r = k, 0 <= x < q**r, os beta^(x*k) = alpha^k mod p
        static void getXModP(BigInteger beta, BigInteger alpha, BigInteger p, BigInteger q, BigInteger r,
                                 out BigInteger xFinal, out BigInteger qPow)
        {
            xFinal = 0; // returns x=x0+x1q+x2q^2+...+xiq^i with 0<=xi<q-1
            qPow = 1;
            BigInteger oDiv = (p - 1) / q; // first divided group order
            BigInteger bCurrent = beta;
            BigInteger alphaRaisedModp = BigInteger.ModPow(alpha, oDiv, p);
            BigInteger alphaInv = 0; // просто обратный 
            BigInteger gcd = 0;
            BigInteger prevy = 0;
            xgcd(alpha, p, out gcd, out alphaInv, out prevy);
                
            for (int i = 0; i < r; i++)
            {
                BigInteger betaRaisedModp = BigInteger.ModPow(bCurrent, oDiv, p);
                BigInteger xCurrent = discreteLogModP(alphaRaisedModp, betaRaisedModp, p);
                xFinal = xFinal + xCurrent * qPow;
                // now we calculate the next beta, power of q, order factor
                bCurrent = (bCurrent * BigInteger.ModPow(alphaInv, xCurrent * qPow, p)) % p;           
                qPow = qPow * q;
                oDiv = oDiv / q;
            }
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
                    x = PohligHellmanModP_linear(a, b, p);
                    date2 = DateTime.Now;
                }
                else
                {
                    date1 = DateTime.Now;
                    x = PohligHellmanModP_parallel(a, b, p);
                    date2 = DateTime.Now;
                }

                Console.WriteLine("p =" + p + " time = " + (date2 - date1).TotalMilliseconds + " x = "+ x);
                
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
        static BigInteger PohligHellmanModP_linear(BigInteger beta, BigInteger alpha,BigInteger p, bool verbose = true)
        {
            List<Tuple<BigInteger, BigInteger>> congruenceList = new List<Tuple<BigInteger, BigInteger>>();
            var temp = factorMultiplicity(p-1);
            Parallel.For(0, temp.Count, i =>
            //foreach (Tuple<BigInteger, BigInteger> pair in temp)
            {
                BigInteger xFinal, qPow;
                getXModP(beta, alpha, p, temp[i].Item1, temp[i].Item2, out xFinal, out qPow);
                congruenceList.Add(Tuple.Create(xFinal, qPow));
            });
            Tuple<BigInteger, BigInteger> result = ChineseRemainder(congruenceList);
            BigInteger x = result.Item1;
            BigInteger m = result.Item2;
            if (x < 0) x = x + m;          
            try
            {
                //Console.WriteLine("Given" + beta + "=" + alpha + "^x mod" + p + "\n" + "x=" + x);
                if (BigInteger.ModPow(alpha, x, p) == (beta % p))
                    return x;
            }
            catch
            {
                return -1;
            }               
            return -1;
        }
        static BigInteger PohligHellmanModP_parallel(BigInteger beta, BigInteger alpha, BigInteger p, bool verbose = true)
        {
            List<Tuple<BigInteger, BigInteger>> congruenceList = new List<Tuple<BigInteger, BigInteger>>();
            var temp = factorMultiplicity(p - 1);
            const int N = 2, W = 4;
            //var tasks = new List<Task>();
            var tasks = new List<Thread>();
            foreach (Tuple<BigInteger, BigInteger> pair in temp)
            { 
                tasks.Add( new Thread(new ThreadStart(() => {
                   BigInteger xFinal, qPow;
                   getXModP(beta, alpha, p, pair.Item1, pair.Item2, out xFinal, out qPow);
                   congruenceList.Add(Tuple.Create(xFinal, qPow));
                })));
                tasks.Last().Start();          
            }
            
            foreach( Thread task in tasks) 
            {
                task.Join();
            }
            
            Tuple<BigInteger, BigInteger> result = ChineseRemainder(congruenceList);
            BigInteger x = result.Item1;
            BigInteger m = result.Item2;
            if (x < 0) x = x + m;
            //Console.WriteLine("Given" + beta + "=" + alpha + "^x mod" + p + "\n" + "x=" + x);
            if (BigInteger.ModPow(alpha, x, p) == (beta % p))
                return x;
            return -1;
        }
         
        static void Main(string[] args)
        {
            StreamReader file = new StreamReader(@"../../../../data/example.txt");
            read_file(file, "linear");
            file = new StreamReader(@"../../../../data/example.txt");
            read_file(file, "parallel");

            PohligHellmanModP_linear(2, 4, 7);
            Console.ReadKey();
        }
    }
}
