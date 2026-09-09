using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace LeastSquaresApprox
{
    public class Parser //шедевро парсер, который дважды клал компилятор
    {
        string s;
        int p;

        public Parser(string e) { s = e.Replace(" ", ""); }

        public double Eval(double x)
        {
            p = 0;
            var r = Expr(x);
            if (p != s.Length) throw new FormatException($"Ошибка в {p}: {s}");
            return r;
        }

        char C => p < s.Length ? s[p] : '\0';

        double Expr(double x)
        {
            var v = Term(x);
            while (C == '+' || C == '-')
            {
                var op = C; p++;
                var r = Term(x);
                v = op == '+' ? v + r : v - r;
            }
            return v;
        }

        double Term(double x)
        {
            var v = Fact(x);
            while (C == '*' || C == '/')
            {
                var op = C; p++;
                var r = Fact(x);
                v = op == '*' ? v * r : v / r;
            }
            return v;
        }

        double Fact(double x)
        {
            var b = Unary(x);
            if (C == '^') { p++; return Math.Pow(b, Fact(x)); }
            return b;
        }

        double Unary(double x)
        {
            if (C == '+') { p++; return Unary(x); }
            if (C == '-') { p++; return -Unary(x); }
            return Prim(x);
        }

        double Prim(double x)
        {
            if (C == '(') { p++; var v = Expr(x); if (C != ')') throw new FormatException(") expected"); p++; return v; }
            if (char.IsDigit(C) || C == '.') return Num();
            if (char.IsLetter(C)) return Id(x);
            throw new FormatException($"Bad char {C} at {p}");
        }

        double Num()
        {
            var start = p;
            while (p < s.Length && (char.IsDigit(s[p]) || s[p] == '.')) p++;
            return double.Parse(s.Substring(start, p - start), CultureInfo.InvariantCulture);
        }

        double Id(double x)
        {
            var start = p;
            while (p < s.Length && char.IsLetter(s[p])) p++;
            var name = s.Substring(start, p - start).ToLower();

            if (name == "pi") return Math.PI;
            if (name == "e" && C != '(') return Math.E;
            if (name == "x" && C != '(') return x;

            if (C == '(')
            {
                p++;
                var args = new List<double> { Expr(x) };
                while (C == ',') { p++; args.Add(Expr(x)); }
                if (C != ')') throw new FormatException(") expected after args");
                p++;
                return Func(name, args);
            }
            throw new FormatException($"Unknown {name}");
        }

        double Func(string name, List<double> a)
        {
            var v = a[0];
            switch (name)
            {
                case "sin": return Math.Sin(v);
                case "cos": return Math.Cos(v);
                case "tan": return Math.Tan(v);
                case "cot": return 1 / Math.Tan(v);
                case "exp": return Math.Exp(v);
                case "ln": return Math.Log(v);
                case "log": return a.Count > 1 ? Math.Log(v, a[1]) : Math.Log10(v);
                case "sqrt": return Math.Sqrt(v);
                case "abs": return Math.Abs(v);
                case "pow": return Math.Pow(v, a[1]);
                default: throw new FormatException($"Unknown func {name}");
            }
        }
    }

    public enum BType { Closed, Open }

    public class Bound
    {
        public double Val;
        public BType Type;
        public Bound(double v, BType t) { Val = v; Type = t; }
        public override string ToString()
        {
            return Type == BType.Closed ? Val.ToString("0.####") : Val.ToString("0.####") + " (open)";
        }
    }

    public class Approx
    {
        public double[] c;
        public double lo, hi;

        public Approx(double[] c, double lo, double hi) { this.c = c; this.lo = lo; this.hi = hi; }

        double Norm(double x) { return (2 * x - (lo + hi)) / (hi - lo); }

        public double Eval(double x)
        {
            var t = Norm(x);
            double r = 0, pow = 1;
            for (int k = 0; k < c.Length; k++) { r += c[k] * pow; pow *= t; }
            return r;
        }
    }

    public static class Builder
    {
        public static Approx Build(double[] xs, double[] ys, int deg)
        {
            if (xs.Length != ys.Length) throw new ArgumentException("xs != ys length");
            if (xs.Length < deg + 1) throw new ArgumentException($"Need at least {deg + 1} points, got {xs.Length}");

            var lo = xs.Min();
            var hi = xs.Max();
            var t = xs.Select(x => (2 * x - (lo + hi)) / (hi - lo)).ToArray();

            int m = deg + 1, n = xs.Length;
            var A = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                var v = 1.0;
                for (int k = 0; k < m; k++) { A[i, k] = v; v *= t[i]; }
            }

            var ATA = new double[m, m];
            var ATy = new double[m];
            for (int pp = 0; pp < m; pp++)
            {
                for (int q = 0; q < m; q++)
                {
                    double sum = 0;
                    for (int i = 0; i < n; i++) sum += A[i, pp] * A[i, q];
                    ATA[pp, q] = sum;
                }
                double sumY = 0;
                for (int i = 0; i < n; i++) sumY += A[i, pp] * ys[i];
                ATy[pp] = sumY;
            }

            return new Approx(Solve(ATA, ATy), lo, hi);
        }

        static double[] Solve(double[,] M, double[] b)
        {
            int n = b.Length;
            for (int col = 0; col < n; col++)
            {
                int pivot = col;
                double max = Math.Abs(M[col, col]);
                for (int r = col + 1; r < n; r++)
                    if (Math.Abs(M[r, col]) > max) { max = Math.Abs(M[r, col]); pivot = r; }

                if (max < 1e-14) throw new InvalidOperationException("Singular matrix");

                if (pivot != col)
                {
                    for (int c = 0; c < n; c++) { double tmp = M[col, c]; M[col, c] = M[pivot, c]; M[pivot, c] = tmp; }
                    double tmpB = b[col]; b[col] = b[pivot]; b[pivot] = tmpB;
                }

                for (int r = col + 1; r < n; r++)
                {
                    var f = M[r, col] / M[col, col];
                    for (int c = col; c < n; c++) M[r, c] -= f * M[col, c];
                    b[r] -= f * b[col];
                }
            }

            var x = new double[n];
            for (int r = n - 1; r >= 0; r--)
            {
                double sum = b[r];
                for (int c = r + 1; c < n; c++) sum -= M[r, c] * x[c];
                x[r] = sum / M[r, r];
            }
            return x;
        }
    }

    internal static class Program
    {
        static void Main()
        {
            Console.OutputEncoding = System.Text.Encoding.UTF8;
            if (Ask("Ручной ввод? yes nou: ")) Interactive(); else Test();
        }

        static bool Ask(string q)
        {
            while (true)
            {
                Console.Write(q);
                var a = (Console.ReadLine() ?? "").Trim().ToLower();
                if (a == "y" || a == "yes" || a == "да") return true;
                if (a == "n" || a == "no" || a == "нет") return false;
                Console.WriteLine("y/n");
            }
        }

        static void Interactive()
        {
            Console.Write("f(x) = ");
            var f = Console.ReadLine() ?? "x";
            var p = new Parser(f);
            try { p.Eval(1); }
            catch (Exception e) { Console.WriteLine($"Error: {e.Message}"); return; }

            var L = ReadBound("левой");
            var R = ReadBound("правой");
            if (L.Val >= R.Val) { Console.WriteLine("left < right"); return; }

            int n = ReadInt("Число точек: ", 2, 100000);
            int deg = ReadInt("Степень: ", 1, Math.Min(n - 1, 20));
            Run(f, p, L, R, n, deg);
        }

        static Bound ReadBound(string label)
        {
            Console.Write($"{label} граница: ");
            var s = Console.ReadLine() ?? "0";
            double v;
            try { v = new Parser(s).Eval(0); }
            catch { v = double.Parse(s, CultureInfo.InvariantCulture); }
            Console.Write($"Включена? [1=да, 2=нет]: ");
            var t = (Console.ReadLine() ?? "1").Trim();
            return new Bound(v, t == "2" ? BType.Open : BType.Closed);
        }

        static int ReadInt(string q, int min, int max)
        {
            while (true)
            {
                Console.Write(q);
                if (int.TryParse(Console.ReadLine(), out int v) && v >= min && v <= max) return v;
                Console.WriteLine($"От {min} до {max}");
            }
        }

        static void Test()
        {
            Console.WriteLine("\n--- Тест: x/sin(x) на (pi, 3pi/2] ---");
            Console.WriteLine("pi — особая точка (sin=0), исключена.\n");
            var f = "x/sin(x)";
            var p = new Parser(f);
            var L = new Bound(Math.PI, BType.Open);
            var R = new Bound(3 * Math.PI / 2, BType.Closed);
            int n = 15, deg = 4;
            Console.WriteLine($"Точек: {n}, степень: {deg}\n");
            Run(f, p, L, R, n, deg);
        }

        static void Run(string f, Parser p, Bound L, Bound R, int n, int deg)
        {
            double eps = (R.Val - L.Val) * 0.1;
            double lo = L.Type == BType.Open ? L.Val + eps : L.Val;
            double hi = R.Type == BType.Open ? R.Val - eps : R.Val;
            if (lo >= hi) { Console.WriteLine("Пустой отрезок"); return; }

            var xs = new double[n];
            var ys = new double[n];
            int cnt = 0;
            var bad = new List<double>();

            for (int i = 0; i < n; i++)
            {
                double x = n == 1 ? lo : lo + (hi - lo) * i / (n - 1);
                double y;
                try
                {
                    y = p.Eval(x);
                    if (double.IsNaN(y) || double.IsInfinity(y)) { bad.Add(x); continue; }
                }
                catch { bad.Add(x); continue; }
                xs[cnt] = x; ys[cnt] = y; cnt++;
            }

            if (bad.Count > 0)
            {
                Console.WriteLine($"Исключено особых точек: {bad.Count}");
                foreach (var v in bad) Console.WriteLine($"  x={v:0.######}");
                Console.WriteLine();
            }

            if (cnt < deg + 1) { Console.WriteLine($"Точек {cnt}, надо {deg + 1}"); return; }
            Array.Resize(ref xs, cnt);
            Array.Resize(ref ys, cnt);

            Approx res;
            try { res = Builder.Build(xs, ys, deg); }
            catch (Exception e) { Console.WriteLine($"Ошибка: {e.Message}"); return; }

            Console.WriteLine($"f(x) = {f}");
            Console.WriteLine($"Отрезок [{L.Val}; {R.Val}]");
            Console.WriteLine($"Рабочий [{lo:0.######}; {hi:0.######}]");
            Console.WriteLine($"Точек {cnt}, степень: {deg}\n");
            Console.WriteLine("Коэффы P(t)");
            for (int k = 0; k < res.c.Length; k++) Console.WriteLine($"  c{k} = {res.c[k]:0.######}");
            Console.WriteLine();

            Console.WriteLine($"{"x",10} {"f(x)",14} {"P(x)",14} {"|f-P|",12}");
            double maxErr = 0, sqSum = 0;
            for (int i = 0; i < xs.Length; i++)
            {
                var a = res.Eval(xs[i]);
                var e = Math.Abs(ys[i] - a);
                maxErr = Math.Max(maxErr, e);
                sqSum += e * e;
                Console.WriteLine($"{xs[i],10:0.####} {ys[i],14:0.######} {a,14:0.######} {e,12:0.######}");
            }
            Console.WriteLine();
            Console.WriteLine($"Max error: {maxErr:0.######}");
            Console.WriteLine($"RMS:       {Math.Sqrt(sqSum / xs.Length):0.######}");
        }
    }
}
