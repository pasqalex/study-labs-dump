using System;
using System.Globalization;
using System.Linq;

namespace LeastSquares
{
    public class ExpressionParser
    {
        private readonly string _expr;
        private int _pos;

        public ExpressionParser(string expression)
        {
            _expr = expression.Replace(" ", "");
        }

        public Func<double, double> Parse()
        {
            _pos = 0;
            var node = ParseExpression();
            if (_pos != _expr.Length) throw new FormatException("bad expression");
            return node;
        }

        private char Current => _pos < _expr.Length ? _expr[_pos] : '\0';

        private Func<double, double> ParseExpression()
        {
            var left = ParseTerm();
            while (Current == '+' || Current == '-')
            {
                char op = Current; _pos++;
                var right = ParseTerm();
                var l = left;
                left = op == '+' ? (x => l(x) + right(x)) : (x => l(x) - right(x));
            }
            return left;
        }

        private Func<double, double> ParseTerm()
        {
            var left = ParseUnary();
            while (Current == '*' || Current == '/')
            {
                char op = Current; _pos++;
                var right = ParseUnary();
                var l = left;
                left = op == '*' ? (x => l(x) * right(x)) : (x => l(x) / right(x));
            }
            return left;
        }

        private Func<double, double> ParseUnary()
        {
            if (Current == '-') { _pos++; var o = ParseUnary(); return x => -o(x); }
            if (Current == '+') { _pos++; return ParseUnary(); }
            return ParsePower();
        }

        private Func<double, double> ParsePower()
        {
            var left = ParseAtom();
            if (Current == '^')
            {
                _pos++;
                var right = ParseUnary();
                var l = left;
                return x => Math.Pow(l(x), right(x));
            }
            return left;
        }

        private Func<double, double> ParseAtom()
        {
            if (Current == '(')
            {
                _pos++;
                var node = ParseExpression();
                if (Current != ')') throw new FormatException("no )");
                _pos++;
                return node;
            }

            if (char.IsDigit(Current) || Current == '.')
            {
                int start = _pos;
                while (_pos < _expr.Length && (char.IsDigit(_expr[_pos]) || _expr[_pos] == '.')) _pos++;
                double val = double.Parse(_expr.Substring(start, _pos - start), CultureInfo.InvariantCulture);
                return x => val;
            }

            if (char.IsLetter(Current))
            {
                int start = _pos;
                while (_pos < _expr.Length && char.IsLetter(_expr[_pos])) _pos++;
                string name = _expr.Substring(start, _pos - start);

                if (name == "x") return x => x;
                if (name == "e") return x => Math.E;
                if (name == "pi") return x => Math.PI;

                if (Current == '(')
                {
                    _pos++;
                    var arg1 = ParseExpression();

                    if (name.Equals("pow", StringComparison.OrdinalIgnoreCase))
                    {
                        if (Current != ',') throw new FormatException("pow needs 2 args");
                        _pos++;
                        var arg2 = ParseExpression();
                        if (Current != ')') throw new FormatException("no )");
                        _pos++;
                        return x => Math.Pow(arg1(x), arg2(x));
                    }

                    if (Current != ')') throw new FormatException("no )");
                    _pos++;

                    Func<double, double> a1 = arg1;
                    switch (name.ToLowerInvariant())
                    {
                        case "sin": return x => Math.Sin(a1(x));
                        case "cos": return x => Math.Cos(a1(x));
                        case "tan": return x => Math.Tan(a1(x));
                        case "exp": return x => Math.Exp(a1(x));
                        case "ln": return x => Math.Log(a1(x));
                        case "log": return x => Math.Log10(a1(x));
                        case "sqrt": return x => Math.Sqrt(a1(x));
                        case "abs": return x => Math.Abs(a1(x));
                        default: throw new FormatException("unknown func");
                    }
                }
                throw new FormatException("unknown id");
            }
            throw new FormatException("bad char");
        }
    }

    public enum BoundaryType { Closed, Open }

    public static class LinearSolver
    {
        public static double[] Solve(double[,] A, double[] b)
        {
            int n = b.Length;
            var M = new double[n, n + 1];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++) M[i, j] = A[i, j];
                M[i, n] = b[i];
            }
            for (int c = 0; c < n; c++)
            {
                int p = c;
                for (int r = c + 1; r < n; r++)
                    if (Math.Abs(M[r, c]) > Math.Abs(M[p, c])) p = r;
                if (Math.Abs(M[p, c]) < 1e-300) throw new InvalidOperationException("singular");
                if (p != c)
                    for (int j = 0; j <= n; j++) { double t = M[c, j]; M[c, j] = M[p, j]; M[p, j] = t; }
                for (int r = c + 1; r < n; r++)
                {
                    double f = M[r, c] / M[c, c];
                    for (int j = c; j <= n; j++) M[r, j] -= f * M[c, j];
                }
            }
            var x = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                double s = M[i, n];
                for (int j = i + 1; j < n; j++) s -= M[i, j] * x[j];
                x[i] = s / M[i, i];
            }
            return x;
        }
    }

    public class LeastSquaresPolynomial
    {
        private readonly double _a, _b;
        public double[] Coefficients { get; private set; }
        public double Residual { get; private set; } // R = sum (y_i - P(x_i))^2

        public LeastSquaresPolynomial(double[] xs, double[] ys, int m, double a, double b)
        {
            _a = a; _b = b;
            int n = xs.Length;
            var ts = xs.Select(ToT).ToArray();

            var c = new double[2 * m + 1];
            var bv = new double[m + 1];
            for (int s = 0; s < n; s++)
            {
                double p = 1.0;
                for (int k = 0; k <= 2 * m; k++)
                {
                    c[k] += p;
                    if (k <= m) bv[k] += ys[s] * p;
                    p *= ts[s];
                }
            }

            var A = new double[m + 1, m + 1];
            for (int i = 0; i <= m; i++)
                for (int j = 0; j <= m; j++)
                    A[i, j] = c[i + j];

            Coefficients = LinearSolver.Solve(A, bv);

            double R = 0;
            for (int s = 0; s < n; s++)
            {
                double d = ys[s] - Evaluate(xs[s]);
                R += d * d;
            }
            Residual = R;
        }

        public double ToT(double x) => (2 * x - _a - _b) / (_b - _a);

        public double Evaluate(double x)
        {
            double t = ToT(x);
            double s = 0;
            for (int k = Coefficients.Length - 1; k >= 0; k--) s = s * t + Coefficients[k];
            return s;
        }
    }

    internal static class Program
    {
        private const int MinControlPoints = 20;

        private static void Main()
        {
            Console.OutputEncoding = System.Text.Encoding.UTF8;

            Console.Write("f(x) = ");
            string exprStr = Console.ReadLine() ?? "";

            Func<double, double> f;
            try { f = new ExpressionParser(exprStr).Parse(); }
            catch { Console.WriteLine("parse error"); return; }

            double a = ReadDouble("a = ", Math.PI);
            BoundaryType aType = ReadBoundaryType("left (1=closed, 2=open) = ");
            double b = ReadDouble("b = ", 1.5 * Math.PI);
            BoundaryType bType = ReadBoundaryType("right (1=closed, 2=open) = ");

            if (b <= a) { Console.WriteLine("b <= a"); return; }

            double eps = 0;
            if (aType == BoundaryType.Open || bType == BoundaryType.Open)
            {
                eps = ReadDouble("eps = ", 1e-2);
                if (eps <= 0) eps = 1e-2;
            }

            double aCalc = aType == BoundaryType.Open ? a + eps : a;
            double bCalc = bType == BoundaryType.Open ? b - eps : b;
            if (bCalc <= aCalc) { Console.WriteLine("empty interval"); return; }

            int m = ReadInt("m = ", 4);
            if (m < 0) m = 0;
            int n = ReadInt("n = ", 40);
            if (n < m) { n = m; }
            int K = ReadInt("control points = ", 25);
            if (K < MinControlPoints) K = MinControlPoints;

            var xs = new double[n + 1];
            var ys = new double[n + 1];
            for (int i = 0; i <= n; i++)
            {
                xs[i] = aCalc + (bCalc - aCalc) * i / n;
                ys[i] = f(xs[i]);
                if (double.IsNaN(ys[i]) || double.IsInfinity(ys[i]))
                {
                    Console.WriteLine("f undefined");
                    return;
                }
            }

            LeastSquaresPolynomial poly;
            try { poly = new LeastSquaresPolynomial(xs, ys, m, aCalc, bCalc); }
            catch { Console.WriteLine("solve error"); return; }

            Console.WriteLine();
            Console.WriteLine($"interval [{aCalc:F6}, {bCalc:F6}]");
            Console.WriteLine($"t = (2x - {aCalc + bCalc:G10}) / {bCalc - aCalc:G10}");
            for (int k = 0; k <= m; k++)
                Console.WriteLine($"a_{k,-2} = {poly.Coefficients[k],15:E6}");
            Console.WriteLine($"R = {poly.Residual:E6}");

            Console.WriteLine();
            Console.WriteLine($"{"j",3}  {"x",12}  {"f(x)",15}  {"P(x)",15}  {"|f-P|",12}");
            double maxErr = 0, sumSq = 0;
            int cnt = 0;
            for (int j = 1; j <= K; j++)
            {
                double x = aCalc + (bCalc - aCalc) * j / (K + 1);
                double fx = f(x);
                double px = poly.Evaluate(x);
                double err = Math.Abs(fx - px);
                if (!double.IsNaN(err) && !double.IsInfinity(err))
                {
                    maxErr = Math.Max(maxErr, err);
                    sumSq += err * err;
                    cnt++;
                }
                Console.WriteLine($"{j,3}  {x,12:F6}  {fx,15:E6}  {px,15:E6}  {err,12:E4}");
            }

            Console.WriteLine();
            Console.WriteLine($"max = {maxErr:E6}");
            Console.WriteLine($"rms = {(cnt > 0 ? Math.Sqrt(sumSq / cnt) : double.NaN):E6}");
        }

        private static double ReadDouble(string prompt, double def)
        {
            Console.Write(prompt);
            string s = Console.ReadLine();
            if (string.IsNullOrWhiteSpace(s)) return def;
            s = s.Trim();
            try { return new ExpressionParser(s).Parse()(0); }
            catch { return double.Parse(s, CultureInfo.InvariantCulture); }
        }

        private static int ReadInt(string prompt, int def)
        {
            Console.Write(prompt);
            string s = Console.ReadLine();
            if (string.IsNullOrWhiteSpace(s)) return def;
            return int.Parse(s.Trim(), CultureInfo.InvariantCulture);
        }

        private static BoundaryType ReadBoundaryType(string prompt)
        {
            Console.Write(prompt);
            string s = (Console.ReadLine() ?? "1").Trim();
            return s == "2" ? BoundaryType.Open : BoundaryType.Closed;
        }
    }
}
