using System;
using System.Globalization;
using System.Linq;

namespace LocalSplines
{
    public class ExpressionParser
    {
        private readonly string _expr;
        private int _pos;

        public ExpressionParser(string expression) { _expr = expression.Replace(" ", ""); _pos = 0; }

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

    public class BarycentricInterpolant
    {
        private readonly double[] _xs;
        private readonly double[] _ys;
        private readonly double[] _w;

        public BarycentricInterpolant(double[] xs, double[] ys)
        {
            _xs = xs; _ys = ys;
            int n = xs.Length;
            _w = new double[n];
            for (int j = 0; j < n; j++)
            {
                double prod = 1.0;
                for (int k = 0; k < n; k++) if (k != j) prod *= (xs[j] - xs[k]);
                _w[j] = 1.0 / prod;
            }
        }

        public double Evaluate(double x)
        {
            int n = _xs.Length;
            for (int j = 0; j < n; j++) if (Math.Abs(x - _xs[j]) < 1e-14) return _ys[j];
            double num = 0, den = 0;
            for (int j = 0; j < n; j++)
            {
                double t = _w[j] / (x - _xs[j]);
                num += t * _ys[j];
                den += t;
            }
            return num / den;
        }

        public double WeightSpread()
        {
            double maxW = _w.Max(v => Math.Abs(v));
            double minW = _w.Min(v => Math.Abs(v));
            return minW > 0 ? maxW / minW : double.PositiveInfinity;
        }
    }

    public class LocalSpline
    {
        private readonly Func<double, double> _f;
        private readonly double[] _meshNodes;
        private readonly int[] _orders;
        private readonly BarycentricInterpolant[] _pieces;

        public double[] MeshNodes => _meshNodes;
        public int[] Orders => _orders;

        public LocalSpline(Func<double, double> f, double[] meshNodes, int[] orders)
        {
            _f = f; _meshNodes = meshNodes; _orders = orders;
            int nSub = meshNodes.Length - 1;
            _pieces = new BarycentricInterpolant[nSub];
            for (int i = 0; i < nSub; i++)
            {
                int m = orders[i];
                var xs = ChebyshevLobattoNodes(m, meshNodes[i], meshNodes[i + 1]);
                var ys = xs.Select(f).ToArray();
                _pieces[i] = new BarycentricInterpolant(xs, ys);
            }
        }

        public static double[] ChebyshevLobattoNodes(int m, double lo, double hi)
        {
            var xi = new double[m + 1];
            for (int k = 0; k <= m; k++) xi[k] = Math.Cos(k * Math.PI / m);
            Array.Sort(xi);
            var t = new double[m + 1];
            for (int k = 0; k <= m; k++) t[k] = 0.5 * (hi - lo) * xi[k] + 0.5 * (hi + lo);
            return t;
        }

        private int FindSubinterval(double x)
        {
            int nSub = _meshNodes.Length - 1;
            for (int i = 0; i < nSub; i++)
                if (x >= _meshNodes[i] - 1e-12 && x <= _meshNodes[i + 1] + 1e-12) return i;
            return x < _meshNodes[0] ? 0 : nSub - 1;
        }

        public double Evaluate(double x) => _pieces[FindSubinterval(x)].Evaluate(x);

        public double[] WeightSpreads() => _pieces.Select(p => p.WeightSpread()).ToArray();

        public (double maxErr, double rmsErr)[] EstimateErrorPerSubinterval(int controlPointsPerSub)
        {
            int nSub = _meshNodes.Length - 1;
            var result = new (double, double)[nSub];
            for (int i = 0; i < nSub; i++)
            {
                double lo = _meshNodes[i], hi = _meshNodes[i + 1];
                double maxE = 0, sumSq = 0;
                int cnt = 0;
                for (int k = 1; k <= controlPointsPerSub; k++)
                {
                    double x = lo + (hi - lo) * k / (controlPointsPerSub + 1);
                    double exact = _f(x);
                    if (double.IsNaN(exact) || double.IsInfinity(exact)) continue;
                    double err = Math.Abs(exact - _pieces[i].Evaluate(x));
                    maxE = Math.Max(maxE, err);
                    sumSq += err * err;
                    cnt++;
                }
                result[i] = cnt > 0 ? (maxE, Math.Sqrt(sumSq / cnt)) : (double.NaN, double.NaN);
            }
            return result;
        }
    }

    internal static class Program
    {
        private const int N = 8;
        private const int ErrorTableProbes = 5;

        private static void Main()
        {
            Console.OutputEncoding = System.Text.Encoding.UTF8;

            Console.Write("f(x) = ");
            string exprStr = Console.ReadLine() ?? "1/cos(x)";
            Func<double, double> f;
            try { f = new ExpressionParser(exprStr).Parse(); }
            catch { Console.WriteLine("parse error"); return; }

            double a = ReadDouble("a = ", Math.PI / 2);
            BoundaryType aType = ReadBoundaryType("left type (1=closed, 2=open) = ");
            double b = ReadDouble("b = ", Math.PI);
            BoundaryType bType = ReadBoundaryType("right type (1=closed, 2=open) = ");
            if (b <= a) { Console.WriteLine("b <= a"); return; }

            double eps = ReadDouble("eps = ", 1e-4);
            if (eps <= 0) eps = 1e-4;

            double aCalc = aType == BoundaryType.Open ? a + eps : a;
            double bCalc = bType == BoundaryType.Open ? b - eps : b;
            if (bCalc <= aCalc) { Console.WriteLine("empty interval"); return; }

            int mBase = ReadInt("m = ", 4);
            if (mBase < 1) mBase = 1;

            int nSub = N + 1;
            double[] meshNodes = new double[nSub + 1];
            for (int i = 0; i <= nSub; i++)
                meshNodes[i] = aCalc + (bCalc - aCalc) * i / nSub;

            int[] orders = Enumerable.Repeat(mBase, nSub).ToArray();

            if (aType == BoundaryType.Open || bType == BoundaryType.Open)
            {
                int mSpec = ReadInt("m_spec = ", mBase * 2);
                if (mSpec < 1) mSpec = 1;
                if (aType == BoundaryType.Open) orders[0] = mSpec;
                if (bType == BoundaryType.Open) orders[nSub - 1] = mSpec;
            }

            int controlPoints = ReadInt("control points = ", 10);
            if (controlPoints < 1) controlPoints = 1;

            var spline = new LocalSpline(f, meshNodes, orders);
            var errs = spline.EstimateErrorPerSubinterval(controlPoints);
            var spreads = spline.WeightSpreads();

            Console.WriteLine();
            Console.WriteLine($"{"i",3} {"range",22} {"m",3} {"max",12} {"rms",12} {"w/w",8}");
            for (int i = 0; i < nSub; i++)
            {
                string range = $"[{meshNodes[i]:F4},{meshNodes[i + 1]:F4}]";
                string maxS = double.IsNaN(errs[i].maxErr) ? "NaN" : errs[i].maxErr.ToString("E3");
                string rmsS = double.IsNaN(errs[i].rmsErr) ? "NaN" : errs[i].rmsErr.ToString("E3");
                Console.WriteLine($"{i,3} {range,22} {orders[i],3} {maxS,12} {rmsS,12} {spreads[i],8:F2}");
            }

            double globalMax = errs.Where(e => !double.IsNaN(e.maxErr)).Select(e => e.maxErr).DefaultIfEmpty(double.NaN).Max();
            Console.WriteLine();
            Console.WriteLine($"spline max error = {globalMax:E3}");

            int globalM = orders.Max();
            var gxs = LocalSpline.ChebyshevLobattoNodes(globalM, aCalc, bCalc);
            var gys = gxs.Select(f).ToArray();
            var globalPoly = new BarycentricInterpolant(gxs, gys);

            double globalPolyMax = 0;
            int totalProbes = nSub * ErrorTableProbes;
            for (int k = 0; k <= totalProbes; k++)
            {
                double x = aCalc + (bCalc - aCalc) * k / totalProbes;
                double exact = f(x);
                if (double.IsNaN(exact) || double.IsInfinity(exact)) continue;
                double approx = globalPoly.Evaluate(x);
                if (double.IsNaN(approx) || double.IsInfinity(approx)) continue;
                globalPolyMax = Math.Max(globalPolyMax, Math.Abs(exact - approx));
            }
            Console.WriteLine($"global poly m = {globalM}, max error = {globalPolyMax:E3}");

            Console.WriteLine();
            Console.WriteLine($"{"i",3} {"x",10} {"f",16} {"S",16} {"|f-S|",10} {"|f-P|",10}");
            for (int i = 0; i < nSub; i++)
            {
                double lo = meshNodes[i], hi = meshNodes[i + 1];
                for (int k = 1; k <= ErrorTableProbes; k++)
                {
                    double x = lo + (hi - lo) * k / (ErrorTableProbes + 1);
                    double fx = f(x);
                    double sx = spline.Evaluate(x);
                    double gx = globalPoly.Evaluate(x);
                    if (double.IsNaN(fx) || double.IsInfinity(fx)) continue;
                    Console.WriteLine($"{i,3} {x,10:F5} {fx,16:E6} {sx,16:E6} {Math.Abs(fx - sx),10:E2} {Math.Abs(fx - gx),10:E2}");
                }
            }
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