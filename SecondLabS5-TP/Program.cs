using System;
using System.Collections.Generic;
using System.Globalization;

namespace FourierOrthogonalPolynomials
{
    class Program
    {
        static double A, B; //отрезок, на котором задана функция пользователя

        static void Main()
        {

            Console.Write("функция ");
            string funcStr = Console.ReadLine();

            Console.Write("левая граница отрезка ");
            A = ReadDouble();
// эту фигню можно снести и пихнуть явно свои значения с желаемым отступом, щас он по дефолту тащит границы отрезка
            Console.Write("правая граница отрезка ");
            B = ReadDouble();

            Console.Write("порядок частной суммы n ");
            int n = ReadInt();

            Func<double, double> f = x => EvalFunction(funcStr, x);

            //ароверим, что функция вычисляется в нескольких точках отрезка
            try
            {
                double mid = (A + B) / 2.0;
                double t = f(mid);
                if (double.IsNaN(t) || double.IsInfinity(t))
                    Console.WriteLine("Проверьте область определения.");
            }
            catch (Exception ex)
            {
                Console.WriteLine("Ошибка при вычислении функции: " + ex.Message);
                return;
            }
            
            double[] aCheb1 = ComputeChebyshev1Coeffs(f, n); //могочлены Чебышева 1го рода
            double err1 = Chebyshev1L2Error(f, aCheb1, n);

            
            double[] aCheb2 = ComputeChebyshev2Coeffs(f, n); //многочлены Чебышева 2го рода
            double err2 = Chebyshev2L2Error(f, aCheb2, n);

            
            double[] aLeg = ComputeLegendreCoeffs(f, n); //многочлены Лежандра
            double err3 = LegendreL2Error(f, aLeg, n);

            PrintResults("Многочлены Чебышева 1 рода", aCheb1, err1);
            PrintResults("Многочлены Чебышева 2 рода", aCheb2, err2);
            PrintResults("Многочлены Лежандра", aLeg, err3);

            Console.WriteLine("\nСравнение погрешностей (L2-норма):");
            Console.WriteLine($"  Чебышев 1 рода : {err1:F10}");
            Console.WriteLine($"  Чебышев 2 рода: {err2:F10}");
            Console.WriteLine($"  Лежандр        : {err3:F10}");

            Console.WriteLine("\nПостроить таблицу значений f(x)?: ");
            string ans = Console.ReadLine();
            if (ans != null && ans.Trim().ToLower() == "y")
            {
                PrintTable(f, aCheb1, aCheb2, aLeg, n);
            }
        }

        static double ReadDouble()
        {
            double val;
            while (!double.TryParse(Console.ReadLine(), System.Globalization.NumberStyles.Any,
                     System.Globalization.CultureInfo.InvariantCulture, out val))
            {
                Console.Write("некорректный ввод ");
            }
            return val;
        }

        static int ReadInt()
        {
            int val;
            while (!int.TryParse(Console.ReadLine(), out val) || val < 0)
            {
                Console.Write("Некорректный ввод (целое неотрицательное число ");
            }
            return val;
        }

        static double EvalFunction(string expr, double x)
        {
            var parser = new ExprParser(expr, x);
            return parser.Parse();
        }

        static double ToX(double t) => (B - A) / 2.0 * t + (A + B) / 2.0;

        static Func<double, double> MapToUnit(Func<double, double> f) => t => f(ToX(t));

        static double LegendreP(int n, double x)
        {
            if (n == 0) return 1.0;
            if (n == 1) return x;
            double Pnm1 = 1.0, Pn = x, Pnp1 = 0.0;
            for (int k = 1; k < n; k++)
            {
                Pnp1 = ((2 * k + 1) * x * Pn - k * Pnm1) / (k + 1);
                Pnm1 = Pn;
                Pn = Pnp1;
            }
            return Pn;
        }

        static double LegendreHat(int n, double x) => Math.Sqrt((2.0 * n + 1) / 2.0) * LegendreP(n, x);

        static double[] ComputeLegendreCoeffs(Func<double, double> fOrig, int n)
        {
            var f = MapToUnit(fOrig);
            double[] coeffs = new double[n + 1];
            for (int k = 0; k <= n; k++)
            {
                int kk = k;
                coeffs[k] = SimpsonComposite(t => f(t) * LegendreHat(kk, t), -1, 1, 2000);
            }
            return coeffs;
        }

        static double LegendreL2Error(Func<double, double> fOrig, double[] coeffs, int n)
        {
            var f = MapToUnit(fOrig);
            Func<double, double> S = t =>
            {
                double s = 0;
                for (int k = 0; k <= n; k++) s += coeffs[k] * LegendreHat(k, t);
                return s;
            };
            double integral = SimpsonComposite(t => Math.Pow(f(t) - S(t), 2), -1, 1, 4000);
            double jac = (B - A) / 2.0;
            return Math.Sqrt(Math.Abs(integral) * jac);
        }

        static double ChebT_hat_theta(int n, double theta)
        {
            if (n == 0) return 1.0 / Math.Sqrt(Math.PI);
            return Math.Sqrt(2.0 / Math.PI) * Math.Cos(n * theta);
        }

        static double[] ComputeChebyshev1Coeffs(Func<double, double> fOrig, int n)
        {
            var f = MapToUnit(fOrig);
            double[] coeffs = new double[n + 1];
            for (int k = 0; k <= n; k++)
            {
                int kk = k;
                coeffs[k] = SimpsonComposite(
                    th => f(Math.Cos(th)) * ChebT_hat_theta(kk, th), 0, Math.PI, 2000);
            }
            return coeffs;
        }

        static double Chebyshev1L2Error(Func<double, double> fOrig, double[] coeffs, int n)
        {
            var f = MapToUnit(fOrig);
            Func<double, double> Stheta = th =>
            {
                double s = 0;
                for (int k = 0; k <= n; k++) s += coeffs[k] * ChebT_hat_theta(k, th);
                return s;
            };
            double integral = SimpsonComposite(th => Math.Pow(f(Math.Cos(th)) - Stheta(th), 2), 0, Math.PI, 4000);
            double jac = (B - A) / 2.0;
            return Math.Sqrt(Math.Abs(integral) * jac);
        }

        static double[] ComputeChebyshev2Coeffs(Func<double, double> fOrig, int n)
        {
            var f = MapToUnit(fOrig);
            double[] coeffs = new double[n + 1];
            for (int k = 0; k <= n; k++)
            {
                int kk = k;
                coeffs[k] = SimpsonComposite(
                    th => f(Math.Cos(th)) * Math.Sqrt(2.0 / Math.PI) * Math.Sin((kk + 1) * th) * Math.Sin(th),
                    0, Math.PI, 2000);
            }
            return coeffs;
        }

        static double ChebU_hat(int n, double x)
        {
            double theta = Math.Acos(x);
            double s = Math.Sin(theta);
            if (Math.Abs(s) < 1e-12)
            {
                double sign = (x > 0) ? 1.0 : -1.0;
                return Math.Sqrt(2.0 / Math.PI) * Math.Pow(sign, n) * (n + 1);
            }
            return Math.Sqrt(2.0 / Math.PI) * Math.Sin((n + 1) * theta) / s;
        }

        static double Chebyshev2L2Error(Func<double, double> fOrig, double[] coeffs, int n)
        {
            var f = MapToUnit(fOrig);
            Func<double, double> Stheta = th =>
            {
                double s = 0;
                for (int k = 0; k <= n; k++) s += coeffs[k] * Math.Sqrt(2.0 / Math.PI) * Math.Sin((k + 1) * th) / Math.Sin(th);
                return s;
            };
            double integral = SimpsonComposite(th =>
            {
                double diff;
                if (Math.Sin(th) < 1e-10)
                {
                    // на концах сингулярности нет (домножаем на sin^2) иберём предел через соседнюю точку
                    diff = Math.Pow(f(Math.Cos(th)) - EvalSTheta(coeffs, th, n), 2) * Math.Pow(Math.Sin(th), 2);
                }
                else
                {
                    diff = Math.Pow(f(Math.Cos(th)) - Stheta(th), 2) * Math.Pow(Math.Sin(th), 2);
                }
                return diff;
            }, 0, Math.PI, 4000);
            double jac = (B - A) / 2.0;
            return Math.Sqrt(Math.Abs(integral) * jac);
        }

        static double EvalSTheta(double[] coeffs, double th, int n)
        {
            double x = Math.Cos(th);
            double s = 0;
            for (int k = 0; k <= n; k++) s += coeffs[k] * ChebU_hat(k, x);
            return s;
        }

        static double SimpsonComposite(Func<double, double> g, double lo, double hi, int subintervals)
        {
            if (subintervals % 2 != 0) subintervals++;
            double h = (hi - lo) / subintervals;
            double sum = g(lo) + g(hi);
            for (int i = 1; i < subintervals; i++)
            {
                double x = lo + i * h;
                sum += (i % 2 == 0 ? 2.0 : 4.0) * g(x);
            }
            return sum * h / 3.0;
        }

        static void PrintResults(string title, double[] coeffs, double l2error)
        {
            Console.WriteLine($"--- {title} ---");
            for (int k = 0; k < coeffs.Length; k++)
                Console.WriteLine($"  a_{k} = {coeffs[k]:F8}");
            Console.WriteLine($"  Погрешность ||f - S_n||_L2 = {l2error:F10}\n");
        }

        static void PrintTable(Func<double, double> fOrig, double[] aCheb1, double[] aCheb2, double[] aLeg, int n)
        {
            Console.WriteLine("\n     x        f(x)        S_Cheb1      S_Cheb2      S_Legendre");
            int steps = 10;
            for (int i = 0; i <= steps; i++)
            {
                double x = A + (B - A) * i / steps;
                double t = 2 * (x - A) / (B - A) - 1; // обратная замена x к t
                if (t < -1) t = -1; if (t > 1) t = 1;

                double fx = fOrig(x);

                double sCheb1 = 0;
                for (int k = 0; k < aCheb1.Length; k++)
                {
                    double th = Math.Acos(Math.Max(-1, Math.Min(1, t)));
                    sCheb1 += aCheb1[k] * ChebT_hat_theta(k, th);
                }

                double sCheb2 = 0;
                for (int k = 0; k < aCheb2.Length; k++)
                    sCheb2 += aCheb2[k] * ChebU_hat(k, t);

                double sLeg = 0;
                for (int k = 0; k < aLeg.Length; k++)
                    sLeg += aLeg[k] * LegendreHat(k, t);

                Console.WriteLine($"{x,8:F4}  {fx,10:F6}  {sCheb1,10:F6}  {sCheb2,10:F6}  {sLeg,10:F6}");
            }
        }
    }

    class ExprParser
    {
        private readonly string _s;
        private int _pos;
        private readonly double _x;

        public ExprParser(string s, double x)
        {
            _s = s.Replace(" ", "").ToLowerInvariant();
            _pos = 0;
            _x = x;
        }

        public double Parse()
        {
            double result = ParseExpr();
            if (_pos != _s.Length)
                throw new Exception($"некорректное выражение около позиции {_pos}: '{_s.Substring(_pos)}'");
            return result;
        }

        private char Peek() => _pos < _s.Length ? _s[_pos] : '\0';

        private double ParseExpr()
        {
            double val = ParseTerm();
            while (Peek() == '+' || Peek() == '-')
            {
                char op = Peek();
                _pos++;
                double rhs = ParseTerm();
                val = (op == '+') ? val + rhs : val - rhs;
            }
            return val;
        }

        private double ParseTerm()
        {
            double val = ParseFactor();
            while (Peek() == '*' || Peek() == '/')
            {
                char op = Peek();
                _pos++;
                double rhs = ParseFactor();
                val = (op == '*') ? val * rhs : val / rhs;
            }
            return val;
        }

        private double ParseFactor()
        {
            double baseVal = ParseUnary();
            if (Peek() == '^')
            {
                _pos++;
                double exponent = ParseFactor(); //право-ассоциативность
                return Math.Pow(baseVal, exponent);
            }
            return baseVal;
        }

        private double ParseUnary()
        {
            if (Peek() == '-')
            {
                _pos++;
                return -ParseUnary();
            }
            if (Peek() == '+')
            {
                _pos++;
                return ParseUnary();
            }
            return ParsePrimary();
        }

        private double ParsePrimary()
        {
            if (Peek() == '(')
            {
                _pos++;
                double val = ParseExpr();
                if (Peek() != ')') throw new Exception("Ожидалась ')'");
                _pos++;
                return val;
            }

            // число
            if (char.IsDigit(Peek()) || Peek() == '.')
            {
                int start = _pos;
                while (_pos < _s.Length && (char.IsDigit(_s[_pos]) || _s[_pos] == '.')) _pos++;
                return double.Parse(_s.Substring(start, _pos - start), CultureInfo.InvariantCulture);
            }

            if (char.IsLetter(Peek()))
            {
                int start = _pos;
                while (_pos < _s.Length && char.IsLetter(_s[_pos])) _pos++;
                string ident = _s.Substring(start, _pos - start);

                if (ident == "x") return _x;
                if (ident == "pi") return Math.PI;
                if (ident == "e" && Peek() != '(') return Math.E;

                // функция: ожидаем '('
                if (Peek() == '(')
                {
                    _pos++;
                    double arg = ParseExpr();
                    if (Peek() != ')') throw new Exception("Ожидалась ')' после аргумента");
                    _pos++;
                    switch (ident)
                    {
                        case "sin": return Math.Sin(arg);
                        case "cos": return Math.Cos(arg);
                        case "tan": return Math.Tan(arg);
                        case "exp": return Math.Exp(arg);
                        case "ln": return Math.Log(arg);
                        case "log10": return Math.Log10(arg);
                        case "sqrt": return Math.Sqrt(arg);
                        case "abs": return Math.Abs(arg);
                        default: throw new Exception($"Неизвестная функция '{ident}'");
                    }
                }
                throw new Exception($"Неизвестный идентификатор '{ident}'");
            }

            throw new Exception($"Неожиданный символ '{Peek()}' на позиции {_pos}");
        }
    }
}