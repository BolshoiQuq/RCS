#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <complex>
#include <limits>

const double inf = std::numeric_limits<double>::infinity();
// const double inf=1.7976931348623157E+308;
const double pi=3.14159265358979;
int ws=0;

template <typename T> int sign(T val)
{
    return (T(0) < val) - (val < T(0));
}

template <typename T> T abv(T val)
{
    return val*sign(val);
}

bool iszero(double x)
{
    if (x>-0.0001 && x<0.0001)
        return true;
    else
        return false;
}

struct Point
{
    double x;
    double y;
    double z;
    Point(double X, double Y, double Z)
    {
        this->x = X;
        this->y = Y;
        this->z = Z;
    }
    Point(double X, double Y)
    {
        this->x = X;
        this->y = Y;
        this->z = 0;
    }
    Point(double X)
    {
        this->x = X;
        this->y = 0;
        this->z = 0;
    }
    Point()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    void revolveR()
    {
        double X, Y, Z;
        X=this->x;
        Y=this->y;
        Z=this->z;
        this->x = Y;
        this->y = Z;
        this->z = X;
    }
    void revolveL()
    {
        double X, Y, Z;
        X=this->x;
        Y=this->y;
        Z=this->z;
        this->x = Z;
        this->y = X;
        this->z = Y;
    }
    double length()
    {
        return sqrt(this->x*this->x+this->y*this->y+this->z*this->z);
    }
    void normalize()
    {
        double L=this->length();
        this->x /= L;
        this->y /= L;
        this->z /= L;
        if (iszero(L))
            std::cerr << "[Warning] Zero length in normalize()\n";
    }
};

std::istream& operator>>(std::istream& is, Point& p)
{
    double x, y, z;
    char c1, c2, c3, c4;
    is >> c1 >> x >> c2 >> y >> c3 >> z >> c4;
    if (!is) return is;
    if (c1!='(' || c2!=';' || c3!=';' || c4!=')')
    {
        is.clear(std::ios_base::failbit);
        return is;
    }
    p = Point(x, y, z);
    return is;
}

std::ostream& operator<<(std::ostream& os, const Point& p)
{
    return os << "(" << p.x << "; " << p.y << "; " << p.z << ")";
}

Point operator+(Point a, Point b)
{
    return Point(a.x+b.x, a.y+b.y, a.z+b.z);
}

Point operator-(Point a, Point b)
{
    return Point(a.x-b.x, a.y-b.y, a.z-b.z);
}

Point operator*(double k, Point a)
{
    return Point(k*a.x, k*a.y, k*a.z);
}

Point operator*(Point a, Point b)
{
    return Point(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

Point operator/(Point a, double q)
{
    if (iszero(q))
        std::cerr << "[Warning] Zero q in Point division\n";
    return Point(a.x/q, a.y/q, a.z/q);
}

double operator/(Point a, Point b)
{
    return (a.x*b.x+a.y*b.y+a.z*b.z);
}

double operator%(Point a, Point b)
{
    if (iszero(b.length()))
        std::cerr << "[Warning] Zero b.length() in Point projection\n";
    return (a/b)/b.length();
}

bool operator<(Point a, Point b)
{
    return a.length()<b.length();
}

bool operator>(Point a, Point b)
{
    return a.length()>b.length();
}

bool operator==(Point a, Point b)
{
    return iszero((a-b).length());
}

bool operator<=(Point a, Point b)
{
    return a.length()<=b.length();
}

bool operator>=(Point a, Point b) {
    return a.length()>=b.length();
}

void operator*=(Point& a, double k)
{
    a = k*a;
}

Point spher_to_rect(double r, double theta, double phi)
{
    return Point(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));
}

struct Line
{
    Point A;
    Point B;

    Line(Point a, Point b)
    {
        this->A = a;
        this->B = b;
    }
    Line(Point a)
    {
        this->A = a;
        this->B = Point();
    }
    Line()
    {
        this->A = Point();
        this->B = Point();
    }
    Point direct()
    {
        return (this->B-this->A);
    }
    Point norm_direct()
    {
        Point d=this->direct();
        d.normalize();
        return d;
    }
    bool belongs(Point C)
    {
        Point m=this->direct()*(C-this->A);
        if (iszero(m.length()))
            return true;
        else
            return false;
    }
};

std::istream& operator>>(std::istream& is, Line& l)
{
    Point A, B;
    char c1, c2, c3;
    is >> c1 >> A >> c2 >> B >> c3;
    if (!is) return is;
    if (c1!='{' || c2!=',' || c3!='}')
    {
        is.clear(std::ios_base::failbit);
        return is;
    }
    l = Line(A, B);
    return is;
}

std::ostream& operator<<(std::ostream& os, const Line& l)
{
    return os << "{" << l.A << ", " << l.B << "}";
}

Point intersection(Line f, Line g, char* flag)
{
    if (iszero((f.direct()*g.direct()).length()))
    {
        *flag='p'; // parallel
        return Point(inf, inf, inf);
    }
    double alpha, beta, gamma, mu, nu;
    alpha=f.direct()/f.direct();
    beta=g.direct()/g.direct();
    gamma=f.direct()/g.direct();
    mu=(f.A-g.A)/f.direct();
    nu=(g.A-f.A)/g.direct();
    double s, t;
    if (iszero(alpha*beta-gamma*gamma))
        std::cerr << "[Warning: " << ++ws << "] Zero alpha*beta-gamma*gamma in intersection\n";
    t=-(mu*beta+nu*gamma)/(alpha*beta-gamma*gamma);
    s=-(mu*gamma+nu*alpha)/(alpha*beta-gamma*gamma);
    if (iszero(sqrt((f.A-g.A+t*f.direct()-s*g.direct())/(f.A-g.A+t*f.direct()-s*g.direct()))))
    {
        if (t>=0 && t<=1)
            if (s>=0 && s<=1)
                *flag='s'; // segment-to-segment intersection
            else
                *flag='z'; // segment-to-line intersection
        else
            *flag='i'; // line-to-line intersection
        return f.A+t*f.direct();
    }
    else
    {
        *flag='n'; // no intersection
        return Point(inf, inf, inf);
    }
}

struct Triangle
{
    Point A;
    Point B;
    Point C;
    Triangle(Point a, Point b, Point c)
    {
        this->A = a;
        this->B = b;
        this->C = c;
    }
    Triangle(Point a, Point b)
    {
        this->A = a;
        this->B = b;
        this->C = Point();
    }
    Triangle(Point a)
    {
        this->A = a;
        this->B = Point();
        this->C = Point();
    }
    Triangle()
    {
        this->A = Point();
        this->B = Point();
        this->C = Point();
    }
    void revolveR()
    {
        Point a, b, c;
        a=this->A;
        b=this->B;
        c=this->C;
        this->A = b;
        this->B = c;
        this->C = a;
    }
    void revolveL()
    {
        Point a, b, c;
        a=this->A;
        b=this->B;
        c=this->C;
        this->A = c;
        this->B = a;
        this->C = b;
    }
    Point ori_area()
    {
        return (this->B-this->A)*(this->C-this->A)/2;
    }
    Point normal()
    {
        Point N=this->ori_area();
        N.normalize();
        return N;
    }
    double area()
    {
        Point N=this->ori_area();
        return N.length();
    }
    Point hinter()
    {
        char g;
        return intersection(Line(this->A, this->A+(this->C-this->B)*this->ori_area()), Line(this->B, this->B+(this->A-this->C)*this->ori_area()), &g);
    }
    Point vsplit(Point v)
    {
        char f;
        Point D, u=this->ori_area()*v;
        D=intersection(Line(this->B, this->C), Line(this->A, this->A+u), &f);
        if (f=='s' || f=='z')
            return D; // splitting line is through A
        else
        {
            D=intersection(Line(this->C, this->A), Line(this->B, this->B+u), &f);
            if (f=='s' || f=='z')
            {
                this->revolveR();
                return D; // splitting line is through B, rotate B to A
            }
            else
            {
                D=intersection(Line(this->A, this->B), Line(this->C, this->C+u), &f);
                if (f=='s' || f=='z')
                {
                    this->revolveL();
                    return D; // splitting line is through C, rotate C to A
                }
                else
                {
                    return Point(inf, inf, inf); // no splitting line
                }
            }
        }
    }
    Point hsplit(char* f)
    {
        char g;
        Point H=intersection(Line(this->B, this->C), Line(this->A, this->A+(this->C-this->B)*this->ori_area()), &g);
        if (g=='s' || g=='z')
            *f='i'; // internal height
        else
            if (g=='i')
                *f='e'; // external height
            else
                *f='x'; // error
        return H;
    }
    double rect_calc(Point v, double k, int sgn)
    {
        Point a=this->A-this->C, b=this->B-this->C, n=sgn*this->ori_area();
        if (!iszero(a/b))
            return inf;
        if (!iszero(a/v))
            if (iszero(b/v))
                std::swap(a, b);
            else
                return inf;
        double al=a.length(), bl=b.length(), lambda=2*pi/k;
        double sm=pi*al*al*bl*bl/(lambda*lambda), cosphi=(n/v)/(n.length()*v.length());
        if (cosphi<0)
            return 0;
        else
        {
            double sinphi=sqrt(1-cosphi*cosphi);
            double x=sin(k*bl*sinphi)/(k*bl*sinphi);
            //return sm*cosphi*cosphi*x*x*x*x;
            return sm*cosphi*x*x;
            //return cosphi*al*bl*sin(k*bl*sinphi)/(2*k*bl*sinphi);
        }
    }
    double calc(Point v, double k)
    {
        Point D=this->vsplit(v);
        Triangle ta(this->B, D, this->A), tb(this->C, this->A, D);
        char f1, f2;
        Point H1=ta.hsplit(&f1), H2=tb.hsplit(&f2);
        Triangle t1(this->A, this->B, H1), t2(this->B, D, H1), t3(this->C, this->A, H2), t4(D, this->C, H2);
        int s1=sign(t1.ori_area()/this->ori_area()), s2=sign(t2.ori_area()/this->ori_area()), s3=sign(t3.ori_area()/this->ori_area()), s4=sign(t4.ori_area()/this->ori_area());
        std::cout << s1*t1.rect_calc(v, k, s1)+s2*t2.rect_calc(v, k, s2)+s3*t3.rect_calc(v, k, s3)+s4*t4.rect_calc(v, k, s4) << "\n";
        return abv(s1*t1.rect_calc(v, k, s1)+s2*t2.rect_calc(v, k, s2)+s3*t3.rect_calc(v, k, s3)+s4*t4.rect_calc(v, k, s4));
    }
};

std::istream& operator>>(std::istream& is, Triangle& t)
{
    Point A, B, C;
    char c1, c2, c3, c4;
    is >> c1 >> A >> c2 >> B >> c3 >> C >> c4;
    if (!is) return is;
    if (c1!='[' || c2!=',' || c3!=',' || c4!=']')
    {
        is.clear(std::ios_base::failbit);
        return is;
    }
    t = Triangle(A, B, C);
    return is;
}

std::ostream& operator<<(std::ostream& os, const Triangle& t)
{
    return os << "[" << t.A << ", " << t.B << ", " << t.C << "]";
}
/*
int main()
{
    Triangle t(Point(0,0,7), Point(3,0,0), Point(0,5,0));
    Point v=spher_to_rect(1, pi/6, 0);
    Point D=t.vsplit(v);
    Triangle ta(t.B, D, t.A), tb(t.C, t.A, D);
    char f1, f2;
    Point H1=ta.hsplit(&f1), H2=tb.hsplit(&f2);
    Triangle t1(t.A, t.B, H1), t2(t.B, D, H1), t3(t.C, t.A, H2), t4(D, t.C, H2);
    std::cout << t.calc(v, 20);
    return 0;
}
*/
