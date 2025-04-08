#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <complex>

#include "3d_modelling3.cpp"

//const double pi=3.14159265358979;

const int freq=9; // frequency of curve segment discretization, 2000, 9
const int arr_size=100000; // maximum point array size
const int greq=18; // frequency of circular discretization, 360, 18
const int std_lng=2000; // standard length of surface
const int graph_size=500; // graph resolution
const double stretch_coeff=1.5; // y/x stretch coefficient, 1
const double srel_pwr=3.0; // superelliptic power, 2

void quad_cur_discr(double X0, double Y0, double X1, double Y1, double X2, double Y2, int freq, double* x, double* y)
{
    for (int i=0; i<=freq; ++i)
    {
        double t=double(i)/double(freq);
        *(x+i)=(1-t)*(1-t)*X0+2*(1-t)*t*X1+t*t*X2;
        *(y+i)=(1-t)*(1-t)*Y0+2*(1-t)*t*Y1+t*t*Y2;
    }
}

void cubic_cur_discr(double X0, double Y0, double X1, double Y1, double X2, double Y2, double X3, double Y3, int freq, double* x, double* y)
{
    for (int i=0; i<=freq; ++i)
    {
        double t=double(i)/double(freq);
        *(x+i)=(1-t)*(1-t)*(1-t)*X0+3*(1-t)*(1-t)*t*X1+3*(1-t)*t*t*X2+t*t*t*X3;
        *(y+i)=(1-t)*(1-t)*(1-t)*Y0+3*(1-t)*(1-t)*t*Y1+3*(1-t)*t*t*Y2+t*t*t*Y3;
    }
}

void curve_offset(double* x, double* y, int n)
{
    double x_min=1000000, y_min=1000000;
    for (int i=0; i<n; ++i)
    {
        if (x[i]<x_min) x_min=x[i];
        if (y[i]<y_min) y_min=y[i];
    }
    for (int i=0; i<n; ++i)
    {
        x[i]-=x_min;
        y[i]-=y_min;
    }
}

void curve_maszstab(double* x, double* y, int n, int sz)
{
    double x_max=-1000000;
    for (int i=0; i<n; ++i)
        if (x[i]>x_max) x_max=x[i];
    for (int i=0; i<n; ++i)
    {
        x[i]=x[i]*sz/x_max;
        y[i]=y[i]*sz/x_max;
    }
}

void curve_reverse(double* x, double* y, int n)
{
    double x_max=0;
    for (int i=0; i<n; ++i)
        if (x[i]>x_max) x_max=x[i];

    for (int i=0; i<n; ++i)
        x[i]=x_max-x[i];

    for (int i=0; i<n/2; ++i)
    {
        std::swap(x[i], x[n-i-1]);
        std::swap(y[i], y[n-i-1]);
    }
}

void curve_invert(double* x, double* y, int n)
{
    double y_max=0;
    for (int i=0; i<n; ++i)
        if (y[i]>y_max) y_max=y[i];

    for (int i=0; i<n; ++i)
        y[i]=y_max-y[i];
}

void volumize(double* z, double* r, int n, int gr, double stcf, double pwr, double* X, double* Y, double* Z)
{
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<=gr; ++j)
        {
            if (iszero(pwr-2))
            {
                X[i*(gr+1)+j]=r[i]*cos(2*pi*j/gr);
                Y[i*(gr+1)+j]=stcf*r[i]*sin(2*pi*j/gr);
                Z[i*(gr+1)+j]=z[i];
            }
            else
            {
                X[i*(gr+1)+j]=r[i]*pow(abv(cos(2*pi*j/gr)), 2/pwr)*sign(cos(2*pi*j/gr));
                Y[i*(gr+1)+j]=stcf*r[i]*pow(abv(sin(2*pi*j/gr)), 2/pwr)*sign(sin(2*pi*j/gr));
                Z[i*(gr+1)+j]=z[i];
            }
        }
    }
}

void objectify(double* X, double* Y, double* Z, int n, int gr, std::string objf)
{
    std::fstream f;
    int k=1;
    f.open(objf, std::ios::out);
    f << std::setprecision(18);
    for (int i=1; i<n; ++i)
    {
        for (int j=0; j<gr; ++j)
        {
            f << "v " << X[(i-1)*(gr+1)+j] << " " << Y[(i-1)*(gr+1)+j] << " " << Z[(i-1)*(gr+1)+j] << "\n";
            f << "v " << X[i*(gr+1)+j] << " " << Y[i*(gr+1)+j] << " " << Z[i*(gr+1)+j] << "\n";
            f << "v " << X[i*(gr+1)+j+1] << " " << Y[i*(gr+1)+j+1] << " " << Z[i*(gr+1)+j+1] << "\n";
            f << "v " << X[(i-1)*(gr+1)+j+1] << " " << Y[(i-1)*(gr+1)+j+1] << " " << Z[(i-1)*(gr+1)+j+1] << "\n";
        }
    }
    for (int i=1; i<n; ++i)
    {
        for (int j=0; j<gr; ++j)
        {
            f << "f " << k << " " << k+1 << " " << k+2 << " " << k+3 << "\n";
            k+=4;
        }
    }
    f.close();
}

double curve_calculate(double* z, double* r, int n, Point v, double k)
{
    double s=0;
    std::cout << "allright so far\n";
    for (int i=1; i<n; ++i)
    {
        for (int j=0; j<greq; ++j)
        {
            Point A(r[i-1]*cos(2*pi*j/greq), r[i-1]*sin(2*pi*j/greq), z[i-1]);
            Point B(r[i]*cos(2*pi*j/greq), r[i]*sin(2*pi*j/greq), z[i]);
            Point C(r[i]*cos(2*pi*(j+1)/greq), r[i]*sin(2*pi*(j+1)/greq), z[i]);
            Point D(r[i-1]*cos(2*pi*(j+1)/greq), r[i-1]*sin(2*pi*(j+1)/greq), z[i-1]);
            //std::cout << A << B << C << D << '\n';
            Triangle t1(A, B, C), t2(C, D, A);
            double sq=t1.calc(v, k)+t2.calc(v, k);
            if (sq==sq)
            {
                s+=abv(sq);
                //std::cout << "sq: " << sq << "\n";
            }
            std::cout << s << "\n";
        }
    }
    return s;
}

void draw_calc(double theta_s, double theta_e, double phi_s, double phi_e)
{
    std::fstream out;
    out.open("calc.data", std::ios::out);
    Triangle t(Point(0, 0, 1), Point(0, 1, 0), Point(1, 0, 0));
    for (double phi=phi_s; phi < phi_e; phi+=pi/180)
    {
        for (double theta=theta_s; theta < theta_e; theta+=pi/180)
        {
            Point v=spher_to_rect(2, theta, phi);
            double s=t.calc(v, 20.0);
            out << s << " ";
        }
        out << "\n";
    }
    out.close();
}

int curve_translate(std::string in, int freq, double* x, double* y)
{
    int i=0, n;
    bool flag=false;
    char c, C, status='m', a[4];
    std::fstream f;

    f.open(in, std::ios::in);

    f >> a[0] >> a[1] >> a[2] >> a[3];
    while (!flag)
    {
        if (a[0]=='p' && a[1]=='a' && a[2]=='t' && a[3]=='h')
            flag=true;
        else
        {
            a[0]=a[1];
            a[1]=a[2];
            a[2]=a[3];
            f >> a[3];
        }
    }

    flag=false;

    while (!flag)
    {
        if (a[0]=='d' && a[1]=='=' && a[2]=='"')
            flag=true;
        else
        {
            a[0]=a[1];
            a[1]=a[2];
            f >> a[2];
        }
    }

    while (c!='"')
    {
        f >> c;
        if (c>='0' && c<='9')
            f.putback(c);
        else
            status=c;

        double dX1, dY1, dX2, dY2, dX3, dY3, X1, Y1, X2, Y2, X3, Y3;
        switch (status)
        {
        case 'M':
        case 'L':
            f >> X1 >> C >> Y1;
            x[i]=X1;
            y[i]=Y1;
            break;
        case 'm':
        case 'l':
            f >> dX1 >> C >> dY1;
            if (i>0)
            {
                x[i]=x[i-1]+dX1;
                y[i]=y[i-1]+dY1;
            }
            else
            {
                x[i]=dX1;
                y[i]=dY1;
            }
            break;
        case 'H':
            f >> X1;
            x[i]=X1;
            y[i]=y[i-1];
            break;
        case 'h':
            f >> dX1;
            x[i]=x[i-1]+dX1;
            y[i]=y[i-1];
            break;
        case 'V':
            f >> Y1;
            x[i]=x[i-1];
            y[i]=Y1;
            break;
        case 'v':
            f >> dY1;
            x[i]=x[i-1];
            y[i]=y[i-1]+dY1;
            break;
        case 'q':
            f >> dX1 >> C >> dY1 >> dX2 >> C >> dY2;
            if (i>0)
            {
                X1=x[i-1]+dX1;
                Y1=y[i-1]+dY1;
                X2=x[i-1]+dX2;
                Y2=y[i-1]+dY2;
            }
            else
            {
                X1=dX1;
                Y1=dY1;
                X2=X1+dX2;
                Y2=X2+dY2;
            }
            quad_cur_discr(x[i-1], y[i-1], X1, Y1, X2, Y2, freq, x+i-1, y+i-1);
            i+=(freq-1);
            break;
        case 'c':
            f >> dX1 >> C >> dY1 >> dX2 >> C >> dY2 >> dX3 >> C >> dY3;
            if (i>0)
            {
                X1=x[i-1]+dX1;
                Y1=y[i-1]+dY1;
                X2=x[i-1]+dX2;
                Y2=y[i-1]+dY2;
                X3=x[i-1]+dX3;
                Y3=y[i-1]+dY3;
            }
            else
            {
                X1=dX1;
                Y1=dY1;
                X2=X1+dX2;
                Y2=X2+dY2;
                X3=X1+dX3;
                Y3=X2+dY3;
            }
            cubic_cur_discr(x[i-1], y[i-1], X1, Y1, X2, Y2, X3, Y3, freq, x+i-1, y+i-1);
            i+=(freq-1);
            break;
        default:
            f >> dX1 >> C >> dY1;
            if (i>0)
            {
                x[i]=x[i-1]+dX1;
                y[i]=y[i-1]+dY1;
            }
            else
            {
                x[i]=dX1;
                y[i]=dY1;
            }
        }
        ++i;
    }
    f.close();

    n=i-1;
    return n;
}

void curve_output(std::string out, int n, double* x, double* y)
{
    std::fstream f;
    f.open(out, std::ios::out);
    for (int i=0; i<n; ++i)
    {
        f << x[i] << " " << y[i] << "\n";
    }
    f.close();
}

void curve_graph(std::string graph, int n, double* x, double* y)
{
    int ms=graph_size/25;
    bool b[graph_size][graph_size];
    std::fstream f;

    for (int i=0; i<graph_size; ++i)
        for (int j=0; j<graph_size; ++j)
            b[i][j]=0;

    for (int i=0; i<n; ++i)
        b[int(ms*y[i]/std_lng)][int(ms*x[i]/std_lng)]=1;

    f.open(graph, std::ios::out);
    for (int i=0; i<graph_size; ++i)
    {
        for (int j=0; j<graph_size; ++j)
            if (b[i][j]==0)
                f << "  ";
            else
                f << "█ ";
        f << "\n";
    }
    f.close();
}

int main()
{
    int n;
    double z[arr_size], r[arr_size];
    std::string in="in2.svg", outa="out_A.txt", graph="graph.txt";
    Point v=spher_to_rect(1, pi/6, 0);
    //Point v=spher_to_rect(1, 0, 0);
    std::cout << "Vector: " << v << "\n";

    n=curve_translate(in, freq, z, r);
    double X[n*(greq+1)], Y[n*(greq+1)], Z[n*(greq+1)];

    curve_invert(z, r, n);
    curve_offset(z, r, n);
    curve_maszstab(z, r, n, std_lng);
    volumize(z, r, n, greq, stretch_coeff, srel_pwr, X, Y, Z);
    //curve_output(outd, n, x, y);
    objectify(X, Y, Z, n, greq, "sed.obj");
    //curve_graph(graph, n, x, y);
    //std::cout << "No problems so far\n";
    std::cout << curve_calculate(z, r, n, v, 20.0);
    //draw_calc(-pi, pi, -pi, pi);
    return 0;
}
