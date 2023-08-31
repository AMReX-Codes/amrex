#include <AMReX.H>
#include <cmath>

// This is intentional. Cannot have std:: in amrex::Parser expressions.
using std::sin;
using std::atan2;
using std::pow;

double f (int icase, double x, double y, double z)
{
    switch (icase)
    {
    case 0:
        return 0 - sin(x);
    case 1:
        return 3 + 4;
    case 2:
        return y + y;
    case 3:
        return 3*sin(z) + sin(z);
    case 4:
        return x + 3*x;
    case 5:
        return 3*(sin(y)+sin(z)) + 4*(sin(y)+sin(z));
    case 6:
        return 3/z + 4/z;
    case 7:
        return x + x*sin(z);
    case 8:
        return sin(x) + y*sin(x);
    case 9:
        return sin(x)*sin(y) + sin(x);
    case 10:
        return sin(x)*sin(y) + sin(y);
    case 11:
        return x*y + x*z;
    case 12:
        return sin(x)*sin(y) + x*sin(x);
    case 13:
        return x*z + y*z;
    case 14:
        return x + (3+y);
    case 15:
        return x + (3-y);
    case 16:
        return (3+x) + y;
    case 17:
        return (3-y) + z;
    case 18:
        return 0.0*y;
    case 19:
        return 1.0*(z*sin(x));
    case 20:
        return 3*4;
    case 21:
        return 3*(2*x);
    case 22:
        return 3*(4/x);
    case 23:
        return (3*sin(x)) * (4*sin(y));
    case 24:
        return z * (x/z);
    case 25:
        return (x*y) * (z/x);
    case 26:
        return (x*y) * (z/y);
    case 27:
        return (sin(x)/sin(y)) * sin(y);
    case 28:
        return x * pow(x,3);
    case 29:
        return pow(atan2(x,y),2) * atan2(x,y);
    case 30:
        return (x/y) * pow(y,3);
    case 31:
        return x*(3*y);
    case 32:
        return x*(3/y);
    case 33:
        return 3*(4-z);
    case 34:
        return 2*(4*x+sin(y));
    case 35:
        return 3*(x+4*y);
    case 36:
        return (4*x) * sin(y);
    case 37:
        return (4/x) * sin(y);
    case 38:
        return (4*x) * y;
    case 39:
        return (4/x) * y;
    case 40:
        return (x/y) * (z/sin(x));
    case 41:
        return pow(x,3) * pow(x,-2.5);
    case 42:
        return 0/x;
    case 43:
        return 4./8.;
    case 44:
        return (sin(x)+sin(y))/(sin(x)+sin(y));
    case 45:
        return (4*x)/3;
    case 46:
        return (4/x)/3;
    case 47:
        return (4*x) / (3*y);
    case 48:
        return (4/x) / (3*y);
    case 49:
        return (4*x) / (3/y);
    case 50:
        return (4/x) / (3/y);
    case 51:
        return x/4;
    case 52:
        return x / (4*y);
    case 53:
        return x / (4/y);
    case 54:
        return x / (x*y);
    case 55:
        return sin(x) / (x*sin(x));
    case 56:
        return x / (x/y);
    case 57:
        return (x*y) / x;
    case 58:
        return (x*sin(x)) / sin(x);
    case 59:
        return (x/y)/x;
    case 60:
        return (x/y) / (sin(y)/sin(x));
    case 61:
        return x / (y/z);
    case 62:
        return (x/y) / z;
    case 63:
        return pow(x,3.3) / pow(x,1.2);
    case 64:
        return x / pow(y,3);
    case 65:
        return (x*y*z*sin(x)*sin(y)) / z;
    case 66:
        return (x*y*z*sin(x)*sin(y)) / ((x+1)*(y+1)*z*(sin(x)+1)*(sin(y)+1));
    case 67:
        return pow(sin(x),0);
    case 68:
        return pow(sin(x),1.0);
    case 69:
        return pow(0.0, x);
    case 70:
        return pow(sin(x), -1.0);
    case 71:
        return 3 + (-1.0)*x;
    case 72:
        return 3 - sin(x);
    case 73:
        return 3 + sin(x);
    case 74:
        return x - y;
    case 75:
        return x + y;
    case 76:
        return x - sin(y);
    case 77:
        return x + sin(y);
    case 78:
        return -x + sin(y);
    case 79:
        return (x+y) + (-1)*z;
    case 80:
        return -sin(x) + sin(y);
    case 81:
        return -sin(x) + sin(y+1);
    case 82:
        return 3*x;
    case 83:
        return 3*sin(x);
    case 84:
        return x*z;
    case 85:
        return x*sin(y);
    case 86:
        return sin(x)*sin(x);
    case 87:
        return 3/x;
    case 88:
        return 2/sin(x);
    case 89:
        return y/z;
    case 90:
        return x/sin(x);
    case 91:
        return sin(x)/x;
    case 92:
        return (sin(x)*sin(y)+sin(z)) / (sin(x)+sin(z));
    case 93:
        return (sin(x)+sin(z)) / (sin(x)*sin(y)+sin(z));
    case 94:
        return pow(sin(x),2.5);
    case 95:
        return pow(sin(y),3);
    case 96:
        return pow(x,y+1);
    case 97:
        return pow(x+1,y);
    default:
        amrex::Abort("Unknown case "+std::to_string(icase));
        return 0.0;
    }
}
