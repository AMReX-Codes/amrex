#include <AMReX.H>
#include <cmath>

// This is intentional. Cannot have std:: in amrex::Parser expressions.
using std::atan2;
using std::pow;
using std::sin;
using std::sqrt;

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
    case 98:
        return x+(y+z);
    case 99:
        return x+(y-z);
    case 100:
        return x+(y+-z);
    case 101:
        return x+(y*z);
    case 102:
        return x+(y*-z);
    case 103:
        return x+-(y*z);
    case 104:
        return x+(y/z);
    case 105:
        return x+-(y/z);
    case 106:
        return x+(-y/z);
    case 107:
        return x+(y/-z);
    case 108:
        return x+(-y/-z);
    case 109:
        return -x+(y*z);
    case 110:
        return -x+(y/z);
    case 111:
        return (x+y) + -z;
    case 112:
        return (x-y) + -z;
    case 113:
        return (x+-y) + -z;
    case 114:
        return (x-(-y)) + -z;
    case 115:
        return (-x+y) + -z;
    case 116:
        return (x*y) + -z;
    case 117:
        return x * (y+z);
    case 118:
        return x * (y-z);
    case 119:
        return x * (-y+z);
    case 120:
        return x * (y*z);
    case 121:
        return x * (-y*-z);
    case 122:
        return x * (y/z);
    case 123:
        return x * (-y/-z);
    case 124:
        return x / (y+z);
    case 125:
        return x / (y-z);
    case 126:
        return x / (-y+z);
    case 127:
        return x / (y*z);
    case 128:
        return x / (-y*-z);
    case 129:
        return x / (y/z);
    case 130:
        return x / (-y/-z);
    case 131:
        return (x*y) / z;
    case 132:
        return (x+y) / z;
    case 133:
        return (x-y) / z;
    case 134:
        return (x/y) / z;
    case 135:
        return (-x/y) / -z;
    case 136:
        return -(3*x)*(y*z);
    case 137:
        return -(3*x)/(y*z);
    case 138:
        return x - (y+z);
    case 139:
        return x - (y-z);
    case 140:
        return x - (y*z);
    case 141:
        return x - (y/z);
    case 142:
        return x*x;
    case 143:
        return pow(x,2);
    case 144:
        return x*x*x;
    case 145:
        return pow(x,3);
    case 146:
        return sqrt(x*x+y*y);
    case 147:
        return sqrt(pow(x,2)+pow(y,2));
    case 148:
        return sqrt(pow(x+y,2)+pow(x+z,2));
    case 149:
        return sqrt(x*x+y*y+z*z);
    case 150:
        return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    case 151:
        return sqrt(pow(x+y,2)+pow(x+z,2)+pow(y+z,2));
    case 152:
        return pow(x,0.5);
    case 153:
        return pow(pow(pow(x,0.3),0.4),0.5);
    case 154:
        return pow(pow(pow(x,0.3),y),0.5);
    case 155:
        return pow(x,-1);
    case 156:
        return (10.1+x) + y;
    case 157:
        return (10.1+x) - y;
    case 158:
        return (10.1+x) * y;
    case 159:
        return (x+y) * 10.1;
    case 160:
        return (10.1+x) * 11.2;
    case 161:
        return (10.1+x) / y;
    case 162:
        return (10.1-y) - z;
    case 163:
        return (9.1-y) * z;
    case 164:
        return (x-y) * 9.1;
    case 165:
        return (7.1-y) / z;
    case 166:
        return (7.1*y) + z;
    case 167:
        return (x*y) + 7.3;
    case 168:
        return (7.1*y) - z;
    case 169:
        return (7.1*y) * z;
    case 170:
        return (7.1*y) / z;
    case 171:
        return (7.1/y) + z;
    case 172:
        return (x/y) + 7.3;
    case 173:
        return (7.1/y) + 7.3;
    case 174:
        return (7.1/y) - z;
    case 175:
        return (7.1/y) * z;
    case 176:
        return 7.1 - (y + z);
    case 177:
        return 7.1 - ( y * z );
    case 178:
        return 7.1 - ( y / z );
    case 179:
        return 7.1 / ( y + z );
    case 180:
        return x / ( 7.2 + z );
    case 181:
        return 7.1 / ( 7.2 + z );
    case 182:
        return 7.1 / ( y - z );
    case 183:
        return x / ( 7.2 - z );
    case 184:
        return 7.1 / ( y - 7.3 );
    case 185:
        return 7.1 / ( 7.2 - z );
    default:
        amrex::Abort("Unknown case "+std::to_string(icase));
        return 0.0;
    }
}
