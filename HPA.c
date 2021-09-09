#include <stdio.h>
#include <math.h>

struct num {
    float val;
    float err;
};
typedef struct num NUM;

void scalar_init(float value, NUM *x);
void scalar_sum(NUM x, NUM y, NUM *z);
void scalar_sub(NUM x, NUM y, NUM *z);
void scalar_mult(NUM x, NUM y, NUM *z);
void scalar_div(NUM x, NUM y, NUM *z);

int main(){

    long double a = -0.000587, b = 26584538;
    NUM x, y;
    scalar_init((float)a, &x);
    scalar_init((float)b, &y);
    
    printf("SUM = %.10Lf\n\n", a + 10*b);

    for (int i = 0; i<10; i++){
        scalar_sum(x, y, &x);
    }
    
    //scalar_sub(z, x2, &z);

    printf("RES = %.10f\n", x.val);
    printf("ERR = %.10f\n\n", x.err);
    printf("SUM2 = %.10f\n\n", x.val + x.err);


    float a1 = -0.000587, b1 = 26584538;
    printf("FLOAT = %.10f\n\n", a1 + 10*b1);

    return 0;
}

void scalar_init(float value, NUM *x){
    float p = value * (1 + pow(2, 12) );
    (*x).val = value - p + p;
    (*x).err = value - (*x).val;
}

void scalar_sum(NUM x, NUM y, NUM *z){
    float r, s;
    r = x.val + y.val;
    if (fabs(x.val) > fabs(y.val)) s = x.val - r + y.val + y.err + x.err;
    else s = y.val - r + x.val + x.err + y.err;
    (*z).val = r + s;
    (*z).err = (r - (*z).val) + s;
}

void scalar_sub(NUM x, NUM y, NUM *z){
    float r, s;
    r = x.val - y.val;
    if (fabs(x.val) > fabs(y.val)) s = x.val - r - y.val - y.err + x.err;
    else s = -y.val - r + x.val + x.err - y.err;
    (*z).val = r + s;
    (*z).err = r - (*z).val + s;
}

void scalar_mult(NUM x, NUM y, NUM *z){
    float p, q, hx, tx, hy, ty;
    p = x.val * (1 + pow(2, 12));
    hx = p - (p - x.val); tx = x.val - hx;
    p = y.val * (1 + pow(2, 12));
    hy = p - (p - y.val); ty = y.val - hy;
    p = hx * hy;
    q = hx * ty + tx * hy;
    float c = p + q;
    float cc = p - c + q + tx * ty;

    cc = x.val * y.err + x.err * y.val + cc;
    (*z).val = c + cc;
    (*z).err = c - (*z).val + cc;
}

void scalar_div(NUM x, NUM y, NUM *z){
    float c = x.val / y.val;

    float p, q, hc, tc, hy, ty;
    p = c * (1 + pow(2, 12));
    hc = p - (p - c); tc = c - hc;
    p = y.val * (1 + pow(2, 12));
    hy = p - (p - y.val); ty = y.val - hy;
    p = hc * hy;
    q = hc * ty + tc * hy;
    float u = p + q;
    float uu = p - u + q + tc * ty;

    float cc = (x.val - u - uu + x.err - c * y.err) / y.val;
    (*z).val = c + cc;
    (*z).err = c - (*z).val + cc;
}