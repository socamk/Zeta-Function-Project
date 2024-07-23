#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <stdio.h>


void main(int argc, char *argv[])
{
    arb_t three;
    arb_t two;
    arb_t x;
    arb_t y;
    acb_t difference;
    acb_t quotient;
    arb_t a;
    arb_t b;
    acb_t z;
    acb_t res;
    acb_init(res);
    arb_set_str(x, argv[1], 100);
    arb_set_str(y, argv[2], 100);
    acb_set_arb_arb(z, x, y);
    arb_set_str(three, "3", 100);
    arb_set_str(two, "-2", 100);
    acb_sub_arb(difference, z, three, 100);
    acb_div_arb(quotient, difference, two, 100);
    acb_digamma(res, quotient, 100);
    arb_set(a, acb_realref(res));
    arb_set(b, acb_imagref(res));
    arb_printn(a, 50, 0);
    printf("\n");
    arb_printn(b, 50, 0);
}