#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <stdio.h>

void main(int argc, char *argv[])
{
    arb_t d;
    arb_t m;
    arb_t one;
    arb_t two;
    arb_t difference;
    arb_t sum;
    arb_t quotient;
    arb_t res;
    arb_set_str(d, argv[1], 100);
    arb_set_str(m, argv[2], 100);
    arb_set_str(one, "1", 100);
    arb_set_str(two, "2", 100);
    arb_sub(difference, one, d, 100);
    arb_add(sum, difference, m, 100);
    arb_div(quotient, sum, two, 100);
    arb_digamma(res, quotient, 100);
    arb_printn(res, 20, 0);
}