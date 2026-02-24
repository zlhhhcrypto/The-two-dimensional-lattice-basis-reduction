#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <mpfr.h>
#include <math.h>

mpz_t q0, q1, q2;
long c = 100;

void Half_Gaussian(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);
bool check_strong_coherent(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);
long get_bit_length(mpz_t ax, mpz_t ay);
double get_double_length(mpz_t ax, mpz_t ay);
int detdim2mat(mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2);
bool is_coherent(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);
void fixup(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t ax0, mpz_t ay0, mpz_t bx0, mpz_t by0, mpz_t ta1, mpz_t ta2, mpz_t tb1, mpz_t tb2, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2, long m);
bool is_terminal(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);
bool is_straddles(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, double length);
void half_gaussian(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2);
void make_addmissible(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);

int main(int argc, char *argv[]){
	mpz_t a[2], b[2];
	mpz_t tmp1, tmp2, tmp3, tmp4, q, tmp5, tmp6;
	
	mpz_init(a[0]); mpz_init(a[1]);mpz_init(b[0]); mpz_init(b[1]);
	mpz_init(tmp1);	mpz_init(tmp2);	mpz_init(tmp3);	mpz_init(tmp4);
	mpz_init(tmp5);	mpz_init(tmp6);	mpz_init(q);
	
	mpz_init(q0);	mpz_init(q1);	mpz_init(q2);
	mpz_set_ui(q0, 0);
	mpz_set_ui(q1, 0);
	mpz_set_ui(q2, 0);
	
	clock_t start, finish;
	
	printf("\nHalf Gaussian algorithm for all kinds of lattices.\n");
	
	FILE *fp = fopen("Input.txt", "r");
	gmp_fscanf(fp, "%Zd %Zd %Zd %Zd", a[0], a[1], b[0], b[1]);
	fclose(fp);
	
	if(mpz_cmp(a[0], b[0]) < 0){
		mpz_swap(a[0], b[0]);
	}
	mpz_set_si(a[1], 0);
    	
	int euctime = 0;
	start = clock();
	Half_Gaussian(a[0], a[1], b[0], b[1]);
	
	mpz_mul(tmp1, b[0], b[0]);
	mpz_mul(tmp2, b[1], b[1]);
	mpz_add(tmp1, tmp1, tmp2);
	
	mpz_sub(tmp3, a[0], b[0]);
	mpz_mul(tmp3, tmp3, tmp3);
	mpz_sub(tmp4, a[1], b[1]);
	mpz_mul(tmp4, tmp4, tmp4);
	mpz_add(tmp3, tmp3, tmp4);
	
	if(mpz_cmp(tmp3, tmp1) < 0){
		mpz_sub(a[0], a[0], b[0]);
		mpz_sub(a[1], a[1], b[1]);
		mpz_swap(a[0], b[0]);
		mpz_swap(a[1], b[1]);
	}
	
	mpz_mul(tmp1, a[0], b[0]);
	mpz_mul(tmp2, a[1], b[1]);
	mpz_add(tmp1, tmp1, tmp2);
	mpz_mul(tmp3, b[0], b[0]);
	mpz_mul(tmp4, b[1], b[1]);
	mpz_add(tmp3, tmp3, tmp4);
	mpz_fdiv_q(q, tmp1, tmp3);
	
	mpz_submul(a[0], q, b[0]);
	mpz_submul(a[1], q, b[1]);
	mpz_sub(tmp1, a[0], b[0]);
	mpz_sub(tmp2, a[1], b[1]);
	
	mpz_mul(tmp3, a[0], a[0]);
	mpz_mul(tmp4, a[1], a[1]);
	mpz_add(tmp3, tmp3, tmp4);
	mpz_mul(tmp5, tmp1, tmp1);
	mpz_mul(tmp6, tmp2, tmp2);
	mpz_add(tmp5, tmp5, tmp6);
	if(mpz_cmp(tmp3, tmp5) > 0){
		mpz_set(a[0], tmp1);
		mpz_set(a[1], tmp2);
	}
	
    finish = clock();	
	
	FILE *fo;
	fo = fopen("Output.txt", "w");
	gmp_fprintf(fo, "%Zd\n",a[0]);
	gmp_fprintf(fo, "%Zd\n",a[1]);
	gmp_fprintf(fo, "%Zd\n",b[0]);
	gmp_fprintf(fo, "%Zd\n",b[1]);
	fclose(fo);
	
	printf("Partial Euclidean cycle times = %d\n",euctime);
	printf("Prog time: %lf\n", (double)(finish-start)/CLOCKS_PER_SEC);
	
	mpz_clear(a[0]); mpz_clear(a[1]); mpz_clear(b[0]); mpz_clear(b[1]);
	mpz_clear(q0);	mpz_clear(q1);	mpz_clear(q2);
	mpz_init(tmp1);	mpz_init(tmp2);	mpz_init(tmp3);	mpz_init(tmp4);
	mpz_init(tmp5);	mpz_init(tmp6);	mpz_init(q);

    return 0;
}

void make_addmissible(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by){
	mpz_t tmp1, tmp2, tmp3, tmp4;
	mpz_init(tmp1);	mpz_init(tmp2);	mpz_init(tmp3);	mpz_init(tmp4);
	
	mpz_mul(tmp1, ax, ax);	mpz_mul(tmp2, ay, ay);	mpz_add(tmp1, tmp1, tmp2);
	mpz_mul(tmp3, bx, bx);	mpz_mul(tmp4, by, by);	mpz_add(tmp3, tmp3, tmp4);
	if(mpz_cmp(tmp1, tmp3) < 0){
		mpz_swap(ax, bx);
		mpz_swap(ay, by);
	}
	if(is_coherent(ax, ay, bx, by) == false){
		mpz_neg(bx, bx);
		mpz_neg(by, by);
	}
	return;
}


void Half_Gaussian(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by){
	mpz_t a1, a2, b1, b2, tax, tay, tbx, tby, tmp1, tmp2, tmp3, tmp4, q;
	mpz_init(a1);	mpz_init(a2);	mpz_init(b1);	mpz_init(b2);
	mpz_init(tax);	mpz_init(tay);	mpz_init(tbx);	mpz_init(tby);
	mpz_init(tmp1);	mpz_init(tmp2);	mpz_init(tmp3);	mpz_init(tmp4);
	mpz_init(q);
	while(1){
		make_addmissible(ax, ay, bx, by);
		
		mpz_set(tax, ax);	mpz_set(tay, ay);
		mpz_set(tbx, bx);	mpz_set(tby, by);
		mpz_set_si(a1, 1);	mpz_set_si(a2, 0);	mpz_set_si(b1, 0);	mpz_set_si(b2, 1);

    	half_gaussian(ax, ay, bx, by, a1, a2, b1, b2);	
    		
    		
    	int det = detdim2mat(a1, a2, b1, b2);
    	mpz_mul(tmp1, tax, b2);	mpz_mul(tmp2, tbx, b1);	mpz_sub(ax, tmp1, tmp2);
    	mpz_mul(tmp1, tay, b2);	mpz_mul(tmp2, tby, b1);	mpz_sub(ay, tmp1, tmp2);
    	mpz_mul(tmp1, tbx, a1);	mpz_mul(tmp2, tax, a2);	mpz_sub(bx, tmp1, tmp2);
    	mpz_mul(tmp1, tby, a1);	mpz_mul(tmp2, tay, a2);	mpz_sub(by, tmp1, tmp2);
    	if(det < 0){
    		mpz_neg(ax, ax);	mpz_neg(ay, ay);
    		mpz_neg(bx, bx);	mpz_neg(by, by);
    	}
    		
    	if(is_terminal(ax, ay, bx, by) == true){
    		break;
    	}
    		
    	mpz_mul(tmp1, ax, bx);	mpz_mul(tmp2, ay, by);	mpz_add(tmp1, tmp1, tmp2);
    	mpz_mul(tmp3, bx, bx);	mpz_mul(tmp4, by, by);	mpz_add(tmp3, tmp3, tmp4);
    		
    	if(mpz_sgn(tmp1) == 0 || mpz_sgn(tmp3) == 0){
    		break;
    	}
    		
    	mpz_fdiv_q(q, tmp1, tmp3);
    	mpz_submul(ax, q, bx);	mpz_submul(ay, q, by);
    	mpz_swap(ax, bx);	mpz_swap(ay, by);
	}
    	
	mpz_clear(a1);	mpz_clear(a2);	mpz_clear(b1);	mpz_clear(b2);
	mpz_clear(tax);	mpz_clear(tay);	mpz_clear(tbx);	mpz_clear(tby);
	mpz_clear(tmp1);	mpz_clear(tmp2);
	mpz_clear(tmp3);	mpz_clear(tmp4);
	mpz_clear(q);
	return;
}

bool check_strong_coherent(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by){
	mpf_t fax, fay, fbx, fby, tmp1, tmp2, tmp3, tmp4, numerator, denominator, cos, strong_coherent, zero;
	
	mpf_init(fax);	mpf_init(fay);	mpf_init(fbx);	mpf_init(fby);
	mpf_init(tmp1);	mpf_init(tmp2);	mpf_init(tmp3);	mpf_init(tmp4);
	mpf_init(numerator);	mpf_init(denominator);
	mpf_init(cos);	mpf_init(strong_coherent);
	mpf_init(zero);

	mpf_set_d(zero, 0.0);
	
	mpf_set_d(strong_coherent, 0.5);
	
	mpf_set_z(fax, ax);	mpf_set_z(fay, ay);
	mpf_set_z(fbx, bx);	mpf_set_z(fby, by);
	
	if(mpf_cmp(fbx, zero) == 0 && mpf_cmp(fby, zero) == 0){
		return false;
	}
	mpf_mul(tmp1, fax, fbx);	mpf_mul(tmp2, fay, fby);	mpf_add(numerator, tmp1, tmp2);
	mpf_mul(tmp1, fax, fax);	mpf_mul(tmp2, fay, fay);	mpf_add(tmp1, tmp1, tmp2);	mpf_sqrt(tmp1, tmp1);
	mpf_mul(tmp3, fbx, fbx);	mpf_mul(tmp4, fby, fby);	mpf_add(tmp3, tmp3, tmp4);	mpf_sqrt(tmp3, tmp3);
	mpf_mul(denominator, tmp1, tmp3);
	mpf_div(cos, numerator, denominator);
	
	if(mpf_cmp(cos, strong_coherent) <= 0){
		
		mpf_clear(fax);	mpf_clear(fay);	mpf_clear(fbx);	mpf_clear(fby);
		mpf_clear(tmp1);	mpf_clear(tmp2);	mpf_clear(tmp3);	mpf_clear(tmp4);
		mpf_clear(numerator);	mpf_clear(denominator);
		mpf_clear(cos);	mpf_clear(strong_coherent);
		
		return false;
	}
	
	mpf_clear(fax);	mpf_clear(fay);	mpf_clear(fbx);	mpf_clear(fby);
	mpf_clear(tmp1);	mpf_clear(tmp2);	mpf_clear(tmp3);	mpf_clear(tmp4);
	mpf_clear(numerator);	mpf_clear(denominator);
	mpf_clear(cos);	mpf_clear(strong_coherent);
	
	return true;
}

double get_double_length(mpz_t ax, mpz_t ay){
	mpz_t tmp1, tmp2;
	
	mpz_init(tmp1);	mpz_init(tmp2);
	
	mpz_mul(tmp1, ax, ax);	mpz_mul(tmp2, ay, ay);	mpz_add(tmp1, tmp1, tmp2);
	
	mpfr_t f_sum, f_sqrt, f_log2;
    mpfr_inits2(256, f_sum, f_sqrt, f_log2, NULL);
    	
    mpfr_set_z(f_sum, tmp1, MPFR_RNDN);
    mpfr_sqrt(f_sqrt, f_sum, MPFR_RNDN);
    mpfr_log2(f_log2, f_sqrt, MPFR_RNDN);
    	
    double result = mpfr_get_d(f_log2, MPFR_RNDN);

    mpfr_clears(f_sum, f_sqrt, f_log2, NULL);
	
	mpz_clear(tmp1);	mpz_clear(tmp2);
	
	return result;
}

long get_bit_length(mpz_t ax, mpz_t ay){
	mpz_t tmp1, tmp2;
	
	mpz_init(tmp1);	mpz_init(tmp2);
	
	mpz_mul(tmp1, ax, ax);	mpz_mul(tmp2, ay, ay);	mpz_add(tmp1, tmp1, tmp2);
	mpfr_t f_sum, f_sqrt, f_log2;
    mpfr_inits2(256, f_sum, f_sqrt, f_log2, NULL);
    	
    mpfr_set_z(f_sum, tmp1, MPFR_RNDN);
    mpfr_sqrt(f_sqrt, f_sum, MPFR_RNDN);
    mpfr_log2(f_log2, f_sqrt, MPFR_RNDN);
    	
    double result = mpfr_get_d(f_log2, MPFR_RNDN);
    long ceil_result = (long)ceil(result);
    mpfr_clears(f_sum, f_sqrt, f_log2, NULL);

	mpz_clear(tmp1);	mpz_clear(tmp2);
	
	return ceil_result;
}

void direct_compute(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2){
	mpz_t tmp1, tmp2, tmp3, tmp4, q;
	
	mpz_init(tmp1);	mpz_init(tmp2);	mpz_init(tmp3);	mpz_init(tmp4);
	mpz_init(q);
	
	mpz_mul(tmp1, ax, bx);	mpz_mul(tmp2, ay, by);	mpz_add(tmp1, tmp1, tmp2);
	mpz_mul(tmp3, bx, bx);	mpz_mul(tmp4, by, by);	mpz_add(tmp3, tmp3, tmp4);
	
	mpz_fdiv_q(q, tmp1, tmp3);
	
	mpz_set(q0, q1);
	mpz_set(q1, q2);
	mpz_set(q2, q);
	
	mpz_submul(ax, q, bx);	
	mpz_submul(ay, q, by);
	mpz_swap(ax, bx);
	mpz_swap(ay, by);
	
	mpz_set(a1, q);
	mpz_set_ui(a2, 1);
	mpz_set_ui(b1, 1);
	mpz_set_ui(b2, 0);
	
	mpz_clear(tmp1);	mpz_clear(tmp2);
	mpz_clear(tmp3);	mpz_clear(tmp4);
	mpz_clear(q);
	
	return;
}

int detdim2mat(mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2){
	mpz_t multp1, multp2;
	mpz_init(multp1); mpz_init(multp2);
	mpz_mul(multp1, a1, b2);
	mpz_mul(multp2, a2, b1);
	if(mpz_cmp(multp1, multp2) < 0){
		mpz_clear(multp1); mpz_clear(multp2);
		return -1;
	}else{
		mpz_clear(multp1); mpz_clear(multp2);
		return 1;
	} 
}

bool is_coherent(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by){
	mpf_t fax, fay, fbx, fby, tmp1, tmp2, tmp3, tmp4, numerator, denominator, cos, strong_coherent, zero;
	
	mpf_init(fax);	mpf_init(fay);	mpf_init(fbx);	mpf_init(fby);
	mpf_init(tmp1);	mpf_init(tmp2);	mpf_init(tmp3);	mpf_init(tmp4);
	mpf_init(numerator);	mpf_init(zero);
	
	mpf_set_z(fax, ax);	mpf_set_z(fay, ay);
	mpf_set_z(fbx, bx);	mpf_set_z(fby, by);
	mpf_set_d(zero, 0.0);
	
	if(mpf_cmp(fbx, zero) == 0 && mpf_cmp(fby, zero) == 0){
		return true;
	}
	
	mpf_mul(tmp1, fax, fbx);	mpf_mul(tmp2, fay, fby);	mpf_add(numerator, tmp1, tmp2);
	
	if(mpf_cmp(numerator, zero) < 0){
		
		mpf_clear(fax);	mpf_clear(fay);	mpf_clear(fbx);	mpf_clear(fby);
		mpf_clear(tmp1);	mpf_clear(tmp2);	mpf_clear(tmp3);	mpf_clear(tmp4);
		mpf_clear(numerator);	mpf_clear(zero);
		
		return false;
	}
	
	mpf_clear(fax);	mpf_clear(fay);	mpf_clear(fbx);	mpf_clear(fby);
	mpf_clear(tmp1);	mpf_clear(tmp2);	mpf_clear(tmp3);	mpf_clear(tmp4);
	mpf_clear(numerator);	mpf_clear(zero);
	
	return true;
}

void fixup(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t ax0, mpz_t ay0, mpz_t bx0, mpz_t by0, mpz_t ta1, mpz_t ta2, mpz_t tb1, mpz_t tb2, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2, long m){

	mpz_t tmp1, tmp2, tmp3, tmp4, axp, ayp, bxp, byp, ax0p, ay0p, bx0p, by0p, one, uxpp, uypp, q;
	
	mpz_init(tmp1);	mpz_init(tmp2);	mpz_init(tmp3);	mpz_init(tmp4);
	mpz_init(axp);	mpz_init(ayp);	mpz_init(bxp);	mpz_init(byp);
	mpz_init(ax0p);	mpz_init(ay0p);	mpz_init(bx0p);	mpz_init(by0p);
	mpz_init(one);	mpz_init(uxpp);	mpz_init(uypp);	mpz_init(q);
	
	mpz_set_ui(one, 1);
	
	int det = detdim2mat(ta1, ta2, tb1, tb2);
	
	mpz_mul(tmp1, ax, tb2);	mpz_mul(tmp2, bx, tb1);	mpz_sub(axp, tmp1, tmp2);
	mpz_mul(tmp1, ay, tb2);	mpz_mul(tmp2, by, tb1);	mpz_sub(ayp, tmp1, tmp2);
	mpz_mul(tmp1, bx, ta1);	mpz_mul(tmp2, ax, ta2);	mpz_sub(bxp, tmp1, tmp2);
	mpz_mul(tmp1, by, ta1);	mpz_mul(tmp2, ay, ta2);	mpz_sub(byp, tmp1, tmp2);
	if(det < 0){
		mpz_neg(axp, axp);	mpz_neg(ayp, ayp);
		mpz_neg(bxp, bxp);	mpz_neg(byp, byp);
	}
	
	mpz_mul(tmp1, ax0, tb2);	mpz_mul(tmp2, bx0, tb1);	mpz_sub(ax0p, tmp1, tmp2);
	mpz_mul(tmp1, ay0, tb2);	mpz_mul(tmp2, by0, tb1);	mpz_sub(ay0p, tmp1, tmp2);
	mpz_mul(tmp1, bx0, ta1);	mpz_mul(tmp2, ax0, ta2);	mpz_sub(bx0p, tmp1, tmp2);
	mpz_mul(tmp1, by0, ta1);	mpz_mul(tmp2, ay0, ta2);	mpz_sub(by0p, tmp1, tmp2);
	if(det < 0){
		mpz_neg(ax0p, ax0p);	mpz_neg(ay0p, ay0p);
		mpz_neg(bx0p, bx0p);	mpz_neg(by0p, by0p);
	}
	
	double u0p_double_length = get_double_length(ax0p, ay0p);
	double v0p_double_length = get_double_length(bx0p, by0p);
	double u0_double_length = get_double_length(ax0, ay0);
	
	mpz_set(a1, ta1);
	mpz_set(a2, ta2);
	mpz_set(b1, tb1);
	mpz_set(b2, tb2);
	
	if((u0p_double_length >= (double)(u0_double_length-(double)m/2+c)) && ((double)(u0_double_length-(double)m/2+c) > v0p_double_length)){
		double up_double_length = get_double_length(axp, ayp);
		double vp_double_length = get_double_length(bxp, byp);
		if(up_double_length <= vp_double_length){
			
			mpz_submul(a1, q2, b1);
			mpz_submul(a2, q2, b2);
			mpz_swap(a1, b1);
			mpz_swap(a2, b2);
			
			mpz_set(q2, q1);
			mpz_set(q1, q0);
		}
		if(is_coherent(axp, ayp, bxp, byp) == false){
			mpz_add(uxpp, axp, bxp);
			mpz_add(uypp, ayp, byp);
			double upp_double_length = get_double_length(uxpp, uypp);
			if(mpz_cmp(q2, one) > 0){
				if(up_double_length > upp_double_length){
					mpz_submul(a1, q2, b1);
					mpz_submul(a2, q2, b2);
					mpz_swap(a1, b1);
					mpz_swap(a2, b2);
			
					mpz_sub(q2, q2, one);
					
					mpz_swap(a1, b1);
					mpz_swap(a2, b2);
					mpz_addmul(a1, q2, b1);
					mpz_addmul(a2, q2, b2);
				
				}else{
					mpz_submul(a1, q2, b1);
					mpz_submul(a2, q2, b2);
					mpz_swap(a1, b1);
					mpz_swap(a2, b2);
			
					mpz_sub(q2, q2, one);
					
					mpz_swap(a1, b1);
					mpz_swap(a2, b2);
					mpz_addmul(a1, q2, b1);
					mpz_addmul(a2, q2, b2);
					
					mpz_set(q0, q1);
					mpz_set(q1, q2);
					mpz_set(q2, one);
					
					mpz_swap(a1, b1);
					mpz_swap(a2, b2);
					mpz_add(a1, a1, b1);
					mpz_add(a2, a2, b2);
				}
			}else{
				if(upp_double_length > up_double_length){
					mpz_submul(a1, q2, b1);
					mpz_submul(a2, q2, b2);
					mpz_swap(a1, b1);
					mpz_swap(a2, b2);
			
					mpz_set(q2, q1);
					mpz_set(q1, q0);
					
				}else{
					mpz_submul(a1, q2, b1);
					mpz_submul(a2, q2, b2);
					mpz_swap(a1, b1);
					mpz_swap(a2, b2);
			
					mpz_set(q2, q1);
					mpz_set(q1, q0);
					
					mpz_submul(a1, q2, b1);
					mpz_submul(a2, q2, b2);
					mpz_swap(a1, b1);
					mpz_swap(a2, b2);
			
					mpz_set(q2, q1);
					mpz_set(q1, q0);
				}
			}
		}
	}
	
	mpz_clear(tmp1);	mpz_clear(tmp2);	mpz_clear(tmp3);	mpz_clear(tmp4);
	mpz_clear(axp);	mpz_clear(ayp);	mpz_clear(bxp);	mpz_clear(byp);
	mpz_clear(ax0p);	mpz_clear(ay0p);	mpz_clear(bx0p);	mpz_clear(by0p);
	mpz_clear(one);	mpz_clear(uxpp);	mpz_clear(uypp);	mpz_clear(q);
	
	return;
}

bool is_terminal(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by){
	mpz_t tmp1, tmp2, tmp3, tmp4;
	mpz_init(tmp1);	mpz_init(tmp2);	mpz_init(tmp3);	mpz_init(tmp4);
	
	double u_length = get_double_length(ax, ay);
	double v_length = get_double_length(bx, by);
	
	if(u_length >= v_length){
		mpz_mul(tmp1, ax, bx);	mpz_mul(tmp2, ay, by);	mpz_add(tmp1, tmp1, tmp2);
		mpz_mul(tmp3, bx, bx);	mpz_mul(tmp4, by, by);	mpz_add(tmp3, tmp3, tmp4);
		if(mpz_cmpabs(tmp1, tmp3) < 0){
			return true;
		}
		return false;
	}else{
		mpz_mul(tmp1, ax, bx);	mpz_mul(tmp2, ay, by);	mpz_add(tmp1, tmp1, tmp2);
		mpz_mul(tmp3, ax, ax);	mpz_mul(tmp4, ay, ay);	mpz_add(tmp3, tmp3, tmp4);
		if(mpz_cmpabs(tmp1, tmp3) < 0){
			return true;
		}
		return false;
	}
}

bool is_straddles(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, double length){
	double u_length = get_double_length(ax, ay);
	double v_length = get_double_length(bx, by);
	if(u_length >= length && length > v_length){
		return true;
	}
	return false;
}

void half_gaussian(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2){
	
	mpz_t ta1, ta2, tb1, tb2, ax0, ay0, bx0, by0;
	mpz_t a11, a12, b11, b12, a21, a22, b21, b22;
	mpz_t tmp1, tmp2, tmp3, tmp4;
	mpz_t q, zero;
	mpz_t tax0, tay0, tbx0, tby0;
	
	mpz_init(ta1);	mpz_init(ta2);	mpz_init(tb1);	mpz_init(tb2);
	mpz_init(ax0);	mpz_init(ay0);	mpz_init(bx0);	mpz_init(by0);
	mpz_init(a11);	mpz_init(a12);	mpz_init(b11);	mpz_init(b12);
	mpz_init(a21);	mpz_init(a22);	mpz_init(b21);	mpz_init(b22);
	mpz_init(tmp1);	mpz_init(tmp2);	mpz_init(tmp3);	mpz_init(tmp4);
	mpz_init(q);	mpz_init(zero);
	mpz_init(tax0);	mpz_init(tay0);	mpz_init(tbx0);	mpz_init(tby0);
	
	mpz_set_ui(zero, 0);
	
	if(check_strong_coherent(ax, ay, bx, by) == false){
		mpz_set_ui(a1, 1);
		mpz_set_ui(a2, 0);
		mpz_set_ui(b1, 0);
		mpz_set_ui(b2, 1);
		return;
	}	
	long n = get_bit_length(ax, ay);
	double u_double_length = get_double_length(ax, ay);
	long m = (long)ceil((double)n / 2);
	
	mpz_set_ui(a11, 1);	mpz_set_ui(a12, 0);	mpz_set_ui(b11, 0);	mpz_set_ui(b12, 1);
	int flag = 0;
	double v_length = get_double_length(bx, by);
	if((v_length < (double)(n - (double)m/2 + c)) || (m <= 8)){
		mpz_set_ui(a11, 1);
		mpz_set_ui(a12, 0);
		mpz_set_ui(b11, 0);
		mpz_set_ui(b12, 1);
		flag = 1;
	}
	if(flag == 0){
		if(mpz_sgn(ax) >= 0){
			mpz_cdiv_q_2exp(ax0, ax, m);
		}else{
			mpz_fdiv_q_2exp(ax0, ax, m);
		}
		if(mpz_sgn(ay) >= 0){
			mpz_cdiv_q_2exp(ay0, ay, m);
		}else{
			mpz_fdiv_q_2exp(ay0, ay, m);
		}
		if(mpz_sgn(bx) >= 0){
			mpz_fdiv_q_2exp(bx0, bx, m);
		}else{
			mpz_cdiv_q_2exp(bx0, bx, m);
		}
		if(mpz_sgn(by) >= 0){
			mpz_fdiv_q_2exp(by0, by, m);
		}else{
			mpz_cdiv_q_2exp(by0, by, m);
		}
		
		mpz_mul(tmp1, ax0, ax0);	mpz_mul(tmp2, ay0, ay0);	mpz_add(tmp1, tmp1, tmp2);
		mpz_mul(tmp3, bx0, bx0);	mpz_mul(tmp4, by0, by0);	mpz_add(tmp3, tmp3, tmp4);
		
		if(mpz_cmp(tmp1, tmp3) == 0){
			direct_compute(ax0, ay0, bx0, by0, a11, a12, b11, b12);
			flag = 1;
		}
		if(flag == 0){
			mpz_set(tax0, ax0);	mpz_set(tay0, ay0);
			mpz_set(tbx0, bx0);	mpz_set(tby0, by0);
			
			half_gaussian(ax0, ay0, bx0, by0, ta1, ta2, tb1, tb2);
			
			fixup(ax, ay, bx, by, tax0, tay0, tbx0, tby0, ta1, ta2, tb1, tb2, a11, a12, b11, b12, m);
		}
	}
	
	int det = detdim2mat(a11, a12, b11, b12);
	mpz_mul(tmp1, ax, b12);	mpz_mul(tmp2, bx, b11);	mpz_sub(ax0, tmp1, tmp2);
	mpz_mul(tmp1, ay, b12);	mpz_mul(tmp2, by, b11);	mpz_sub(ay0, tmp1, tmp2);
	mpz_mul(tmp1, bx, a11);	mpz_mul(tmp2, ax, a12);	mpz_sub(bx0, tmp1, tmp2);
	mpz_mul(tmp1, by, a11);	mpz_mul(tmp2, ay, a12);	mpz_sub(by0, tmp1, tmp2);
	
	if(det == -1){
		mpz_neg(ax0, ax0);
		mpz_neg(ay0, ay0);
		mpz_neg(bx0, bx0);
		mpz_neg(by0, by0);
	}
	
	if((is_terminal(ax0, ay0, bx0, by0) == true) || (is_straddles(ax0, ay0, bx0, by0, (double)((double)n/2+c)) == true)){
		mpz_set(a1, a11);	mpz_set(a2, a12);
		mpz_set(b1, b11);	mpz_set(b2, b12);
		
		mpz_clear(ta1);	mpz_clear(ta2);	mpz_clear(tb1);	mpz_clear(tb2);
		mpz_clear(ax0);	mpz_clear(ay0);	mpz_clear(bx0);	mpz_clear(by0);
		mpz_clear(a11);	mpz_clear(a12);	mpz_clear(b11);	mpz_clear(b12);
		mpz_clear(a21);	mpz_clear(a22);	mpz_clear(b21);	mpz_clear(b22);
		mpz_clear(tmp1);	mpz_clear(tmp2);
		mpz_clear(tmp3);	mpz_clear(tmp4);
		mpz_clear(q);
		mpz_clear(tax0);	mpz_clear(tay0);
		mpz_clear(tbx0);	mpz_clear(tby0);
		
		return;
	}
	
	mpz_set(ax, ax0);	mpz_set(ay, ay0);
	mpz_set(bx, bx0);	mpz_set(by, by0);
	
	while(get_double_length(ax, ay) > (n*3)/4+c){
		
		mpz_mul(tmp1, ax, bx);	mpz_mul(tmp2, ay, by);	mpz_add(tmp1, tmp1, tmp2);
		mpz_mul(tmp3, bx, bx);	mpz_mul(tmp4, by, by);	mpz_add(tmp3, tmp3, tmp4);
		mpz_fdiv_q(q, tmp1, tmp3);
		mpz_submul(ax, q, bx);
		mpz_submul(ay, q, by);
		mpz_swap(ax, bx);
		mpz_swap(ay, by);

		mpz_set(q0, q1);
		mpz_set(q1, q2);
		mpz_set(q2, q);

		mpz_swap(a11, b11);
		mpz_swap(a12, b12);
		mpz_addmul(a11, q, b11); 
		mpz_addmul(a12, q, b12);
		
		if((is_terminal(ax, ay, bx, by) == true) || (is_straddles(ax, ay, bx, by, (double)((double)n/2+c)) == true)){
			mpz_set(a1, a11);	mpz_set(a2, a12);
			mpz_set(b1, b11);	mpz_set(b2, b12);
			
			mpz_clear(ta1);	mpz_clear(ta2);	mpz_clear(tb1);	mpz_clear(tb2);
			mpz_clear(ax0);	mpz_clear(ay0);	mpz_clear(bx0);	mpz_clear(by0);
			mpz_clear(a11);	mpz_clear(a12);	mpz_clear(b11);	mpz_clear(b12);
			mpz_clear(a21);	mpz_clear(a22);	mpz_clear(b21);	mpz_clear(b22);
			mpz_clear(tmp1);	mpz_clear(tmp2);
			mpz_clear(tmp3);	mpz_clear(tmp4);
			mpz_clear(q);
			mpz_clear(tax0);	mpz_clear(tay0);
			mpz_clear(tbx0);	mpz_clear(tby0);
			return;
		}
	}
	
	long np = get_bit_length(ax, ay);
	long mp = 2*m - np + 1;
	
	flag = 0;
	mpz_set_ui(a21, 1);	mpz_set_ui(a22, 0);	mpz_set_ui(b21, 0);	mpz_set_ui(b22, 1);
	v_length = get_double_length(bx, by);
	if((v_length < (double)(np - (double)mp/2 + c)) || (mp <= 8)){
		mpz_set_ui(a21, 1);
		mpz_set_ui(a22, 0);
		mpz_set_ui(b21, 0);
		mpz_set_ui(b22, 1);
		flag = 1;
	}
	if(flag == 0){
		if(mpz_sgn(ax) >= 0){
			mpz_cdiv_q_2exp(ax0, ax, mp);
		}else{
			mpz_fdiv_q_2exp(ax0, ax, mp);
		}
		if(mpz_sgn(ay) >= 0){
			mpz_cdiv_q_2exp(ay0, ay, mp);
		}else{
			mpz_fdiv_q_2exp(ay0, ay, mp);
		}
		if(mpz_sgn(bx) >= 0){
			mpz_fdiv_q_2exp(bx0, bx, mp);
		}else{
			mpz_cdiv_q_2exp(bx0, bx, mp);
		}
		if(mpz_sgn(by) >= 0){
			mpz_fdiv_q_2exp(by0, by, mp);
		}else{
			mpz_cdiv_q_2exp(by0, by, mp);
		}
		mpz_mul(tmp1, ax0, ax0);	mpz_mul(tmp2, ay0, ay0);	mpz_add(tmp1, tmp1, tmp2);
		mpz_mul(tmp3, bx0, bx0);	mpz_mul(tmp4, by0, by0);	mpz_add(tmp3, tmp3, tmp4);
		if(mpz_cmp(tmp1, tmp3) == 0){
			direct_compute(ax0, ay0, bx0, by0, a21, a22, b21, b22);
			flag = 1;
		}

		if(flag == 0){
			mpz_set(tax0, ax0);	mpz_set(tay0, ay0);
			mpz_set(tbx0, bx0);	mpz_set(tby0, by0);
			
			half_gaussian(ax0, ay0, bx0, by0, ta1, ta2, tb1, tb2);

			fixup(ax, ay, bx, by, tax0, tay0, tbx0, tby0, ta1, ta2, tb1, tb2, a21, a22, b21, b22, mp);
		}
	}
	
	det = detdim2mat(a21, a22, b21, b22);
	mpz_mul(tmp1, ax, b22);	mpz_mul(tmp2, bx, b21);	mpz_sub(ax0, tmp1, tmp2);
	mpz_mul(tmp1, ay, b22);	mpz_mul(tmp2, by, b21);	mpz_sub(ay0, tmp1, tmp2);
	mpz_mul(tmp1, bx, a21);	mpz_mul(tmp2, ax, a22);	mpz_sub(bx0, tmp1, tmp2);
	mpz_mul(tmp1, by, a21);	mpz_mul(tmp2, ay, a22);	mpz_sub(by0, tmp1, tmp2);
	
	if(det == -1){
		mpz_neg(ax0, ax0);
		mpz_neg(ay0, ay0);
		mpz_neg(bx0, bx0);
		mpz_neg(by0, by0);
	}
	
	mpz_set(ax, ax0);	mpz_set(ay, ay0);
	mpz_set(bx, bx0);	mpz_set(by, by0);
	
	mpz_mul(tmp1, a21, a11);	mpz_mul(tmp2, a22, b11);	mpz_add(a1, tmp1, tmp2);
	mpz_mul(tmp1, a21, a12);	mpz_mul(tmp2, a22, b12);	mpz_add(a2, tmp1, tmp2);
	mpz_mul(tmp1, b21, a11);	mpz_mul(tmp2, b22, b11);	mpz_add(b1, tmp1, tmp2);
	mpz_mul(tmp1, b21, a12);	mpz_mul(tmp2, b22, b12);	mpz_add(b2, tmp1, tmp2);
	
	
	while(get_double_length(bx, by) > n/2+c){
		
		mpz_mul(tmp1, ax, bx);	mpz_mul(tmp2, ay, by);	mpz_add(tmp1, tmp1, tmp2);
		mpz_mul(tmp3, bx, bx);	mpz_mul(tmp4, by, by);	mpz_add(tmp3, tmp3, tmp4);
		mpz_fdiv_q(q, tmp1, tmp3);
		mpz_submul(ax, q, bx);
		mpz_submul(ay, q, by);
		mpz_swap(ax, bx);
		mpz_swap(ay, by);

		mpz_set(q0, q1);
		mpz_set(q1, q2);
		mpz_set(q2, q);

		mpz_swap(a1, b1);
		mpz_swap(a2, b2);
		mpz_addmul(a1, q, b1); 
		mpz_addmul(a2, q, b2);
		
		if((is_terminal(ax, ay, bx, by) == true) || (is_straddles(ax, ay, bx, by, (double)((double)n/2+c)) == true)){
			mpz_clear(ta1);	mpz_clear(ta2);	mpz_clear(tb1);	mpz_clear(tb2);
			mpz_clear(ax0);	mpz_clear(ay0);	mpz_clear(bx0);	mpz_clear(by0);
			mpz_clear(a11);	mpz_clear(a12);	mpz_clear(b11);	mpz_clear(b12);
			mpz_clear(a21);	mpz_clear(a22);	mpz_clear(b21);	mpz_clear(b22);
			mpz_clear(tmp1);	mpz_clear(tmp2);
			mpz_clear(tmp3);	mpz_clear(tmp4);
			mpz_clear(q);
			mpz_clear(tax0);	mpz_clear(tay0);
			mpz_clear(tbx0);	mpz_clear(tby0);
			return;
		}
	}
	
	
	mpz_clear(ta1);	mpz_clear(ta2);	mpz_clear(tb1);	mpz_clear(tb2);
	mpz_clear(ax0);	mpz_clear(ay0);	mpz_clear(bx0);	mpz_clear(by0);
	mpz_clear(a11);	mpz_clear(a12);	mpz_clear(b11);	mpz_clear(b12);
	mpz_clear(a21);	mpz_clear(a22);	mpz_clear(b21);	mpz_clear(b22);
	mpz_clear(tmp1);	mpz_clear(tmp2);
	mpz_clear(tmp3);	mpz_clear(tmp4);
	mpz_clear(q);
	mpz_clear(tax0);	mpz_clear(tay0);
	mpz_clear(tbx0);	mpz_clear(tby0);

	return;
}

