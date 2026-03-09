#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>

void getSnumber(mpz_t p, mpz_t nv, long *s);
 
int main(int argc, char *argv[]){
	clock_t start, finish;  //time
	
	mpz_t ax, ay, bx, by;
	mpz_t nu, nv, p, tmp, temp1, temp2, two, zero;
	long s, t1s, t2s;
	
	mpz_init(ax);	mpz_init(ay);	mpz_init(bx); mpz_init(by);
	mpz_init(nu);	mpz_init(nv);	mpz_init(p);
	mpz_init(tmp);	mpz_init(temp1);	mpz_init(temp2);
	mpz_init(two);	mpz_init(zero);
	
	mpz_set_si(two, 2);	mpz_set_si(zero, 0);

	FILE *fp = fopen("Input.txt", "r");
	gmp_fscanf(fp, "%Zd %Zd %Zd %Zd", ax, ay, bx, by);
	fclose(fp);
	
	start = clock();
	
	mpz_mul(temp1, ax, ax);
	mpz_mul(temp2, ay, ay);
	mpz_add(nu, temp1, temp2);
	
	mpz_mul(temp1, bx, bx);
	mpz_mul(temp2, by, by);
	mpz_add(nv, temp1, temp2);
	
	mpz_mul(temp1, ax, bx);
	mpz_mul(temp2, ay, by);
	mpz_add(p, temp1, temp2);
	
	while(1){
		if(mpz_cmp(nu, nv) < 0){
			mpz_swap(ax, bx);
			mpz_swap(ay, by);
			mpz_swap(nu, nv);
		}
		mpz_abs(tmp, p);
		mpz_mul(tmp, tmp, two);
		if(mpz_cmp(tmp, nv) <= 0){
			mpz_swap(ax, bx);
			mpz_swap(ay, by);
			break;
		}
		getSnumber(p, nv, &s);
		if(mpz_cmp(p, zero) > 0){
			mpz_mul_2exp(temp1, bx, s);
			mpz_mul_2exp(temp2, by, s);
			mpz_sub(ax, ax, temp1);
			mpz_sub(ay, ay, temp2);
			t1s = s*2;
			t2s = s+1;
			mpz_mul_2exp(temp1, nv, t1s);
			mpz_mul_2exp(temp2, p, t2s);
			mpz_add(nu, nu, temp1);
			mpz_sub(nu, nu, temp2);
			mpz_mul_2exp(tmp, nv, s);
			mpz_sub(p, p, tmp);
		}else{
			mpz_mul_2exp(temp1, bx, s);
			mpz_mul_2exp(temp2, by, s);
			mpz_add(ax, ax, temp1);
			mpz_add(ay, ay, temp2);
			t1s = s*2;
			t2s = s+1;
			mpz_mul_2exp(temp1, nv, t1s);
			mpz_mul_2exp(temp2, p, t2s);
			mpz_add(nu, nu, temp1);
			mpz_add(nu, nu, temp2);
			mpz_mul_2exp(tmp, nv, s);
			mpz_add(p, p, tmp);
		}
	}
    finish = clock();	
	
	FILE *fo;
	fo = fopen("Output.txt", "w");
	gmp_fprintf(fo, "%Zd\n",ax);
	gmp_fprintf(fo, "%Zd\n",ay);
	gmp_fprintf(fo, "%Zd\n",bx);
	gmp_fprintf(fo, "%Zd\n",by);
	fclose(fo);
	
	printf("\nProg time: %lf\n", (double)(finish-start)/CLOCKS_PER_SEC);
	
	mpz_clear(ax); mpz_clear(ay); mpz_clear(bx); mpz_clear(by);
	mpz_clear(p); mpz_clear(nu);	mpz_clear(nv);	mpz_clear(tmp);
	mpz_clear(temp1);	mpz_clear(temp2);	mpz_clear(two);
	mpz_clear(zero);
    return 0;
}

void getSnumber(mpz_t p, mpz_t nv, long *s){
	mpz_t tmp;
	mpz_init(tmp);
	long p_len, nv_len, result;
	//get p_len
	if(mpz_sgn(p) == 0){
		p_len = 0;
	}else if(mpz_sgn(p) > 0){
		p_len = mpz_sizeinbase(p, 2);
	}else{
		mpz_neg(tmp, p);
		mpz_sub_ui(tmp, tmp, 1);
		if(mpz_sgn(tmp) == 0){
			p_len = 0;
		}else{
			p_len = mpz_sizeinbase(tmp, 2);
		}
	}
	//get nv_len
	if(mpz_sgn(nv) == 0){
		nv_len = 0;
	}else if(mpz_sgn(nv) > 0){
		nv_len = mpz_sizeinbase(nv, 2);
	}else{
		mpz_neg(tmp, nv);
		mpz_sub_ui(tmp, tmp, 1);
		if(mpz_sgn(tmp) == 0){
			nv_len = 0;
		}else{
			nv_len = mpz_sizeinbase(tmp, 2);
		}
	}
	
	result = p_len-nv_len;
	if(result > 0){
		(*s) = result;
	}else{
		(*s) = 0;
	}
	mpz_clear(tmp);
}
