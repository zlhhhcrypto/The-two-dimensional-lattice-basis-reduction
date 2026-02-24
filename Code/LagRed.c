#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

mpz_t a[1001], b[1001], c[1001], d[1001];
int dim;
int is_finish = 0;

void findanswer(){
	mpz_t first, second, temp, zero;
	mpf_t mid, ftmp1, ftmp2, fq;
	
	mpz_init(temp);
	mpz_init(zero);
	mpz_set_ui(zero, 0);
	mpz_init(first);
	mpz_set_ui(first, 0);
	mpz_init(second);
	mpz_set_ui(second, 0);
	
	mpf_init(mid);	mpf_init(ftmp1);	
	mpf_init(ftmp2);mpf_init(fq);
	
	mpf_set_d(mid, 0.5);
	
	for(int i = 0; i < dim; i++){
		mpz_init(c[i]);
		mpz_init(d[i]);
	}

	for(int i = 0; i < dim; i++){
		mpz_mul(temp, a[i], a[i]);
		mpz_add(first, first, temp);
	}
	for(int i = 0; i < dim; i++){
		mpz_mul(temp, b[i], b[i]);
		mpz_add(second, second, temp);
	}
	if(mpz_cmp(first, second) > 0){
		for(int i = 0; i < dim; i++){
			mpz_swap(a[i], b[i]);
		}
	}
	mpz_set(first, zero);
	mpz_set(second, zero);
	for(int i = 0; i < dim; i++){
		mpz_sub(c[i], a[i], b[i]);
	}
	for(int i = 0; i < dim; i++){
		mpz_add(d[i], a[i], b[i]);
	}
	for(int i = 0; i < dim; i++){
		mpz_mul(temp, c[i], c[i]);
		mpz_add(first, first, temp);
	}
	for(int i = 0; i < dim; i++){
		mpz_mul(temp, d[i], d[i]);
		mpz_add(second, second, temp);
	}
	if(mpz_cmp(first, second) > 0){
		for(int i = 0; i < dim; i++){
			mpz_neg(b[i], b[i]);
		}
	}
	mpz_set(first, zero);
	mpz_set(second, zero);
	for(int i = 0; i < dim; i++){
		mpz_sub(c[i], a[i], b[i]);
	}
	for(int i = 0; i < dim; i++){
		mpz_mul(temp, c[i], c[i]);
		mpz_add(first, first, temp);
	}
	for(int i = 0; i < dim; i++){
		mpz_mul(temp, b[i], b[i]);
		mpz_add(second, second, temp);
	}
	if(mpz_cmp(first, second) >= 0){
		return;
	}
	
	int isloop = 0;
	mpz_set(first, zero);
	mpz_set(second, zero);
	for(int i = 0; i < dim; i++){
		mpz_sub(c[i], a[i], b[i]);
	}
	for(int i = 0; i < dim; i++){
		mpz_mul(temp, c[i], c[i]);
		mpz_add(first, first, temp);
	}
	for(int i = 0; i < dim; i++){
		mpz_mul(temp, a[i], a[i]);
		mpz_add(second, second, temp);
	}
	if(mpz_cmp(first, second) >= 0){
		isloop = 1;
	}
	
	if(isloop == 0){
		mpz_set(first, zero);
		mpz_set(second, zero);
		for(int i = 0; i < dim; i++){
			mpz_mul(temp, a[i], a[i]);
			mpz_add(first, first, temp);
		}
		for(int i = 0; i < dim; i++){
			mpz_mul(temp, b[i], b[i]);
			mpz_add(second, second, temp);
		}
		if(mpz_cmp(first, second) == 0){
			for(int i = 0; i < dim; i++){
				mpz_sub(b[i], a[i], b[i]);
				mpz_swap(a[i], b[i]);
			}
			return;
		}
		for(int i = 0; i < dim; i++){
			mpz_sub(a[i], b[i], a[i]);
		}
	}
	mpz_t u1, u2, one, curfirst, cursecond;
	
	mpz_init(u1);	
	mpz_init(u2);
	mpz_init(curfirst);
	mpz_init(cursecond);
	mpz_init(one);
	mpz_set_ui(one, 1);

	while(1){
		mpz_set(first, zero);
		mpz_set(second, zero);
		for(int i = 0; i < dim; i++){
			mpz_mul(temp, a[i], b[i]);
			mpz_add(first, first, temp);
		}
		for(int i = 0; i < dim; i++){
			mpz_mul(temp, a[i], a[i]);
			mpz_add(second, second, temp);
		}
		
		if(mpz_sgn(first) == 0 || mpz_sgn(second) == 0){
			return;
		}
			
		mpf_set_z(ftmp1, first);
		mpf_set_z(ftmp2, second);
		mpf_div(fq, ftmp1, ftmp2);
		
		if(mpf_sgn(fq) >= 0){
			mpf_add(fq, fq, mid);
		}else{
			mpf_sub(fq, fq, mid);
		}
		mpz_set_f(u1, fq);
		
		if(mpz_sgn(u1) == 0){
			return;
		}
		
		for(int i = 0; i < dim; i++){
			mpz_mul(temp, u1, a[i]);
			mpz_sub(b[i], b[i], temp);
		}
		mpz_set(first, zero);
		mpz_set(second, zero);
		for(int i = 0; i < dim; i++){
			mpz_sub(c[i], a[i], b[i]);
		}
		for(int i = 0; i < dim; i++){
			mpz_add(d[i], a[i], b[i]);
		}
		for(int i = 0; i < dim; i++){
			mpz_abs(c[i], c[i]);
			mpz_add(first, first, c[i]);
		}
		for(int i = 0; i < dim; i++){
			mpz_abs(d[i], d[i]);
			mpz_add(second, second, d[i]);
		}
		if(mpz_cmp(first, second) > 0){
			for(int i = 0; i < dim; i++){
				mpz_neg(b[i], b[i]);
			}
		}
		for(int i = 0; i < dim; i++){
			mpz_swap(a[i], b[i]);
		}
	}
}

int main(){
	
	printf("\nLag Algorithm!\n");
	
	clock_t start, finish;
	dim = 2;
	
	for(int i = 0; i < dim; i++){
		mpz_init(a[i]);
		mpz_init(b[i]);
	}
	
	FILE *fp = fopen("Input.txt", "r");
	for(int i = 0; i < 2*dim; i++){
		if(i < dim){
			gmp_fscanf(fp, "%Zd", a[i]);
		}else{
			gmp_fscanf(fp, "%Zd", b[i-dim]);
		}
	}
	fclose(fp);	
	
	if(mpz_cmp(a[0], b[0]) < 0){
		mpz_swap(a[0], b[0]);
	}
	mpz_set_si(a[1], 0);
	start = clock();
	findanswer();
	finish = clock();
	
	FILE *fo;
	fo = fopen("Output.txt", "w");
	for(int i = 0; i < 2*dim; i++){
		if(i < dim){
			gmp_fprintf(fo, "%Zd\n",a[i]);
		}else{
			gmp_fprintf(fo, "%Zd\n",b[i-dim]);
		}
	}
	fclose(fo);
	
	printf("Lag time: %lf\n", (double)(finish-start)/CLOCKS_PER_SEC);
	
	for(int i = 0; i < dim; i++){
		mpz_clear(a[i]);
		mpz_clear(b[i]);
	}
	
	return 0;
}
