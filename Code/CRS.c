#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>

int main(int argc, char *argv[]){
	mpz_t a[2], b[2], c[2];
	mpz_t tmp1, tmp2, tmp3, tmp4, q, zero, tmp5, tmp6;
	
	mpz_init(a[0]); mpz_init(a[1]);mpz_init(b[0]); mpz_init(b[1]);
	mpz_init(c[0]);	mpz_init(c[1]);
	mpz_init(tmp1);	mpz_init(tmp2);	mpz_init(tmp3);	mpz_init(tmp4);
	mpz_init(q);	mpz_init(zero);
	mpz_init(tmp5);	mpz_init(tmp6);
	
	clock_t start, finish;
	
	printf("\nCRS Algorithm.\n");
	
	FILE *fp = fopen("Input.txt", "r");
	gmp_fscanf(fp, "%Zd %Zd %Zd %Zd", a[0], a[1], b[0], b[1]);
	fclose(fp);
	
	if(mpz_cmp(a[0], b[0]) < 0){
		mpz_swap(a[0], b[0]);
	}
	//mpz_set_si(a[1], 0);
    	
    mpz_set_ui(zero, 0);
	mpz_mul(tmp1, a[0], b[0]);
	mpz_mul(tmp2, a[1], b[1]);
	mpz_add(tmp1, tmp1, tmp2);
	
	mpz_mul(tmp3, b[0], b[0]);
	mpz_mul(tmp4, b[1], b[1]);
	mpz_add(tmp3, tmp3, tmp4);
	if(mpz_cmpabs(tmp1, tmp3) < 0){
		mpz_swap(a[0], b[0]);
		mpz_swap(a[1], b[1]);
	}
    	
	int euctime = 0;
	int cnt = 0;
	
	start = clock();
	
	while(1){
		cnt = cnt + 1;

		mpz_mul(tmp1, a[0], b[0]);
		mpz_mul(tmp2, a[1], b[1]);
		mpz_add(tmp1, tmp1, tmp2);
		
		mpz_mul(tmp3, b[0], b[0]);
		mpz_mul(tmp4, b[1], b[1]);
		mpz_add(tmp3, tmp3, tmp4);
		
		if(mpz_sgn(tmp1) == 0 || mpz_sgn(tmp3) == 0){
			break;
		}
		
		mpz_fdiv_q(q, tmp1, tmp3);
		
		if(mpz_sgn(q) == 0){
			break;
		}
		
		mpz_mul(tmp1, q, b[0]);
		mpz_sub(a[0], a[0], tmp1);
		mpz_mul(tmp1, q, b[1]);
		mpz_sub(a[1], a[1], tmp1);
		mpz_swap(a[0], b[0]);
		mpz_swap(a[1], b[1]);
	}
	
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
	
		
	printf("Partial Euclidean cycle times = %d\n",cnt);
	printf("Prog time: %lf\n", (double)(finish-start)/CLOCKS_PER_SEC);
	
	mpz_clear(a[0]); mpz_clear(a[1]); mpz_clear(b[0]); mpz_clear(b[1]);
	mpz_clear(c[0]); mpz_clear(c[1]);
	mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3); mpz_clear(tmp4);
	mpz_clear(q);
	mpz_clear(tmp5);	mpz_clear(tmp6);
    
	return 0;
}

