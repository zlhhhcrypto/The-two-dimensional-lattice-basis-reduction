#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>

void findanswer(mpz_t min1, mpz_t min2, mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by); 
int ParEuc(mpz_t ansp, mpz_t ansq, mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);
 
int main(int argc, char *argv[]){
	mpz_t a[2], b[2];
	mpz_t p, q, min1, min2;
	
	mpz_init(a[0]); mpz_init(a[1]);mpz_init(b[0]); mpz_init(b[1]);
	mpz_init(p); mpz_init(q);
	mpz_init(min1);	mpz_init(min2);
	
	clock_t start, finish;
	
	printf("\nPartial Euclidean algorithm for all kinds of lattices.\n");
	
	FILE *fp = fopen("Input.txt", "r");
	gmp_fscanf(fp, "%Zd %Zd %Zd %Zd", a[0], a[1], b[0], b[1]);
	fclose(fp);
	
	int euctime = 0;
	
	start = clock();
	
	euctime = ParEuc(p, q, a[0], a[1], b[0], b[1]);	

    finish = clock();	
	
	findanswer(min1, min2, a[0], a[1], b[0], b[1]);
	
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
	mpz_clear(p); mpz_clear(q);
    return 0;
}

void findanswer(mpz_t min1, mpz_t min2, mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by){
	mpz_t tax, tay, tbx, tby;
	
	mpz_init(tax);	mpz_init(tay);
	mpz_init(tbx);	mpz_init(tby);
	
	mpz_abs(tax, ax);
	mpz_abs(tay, ay);
	mpz_abs(tbx, bx);
	mpz_abs(tby, by);
	
	mpz_set(min1, tax);
	if(mpz_cmp(tax, tay) < 0){
		mpz_set(min1, tay);
	}
	mpz_set(min2, tbx);
	if(mpz_cmp(tbx, tby) < 0){
		mpz_set(min2, tby);
	}
}

int ParEuc(mpz_t ansp, mpz_t ansq, mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by){
	
	int sgnax, sgnay, sgnbx, sgnby;
	int recnt = 0;
	
	sgnax = mpz_sgn(ax);
	sgnay = mpz_sgn(ay);
	sgnbx = mpz_sgn(bx);
	sgnby = mpz_sgn(by);
	
	int flag = 0;

	if(sgnax * sgnay * sgnbx * sgnby <= 0){
		flag = 1;
		if(mpz_cmpabs(ax, ay) >= 0 && mpz_cmpabs(bx, by) <= 0){
			return recnt;
		}
		if(mpz_cmpabs(ax, ay) <= 0 && mpz_cmpabs(bx, by) >= 0){
			return recnt;
		}
	}

	int next[3] = {1,2,0};
	int pre = 0, cur = 1, net = 2, index = 0;
	int pre_p = 0, cur_p = 1, net_p = 2;
	
	mpz_t R[3], P[3], q, r, tmp, zero, one, tt;

	mpz_init_set(R[pre], ax);
	mpz_init_set(R[cur], bx);
	mpz_init(R[net]);
	
	mpz_init_set(P[pre_p], ay);
	mpz_init_set(P[cur_p], by);
	mpz_init(P[net]);

	mpz_init(q);
	mpz_init(r);
	mpz_init(tmp);
	mpz_init(zero);
	mpz_init(one);
	mpz_init(tt);
	
	mpz_set_si(zero, 0);
	mpz_set_si(one, 1);
	
	if(flag == 0 && mpz_sgn(R[pre]) != mpz_sgn(R[cur])){
		mpz_neg(R[cur], R[cur]);
		mpz_neg(P[cur_p], P[cur_p]);
	}
	
	while(1){
		if(mpz_cmpabs(R[pre], P[pre_p]) >= 0 && mpz_sgn(R[cur]) != 0){
			recnt++;
			mpz_tdiv_q(q, R[pre], R[cur]);
			mpz_submul(R[pre], q, R[cur]);
			mpz_submul(P[pre_p], q, P[cur_p]);
			if(flag == 0 && mpz_sgn(P[cur_p]) != mpz_sgn(P[pre_p])){
				flag = 1;
			}
			
			if(flag == 0 && mpz_cmpabs(P[cur_p], P[pre_p]) <= 0){
				mpz_sub(R[pre], R[pre], R[cur]);
				mpz_sub(P[pre_p], P[pre_p], P[cur_p]);
				flag = 1;
			}
			mpz_swap(R[pre], R[cur]);
			mpz_swap(P[pre_p], P[cur_p]);
		}else{
			recnt++;
			mpz_tdiv_q(q, P[pre_p], P[cur_p]);
			mpz_submul(P[pre_p], q, P[cur_p]);
			mpz_submul(R[pre], q, R[cur]);
			if(flag == 0 && mpz_sgn(R[cur]) != mpz_sgn(R[pre])){
				flag = 1;
			}
			
			if(flag == 0 && mpz_cmpabs(R[cur], R[pre]) <= 0){
				mpz_sub(R[pre], R[pre], R[cur]);
				mpz_sub(P[pre_p], P[pre_p], P[cur_p]);
				flag = 1;
			}
			mpz_swap(R[pre], R[cur]);
			mpz_swap(P[pre_p], P[cur_p]);
		}
		if(flag == 1){
			if(mpz_cmpabs(R[pre], P[pre_p]) >= 0 && mpz_cmpabs(R[cur], P[cur_p]) <= 0){
				break;
			}
			if(mpz_cmpabs(R[pre], P[pre_p]) <= 0 && mpz_cmpabs(R[cur], P[cur_p]) >= 0){
				break;
			}
		}
	}
	
	sgnax = mpz_sgn(R[pre]);
	sgnay = mpz_sgn(P[pre_p]);
	sgnbx = mpz_sgn(R[cur]);
	sgnby = mpz_sgn(P[cur_p]);
	if(sgnax * sgnbx == 0){
		if(sgnay != sgnby){
			mpz_neg(R[cur], R[cur]);
			mpz_neg(P[cur_p], P[cur_p]);
		}
	}
	if(sgnay * sgnby == 0){
		if(sgnax != sgnbx){
			mpz_neg(R[cur], R[cur]);
			mpz_neg(P[cur_p], P[cur_p]);
		}
	}
	
	mpz_t TP1[2], TP2[2], divd, disr, mina, minb;
	mpz_init(TP1[0]); mpz_init(TP1[1]); mpz_init(TP2[0]); mpz_init(TP2[1]);
	mpz_init(divd); mpz_init(disr);
	mpz_init(mina);	mpz_init(minb);
	
	mpz_set(mina, R[pre]);
	if(mpz_cmpabs(R[pre], P[pre_p]) < 0){
		mpz_set(mina, P[pre_p]);
	}
	mpz_set(minb, R[cur]);
	if(mpz_cmpabs(R[cur], P[cur_p]) < 0){
		mpz_set(minb, P[cur_p]);
	}
	
	if(mpz_cmpabs(mina, minb) > 0){
		mpz_swap(R[pre], R[cur]);
		mpz_swap(P[pre_p], P[cur_p]);
	}
	
	mpz_set(ax, R[pre]); 
	mpz_set(ay, P[pre_p]);

	if((mpz_cmp(ax, zero) == 0) && mpz_cmp(ay, zero) == 0){
		return recnt;
	}	

	sgnay = mpz_sgn(P[pre_p]);
	sgnby = mpz_sgn(P[cur_p]);
	
	if(sgnay * sgnby <= 0){
		mpz_abs(divd, P[cur_p]);
		mpz_abs(tmp, R[cur]);
		mpz_sub(divd, tmp, divd);
	}else{
		mpz_abs(divd, P[cur_p]);
		mpz_abs(tmp, R[cur]);
		mpz_sub(divd, divd, tmp);
	}
	
	mpz_abs(disr, P[pre_p]);
	mpz_abs(tmp, R[pre]);
	mpz_add(disr, disr, tmp);	

	mpz_tdiv_q(q, divd, disr);

	mpz_set(TP1[0], R[cur]);
	mpz_set(TP1[1], P[cur_p]);

	mpz_submul(TP1[0], q, R[pre]);
	mpz_submul(TP1[1], q, P[pre_p]);
	if(mpz_sgn(q) >= 0){
		mpz_sub(TP2[0], TP1[0], R[pre]);
		mpz_sub(TP2[1], TP1[1], P[pre_p]);
	}else{
		mpz_add(TP2[0], TP1[0], R[pre]);
		mpz_add(TP2[1], TP1[1], P[pre_p]);
	}

	int tp1, tp2;
	if(mpz_cmpabs(TP1[0], TP1[1]) >= 0) tp1 = 0;
	else tp1 = 1;
	if(mpz_cmpabs(TP2[0], TP2[1]) >= 0) tp2 = 0;
	else tp2 = 1;

	if(mpz_cmpabs(TP1[tp1], TP2[tp2]) > 0){
		mpz_set(bx, TP2[0]);
		mpz_set(by, TP2[1]);
	}else{
		mpz_set(bx, TP1[0]);
		mpz_set(by, TP1[1]);
	}

	mpz_clear(TP1[0]); mpz_clear(TP1[1]); mpz_clear(TP2[0]); mpz_clear(TP2[1]);
	mpz_clear(divd); mpz_clear(disr);	mpz_clear(tt);
	
	mpz_clear(R[0]); mpz_clear(R[1]); mpz_clear(R[2]);
	mpz_clear(P[0]); mpz_clear(P[1]);
	mpz_clear(q);	 mpz_clear(r);
	return recnt;
}
