#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>

int detdim2mat(mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2);
void updateAB(mpz_t ax0, mpz_t ay0, mpz_t bx0, mpz_t by0, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2);
void bitlen(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, long * maxlen, long * minlen);
void integerbitlen(mpz_t ax, mpz_t bx, long * minlen);
void bitsublen(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, long * sublen);
void bitaddlen(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, long * addlen);
void bitsubqlen(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t q, long * subqlen);
void bitaddqlen(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t q, long * addqlen);
void hgcd(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2);
void HGCD(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);
void FinalRevise(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);
void preprocess(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);
void findanswer(mpz_t min1, mpz_t min2, mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);
 
int main(int argc, char *argv[]){
	mpz_t a[2], b[2];
	mpz_t p, q, min1, min2, tq, maxa, maxb;
	
	mpz_init(a[0]); mpz_init(a[1]);mpz_init(b[0]); mpz_init(b[1]);
	mpz_init(p); mpz_init(q);
	mpz_init(min1);	mpz_init(min2);
	mpz_init(tq);
	mpz_init(maxa);	mpz_init(maxb);
	
	mpz_t presublen, newsublen;
    mpz_init(presublen);	mpz_init(newsublen);
	
	clock_t start, finish;
	int len = 0;
	
	printf("\nHGCD-Partial Euclidean algorithm for all kinds of lattices.\n");
	
	FILE *fp = fopen("Input.txt", "r");
	gmp_fscanf(fp, "%Zd %Zd %Zd %Zd", a[0], a[1], b[0], b[1]);
	fclose(fp);
	
	int sgnax, sgnay, sgnbx, sgnby;
	long maxlen, minlen, addlen, sublen;
	long N, S;
	int euctime = 0;
	int cnt = 0;
	
	start = clock();

    while(1){
    	HGCD(a[0], a[1], b[0], b[1]);
    	cnt++;
    	
    	sgnax = mpz_sgn(a[0]);
		sgnay = mpz_sgn(a[1]);
		sgnbx = mpz_sgn(b[0]);
		sgnby = mpz_sgn(b[1]);
		
		if(sgnax * sgnay * sgnbx * sgnby <= 0){
			if(mpz_cmpabs(a[0], a[1]) >= 0 && mpz_cmpabs(b[0], b[1]) <= 0){
				break;
			}
			if(mpz_cmpabs(a[0], a[1]) <= 0 && mpz_cmpabs(b[0], b[1]) >= 0){
				break;
			}
		}
		
		mpz_abs(maxa, a[0]);
		if(mpz_cmpabs(a[0], a[1]) < 0){
			mpz_abs(maxa, a[1]);
		}
		mpz_abs(maxb, b[0]);
		if(mpz_cmpabs(b[0], b[1]) < 0){
			mpz_abs(maxb, b[1]);
		}
		if(mpz_cmpabs(maxa, maxb) < 0){
			mpz_swap(a[0], b[0]);
			mpz_swap(a[1], b[1]);
		}
    		
    		bitlen(a[0], a[1], b[0], b[1], &maxlen, &minlen);
		N = maxlen;
		S = (N / 2) + 1;
		
		if(minlen <= S){
			if(mpz_cmpabs(a[0], a[1]) >= 0){
				mpz_tdiv_q(q, a[0], b[0]);
				
			}else{
				mpz_tdiv_q(q, a[1], b[1]);
			}
			mpz_submul(a[0], q, b[0]);
			mpz_submul(a[1], q, b[1]);
			mpz_swap(a[0], b[0]);
			mpz_swap(a[1], b[1]);
			continue;
		}
		
		bitaddlen(a[0], a[1], b[0], b[1], &addlen);
		bitsublen(a[0], a[1], b[0], b[1], &sublen);
		if(addlen <= S){
			mpz_add(a[0], a[0], b[0]);
			mpz_add(a[1], a[1], b[1]);
			mpz_swap(a[0], b[0]);
			mpz_swap(a[1], b[1]);
		}
		if(sublen <= S){
			mpz_sub(a[0], a[0], b[0]);
			mpz_sub(a[1], a[1], b[1]);
			mpz_swap(a[0], b[0]);
			mpz_swap(a[1], b[1]);
		}
    }
	
	sgnax = mpz_sgn(a[0]);
	sgnay = mpz_sgn(a[1]);
	sgnbx = mpz_sgn(b[0]);
	sgnby = mpz_sgn(b[1]);
	if(sgnax * sgnbx == 0){
		if(sgnay != sgnby){
			mpz_neg(b[0], b[0]);
			mpz_neg(b[1], b[1]);
		}
	}
	if(sgnay * sgnby == 0){
		if(sgnax != sgnbx){
			mpz_neg(b[0], b[0]);
			mpz_neg(b[1], b[1]);
		}
	}
	
	FinalRevise(a[0], a[1], b[0], b[1]);
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
	printf("HGCD cycle times = %d\n",cnt);
	printf("Prog time: %lf\n", (double)(finish-start)/CLOCKS_PER_SEC);
	
	mpz_clear(a[0]); mpz_clear(a[1]); mpz_clear(b[0]); mpz_clear(b[1]);
	mpz_clear(p); mpz_clear(q);	mpz_clear(presublen);	mpz_clear(newsublen);
	mpz_clear(tq);	mpz_clear(maxa);	mpz_clear(maxb);
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

void HGCD(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by){
	mpz_t a1, a2, b1, b2;
	mpz_init(a1);	mpz_init(a2);	mpz_init(b1);	mpz_init(b2);
	
	mpz_set_si(a1, 1);	mpz_set_si(a2, 0);	mpz_set_si(b1, 0);	mpz_set_si(b2, 1);
    hgcd(ax, ay, bx, by, a1, a2, b1, b2);	
	
	mpz_clear(a1);	mpz_clear(a2);	mpz_clear(b1);	mpz_clear(b2);
	return;
}

void FinalRevise(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by){

	int pre = 0, cur = 1;
	int pre_p = 0, cur_p = 1;
	
	mpz_t R[2], P[2], q, tmp, zero;

	mpz_init_set(R[pre], ax);
	mpz_init_set(R[cur], bx);
	mpz_init_set(P[pre_p], ay);
	mpz_init_set(P[cur_p], by);

	mpz_init(q);
	mpz_init(tmp);
	mpz_init(zero);
	
	mpz_set_si(zero, 0);
	
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
		return;
	}
	int sgnay = mpz_sgn(P[pre_p]);
	int sgnby = mpz_sgn(P[cur_p]);
	
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
	mpz_clear(divd); mpz_clear(disr);	mpz_clear(tmp);
	
	mpz_clear(R[0]); mpz_clear(R[1]);
	mpz_clear(P[0]); mpz_clear(P[1]);
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

void updateAB(mpz_t ax0, mpz_t ay0, mpz_t bx0, mpz_t by0, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2){
	mpz_t atx, aty, btx, bty;
	mpz_init(atx);	mpz_init(aty);	mpz_init(btx);	mpz_init(bty);
	
	mpz_mul(atx, ax0, b2);	mpz_submul(atx, bx0, b1);
	mpz_mul(aty, ay0, b2);	mpz_submul(aty, by0, b1);
	mpz_mul(btx, bx0, a1);	mpz_submul(btx, ax0, a2);
	mpz_mul(bty, by0, a1);	mpz_submul(bty, ay0, a2);
	
	if(detdim2mat(a1, a2, b1, b2) < 0){
		mpz_neg(atx, atx);
		mpz_neg(aty, aty);
		mpz_neg(btx, btx);
		mpz_neg(bty, bty);
	}
	
	mpz_set(ax0, atx);	mpz_set(ay0, aty);
	mpz_set(bx0, btx);	mpz_set(by0, bty);
	
	mpz_clear(atx);	mpz_clear(aty);	mpz_clear(btx);	mpz_clear(bty);
}

void bitlen(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, long * maxlen, long * minlen){
	long l[4];
	l[0] = mpz_sizeinbase(ax, 2);
	l[1] = mpz_sizeinbase(ay, 2);
	l[2] = mpz_sizeinbase(bx, 2);
	l[3] = mpz_sizeinbase(by, 2);
	
	(*maxlen) = l[0];
	(*minlen) = l[0];
	for(int i = 1; i < 4; i++){
		if(l[i] > (*maxlen)){
			(*maxlen) = l[i];
		}
		if(l[i] < (*minlen)){
			(*minlen) = l[i];
		}
	}
}

void integerbitlen(mpz_t ax, mpz_t bx, long * minlen){
	long l[2];
	l[0] = mpz_sizeinbase(ax, 2);
	l[1] = mpz_sizeinbase(bx, 2);
	
	(*minlen) = l[0];
	if(l[1] < l[0]){
		(*minlen = l[1]);
	}
}

void bitsublen(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, long * sublen){
	mpz_t sub1, sub2;
	mpz_init(sub1);	mpz_init(sub2);
	
	mpz_sub(sub1, ax, bx);
	mpz_sub(sub2, ay, by);
	
	long s1 = mpz_sizeinbase(sub1, 2);
	long s2 = mpz_sizeinbase(sub2, 2);
	
	if(s1 >= s2){
		(*sublen) = s1;
	}else{
		(*sublen) = s2;
	}
	
	mpz_clear(sub1);	mpz_clear(sub2);
}

void bitaddlen(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, long * addlen){
	mpz_t add1, add2;
	mpz_init(add1);	mpz_init(add2);
	
	mpz_add(add1, ax, bx);
	mpz_add(add2, ay, by);
	
	long s1 = mpz_sizeinbase(add1, 2);
	long s2 = mpz_sizeinbase(add2, 2);
	
	if(s1 >= s2){
		(*addlen) = s1;
	}else{
		(*addlen) = s2;
	}
	
	mpz_clear(add1);	mpz_clear(add2);
}

void bitsubqlen(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t q, long * subqlen){
	mpz_t sub1, sub2;
	mpz_init(sub1);	mpz_init(sub2);
	
	mpz_set(sub1, ax);	mpz_set(sub2, ay);
	mpz_submul(sub1, q, bx);
	mpz_submul(sub2, q, by);
	
	long s1 = mpz_sizeinbase(sub1, 2);
	long s2 = mpz_sizeinbase(sub2, 2);
	
	if(s1 >= s2){
		(*subqlen) = s1;
	}else{
		(*subqlen) = s2;
	}
	
	mpz_clear(sub1);	mpz_clear(sub2);
}

void bitaddqlen(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t q, long * addqlen){
	mpz_t add1, add2;
	mpz_init(add1);	mpz_init(add2);
	
	mpz_set(add1, ax);	mpz_set(add2, ay);
	mpz_addmul(add1, q, bx);
	mpz_addmul(add2, q, by);
	
	long s1 = mpz_sizeinbase(add1, 2);
	long s2 = mpz_sizeinbase(add2, 2);
	
	if(s1 >= s2){
		(*addqlen) = s1;
	}else{
		(*addqlen) = s2;
	}
	
	mpz_clear(add1);	mpz_clear(add2);
}

//
void hgcd(mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2){
	int sgnax, sgnay, sgnbx, sgnby;
	sgnax = mpz_sgn(ax);
	sgnay = mpz_sgn(ay);
	sgnbx = mpz_sgn(bx);
	sgnby = mpz_sgn(by);
	
	if(sgnax * sgnay * sgnbx * sgnby <= 0){
		if(mpz_cmpabs(ax, ay) >= 0 && mpz_cmpabs(bx, by) <= 0){
			return;
		}
		if(mpz_cmpabs(ax, ay) <= 0 && mpz_cmpabs(bx, by) >= 0){
			return;
		}
	}
	
	mpz_t ax0, ax1, ay0, ay1, bx0, bx1, by0, by1;
	mpz_t M1a1, M1a2, M1b1, M1b2;
	mpz_t M2a1, M2a2, M2b1, M2b2;

	mpz_init(ax0);	mpz_init(ax1);	mpz_init(ay0);	mpz_init(ay1);
	mpz_init(bx0);	mpz_init(bx1);	mpz_init(by0);	mpz_init(by1);
	mpz_init(M1a1);	mpz_init(M1a2);	mpz_init(M1b1);	mpz_init(M1b2);
	mpz_init(M2a1);	mpz_init(M2a2);	mpz_init(M2b1);	mpz_init(M2b2);
	
	long maxlen, minlen, sublen, addlen, subqlen, addqlen, minalen, minblen;
	long N, N2, S, p1, p2;
	
	bitlen(ax, ay, bx, by, &maxlen, &minlen);
	N = maxlen;
	S = (N / 2) + 1;
	
	integerbitlen(ax, bx, &minalen);
	integerbitlen(ay, by, &minblen);
	
	if((minalen >= (((N*3) / 4) + 2)) || (minblen >= (((N*3) / 4) + 2))){
		p1 = (N / 2);
		
		mpz_fdiv_q_2exp(ax1, ax, p1);
		mpz_fdiv_r_2exp(ax0, ax, p1);
		mpz_fdiv_q_2exp(ay1, ay, p1);
		mpz_fdiv_r_2exp(ay0, ay, p1);
		mpz_fdiv_q_2exp(bx1, bx, p1);
		mpz_fdiv_r_2exp(bx0, bx, p1);
		mpz_fdiv_q_2exp(by1, by, p1);
		mpz_fdiv_r_2exp(by0, by, p1);
		
		mpz_set_ui(M1a1, 1); mpz_set_ui(M1a2, 0); mpz_set_ui(M1b1, 0); mpz_set_ui(M1b2, 1);
		
		hgcd(ax1, ay1, bx1, by1, M1a1, M1a2, M1b1, M1b2);	
		
		updateAB(ax0, ay0, bx0, by0, M1a1, M1a2, M1b1, M1b2);
		
		mpz_mul_2exp(ax1, ax1, p1);
		mpz_mul_2exp(ay1, ay1, p1);
		mpz_mul_2exp(bx1, bx1, p1);
		mpz_mul_2exp(by1, by1, p1);
		mpz_add(ax, ax1, ax0);
		mpz_add(ay, ay1, ay0);
		mpz_add(bx, bx1, bx0);
		mpz_add(by, by1, by0);
		
		mpz_mul(M2a1, M1a1, a1);	mpz_addmul(M2a1, M1a2, b1);
		mpz_mul(M2a2, M1a1, a2);	mpz_addmul(M2a2, M1a2, b2);
		mpz_mul(M2b1, M1b1, a1);	mpz_addmul(M2b1, M1b2, b1);
		mpz_mul(M2b2, M1b1, a2);	mpz_addmul(M2b2, M1b2, b2);
		mpz_set(a1, M2a1);	mpz_set(a2, M2a2);
		mpz_set(b1, M2b1);	mpz_set(b2, M2b2);
	}
	
	mpz_t q1, q2, q, one;
	mpz_t tax, tay;
	
	mpz_init(q1);	mpz_init(q2);	mpz_init(q);	mpz_init(one);
	mpz_init(tax);	mpz_init(tay);

	mpz_set_si(one, 1);
	
	bitlen(ax, ay, bx, by, &maxlen, &minlen);
	bitsublen(ax, ay, bx, by, &sublen);
	bitaddlen(ax, ay, bx, by, &addlen);	
	
	bool flag = true;
	sgnax = mpz_sgn(ax);
	sgnay = mpz_sgn(ay);
	sgnbx = mpz_sgn(bx);
	sgnby = mpz_sgn(by);
	if(sgnax * sgnay * sgnbx * sgnby <= 0){
		flag = false;
	}
	if(flag == true && mpz_sgn(ax) != mpz_sgn(bx)){
		mpz_neg(bx, bx);
		mpz_neg(by, by);
		
		mpz_neg(b1, b1);
		mpz_neg(b2, b2);
	}
	
	
	
	while((maxlen > (((N*3) / 4) + 1)) && (sublen > S) && (addlen > S)){
		if(flag == false){
			if((mpz_cmpabs(ax, ay) >= 0 && mpz_cmpabs(bx, by) <= 0) || (mpz_cmpabs(ax, ay) <= 0 && mpz_cmpabs(bx, by) >= 0)){
				return;
			}
		}
		
		if(mpz_cmpabs(ax, ay) >= 0 && mpz_sgn(bx) != 0){
			mpz_tdiv_q(q, ax, bx);
			if(flag == true){
				mpz_mul(tay, q, by);
				mpz_sub(tay, ay, tay);
				if(mpz_sgn(tay) != mpz_sgn(ay)){
					flag = false;
				}
				if(flag == true && mpz_cmpabs(by, tay) <= 0){
					mpz_add(q, q, one);
					mpz_sub(tay, tay, by);
					flag = false;
				}
				bitsubqlen(ax, ay, bx, by, q, &subqlen);
				if(subqlen <= S){
					mpz_sub(q, q, one);
					mpz_add(tay, tay, by);
				}
				mpz_submul(ax, q, bx);
				mpz_set(ay, tay);
				mpz_swap(ax, bx);
				mpz_swap(ay, by);
			}else{
				bitsubqlen(ax, ay, bx, by, q, &subqlen);
				if(subqlen <= S){
					mpz_sub(q, q, one);
				}
				mpz_submul(ax, q, bx);
				mpz_submul(ay, q, by);
				mpz_swap(ax, bx);
				mpz_swap(ay, by);
			}
				
			mpz_swap(a1, b1);
			mpz_swap(a2, b2);
			mpz_addmul(a1, q, b1);
			mpz_addmul(a2, q, b2);
		}else{
			mpz_tdiv_q(q, ay, by);
			if(flag == true){
				mpz_mul(tax, q, bx);
				mpz_sub(tax, ax, tax);
				if(mpz_sgn(tax) != mpz_sgn(ax)){
					flag = false;
				}
				if(flag == true && mpz_cmpabs(bx, tax) <= 0){
					mpz_add(q, q, one);
					mpz_sub(tax, tax, bx);
					flag = false;
				}
				bitsubqlen(ax, ay, bx, by, q, &subqlen);
				if(subqlen <= S){
					mpz_sub(q, q, one);
					mpz_add(tax, tax, bx);
				}
				mpz_submul(ay, q, by);
				mpz_set(ax, tax);
				mpz_swap(ax, bx);
				mpz_swap(ay, by);
			}else{
				bitsubqlen(ax, ay, bx, by, q, &subqlen);
				if(subqlen <= S){
					mpz_sub(q, q, one);
				}
				mpz_submul(ax, q, bx);
				mpz_submul(ay, q, by);
				mpz_swap(ax, bx);
				mpz_swap(ay, by);
			}

				

			mpz_swap(a1, b1);
			mpz_swap(a2, b2);
			mpz_addmul(a1, q, b1);
			mpz_addmul(a2, q, b2);
		}
		
		bitlen(ax, ay, bx, by, &maxlen, &minlen);
		bitsublen(ax, ay, bx, by, &sublen);
		bitaddlen(ax, ay, bx, by, &addlen);
	}
	
	integerbitlen(ax, bx, &minalen);
	integerbitlen(ay, by, &minblen);
	if((minalen > S + 2) || (minblen > S + 2)){
		N2 = maxlen;
		p2 = S*2 - N2 +1 ;
		mpz_fdiv_q_2exp(ax1, ax, p2);
		mpz_fdiv_r_2exp(ax0, ax, p2);
		mpz_fdiv_q_2exp(ay1, ay, p2);
		mpz_fdiv_r_2exp(ay0, ay, p2);
		mpz_fdiv_q_2exp(bx1, bx, p2);
		mpz_fdiv_r_2exp(bx0, bx, p2);
		mpz_fdiv_q_2exp(by1, by, p2);
		mpz_fdiv_r_2exp(by0, by, p2);
		
		mpz_set_ui(M2a1, 1); mpz_set_ui(M2a2, 0); mpz_set_ui(M2b1, 0); mpz_set_ui(M2b2, 1);
		
		hgcd(ax1, ay1, bx1, by1, M2a1, M2a2, M2b1, M2b2);
		updateAB(ax0, ay0, bx0, by0, M2a1, M2a2, M2b1, M2b2);
		mpz_mul_2exp(ax1, ax1, p2);
		mpz_mul_2exp(ay1, ay1, p2);
		mpz_mul_2exp(bx1, bx1, p2);
		mpz_mul_2exp(by1, by1, p2);
		mpz_add(ax, ax1, ax0);
		mpz_add(ay, ay1, ay0);
		mpz_add(bx, bx1, bx0);
		mpz_add(by, by1, by0);
		
		mpz_mul(M1a1, M2a1, a1);	mpz_addmul(M1a1, M2a2, b1);
		mpz_mul(M1a2, M2a1, a2);	mpz_addmul(M1a2, M2a2, b2);
		mpz_mul(M1b1, M2b1, a1);	mpz_addmul(M1b1, M2b2, b1);
		mpz_mul(M1b2, M2b1, a2);	mpz_addmul(M1b2, M2b2, b2);
		mpz_set(a1, M1a1);	mpz_set(a2, M1a2);
		mpz_set(b1, M1b1);	mpz_set(b2, M1b2);
	}

	bitsublen(ax, ay, bx, by, &sublen);
	bitaddlen(ax, ay, bx, by, &addlen);
	
	flag = true;
	sgnax = mpz_sgn(ax);
	sgnay = mpz_sgn(ay);
	sgnbx = mpz_sgn(bx);
	sgnby = mpz_sgn(by);
	if(sgnax * sgnay * sgnbx * sgnby <= 0){
		flag = false;
	}
	
	if(flag == true && mpz_sgn(ax) != mpz_sgn(bx)){
		mpz_neg(bx, bx);
		mpz_neg(by, by);
		
		mpz_neg(b1, b1);
		mpz_neg(b2, b2);
	}

	while((sublen > S) && (addlen > S)){
		if(flag == false){
			if((mpz_cmpabs(ax, ay) >= 0 && mpz_cmpabs(bx, by) <= 0) || (mpz_cmpabs(ax, ay) <= 0 && mpz_cmpabs(bx, by) >= 0)){
				return;
			}
		}
		
		if(mpz_cmpabs(ax, ay) >= 0 && mpz_sgn(bx) != 0){
			mpz_tdiv_q(q, ax, bx);
			if(flag == true){
				mpz_mul(tay, q, by);
				mpz_sub(tay, ay, tay);
				if(mpz_sgn(tay) != mpz_sgn(ay)){
					flag = false;
				}
				if(flag == true && mpz_cmpabs(by, tay) <= 0){
					mpz_add(q, q, one);
					mpz_sub(tay, tay, by);
					flag = false;
				}
				bitsubqlen(ax, ay, bx, by, q, &subqlen);
				if(subqlen <= S){
					mpz_sub(q, q, one);
					mpz_add(tay, tay, by);
				}
				mpz_submul(ax, q, bx);
				mpz_set(ay, tay);
				mpz_swap(ax, bx);
				mpz_swap(ay, by);
			}else{
				bitsubqlen(ax, ay, bx, by, q, &subqlen);
				if(subqlen <= S){
					mpz_sub(q, q, one);
				}
				mpz_submul(ax, q, bx);
				mpz_submul(ay, q, by);
				mpz_swap(ax, bx);
				mpz_swap(ay, by);
			}

			mpz_swap(a1, b1);
			mpz_swap(a2, b2);
			mpz_addmul(a1, q, b1);
			mpz_addmul(a2, q, b2);
		}else{
			mpz_tdiv_q(q, ay, by);
			if(flag == true){
				mpz_mul(tax, q, bx);
				mpz_sub(tax, ax, tax);
				if(mpz_sgn(tax) != mpz_sgn(ax)){
					flag = false;
				}
				if(flag == true && mpz_cmpabs(bx, tax) <= 0){
					mpz_add(q, q, one);
					mpz_sub(tax, tax, bx);
					flag = false;
				}
				bitsubqlen(ax, ay, bx, by, q, &subqlen);
				if(subqlen <= S){
					mpz_sub(q, q, one);
					mpz_add(tax, tax, bx);
				}
				mpz_submul(ay, q, by);
				mpz_set(ax, tax);
				mpz_swap(ax, bx);
				mpz_swap(ay, by);
			}else{
				bitsubqlen(ax, ay, bx, by, q, &subqlen);
				if(subqlen <= S){
					mpz_sub(q, q, one);
				}
				mpz_submul(ax, q, bx);
				mpz_submul(ay, q, by);
				mpz_swap(ax, bx);
				mpz_swap(ay, by);
			}
				
			mpz_swap(a1, b1);
			mpz_swap(a2, b2);
			mpz_addmul(a1, q, b1);
			mpz_addmul(a2, q, b2);
		}
		
		bitsublen(ax, ay, bx, by, &sublen);
		bitaddlen(ax, ay, bx, by, &addlen);
	}
	
	mpz_clear(M1a1); mpz_clear(M1a2); mpz_clear(M1b1); mpz_clear(M1b2);
	mpz_clear(M2a1); mpz_clear(M2a2); mpz_clear(M2b1); mpz_clear(M2b2);
	mpz_clear(q1);	mpz_clear(q2);	mpz_clear(q);	mpz_clear(one);
	mpz_clear(ax0);	mpz_clear(ax1);	mpz_clear(ay0);	mpz_clear(ay1);	
	mpz_clear(bx0);	mpz_clear(bx1);	mpz_clear(by0);	mpz_clear(by1);
	mpz_clear(tax);	mpz_clear(tay);	
	
	return;
}
