#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int Gol_Euc(mpz_t ansp, mpz_t ansq, mpz_t vax, mpz_t vay, mpz_t vbx, mpz_t vby);
void search(mpz_t result, mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);
void findanswer(mpz_t min1, mpz_t min2, mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by);

int main(int argc, char *argv[]){
	mpz_t a[2], b[2];
	mpz_t p, q, min1, min2;
	clock_t start, finish;  
	
	mpz_init(a[0]);	mpz_init(a[1]);	mpz_init(b[0]);	mpz_init(b[1]);
	mpz_init(p); mpz_init(q);
	mpz_init(min1);	mpz_init(min2);
	
	printf("\nGolbal Euclidean algorithm.\n");
	
	FILE *fp = fopen("Input.txt", "r");
	gmp_fscanf(fp, "%Zd %Zd %Zd %Zd", a[0], a[1], b[0], b[1]);
	
	start = clock();
	int rectime = Gol_Euc(p, q, a[0], a[1], b[0], b[1]);
	finish = clock();
	
	findanswer(min1, min2, a[0], a[1], b[0], b[1]);
	
	FILE *fo;
	fo = fopen("Output.txt", "w");
	gmp_fprintf(fo, "%Zd\n",a[0]);
	gmp_fprintf(fo, "%Zd\n",a[1]);
	gmp_fprintf(fo, "%Zd\n",b[0]);
	gmp_fprintf(fo, "%Zd\n",b[1]);
	fclose(fo);
	
	printf("Golbal Euclidean cycle times = %d\n",rectime);	
	printf("Prog time: %lf\n", (double)(finish-start)/CLOCKS_PER_SEC);
	
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


void search(mpz_t res, mpz_t ax, mpz_t ay, mpz_t bx, mpz_t by){
    int f1, f2, f3, f4;
    int n1 = 0, n2 = 1, n3 = 2, n4 = 3;
    mpz_t tp[4];
    mpz_init(tp[0]); mpz_init(tp[1]); mpz_init(tp[2]); mpz_init(tp[3]);

	int m1, m2;
	if(mpz_sgn(ax) == -1)	f1 = 1;	
	else f1 = 0;
	if(mpz_sgn(ay) == -1)	f2 = 1;	
	else f2 = 0;
	if(mpz_sgn(bx) == -1)	f3 = 1;	
	else f3 = 0;
	if(mpz_sgn(by) == -1)	f4 = 1;	
	else f4 = 0;
	
	if(f3)	mpz_neg(bx, bx);
	if(f4)	mpz_neg(by, by);
	if(f1)	mpz_neg(ax, ax);
	if(f2)	mpz_neg(ay, ay);
	
	if((f1^f2^f3^f4) == 1){
		mpz_sub(tp[n1], bx, by);
		mpz_add(tp[n2], ax, ay);
		mpz_abs(tp[n1], tp[n1]);
	}else{
		mpz_add(tp[n1], bx, by);
		mpz_add(tp[n2], ax, ay);
	}
	if(f3)	mpz_neg(bx, bx);
	if(f4)	mpz_neg(by, by);
	if(f1)	mpz_neg(ax, ax);
	if(f2)	mpz_neg(ay, ay);
	
    mpz_fdiv_q(res, tp[n1], tp[n2]);
	mpz_mul(tp[n1], ax, res);
	mpz_mul(tp[n2], ay, res);
	mpz_sub(tp[n1], bx, tp[n1]);
    	mpz_sub(tp[n2], by, tp[n2]);
	mpz_sub(tp[n3], tp[n1], ax);
    	mpz_sub(tp[n4], tp[n2], ay);
    
	if(mpz_cmpabs(tp[n1],tp[n2]) > 0)	m1 = n1;
	else	m1 = n2;
	if(mpz_cmpabs(tp[n3],tp[n4]) > 0)	m2 = n3;
	else	m2 = n4;
	
	if(mpz_cmpabs(tp[m1],tp[m2]) > 0){
		mpz_set(bx, tp[n3]);
		mpz_set(by, tp[n4]);
	}else{
		mpz_set(bx, tp[n1]);
		mpz_set(by, tp[n2]);
	}

    mpz_clear(tp[0]); mpz_clear(tp[1]); mpz_clear(tp[2]); mpz_clear(tp[3]);
    return;
}

int Gol_Euc(mpz_t ansp, mpz_t ansq, mpz_t vax, mpz_t vay, mpz_t vbx, mpz_t vby){
    mpz_t V[4],T[4],u;
    mpz_init_set(V[0], vax); mpz_init_set(V[1], vay); mpz_init_set(V[2], vbx); mpz_init_set(V[3], vby);
	mpz_init(T[0]);	mpz_init(T[1]);	mpz_init(T[2]);	mpz_init(T[3]);
	mpz_init(u);

    int ax = 0, ay = 1, bx = 2, by = 3;
	int tlx = 0, tly = 1, tsx = 2, tsy = 3;
    int nfind = 1;
	int ma,mb;
    int recnt = 0;
    int vch[4] = {2,3,0,1};

	if(mpz_cmpabs(V[ax],V[ay]) > 0)	ma = ax;
	else	ma = ay;
	if(mpz_cmpabs(V[bx],V[by]) > 0)	mb = bx;
	else	mb = by;
	if(mpz_cmpabs(V[ma],V[mb]) > 0){
		ax = vch[ax];
		ay = vch[ay];
		bx = vch[bx];
		by = vch[by];
	}
	mpz_sub(T[tsx], V[ax], V[bx]);
	mpz_sub(T[tsy], V[ay], V[by]);
	mpz_add(T[tlx], V[ax], V[bx]);
	mpz_add(T[tly], V[ay], V[by]);
	if(mpz_cmpabs(T[tsx],T[tsy]) > 0)	ma = tsx;
	else	ma = tsy;
	if(mpz_cmpabs(T[tlx],T[tly]) > 0)	mb = tlx;
	else	mb = tly;
	if(mpz_cmpabs(T[ma],T[mb]) > 0){
		mpz_neg(V[bx],V[bx]);
		mpz_neg(V[by],V[by]);
		tlx = vch[tlx];
		tly = vch[tly];
		tsx = vch[tsx];
		tsy = vch[tsy];
	}
	if(mpz_cmpabs(T[tsx],T[tsy]) > 0)	ma = tsx;
	else	ma = tsy;
	if(mpz_cmpabs(V[bx],V[by]) > 0)	mb = bx;
	else	mb = by;
	if(mpz_cmpabs(V[mb],T[ma]) <= 0)	nfind = 0;
	
	if(nfind){
		int flag = 1;
		if(mpz_cmpabs(T[tsx],T[tsy]) > 0)	ma = tsx;
		else	ma = tsy;
		if(mpz_cmpabs(V[ax],V[ay]) > 0)	mb = ax;
		else	mb = ay;
		if(mpz_cmpabs(V[mb],T[ma]) < 0)	flag = 0;
	
		if(flag){
			if(mpz_cmpabs(V[ax],V[ay]) > 0)	ma = ax;
			else	ma = ay;
			if(mpz_cmpabs(V[bx],V[by]) > 0)	mb = bx;
			else	mb = by;
			
			if(mpz_cmpabs(V[ma],V[mb]) == 0){
				mpz_sub(V[bx], V[ax], V[bx]);
				mpz_sub(V[by], V[ay], V[by]);
				nfind = 0;
			}else{
				mpz_sub(V[ax], V[bx], V[ax]);
				mpz_sub(V[ay], V[by], V[ay]);
			}	
		}
	}
	
	while(nfind){	
		recnt++;
		search(u, V[ax], V[ay], V[bx], V[by]);
		
		ax = vch[ax];
		ay = vch[ay];
		bx = vch[bx];
		by = vch[by];
		
		mpz_sub(T[tsx], V[ax], V[bx]);
		mpz_sub(T[tsy], V[ay], V[by]);
		mpz_add(T[tlx], V[ax], V[bx]);
		mpz_add(T[tly], V[ay], V[by]);
		if(mpz_cmpabs(T[tsx],T[tsy]) > 0)	ma = tsx;
		else	ma = tsy;
		if(mpz_cmpabs(T[tlx],T[tly]) > 0)	mb = tlx;
		else	mb = tly;
		if(mpz_cmpabs(T[ma],T[mb]) > 0){
			mpz_neg(V[bx],V[bx]);
			mpz_neg(V[by],V[by]);
			ma = mb;
		}
		
		if(mpz_cmpabs(V[bx],V[by]) > 0)	mb = bx;
		else	mb = by;
		if(mpz_cmpabs(V[mb],T[ma]) <= 0){
			if(mpz_cmpabs(V[ax],V[ay]) > 0)	ma = ax;
			else	ma = ay;
			if(mpz_cmpabs(V[bx],V[by]) > 0)	mb = bx;
			else	mb = by;
			if(mpz_cmpabs(V[ma],V[mb]) > 0){
				ax = vch[ax];
				ay = vch[ay];
				bx = vch[bx];
				by = vch[by];
			}

			break;
		}
	}

	mpz_set(vax, V[ax]);
	mpz_set(vay, V[ay]);
	mpz_set(vbx, V[bx]);
	mpz_set(vby, V[by]);

	mpz_clear(V[0]); mpz_clear(V[1]); mpz_clear(V[2]); mpz_clear(V[3]);
	mpz_clear(T[0]); mpz_clear(T[1]); mpz_clear(T[2]); mpz_clear(T[3]);
	mpz_clear(u);
    return recnt;
}
