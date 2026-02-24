#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

int main()
{	
	int a;
	FILE *fpWrite1=fopen("Input.txt","w");
	srand(time(NULL));
	
	while(1){
		a = rand()%10;
		if(a!=0)
			break;
	}
	fprintf(fpWrite1,"%d",a);
	for(int i = 1 ; i <= 10; i++){
		fprintf(fpWrite1,"%d",rand()%10);
	}
	fprintf(fpWrite1,"\n");
	
	while(1){
		a = rand()%10;
		if(a!=0)
			break;
	}
	fprintf(fpWrite1,"%d",a);
	for(int i = 1 ; i <= 2; i++){
		fprintf(fpWrite1,"%d",rand()%10);
	}
	fprintf(fpWrite1,"\n");
	
	while(1){
		a = rand()%10;
		if(a!=0)
			break;
	}
	fprintf(fpWrite1,"%d",a);
	for(int i = 1 ; i <= 10; i++){
		fprintf(fpWrite1,"%d",rand()%10);
	}
	fprintf(fpWrite1,"\n");
	
	while(1){
		a = rand()%10;
		if(a!=0)
			break;
	}
	fprintf(fpWrite1,"%d",a);
	for(int i = 1 ; i <= 1; i++){
		fprintf(fpWrite1,"%d",rand()%10);
	}
	fprintf(fpWrite1,"\n");
	fclose(fpWrite1);
	
	return 0;
}

