#include <stdio.h>
#include<Windows.h>
#include<vector>
#include<iostream>

using namespace std;



void main(void)
{
	int n = 10;
	int num = 0;
	for (int i = 2; i < n - 1; i++){
		num += n / i;
	}

	//////////////////creating S array///////////////////////////////////
	int *S = (int*)(calloc(num, sizeof(int)));
	int kk = 0, j;
	for (int i = 2; i < n - 1; i++){
		for (int j = 0; j < n / i; j++){
			S[kk++] = i;
		}

	}
	///////////////////////////////////////////////////////////////////////
	////////////////////////        //////////////////////////////////////
	int X[5];
	int k = 0; int i = 0;
	int sum = 0;
	while (k >= 0 && X[0]<=n/2){
		X[k] = S[i];
		sum += X[k];
		/*for (int m = 0; m <= k; m++)
			sum = sum + X[m];*/
		if (sum < n )
		{
			if (sum + S[i + 1] <= n){
				k++; i++;
			}
			else{
				i++;
				sum-=X[k];
			}
		}
		else if (sum > n) { 
			sum -= X[k];
			k --; 
			sum -= X[k];
			for (int j = 0; j <num; j++){//find the minimum number in S that is greater than X[k](a formula can be replaced)
				if (S[j] > X[k]){
					i = j;
					break;
				}
			}
		}
		else {
			for (int p = 0; p <= k; p++) printf("%d", X[p]); // to be extended for the next step
			printf("\n");
			sum -= X[k];
			k--;
			sum -= X[k];
			for (int j = 0; j <num; j++){//find the minimum number in S that is greater than X[k](a formula can be replaced)
				if (S[j] > X[k]){
					i = j;
					break;
				}
			}
		}
	}
	printf("%d\n", n);
	getchar();
}

