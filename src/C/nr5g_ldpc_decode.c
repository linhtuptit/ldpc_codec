#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Create zeros matrix with 'row' rows and 'col' columns
double** init_zeros_mat(int row, int col)
{
	double **M;
	int rdx, cdx;
	
	M = (double**) malloc(row * sizeof(double*));
	if (!M) 	return 0;
	for (rdx = 0; rdx < row; rdx++){
		M[rdx] = (double*) malloc(col * sizeof(double));
		if(!M[rdx]) 	return 0;
	}
	
	for(rdx = 0; rdx < row; rdx++){
		for(cdx = 0; cdx < col; cdx++){
			M[rdx][cdx] = 0;
		}
	}
	return M;
}

// Free memory
void FreeMem(double** Buff, int row)
{
	int idx;
	for(idx = 0; idx < row; idx++){
		if(Buff[idx]){
			free(Buff[idx]);
			Buff[idx] = NULL;
		}
	}
	if(Buff)	free(Buff);
}

// Define sign of float number using fabsf() funtion
int sign(double x)
{
	int out;
	if(x == 0)	out = 1;
	else	out = int (fabs(x)/x);
	return out;
}

// Define min of two double numbers
double min(double in1, double in2)
{
	return ((in1 > in2) ? in2 : in1);
}

// Define phi function for check-to-variable processing
double ldpc_phi(double in)
{
	double out;
	if(in < 0)	printf("\n Input value of phi function must be positive");
	else if(in == 0)	out = -20.0;
	else	out = log10((1.0 + exp(-in))/(1.0 - exp(-in)));
	return out;
}

// Define approximation function to calculate logarit nepe 
double approx(double in)
{
	double out;
	if(fabs(in) < 2.5)	out = 0.6 - 0.24 * fabs(in);
	else	out = 0;
	return out;
}

// Define box-plus function for check-to-variable processing
double ldpc_boxplus(double in1, double in2)
{
	double out;
	out = sign(in1) * sign(in2) * min(fabs(in1), fabs(in2)) + (approx(in1+in2) - approx(in1-in2));
	return out;
}

// This function to define each bit in the output bit stream is 0 or 1
int* ldpc_hard_decode(double* in, int len)
{
	int* out = (int*) malloc(len * sizeof(int));
	int idx;
	for(idx = 0; idx < len; idx++){
		if(in[idx] >= 0) 	out[idx] = 0;
		else	out[idx] = 1;
	}
	return out; 
}

// This function to check whether output bit stream is ...
// satisfied with parity condition or not
int check_syndrome(int *in, int **H, int row, int col)
{
	int* chk_vec = (int*) malloc(row * sizeof(int));
	int flag = 0;
	int idx, jdx;
	for(idx = 0; idx < row; idx++){
		int chk_sum = 0;
		for(jdx = 0; jdx < col; jdx++){
			chk_sum += in[jdx] * H[idx][jdx];
		}
		chk_vec[idx] = (chk_sum % 2);
		if (chk_vec[idx] == 1){
			flag = 1;
			break;
		} 	
		else	continue;
	}
	free(chk_vec);
	return flag;
}

// LDPC decoder based on Belief Propagation Algorithm ...
// using phi funtion to execute check-to-variable processing
int* ldpc_decode_blfprg(double* LLR_in, int** H, int row, int col, int max_iter)
{
	int* bit_out = (int*) malloc(col * sizeof(int));
	
	double *LLR_out = (double*) malloc(col * sizeof(double));
	double **VC_msg_, **CV_msg_;
	VC_msg_ = (double**) init_zeros_mat(col, row);
	CV_msg_ = (double**) init_zeros_mat(row, col);
	int rdx, cdx, kdx;
	
	/* Initialization phase*/
	for(rdx = 0; rdx < row; rdx++){
		for(cdx = 0; cdx < col; cdx++){
			if(H[rdx][cdx] == 1)	VC_msg_[cdx][rdx] = LLR_in[cdx];
		}
	}
	/* LOOP phase*/
	int iter = 1;
	while (1)
	{
		// Check node processing
		for(rdx = 0; rdx < row; rdx++){
			for(cdx = 0; cdx < col; cdx++){
				if(H[rdx][cdx] == 1){
					double prod_cv = 1; double sum_cv = 0;
					for(kdx = 0; kdx < col; kdx++){
						if((H[rdx][kdx] == 1) && (kdx != cdx)){
							prod_cv *= sign(VC_msg_[kdx][rdx]);
							sum_cv += ldpc_phi(fabs(VC_msg_[kdx][rdx]));
						}
					}
					CV_msg_[rdx][cdx] = prod_cv*ldpc_phi(fabs(sum_cv));
				}
			}
		}
		// Variable node update
		for(cdx = 0; cdx < col; cdx++){
			double sum_vp = 0;
			for(rdx = 0; rdx < row; rdx++){
				sum_vp += CV_msg_[rdx][cdx];
			}
			LLR_out[cdx] = LLR_in[cdx] + sum_vp;
		}
		// Variable node processing
		for(rdx = 0; rdx < row; rdx++){
			for(cdx = 0; cdx < col; cdx++){
				if(H[rdx][cdx] == 1){
					VC_msg_[cdx][rdx] = LLR_out[cdx] - CV_msg_[rdx][cdx];
				}
			}
		}
		// Hard decoder
		bit_out = ldpc_hard_decode(LLR_out, col);
		// check syndrome
		int flag;
		flag = check_syndrome(bit_out, H, row, col);
		if ((flag == 0) || (iter == max_iter)){
			printf("\n Number of iteration is: %d", iter);
			break;
		}
		else if ((flag == 1) && (iter < max_iter)){
			iter += 1;
		}
	}
	// FreeMem(VC_msg_, row); FreeMem(CV_msg_, col);
	return bit_out;
}

// LDPC decoder based on Minimum Summation algorithm ...
// using minimum function to execute check-to-variable processing
int* ldpc_decode_minsum(double* LLR_in, int** H, int row, int col, int max_iter)
{
	int* bit_out = (int*) malloc(col * sizeof(int));
	
	double *LLR_out = (double*) malloc(col * sizeof(double));
	double **VC_msg_, **CV_msg_;
	VC_msg_ = (double**) init_zeros_mat(col, row);
	CV_msg_ = (double**) init_zeros_mat(row, col);
	int rdx, cdx, kdx;
	
	/* Initialization phase*/
	for(rdx = 0; rdx < row; rdx++){
		for(cdx = 0; cdx < col; cdx++){
			if(H[rdx][cdx] == 1)	VC_msg_[cdx][rdx] = LLR_in[cdx];
		}
	}
	/* LOOP phase*/
	int iter = 1;
	while (1)
	{
		// Check node processing
		for(rdx = 0; rdx < row; rdx++){
			for(cdx = 0; cdx < col; cdx++){
				if(H[rdx][cdx] == 1){
					double prod_cv = 1; double min_cv = 1000.0;
					for(kdx = 0; kdx < col; kdx++){
						if((H[rdx][kdx] == 1) && (kdx != cdx)){
							prod_cv *= sign(VC_msg_[kdx][rdx]);
							if(fabs(VC_msg_[kdx][rdx]) < min_cv){
								min_cv = fabs(VC_msg_[kdx][rdx]);
							}
						}
					}
					CV_msg_[rdx][cdx] = prod_cv*min_cv;
				}
			}
		}
		// Variable node update
		for(cdx = 0; cdx < col; cdx++){
			double sum_vp = 0;
			for(rdx = 0; rdx < row; rdx++){
				sum_vp += CV_msg_[rdx][cdx];
			}
			LLR_out[cdx] = LLR_in[cdx] + sum_vp;
		}
		// Variable node processing
		for(rdx = 0; rdx < row; rdx++){
			for(cdx = 0; cdx < col; cdx++){
				if(H[rdx][cdx] == 1){
					VC_msg_[cdx][rdx] = LLR_out[cdx] - CV_msg_[rdx][cdx];
				}
			}
		}
		// Hard decoder
		bit_out = ldpc_hard_decode(LLR_out, col);
		// check syndrome
		int flag;
		flag = check_syndrome(bit_out, H, row, col);
		if ((flag == 0) || (iter == max_iter)){
			printf("\n Number of iteration is: %d", iter);
			break;
		}
		else if ((flag == 1) && (iter < max_iter)){
			iter += 1;
		}
	}
	// FreeMem(VC_msg_, row); FreeMem(CV_msg_, col);
	return bit_out;
}

// LDPC decoder based on Box Plus algorithm ...
// using Box-Plus function to execute check-to-variable processing
int* ldpc_decode_boxpls(double* LLR_in, int** H, int row, int col, int max_iter)
{
	int* bit_out = (int*) malloc(col * sizeof(int));
	
	double *LLR_out = (double*) malloc(col * sizeof(double));
	double **VC_msg_, **CV_msg_;
	VC_msg_ = (double**) init_zeros_mat(col, row);
	CV_msg_ = (double**) init_zeros_mat(row, col);
	int rdx, cdx, kdx;
	
	/* Initialization phase*/
	for(rdx = 0; rdx < row; rdx++){
		for(cdx = 0; cdx < col; cdx++){
			if(H[rdx][cdx] == 1)	VC_msg_[cdx][rdx] = LLR_in[cdx];
		}
	}
	/* LOOP phase*/
	int iter = 1;
	while (1)
	{
		// Check node processing
		for(rdx = 0; rdx < row; rdx++){
			for(cdx = 0; cdx < col; cdx++){
				if(H[rdx][cdx] == 1){
					double box_val_ = 1000.0;
					for(kdx = 0; kdx < col; kdx++){
						if((H[rdx][kdx] == 1) && (kdx != cdx)){
							box_val_ = ldpc_boxplus(VC_msg_[kdx][rdx], box_val_);
						}
					}
					CV_msg_[rdx][cdx] = box_val_;
				}
			}
		}
		// Variable node update
		for(cdx = 0; cdx < col; cdx++){
			double sum_vp = 0;
			for(rdx = 0; rdx < row; rdx++){
				sum_vp += CV_msg_[rdx][cdx];
			}
			LLR_out[cdx] = LLR_in[cdx] + sum_vp;
		}
		// Variable node processing
		for(rdx = 0; rdx < row; rdx++){
			for(cdx = 0; cdx < col; cdx++){
				if(H[rdx][cdx] == 1){
					VC_msg_[cdx][rdx] = LLR_out[cdx] - CV_msg_[rdx][cdx];
				}
			}
		}
		// Hard decoder
		bit_out = ldpc_hard_decode(LLR_out, col);
		// check syndrome
		int flag;
		flag = check_syndrome(bit_out, H, row, col);
		if ((flag == 0) || (iter == max_iter)){
			printf("\n Number of iteration is: %d", iter);
			break;
		}
		else if ((flag == 1) && (iter < max_iter)){
			iter += 1;
		}
	}
	// FreeMem(VC_msg_, row); FreeMem(CV_msg_, col);
	return bit_out;
}

// Testing
int main()
{
	double LLR_in[8] = {1.0,0.5,0.5,2.0,1.0,-1.5,1.5,-1.0};
	int bit_in[8] = {0,0,0,0,1,1,0,1};
	int **H;
	int idx, jdx;
	
	H = (int**) malloc(4 * sizeof(int*));
	for(idx = 0; idx < 4; idx++){
		H[idx] = (int*) malloc(8 * sizeof(int));
	}

	FILE *file;
	file = fopen("H.txt", "r");
	for(idx = 0; idx < 4; idx++){
		for(jdx = 0; jdx < 8; jdx++){
			fscanf(file, "%d", &H[idx][jdx]);
		}
	}
	printf("\n Matrix H is \n");
	for(idx = 0; idx < 4; idx++){
		for(jdx = 0; jdx < 8; jdx++){
			printf("  %d", H[idx][jdx]);
		}
		printf("\n");
	}
	
	//
	int* bit_out = (int*) malloc(8 * sizeof(int));
	bit_out = ldpc_decode_blfprg(LLR_in, H, 4, 8, 3);
	printf("\n Output bit stream of LDPC decode is:");
	for(idx = 0; idx < 8; idx++){
		printf(" %d", bit_out[idx]);
	}
	return 0;
}
