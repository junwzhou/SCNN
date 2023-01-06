#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <sys/time.h>
// Multi-processing
#include <sys/wait.h>
#include <unistd.h>
// tfhe-lib
#include "tfhe.h"
#include "tfhe_garbage_collector.h"
#include "tlwe.h"
#include "tgsw.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "polynomials.h"
using namespace std;
// Defines
#define VERBOSE 1
#define STATISTICS true
#define WRITELATEX true
#define N_PROC 4

// Security constants
#define SECLEVEL 80
#define SECNOISE true
#define SECALPHA pow(2., -20)S
#define SEC_PARAMS_STDDEV    pow(2., -30)
#define SEC_PARAMS_n  400         ///  LweParams
#define SEC_PARAMS_N 1024                   /// TLweParams
#define SEC_PARAMS_k    1                   /// TLweParams
#define SEC_PARAMS_BK_STDDEV pow(2., -36)   /// TLweParams
#define SEC_PARAMS_BK_BASEBITS 10           /// TGswParams
#define SEC_PARAMS_BK_LENGTH    3           /// TGswParams
#define SEC_PARAMS_KS_STDDEV pow(2., -25)   /// Key Switching Params
#define SEC_PARAMS_KS_BASEBITS  1           /// Key Switching Params
#define SEC_PARAMS_KS_LENGTH   18           /// Key Switching Params


#define torus_test 1000
#define num_running  10
const double alpha       = pow(2., -20);

using namespace std;



/// Generate gate bootstrapping parameters for FHE_NN
TFheGateBootstrappingParameterSet *our_default_gate_bootstrapping_parameters(int minimum_lambda)
{
    if (minimum_lambda > 128)
        cerr << "Sorry, for now, the parameters are only implemented for about 128bit of security!\n";

    static const int n = SEC_PARAMS_n;
    static const int N = SEC_PARAMS_N;
    static const int k = SEC_PARAMS_k;
    static const double max_stdev = SEC_PARAMS_STDDEV;

    static const int bk_Bgbit    = SEC_PARAMS_BK_BASEBITS;  //<-- ld, thus: 2^10
    static const int bk_l        = SEC_PARAMS_BK_LENGTH;
    static const double bk_stdev = SEC_PARAMS_BK_STDDEV;

    static const int ks_basebit  = SEC_PARAMS_KS_BASEBITS; //<-- ld, thus: 2^1
    static const int ks_length   = SEC_PARAMS_KS_LENGTH;
    static const double ks_stdev = SEC_PARAMS_KS_STDDEV;


    LweParams  *params_in    = new_LweParams (n,    ks_stdev, max_stdev);
    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
    TGswParams *params_bk    = new_TGswParams(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    return new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
}
 TFheGateBootstrappingParameterSet *params = our_default_gate_bootstrapping_parameters(80);
    //TFheGateBootstrappingParameterSet *params = new_default_gate_bootstrapping_parameters(minimum_lambda);
const LweParams *in_out_params   = params->in_out_params;
TFheGateBootstrappingSecretKeySet *secret = new_random_gate_bootstrapping_secret_keyset(params);
const LweBootstrappingKeyFFT *bs_key = secret->cloud.bkFFT;
struct C
{
LweSample*ciphertext=new_LweSample(in_out_params);
};
void encrypted_data(LweSample*output,float input)
{
LweSample*tmp=new_LweSample(in_out_params);
lweClear(tmp,in_out_params);
Torus32 mu;
mu=floatmodSwitchToTorus32(1,torus_test);
mu=mu*input;
lweSymEncrypt(tmp, mu, alpha, secret->lwe_key);
lweCopy(output,tmp,in_out_params);
delete_LweSample(tmp);
}

void encrypted_int_data(LweSample*output,int input)
{
int lable;
if(input>0.5)
{
lable=1;
}
if(input<0.5)
{
lable=-1;
}
LweSample*tmp=new_LweSample(in_out_params);
lweClear(tmp,in_out_params);
Torus32 mu;
mu=floatmodSwitchToTorus32(lable,8);
//mu=mu*lable;
lweSymEncrypt(tmp, mu, alpha, secret->lwe_key);
lweCopy(output,tmp,in_out_params);
delete_LweSample(tmp);
}
void SIGN(LweSample*output,LweSample*input,int nbit)
{
LweSample*sign=new_LweSample(in_out_params);
Torus32 mu_con=modSwitchToTorus32(1,torus_test);
tfhe_bootstrap_FFT(sign,bs_key,mu_con,input);
LweSample*tmp1=new_LweSample_array(nbit,in_out_params);
LweSample*tmp1_1=new_LweSample(in_out_params);
lweCopy(tmp1_1,input,in_out_params);
LweSample*tmp2=new_LweSample_array(nbit,in_out_params);
LweSample*tmp2_2=new_LweSample(in_out_params);
LweSample*N_input=new_LweSample(in_out_params);
lweNegate(N_input,input,in_out_params);
lweCopy(tmp2_2,N_input,in_out_params);
int L=nbit-2;
LweSample*b=new_LweSample(in_out_params);
encrypted_data(b,1);
LweSample*b_2=new_LweSample(in_out_params);
encrypted_data(b_2,0.5);
LweSample*b_3=new_LweSample(in_out_params);
encrypted_data(b_3,-3);
for(int i=L;i>=0;i--)
{
int p=pow(2,i);
lweSubMulTo(tmp1_1,p,b,in_out_params);
LweSample*s_1=new_LweSample_array(1,in_out_params);
lweCopy(s_1,tmp1_1,in_out_params);
lweAddMulTo(tmp1_1,1,b_2,in_out_params);
LweSample*ls=new_LweSample(in_out_params);
lweCopy(ls,tmp1_1,in_out_params);
lweAddMulTo(tmp1_1,2,ls,in_out_params);
delete_LweSample(ls);
LweSample*s_2=new_LweSample_array(1,in_out_params);
lweCopy(s_2,tmp1_1,in_out_params);
tfhe_bootstrap_FFT(tmp1_1,bs_key,mu_con,s_2);
lweCopy(tmp1+i,tmp1_1,in_out_params);
lweCopy(s_2,tmp1_1,in_out_params);
LweSample*s_3_1=new_LweSample(in_out_params);
lweNegate(s_3_1,s_2,in_out_params);
lweAddTo(s_3_1,b,in_out_params);
LweSample*lss=new_LweSample(in_out_params);
lweClear(lss,in_out_params);
floatlweAddMulTo(lss,0.5,s_3_1,in_out_params);
lweCopy(s_3_1,lss,in_out_params);
lweAddMulTo(s_1,p,s_3_1,in_out_params);
lweCopy(tmp1_1,s_1,in_out_params);
delete_LweSample(s_1);
delete_LweSample(s_2);
delete_LweSample(s_3_1);
delete_LweSample(lss);
}
for(int i=L;i>=0;i--)
{
int p=pow(2,i);
lweSubMulTo(tmp2_2,p,b,in_out_params);
LweSample*s_1=new_LweSample_array(1,in_out_params);
lweCopy(s_1,tmp2_2,in_out_params);
lweAddMulTo(tmp2_2,1,b_2,in_out_params);
LweSample*ls=new_LweSample(in_out_params);
lweCopy(ls,tmp2_2,in_out_params);
lweAddMulTo(tmp2_2,2,ls,in_out_params);
delete_LweSample(ls);
LweSample*s_2=new_LweSample_array(1,in_out_params);
lweCopy(s_2,tmp2_2,in_out_params);
tfhe_bootstrap_FFT(tmp2_2,bs_key,mu_con,s_2);
lweCopy(tmp2+i,tmp2_2,in_out_params);
lweCopy(s_2,tmp2_2,in_out_params);
LweSample*s_3_1=new_LweSample(in_out_params);
lweNegate(s_3_1,s_2,in_out_params);
lweAddTo(s_3_1,b,in_out_params);
LweSample*lss=new_LweSample(in_out_params);
lweClear(lss,in_out_params);
floatlweAddMulTo(lss,0.5,s_3_1,in_out_params);
lweCopy(s_3_1,lss,in_out_params);
lweAddMulTo(s_1,p,s_3_1,in_out_params);
lweCopy(tmp2_2,s_1,in_out_params);
delete_LweSample(s_1);
delete_LweSample(s_2);
delete_LweSample(s_3_1);
delete_LweSample(lss);
}
LweSample*f=new_LweSample(in_out_params);
lweClear(f,in_out_params);
lweAddTo(f,tmp1+L,in_out_params);
lweAddTo(f,tmp2+L,in_out_params);
lweAddTo(f,b_3,in_out_params);
LweSample*boot=new_LweSample(in_out_params);
tfhe_bootstrap_FFT(boot,bs_key,mu_con,f);

LweSample*result=new_LweSample(in_out_params);
lweClear(result,in_out_params);
lweAddTo(result,boot,in_out_params);
lweAddTo(result,sign,in_out_params);
lweAddTo(result,b,in_out_params);
LweSample*C=new_LweSample(in_out_params);
tfhe_bootstrap_FFT(C,bs_key,mu_con,result);
lweCopy(output,C,in_out_params);

delete_LweSample(f);
delete_LweSample(boot);
delete_LweSample(result);
delete_LweSample(C);
delete_LweSample(sign);
delete_LweSample(tmp1_1);
delete_LweSample(tmp2_2);
delete_LweSample_array(nbit,tmp1);
delete_LweSample_array(nbit,tmp2);
delete_LweSample(N_input);
delete_LweSample(b);
delete_LweSample(b_2);
delete_LweSample(b_3);
}
void test_decryption(LweSample*ciphertext)
{
Torus32 mu_con=modSwitchToTorus32(1,torus_test);
Torus32 R_sum=lwePhase(ciphertext,secret->lwe_key);
float a1=R_sum;
float b2=mu_con;
//cout<<"value :"<<a1/b2<<endl;
}
int fully_test_decryption(LweSample*ciphertext)
{
Torus32 mu_con=modSwitchToTorus32(1,8);
Torus32 R_sum=lwePhase(ciphertext,secret->lwe_key);
float a1=R_sum;
float b2=mu_con;
float c3=a1/b2;
if(c3>0)
{
return 1;
}
else{
return 0;
}
}

void test2_decryption(LweSample*ciphertext)
{
Torus32 mu_con=modSwitchToTorus32(1,1000);
Torus32 R_sum=lwePhase(ciphertext,secret->lwe_key);
float a1=R_sum;
float b2=mu_con;
//cout<<"value :"<<a1/b2<<endl;
}
void test3_decryption(LweSample*ciphertext)
{
Torus32 mu_con=modSwitchToTorus32(1,8);
Torus32 R_sum=lwePhase(ciphertext,secret->lwe_key);
float a1=R_sum;
float b2=mu_con;
//cout<<"value :"<<a1<<" "<<b2<<" "<<a1/b2<<endl;
}
float test4_decryption(LweSample*ciphertext)
{
Torus32 mu_con=modSwitchToTorus32(1,8);
Torus32 R_sum=lwePhase(ciphertext,secret->lwe_key);
float a1=R_sum;
float b2=mu_con;
float result=a1/b2;
if(result>0)
{
return 1;
}
else
{
return -1;
}
}
void convolution_operations(C output[3][24][24],C input[28][28],float params[3][5][5],float bias[3])
{
for(int lay=0;lay<3;lay++)
{
    for(int x=0;x<24;x++)
    {
        for(int y=0;y<24;y++)
        {
        LweSample*sum=new_LweSample(in_out_params);
        lweClear(sum,in_out_params);
        //test_decryption(sum);
        LweSample*sum_ex=new_LweSample(in_out_params);
        lweClear(sum_ex,in_out_params);
        for(int x_keral=0;x_keral<5;x_keral++)
        {
            for(int y_keral=0;y_keral<5;y_keral++)
            {
                LweSample*curry_sum=new_LweSample(in_out_params);
                floatlweAddMulTo(curry_sum,params[lay][x_keral][y_keral],input[x+x_keral][y+y_keral].ciphertext,in_out_params);
                lweAddTo(sum,curry_sum,in_out_params);
                //test_decryption(input[x+x_keral][y+y_keral].ciphertext);
                //test_decryption(sum);
            }
        }
        lweCopy(output[lay][x][y].ciphertext,sum,in_out_params);
        //delete_LweSample(sum);
        LweSample*tmp=new_LweSample(in_out_params);
        lweClear(tmp,in_out_params);
        encrypted_data(tmp,bias[lay]*50);
        lweAddTo(output[lay][x][y].ciphertext,tmp,in_out_params);
        lweCopy(sum,output[lay][x][y].ciphertext,in_out_params);
        //lweAddMulTo(sum_ex,100,sum,in_out_params); change
        lweAddMulTo(sum_ex,1,sum,in_out_params);
        Torus32 mu_con=modSwitchToTorus32(1,torus_test);
        //delete_LweSample(tmp);
        //test_decryption(sum_ex);
        //SIGN(output[lay][x][y].ciphertext,sum_ex,8);
        tfhe_bootstrap_FFT(output[lay][x][y].ciphertext,bs_key,mu_con,sum_ex);
        delete_LweSample(sum);
        delete_LweSample(tmp);
        delete_LweSample(sum_ex);
        }
    }
}
}
void H_OR(LweSample*output,LweSample*a,LweSample*b,LweSample*c,LweSample*d)
{
LweSample*tmp=new_LweSample(in_out_params);
lweClear(tmp,in_out_params);
lweAddTo(tmp,a,in_out_params);
lweAddTo(tmp,b,in_out_params);
lweAddTo(tmp,c,in_out_params);
lweAddTo(tmp,d,in_out_params);

LweSample*tmp2=new_LweSample(in_out_params);
lweClear(tmp2,in_out_params);
Torus32 mu;
mu=modSwitchToTorus32(3,torus_test);
lweSymEncrypt(tmp2, mu, alpha, secret->lwe_key);
lweAddTo(tmp,tmp2,in_out_params);
LweSample*output_ex=new_LweSample(in_out_params);
Torus32 mu_con=modSwitchToTorus32(1,torus_test);
tfhe_bootstrap_FFT(output_ex,bs_key,mu_con,tmp);
Torus32 mu_con_f=modSwitchToTorus32(1,8);
tfhe_bootstrap_FFT(output,bs_key,mu_con_f,output_ex);
delete_LweSample(tmp);
delete_LweSample(tmp2);
delete_LweSample(output_ex);
}
void pooling(C output[3][12][12],C input[3][24][24])
{
for(int lay=0;lay<3;lay++)
{
    for(int x=0;x<12;x++)
    {
        for(int y=0;y<12;y++)
        {
        H_OR(output[lay][x][y].ciphertext,input[lay][2*x][2*y].ciphertext,input[lay][2*x+1][2*y].ciphertext,
        input[lay][2*x][2*y+1].ciphertext,input[lay][2*x+1][2*y+1].ciphertext);
        //cout<<"lay : "<<lay<<" x : "<<x<<" y : "<<y<<endl;
        test3_decryption(output[lay][x][y].ciphertext);
        //Torus32 s=lwePhase(output[lay][x][y].ciphertext,secret->lwe_key);

        }
    }
}
}
void h_xnor(LweSample*out,LweSample*c,float p)
{
if(p==-1)
{
LweSample*tmp_xnor=new_LweSample(in_out_params);
lweClear(tmp_xnor,in_out_params);
lweNegate(tmp_xnor,c,in_out_params);
lweCopy(out,tmp_xnor,in_out_params);
delete_LweSample(tmp_xnor);
}
else
{
LweSample*tmp_xnor=new_LweSample(in_out_params);
lweClear(tmp_xnor,in_out_params);
floatlweAddMulTo(tmp_xnor,p,c,in_out_params);
lweCopy(out,tmp_xnor,in_out_params);
delete_LweSample(tmp_xnor);
}
}
void fully_connection(C output[10],C input[432],float params[4320],float bias[10])
{
for(int lay=0;lay<10;lay++)
{
    LweSample*sum=new_LweSample(in_out_params);
    lweClear(sum,in_out_params);
    for(int i=0;i<432;i++)
    {
    LweSample*tmp=new_LweSample(in_out_params);
    lweClear(tmp,in_out_params);
    h_xnor(tmp,input[i].ciphertext,params[lay*432+i]);
    lweAddTo(sum,tmp,in_out_params);
    delete_LweSample(tmp);
    }
    LweSample*tmpf=new_LweSample(in_out_params);
    lweClear(tmpf,in_out_params);
    encrypted_data(tmpf,bias[lay]);
    lweAddTo(sum,tmpf,in_out_params);
    lweCopy(output[lay].ciphertext,sum,in_out_params);
    delete_LweSample(tmpf);
    delete_LweSample(sum);
}
}
void clear_array(LweSample*input,int n)
{
for(int i=0;i<n;i++)
{
    lweClear(input+i,in_out_params);
}
}
void hybrid_and(LweSample*result,LweSample*input,int k)
{
    if(k<0.5)
    {
        lweClear(result,in_out_params);
        LweSample*tmp=new_LweSample(in_out_params);
        encrypted_data(tmp,-1);
        lweAddTo(result,tmp,in_out_params);
        delete_LweSample(tmp);
    }
    else{
    lweCopy(result,input,in_out_params);
    }
}
void adder_calculation(LweSample*result,LweSample*x,LweSample*y,int nbits)
{
LweSample*carry=new_LweSample_array(2,in_out_params);
LweSample*tmp=new_LweSample_array(3,in_out_params);
encrypted_int_data(carry,0);
for(int i=0;i<nbits;i++)
{
    //cout<<"time : "<<i<<endl;
    bootsXOR(tmp,x+i,y+i,&secret->cloud);
    //test3_decryption(tmp);
    bootsXOR(result+i,tmp,carry,&secret->cloud);
    //test3_decryption(result+i);
    bootsAND(tmp+1,x+i,y+i,&secret->cloud);
    //test3_decryption(tmp+1);
    bootsAND(tmp+2,carry,tmp,&secret->cloud);
    //test3_decryption(tmp+2);
    bootsXOR(carry+1,tmp+1,tmp+2,&secret->cloud);
    //test3_decryption(carry+1);
    lweCopy(carry,carry+1,in_out_params);
}
delete_LweSample_array(3,tmp);
delete_LweSample_array(2,carry);
}
void fully_connected(float output[10], C input[432] , float weight[4320], float bias[10])
{
	for(int lay=0; lay<10; lay++){
		for(int ix=0; ix<432; ix++){
            float tmp=test4_decryption(input[ix].ciphertext);
			output[lay] += tmp * weight[lay * 432 + ix];
			}
		output[lay] += bias[lay];
	}
}
void full_connection_bitcalculation(int32_t final_sum[10], float params[4320],float bias[10],C input[432])
{
for(int lay=0;lay<10;lay++)
{
    LweSample*sum=new_LweSample_array(16,in_out_params);
    clear_array(sum,16);
    for(int i=0;i<432;i++)
    {
        float p=params[lay*432+i]*100;
        int int_p=p;
        int int_p_n=-int_p;
        int param_1[16]={0};
        int param_2[16]={0};
        for(int j=0;j<16;j++)
        {
            int message=(int_p>>j)&1;
            param_1[j]=message;
        }

        for(int j=0;j<16;j++)
        {
            int message=(int_p_n>>j)&1;
            param_2[j]=message;
        }

        LweSample*param_1c=new_LweSample_array(16,in_out_params);
        LweSample*param_2c=new_LweSample_array(16,in_out_params);
        LweSample*param_3c=new_LweSample_array(16,in_out_params);
        clear_array(param_1c,16);
        clear_array(param_2c,16);
        clear_array(param_3c,16);

        for(int j=0;j<16;j++)
        {
            hybrid_and(param_1c+j,input[i].ciphertext,param_1[j]);
        }
        for(int j=0;j<16;j++)
        {
            hybrid_and(param_2c+j,input[i].ciphertext,param_2[j]);
        }
        for(int j=0;j<16;j++)
        {
            bootsXOR(param_3c+j,param_1c+j,param_2c+j,&secret->cloud);
        }
        LweSample*cur_sum=new_LweSample_array(16,in_out_params);
        clear_array(cur_sum,16);
        adder_calculation(cur_sum,param_3c,sum,16);
        for(int k_number=0;k_number<16;k_number++)
        {
        lweCopy(sum+k_number,cur_sum+k_number,in_out_params);
        }
        //float M_number=test4_decryption(sum);
        //cout<<"lay : "<<lay<<" i :"<<i<<" number : "<<M_number<<endl;
        delete_LweSample_array(16,cur_sum);
        delete_LweSample_array(16,param_1c);
        delete_LweSample_array(16,param_2c);
        delete_LweSample_array(16,param_3c);
    }
    LweSample*bias_array=new_LweSample_array(16,in_out_params);
    int bias_int=100*bias[lay];
    for(int b=0;b<16;b++)
    {
    encrypted_int_data(bias_array+b,(bias_int>>b)&1);
    }
    LweSample*R=new_LweSample_array(16,in_out_params);
    adder_calculation(R,sum,bias_array,16);
    int32_t ternal=0;
    for(int N_bits=0;N_bits<15;N_bits++)
    {
        int M=fully_test_decryption(R+N_bits);
        //Array_result[lay][N_bits]=M;
        ternal+=(pow(2,N_bits)*M);
    }
    int M_16=fully_test_decryption(R+15);
    ternal+=(-pow(2,15)*M_16);
    delete_LweSample_array(16,bias_array);
    delete_LweSample_array(16,R);
    delete_LweSample_array(16,sum);
}
}


Torus32 Decrypt(LweSample*c)
{
Torus32 s=lwePhase(c,secret->lwe_key);
return s;
}
int lable_result(int32_t R[10])
{
int p=0;
int32_t m=R[0];
for(int i=1;i<10;i++)
{
    if(R[i]>m)
    {
        p=i;
        m=R[i];
    }
}

return p;
}

inline int ReverseInt(int i)
{
	unsigned char ch1, ch2, ch3, ch4;
	ch1 = i & 255;
	ch2 = (i >> 8) & 255;
	ch3 = (i >> 16) & 255;
	ch4 = (i >> 24) & 255;
	return((int)ch1 << 24) + ((int)ch2 << 16) + ((int)ch3 << 8) + ch4;
}

inline void read_Mnist_Label(string filename, vector<float>&labels)
{
	ifstream file(filename, ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int number_of_images = 0;
		file.read((char*)&magic_number, sizeof(magic_number));
		file.read((char*)&number_of_images, sizeof(number_of_images));
		magic_number = ReverseInt(magic_number);
		number_of_images = ReverseInt(number_of_images);
		cout << "magic number = " << magic_number << endl;
		cout << "number of images = " << number_of_images << endl;

		for (int i = 0; i < number_of_images; i++)
		{
			unsigned char label = 0;
			file.read((char*)&label, sizeof(label));
			labels.push_back((double)label);
		}

	}
}

inline void read_Mnist_Images(string filename, vector<vector<float>>&images)
{
	ifstream file(filename, ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int number_of_images = 0;
		int n_rows = 0;
		int n_cols = 0;
		//unsigned char label;
		file.read((char*)&magic_number, sizeof(magic_number));
		file.read((char*)&number_of_images, sizeof(number_of_images));
		file.read((char*)&n_rows, sizeof(n_rows));
		file.read((char*)&n_cols, sizeof(n_cols));
		magic_number = ReverseInt(magic_number);
		number_of_images = ReverseInt(number_of_images);
		n_rows = ReverseInt(n_rows);
		n_cols = ReverseInt(n_cols);

		cout << "magic number = " << magic_number << endl;
		cout << "number of images = " << number_of_images << endl;
		cout << "rows = " << n_rows << endl;
		cout << "cols = " << n_cols << endl;

		for (int i = 0; i < number_of_images; i++)
		{
			vector<float>tp;
			for (int r = 0; r < n_rows; r++)
			{
				for (int c = 0; c < n_cols; c++)
				{
					unsigned char image = 0;
					file.read((char*)&image, sizeof(image));
					tp.push_back(image);
				}
			}
			images.push_back(tp);
		}
	}
}
void calculation(LweSample*result,const float p,LweSample*sample,const LweParams* params)
{

    const int32_t n = params->n;
    cout<<"time : "<<n<<endl;
    for (int32_t i = 0; i < n; ++i) {
    int32_t tmp=p*sample->a[i];
    result->a[i] = tmp;
    }
    int32_t cur=p*sample->b;
    result->b = cur;
    cout<<float(float(result->b)/float(cur))<<endl;
}
int change(int K,float number)
{
float k=K;
int tmp=number*k;
float cur=number*k;
float j=cur-tmp;
float j_t=abs(j);
if(j_t>0.5 && tmp>=0)
{
tmp++;
}
if(j_t>0.5 && tmp<=0)
{
tmp--;
}
return tmp;
}
int main()
{

     vector<float>labels;
	read_Mnist_Label("t10k-labels.idx1-ubyte", labels);
	//for (auto iter = labels.begin(); iter != labels.end(); iter++)
	//{
		///cout << *iter << " ";
	//}

	vector<vector<float>>images;
	read_Mnist_Images("t10k-images.idx3-ubyte", images);
	cout<<"the col size is "<<images[0].size()<<endl;
	cout<<"the total size is "<<images.size()<<endl;
int count_cre=0;
vector<int> comp;
clock_t starttime, endtime;
starttime = clock();
for(int iter=0; iter<num_running; iter++)
{
float image_value[28][28]={0};//this is the image array
int index=0;
for(int i=0; i<28; i++){
     for(int j=0; j<28; j++)
     {
        image_value[i][j] = images[iter][index++];
        image_value[i][j]=image_value[i][j]/255;
     }
}
int lable=labels[iter];//this is lable
C encrypted_image[28][28];
for(int i=0;i<28;i++)
{
    for(int j=0;j<28;j++)
    {
        encrypted_data(encrypted_image[i][j].ciphertext,image_value[i][j]);
    }
}
cout<<"finish encrypt operation"<<endl;
cout<<"start convolution operations"<<endl;
float convolution_params[3][5][5]={0.1};
float c1[25]={-0.0100,  0.1600,  0.0900,  0.0500,  0.6200,0.0300,  0.2300, -0.0200,  0.1000,  0.3100,0.1700,  0.5500,  0.4800,  0.5800, -0.0600,
0.1800,  0.8600,  0.9000,  0.5500, -1.6800,-2.4200, -2.8400, -2.9200, -2.1000, -0.8800};
float c2[25]={0.3700,  0.0600, -0.4100, -1.1000, -0.8600,0.5500, -0.1200,  0.5800, -1.1400, -0.5500,0.4700,  0.0400,  0.5100, -3.3400,  0.2000,
0.2400,  0.4600,  0.1200, -3.1200,  1.2400,-0.0100, -0.4100, -1.8200, -1.3600,  0.2400};
float c3[25]={-0.1000,  0.0800,  0.3600, -0.0600,  1.1500,-1.6900, -1.1300,  0.2800,  0.3100,  0.4500,-1.1700, -1.6600, -2.0400, -1.3400, -0.3200,
0.2100, -0.2100, -1.1100, -1.2400, -2.0300,1.1800,  0.6900,  0.6500,  0.3200, -0.7900};
for(int j=0;j<5;j++)
{
    for(int k=0;k<5;k++)
        {
        convolution_params[0][j][k]=c1[5*j+k];
        }
}
for(int j=0;j<5;j++)
{
    for(int k=0;k<5;k++)
        {
        convolution_params[1][j][k]=c2[5*j+k];
        }
}
for(int j=0;j<5;j++)
{
    for(int k=0;k<5;k++)
        {
        convolution_params[2][j][k]=c3[5*j+k];
        }
}
float convolution_bias[3]={0.0102, 0.0008, 0.0209};
for(int i=0;i<3;i++)
{
convolution_bias[i]=1*convolution_bias[i];
};
C input_pool[3][24][24];
clock_t conv_starttime, conv_endtime;
conv_starttime = clock();
convolution_operations(input_pool,encrypted_image,convolution_params,convolution_bias);
conv_endtime = clock();
cout<<"convolutional layer's time is :"<<double(conv_endtime-conv_starttime)/CLOCKS_PER_SEC<<endl;
for(int i=0;i<28;i++)
{
    for(int j=0;j<28;j++)
    {
    delete_LweSample(encrypted_image[i][j].ciphertext);
    }
}
cout<<"finish convolution operations and start pooling"<<endl;
C input_fully[3][12][12];
clock_t pool_starttime, pool_endtime;
pool_starttime = clock();
pooling(input_fully,input_pool);
pool_endtime = clock();
cout<<"pooling layer's time is :"<<double(pool_endtime-pool_starttime)/CLOCKS_PER_SEC<<endl;
for(int i=0;i<3;i++)
{
    for(int j=0;j<24;j++)
    {
        for(int k=0;k<24;k++)
        {
        delete_LweSample(input_pool[i][j][k].ciphertext);
        }
    }
}
cout<<"finish pooling operations and start fully connection"<<endl;
cout<<"first we should flatten output pooling of ciphertext"<<endl;
C Flatten_ciphertext[432];
for(int i=0;i<3;i++)
{
    for(int j=0;j<12;j++)
    {
        for(int k=0;k<12;k++)
        {
            lweCopy(Flatten_ciphertext[i*144+12*j+k].ciphertext,input_fully[i][j][k].ciphertext,in_out_params);
        }
    }
}
//float check[10]={0};
for(int i=0;i<3;i++)
{
    for(int j=0;j<12;j++)
    {
        for(int k=0;k<12;k++)
        {
            delete_LweSample(input_fully[i][j][k].ciphertext);
        }
    }
}
float fully_params[4320]={0.0};
ifstream infile12;
infile12.open("param.txt",ios::in);
        float *ptr12=&fully_params[0];
        while(!infile12.eof())
        {
            infile12>>*ptr12;
            ptr12++;
        }
        infile12.close();
float fully_bias[10]={-1.9485, -1.1830,  2.0096, -0.3869, -1.4589,  2.3573, -1.2532,  3.1580,
        -0.2837, -1.8581};

int32_t LJJ[10]={0};
clock_t fully_starttime,fully_endtime;
fully_starttime = clock();
full_connection_bitcalculation(LJJ,fully_params,fully_bias,Flatten_ciphertext);
fully_endtime = clock();
cout<<"fully conect layer's time is :"<<double(fully_endtime-fully_starttime)/CLOCKS_PER_SEC<<endl;
for(int i=0;i<432;i++)
{
delete_LweSample(Flatten_ciphertext[i].ciphertext);
}

int F_result=lable_result(LJJ);
cout<<"F_result : "<<F_result<<endl;
comp.push_back(F_result);
if(F_result==lable)
{
count_cre++;
cout<<"inference correct"<<endl;
}
else
{
cout<<"inference error"<<endl;
}
}
endtime = clock();
cout<<"total running number is "<<num_running<<endl;
cout<<"correct number is "<<count_cre<<endl;
cout<<"total running time is "<<(double)(endtime - starttime) /CLOCKS_PER_SEC<<endl;

for(int ix=0; ix<num_running; ix++)
    cout<<"prediction is "<<comp[ix]<<" and true lable is "<<labels[ix]<<endl;

/*for(int ix=0; ix<num_running; ix++){
   for(int jx=0; jx<784; jx++){
       cout<<images[ix][jx]<<" ";
   }
   cout<<endl;
}*/

return 0;

}
