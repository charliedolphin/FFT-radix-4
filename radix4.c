// radix4.c
// 64 point FFT radix 4 Algorithm
// 14.09.2016 by Dmitry Chistyakov 
//
// Updates:
//      11.03.2017 - release.
//*****************************************************************************
//*****************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>

// Arrays for REAL and IM parts of spectrum
float FFT_r[64];    
float FFT_j[64];    

// Arrays for REAL and IM parts of complex exponents
float exp_r[64];
float exp_j[64];


float   rimd0, jimd0,
        rimd1, jimd1,
        rimd2, jimd2,
        rimd3, jimd3;

float pi = 3.1415926535897932384626433832795;
int n;
int step;
int i, kk;
unsigned int N = 64;
float arg;


// REAL coefficients of butterfly
int bat_s1r[4] = {1, 1, 1, 1};
int bat_s2r[4] = {1, 0, -1, 0};
int bat_s3r[4] = {1, -1, 1, -1};
int bat_s4r[4] = {1, 0, -1, 0};

// IM coefficients of butterfly
int bat_s1j[4] = {0, 0, 0, 0};
int bat_s2j[4] = {0, -1, 0, 1};
int bat_s3j[4] = {0, 0, 0, 0};
int bat_s4j[4] = {0, 1, 0, -1};

int index0, index1, index2, index3;

// Return REAL part of multiplying
float m_real (float x1, float y1, float x2, float y2)
{
    float result;
    result = x1*x2-y1*y2;
    return result;
}

// Return IM part of multiplying
float m_im (float x1, float y1, float x2, float y2)
{
    float result;
    result = x1*y2+y1*x2;
    return result;
}

int main()
{

    for (i=0; i<N; i++){
        // input data fills "fr" array
        FFT_r[i]=i; 
        FFT_j[i]=0;
    }

    // Perform first stage
    for (n=0; n<16; n++){

            index0 = 0+n;
            index1 = 16+n;
            index2 = 32+n;
            index3 = 48+n;

            rimd0 = FFT_r[index0];
            jimd0 = FFT_j[index0];
            rimd1 = FFT_r[index1];
            jimd1 = FFT_j[index1];
            rimd2 = FFT_r[index2];
            jimd2 = FFT_j[index2];
            rimd3 = FFT_r[index3];
            jimd3 = FFT_j[index3];

            FFT_r[index0] = m_real(bat_s1r[0], bat_s1j[0],   rimd0, jimd0)+
                    m_real(bat_s1r[1], bat_s1j[1], rimd1, jimd1)+
                    m_real(bat_s1r[2], bat_s1j[2], rimd2, jimd2)+
                    m_real(bat_s1r[3], bat_s1j[3], rimd3, jimd3);
            FFT_j[index0] = m_im(bat_s1r[0], bat_s1j[0], rimd0, jimd0)+
                    m_im(bat_s1r[1], bat_s1j[1], rimd1, jimd1)+
                    m_im(bat_s1r[2], bat_s1j[2], rimd2, jimd2)+
                    m_im(bat_s1r[3], bat_s1j[3], rimd3, jimd3);

            FFT_r[index1] = m_real(bat_s2r[0], bat_s2j[0], rimd0, jimd0)+
                    m_real(bat_s2r[1], bat_s2j[1], rimd1, jimd1)+
                    m_real(bat_s2r[2], bat_s2j[2],  rimd2, jimd2)+
                    m_real(bat_s2r[3], bat_s2j[3], rimd3, jimd3);
            FFT_j[index1] = m_im(bat_s2r[0], bat_s2j[0], rimd0, jimd0)+
                    m_im(bat_s2r[1], bat_s2j[1], rimd1, jimd1)+
                    m_im(bat_s2r[2], bat_s2j[2], rimd2, jimd2)+
                    m_im(bat_s2r[3], bat_s2j[3], rimd3, jimd3);

            FFT_r[index2] = m_real(bat_s3r[0], bat_s3j[0], rimd0, jimd0)+
                    m_real(bat_s3r[1], bat_s3j[1], rimd1, jimd1)+
                    m_real(bat_s3r[2], bat_s3j[2], rimd2, jimd2)+
                    m_real(bat_s3r[3], bat_s3j[3], rimd3, jimd3);
            FFT_j[index2] = m_im(bat_s3r[0], bat_s3j[0], rimd0, jimd0)+
                    m_im(bat_s3r[1], bat_s3j[1], rimd1, jimd1)+
                    m_im(bat_s3r[2], bat_s3j[2], rimd2, jimd2)+
                    m_im(bat_s3r[3], bat_s3j[3], rimd3, jimd3);

            FFT_r[index3] = m_real(bat_s4r[0], bat_s4j[0], rimd0, jimd0)+
                    m_real(bat_s4r[1], bat_s4j[1], rimd1, jimd1)+
                    m_real(bat_s4r[2], bat_s4j[2], rimd2, jimd2)+
                    m_real(bat_s4r[3], bat_s4j[3], rimd3, jimd3);
                    FFT_j[index3] = m_im(bat_s4r[0], bat_s4j[0], rimd0, jimd0)+
                    m_im(bat_s4r[1], bat_s4j[1], rimd1, jimd1)+
                    m_im(bat_s4r[2], bat_s4j[2], rimd2, jimd2)+
                    m_im(bat_s4r[3], bat_s4j[3], rimd3, jimd3);

            // Multiply results on complex exponents
            if (step!=16){
                arg = 0;
                exp_r[index0] = cos(arg);
                exp_j[index0] = sin(-arg);

                arg = (float)(2*pi*n)/(float)(N);
                exp_r[index1] = cos(arg);
                exp_j[index1] = sin(-arg);

                arg = arg * 2;
                exp_r[index2] = cos(arg);
                exp_j[index2] = sin(-arg);

                arg = arg * 3 / 2;
                exp_r[index3] = cos(arg);
                exp_j[index3] = sin(-arg);


                rimd0=FFT_r[index0];
                jimd0=FFT_j[index0];


                rimd1=FFT_r[index1];
                jimd1=FFT_j[index1];


                rimd2=FFT_r[index2];
                jimd2=FFT_j[index2];


                rimd3=FFT_r[index3];
                jimd3=FFT_j[index3];

                FFT_r[index0] = m_real(rimd0, jimd0, exp_r[index0], exp_j[index0]);
                FFT_j[index0] = m_im(rimd0, jimd0, exp_r[index0], exp_j[index0]);

                FFT_r[index1] = m_real(rimd1, jimd1, exp_r[index1], exp_j[index1]);
                FFT_j[index1] = m_im(rimd1, jimd1, exp_r[index1], exp_j[index1]);

                FFT_r[index2] = m_real(rimd2, jimd2, exp_r[index2], exp_j[index2]);
                FFT_j[index2] = m_im(rimd2, jimd2, exp_r[index2], exp_j[index2]);

                FFT_r[index3] = m_real(rimd3, jimd3, exp_r[index3], exp_j[index3]);
                FFT_j[index3] = m_im(rimd3, jimd3, exp_r[index3], exp_j[index3]);
            }
    }

    
//

    // Perform second stage
    for (kk=0; kk<N; kk+=16){
        for (step=1; step<=4; step*=4)
        {
            for (n=0; n<4; n++){

                    index0 = 0/step+n*step+kk;
                    index1 = 4/step+n*step+kk;
                    index2 = 8/step+n*step+kk;
                    index3 = 12/step+n*step+kk;

                    rimd0 = FFT_r[index0];
                    jimd0 = FFT_j[index0];
                    rimd1 = FFT_r[index1];
                    jimd1 = FFT_j[index1];
                    rimd2 = FFT_r[index2];
                    jimd2 = FFT_j[index2];
                    rimd3 = FFT_r[index3];
                    jimd3 = FFT_j[index3];

                    FFT_r[index0] = m_real(bat_s1r[0], bat_s1j[0],   rimd0, jimd0)+
                            m_real(bat_s1r[1], bat_s1j[1], rimd1, jimd1)+
                            m_real(bat_s1r[2], bat_s1j[2], rimd2, jimd2)+
                            m_real(bat_s1r[3], bat_s1j[3], rimd3, jimd3);
                    FFT_j[index0] = m_im(bat_s1r[0], bat_s1j[0], rimd0, jimd0)+
                            m_im(bat_s1r[1], bat_s1j[1], rimd1, jimd1)+
                            m_im(bat_s1r[2], bat_s1j[2], rimd2, jimd2)+
                            m_im(bat_s1r[3], bat_s1j[3], rimd3, jimd3);

                    FFT_r[index1] = m_real(bat_s2r[0], bat_s2j[0], rimd0, jimd0)+
                            m_real(bat_s2r[1], bat_s2j[1], rimd1, jimd1)+
                            m_real(bat_s2r[2], bat_s2j[2], rimd2, jimd2)+
                            m_real(bat_s2r[3], bat_s2j[3],  rimd3, jimd3);
                    FFT_j[index1] = m_im(bat_s2r[0], bat_s2j[0], rimd0, jimd0)+
                            m_im(bat_s2r[1], bat_s2j[1], rimd1, jimd1)+
                            m_im(bat_s2r[2], bat_s2j[2], rimd2, jimd2)+
                            m_im(bat_s2r[3], bat_s2j[3], rimd3, jimd3);

                    FFT_r[index2] = m_real(bat_s3r[0], bat_s3j[0], rimd0, jimd0)+
                            m_real(bat_s3r[1], bat_s3j[1], rimd1, jimd1)+
                            m_real(bat_s3r[2], bat_s3j[2], rimd2, jimd2)+
                            m_real(bat_s3r[3], bat_s3j[3], rimd3, jimd3);
                    FFT_j[index2] = m_im(bat_s3r[0], bat_s3j[0], rimd0, jimd0)+
                            m_im(bat_s3r[1], bat_s3j[1], rimd1, jimd1)+
                            m_im(bat_s3r[2], bat_s3j[2], rimd2, jimd2)+
                            m_im(bat_s3r[3], bat_s3j[3], rimd3, jimd3);

                    FFT_r[index3] = m_real(bat_s4r[0], bat_s4j[0], rimd0, jimd0)+
                            m_real(bat_s4r[1], bat_s4j[1], rimd1, jimd1)+
                            m_real(bat_s4r[2], bat_s4j[2], rimd2, jimd2)+
                            m_real(bat_s4r[3], bat_s4j[3], rimd3, jimd3);
                    FFT_j[index3] = m_im(bat_s4r[0], bat_s4j[0], rimd0, jimd0)+
                            m_im(bat_s4r[1], bat_s4j[1], rimd1, jimd1)+
                            m_im(bat_s4r[2], bat_s4j[2], rimd2, jimd2)+
                            m_im(bat_s4r[3], bat_s4j[3], rimd3, jimd3);

                    // Multiply on complex expornents
                    if (step!=4){
                        arg = 0;
                        exp_r[index0] = cos(arg);
                        exp_j[index0] = sin(-arg);

                        arg = (float)(2*pi*n)/16.0;
                        exp_r[index1] = cos(arg);
                        exp_j[index1] = sin(-arg);

                        arg = arg * 2;
                        exp_r[index2] = cos(arg);
                        exp_j[index2] = sin(-arg);

                        arg = arg * 3 / 2;
                        exp_r[index3] = cos(arg);
                        exp_j[index3] = sin(-arg);


                        rimd0=FFT_r[index0];
                        jimd0=FFT_j[index0];


                        rimd1=FFT_r[index1];
                        jimd1=FFT_j[index1];


                        rimd2=FFT_r[index2];
                        jimd2=FFT_j[index2];


                        rimd3=FFT_r[index3];
                        jimd3=FFT_j[index3];

                        FFT_r[index0] = m_real(rimd0, jimd0, exp_r[index0], exp_j[index0]);
                        FFT_j[index0] = m_im(rimd0, jimd0, exp_r[index0], exp_j[index0]);

                        FFT_r[index1] = m_real(rimd1, jimd1, exp_r[index1], exp_j[index1]);
                        FFT_j[index1] = m_im(rimd1, jimd1, exp_r[index1], exp_j[index1]);

                        FFT_r[index2] = m_real(rimd2, jimd2, exp_r[index2], exp_j[index2]);
                        FFT_j[index2] = m_im(rimd2, jimd2, exp_r[index2], exp_j[index2]);

                        FFT_r[index3] = m_real(rimd3, jimd3, exp_r[index3], exp_j[index3]);
                        FFT_j[index3] = m_im(rimd3, jimd3, exp_r[index3], exp_j[index3]);
                    }
            }
/*
            // Renew source arrays
            for (i=0; i<N/4; i++)
            {
                fr[i+kk]=FFT_r[i+kk];
                fj[i+kk]=FFT_j[i+kk];
            } */
        }
    }


    // See results
    for (i=0; i<N; i++)
    {
        printf(" %f", FFT_r[i]);
        printf("  %f\n", FFT_j[i]);
    }
    getch();
    return 0;
}
