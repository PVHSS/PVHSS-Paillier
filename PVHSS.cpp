#include "PVHSS.h"

void KeyGen(Para &param, EK &ek0, EK &ek1)
{
    HSS_Gen(param.pk, ek0, ek1);
    Biv_ComGen(param.ck);

    ZZ A_ZZ;
    RandomBits(A_ZZ, param.pk.k - param.msg_bits * param.degree_f);

    HSS_Input(param.I_A, param.pk, A_ZZ);

    bn_t A;
    ep_t g1_A;
    bn_new(A);
    ep_new(g1_A);
    ZZ2bn(A, A_ZZ);
    ep_mul_gen(g1_A, A); // g_2^A

    fp12_new(param.gT_A);
    pp_map_oatep_k12(param.ck.gT_A, g1_A, param.ck.g2);
}

void KeyFree(Para &param, EK &ek0, EK &ek1)
{
    HSS_Free(param.pk, ek0, ek1);
    param.I_A.kill();
}


void ProbGen(Mat<ZZ> &Ix, Para param, Vec<ZZ> x)
{
    Ix.kill();
    Ix.SetDims(param.msg_num, param.pk.l + 1);
    int i;
    for (i = 0; i < Ix.NumRows(); ++i)
    {   
        HSS_Input(Ix[i], param.pk, x[i]);
    }
}

void evaluate(bn_t y, PROOF &pi, int b, Para param, EK ekb, Mat<ZZ> Ix, vector<vector<int>> F_TEST)
{
    int prf_key = 0;
    Vec<ZZ> y_b_res;
    HSS_Evaluate(y_b_res, b, Ix, param.pk, ekb, prf_key, F_TEST);

    bn_null(y);
    bn_new(y);
    ZZ2bn(y, y_b_res[0] % param.ck.g1_order_ZZ);

    // Y_b = HSS.HSS_Mul(y_b, I_A)
    Vec<ZZ> Y_ZZ;
    HSS_Mul(Y_ZZ, param.pk, ekb, param.I_A, y_b_res, prf_key);

    BivPE_Prove(pi, b, Y_ZZ[0], param.ck, prf_key);
}

bool verify(bn_t y0, bn_t y1, PROOF pi0, PROOF pi1, Para param, bool outFlag)
{
    return BivPE_Ver(y0, y1, pi0, pi1, param.ck, outFlag);
}