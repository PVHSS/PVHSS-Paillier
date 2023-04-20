#pragma once
#include "CAPSNARK.h"
#include "HSS.h"
typedef struct
{
    PK pk;
    CK ck;

    Vec<ZZ> I_A; // A token for evaluation

    int msg_bits;
    int degree_f;
    int msg_num;
} Para;

void KeyGen(Para &param, EK &ek0, EK &ek1);
void KeyFree(Para &param, EK &ek0, EK &ek1);

void ProbGen(Mat<ZZ> &Ix, Para param, Vec<ZZ> x);
void evaluate(bn_t y, PROOF &proof, int b, Para param, EK ekb, Mat<ZZ> Ix, vector<vector<int>> F_TEST);
bool verify(bn_t y0, bn_t y1, PROOF pi0, PROOF pi1, Para param, bool outFlag = true);