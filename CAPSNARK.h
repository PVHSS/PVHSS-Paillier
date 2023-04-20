#pragma once
#include "helper.h"

typedef struct
{
    // e(g_1,g_2) = g_T
    fp12_t gT_A;
    ZZ g1_order_ZZ; // g2_order_ZZ;
    bn_t g1_order; // g2_order;

    // commitment  C = (c_0,c_1) = (hg,hat_{gh})
    ep_t h; // h = g^z
    ep_t g1;
    ep2_t g2;
    fp12_t gT;

    ep_t hat_h;   // hat_h = h^\alpha
    ep_t hat_g1;  // hat_g1 = g_1^\alpha
    ep2_t hat_g2; // hat_g2 =  g_2^\alpha
} CK;

typedef struct
{
    ep_t c;
    ep_t hat_c;
} COMMITMENT;

typedef struct
{
    bn_t e;   // e
    bn_t tau; // tau_b
    COMMITMENT D;
} PROOF;

void Biv_ComGen(CK &ck);
void Biv_Com(COMMITMENT &C, bn_t rho, CK ck, ZZ x_ZZ);
// bool Biv_ComVer(COMMITMENT C, CK ck);
bool Biv_OpenVer(CK ck, COMMITMENT C, ZZ x_ZZ, bn_t rho);

void BivPE_Gen(CK &ck);
void BivPE_Prove(PROOF &pi, int b, ZZ yb, CK ck, int &prf_key);
bool BivPE_Ver(bn_t y0, bn_t y1, PROOF pi0, PROOF pi1, CK ck, bool outFlag);

void Hash(bn_t e, fp12_t U, CK ck);
