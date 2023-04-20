#include "CAPSNARK.h"

void Biv_ComGen(CK &ck)
{
    core_init();
    ep_curve_init();
    //ep_param_set(B12_P638);
    ep_param_set(B12_P381);
    ep2_curve_init();
    ep2_curve_set_twist(RLC_EP_MTYPE);

    ep_new(ck.g1);
    ep2_new(ck.g2);
    fp12_new(ck.gT);
    ep_curve_get_gen(ck.g1);
    ep2_curve_get_gen(ck.g2);
    pp_map_oatep_k12(ck.gT, ck.g1, ck.g2);

    bn_new(ck.g1_order);
    //bn_new(ck.g2_order);
    ep_curve_get_ord(ck.g1_order);
    //ep2_curve_get_ord(ck.g2_order);
    bn2ZZ(ck.g1_order_ZZ, ck.g1_order);
    //bn2ZZ(ck.g2_order_ZZ, ck.g2_order);

    ep_new(ck.h);
    ep_new(ck.hat_h);
    ep_new(ck.hat_g1);
    ep2_new(ck.hat_g2);

    bn_t alpha;
    bn_t z;
    bn_new(alpha);
    bn_new(z);
    bn_rand_mod(alpha, ck.g1_order);
    bn_rand_mod(z, ck.g1_order);
    ep_mul_gen(ck.h, z);
    ep_mul(ck.hat_h, ck.h, alpha);
    ep_mul_gen(ck.hat_g1, alpha);
    ep2_mul_gen(ck.hat_g2, alpha);
}

void Biv_Com(COMMITMENT &C, bn_t rho, CK ck, ZZ x_ZZ)
{
    bn_null(rho);
    bn_new(rho);

    bn_t x;
    bn_new(x);
    ZZ2bn(x, x_ZZ % ck.g1_order_ZZ);

    bn_rand_mod(rho, ck.g1_order);

    ep_null(C.c);
    ep_null(C.hat_c);
    ep_new(C.c);
    ep_new(C.hat_c);

    ep_t t1, t2;
    ep_new(t1);
    ep_new(t2);

    ep_mul(t1, ck.h, rho); // h^r
    ep_mul(t2, ck.g1, x);  // g1^x
    ep_add(C.c, t1, t2);   // c= h^r g1^x

    ep_mul(t1, ck.hat_h, rho);
    ep_mul(t2, ck.hat_g1, x);
    ep_add(C.hat_c, t1, t2);
}

bool Biv_ComVer(COMMITMENT C, CK ck)
{
    fp12_t left_eq, right_eq;
    fp12_new(left_eq);
    fp12_new(right_eq);

    pp_map_oatep_k12(left_eq, C.c, ck.hat_g2);
    pp_map_oatep_k12(right_eq, C.hat_c, ck.g2);

    return fp12_cmp(left_eq, right_eq) == RLC_EQ;
}

bool Biv_OpenVer(CK ck, COMMITMENT C, ZZ x_ZZ, bn_t rho)
{
    bn_t x;
    bn_new(x);
    ZZ2bn(x, x_ZZ % ck.g1_order_ZZ);

    ep_t t1, t2, cc;
    ep_new(t1);
    ep_new(t2);
    ep_new(cc);

    ep_mul(t1, ck.h, rho);
    ep_mul(t2, ck.g1, x);
    ep_add(cc, t1, t2);

    return (ep_cmp(C.c, cc) == RLC_EQ); //&& Biv_ComVer(C, ck));
}

void BivPE_Gen(CK &ck)
{
    Biv_ComGen(ck);
}

void BivPE_Prove(PROOF &pi, int b, ZZ yb, CK ck, int &prf_key)
{
    bn_t rho;
    Biv_Com(pi.D, rho, ck, yb);

    // k = k_1 - k_0
    bn_t k;
    bn_new(k);
    PRF_bn(k, prf_key++, ck.g1_order_ZZ);

    bn_t k_b;
    bn_new(k_b);
    PRF_bn(k_b, prf_key++, ck.g1_order_ZZ);

    if (b == 1)
    {
        bn_add(k_b, k_b, k);
    }

    // U = (h^k,g)
    ep_t h_k;
    ep_new(h_k);
    ep_mul(h_k, ck.h, k);
    fp12_t u;
    fp12_new(u);
    pp_map_oatep_k12(u, h_k, ck.g2);

    bn_null(pi.e);
    bn_null(pi.tau);
    bn_new(pi.e);
    bn_new(pi.tau);
    // tau = k_b - e * r
    Hash(pi.e, u, ck);
    bn_mul(pi.tau, pi.e, rho);
    bn_sub(pi.tau, k_b, pi.tau);
}

bool BivPE_Ver(bn_t y0, bn_t y1, PROOF pi0, PROOF pi1, CK ck, bool outFlag)
{
    // Start: Check e_0 == e_1
    bool b1, b2, b3, b4;
    b1 = (bn_cmp(pi0.e, pi1.e) == RLC_EQ);
    if (!b1)
    {
        if (outFlag)
        {
            printf("******************** Check 1 ERROR ********************\n");
        }
        return false;
    }
    // End: Check e_0 == e_1
    // Start: Batched ComVer
    ep_t t11;
    ep_new(t11);
    ep_sub(t11, pi1.D.c, pi0.D.c);
    ep_t t22;
    ep_new(t22);
    ep_sub(t22, pi1.D.hat_c, pi0.D.hat_c);
    fp12_t left_eq, right_eq;
    fp12_new(left_eq);
    fp12_new(right_eq);

    pp_map_oatep_k12(left_eq, t11, ck.hat_g2);
    pp_map_oatep_k12(right_eq, t22, ck.g2);
    b2 = (fp12_cmp(left_eq, right_eq) == RLC_EQ);
    if (!b2)
    {
        if (outFlag)
        {
            printf("******************** Check 2 ERROR ********************\n");
        }
        return false;
    }
    // End: Batched ComVer

    // Start: ComVer Check
    // Biv.ComVer(ck,D_0)
    // b2 = Biv_ComVer(pi0.D, ck);
    // if (!b2)
    // {
    //     if (outFlag)
    //     {
    //         printf("******************** Check 2 ERROR ********************\n");
    //     }
    //     return false;
    // }
    // Biv.ComVer(ck,D_1)
    // b2 = Biv_ComVer(pi1.D, ck);
    // if (!b2)
    // {
    //     if (outFlag)
    //     {
    //         printf("******************** Check 2 ERROR ********************\n");
    //     }
    //     return false;
    // }
    // End: ComVer Check

    // Start: Schnorr Proof Verification
    // y = y_1 - y_0
    // tau = tau_1 - tau_0
    bn_t y, tau;
    bn_new(y);
    bn_new(tau);

    bn_sub(y, y1, y0);
    bn_sub(tau, pi1.tau, pi0.tau);
    
    // t11 = d_1 / d_0
    // e(d_1/d_0,g)
    fp12_t t2, t3;
    fp12_new(t2);
    pp_map_oatep_k12(t2, t11, ck.g2);

    // e(g,g)^{-Ay}
    fp12_new(t3);
    fp12_exp(t3, ck.gT_A, y); // g_T^Ay
    fp12_inv(t3, t3);

    // mathbb{A^e}
    fp12_mul(t3, t2, t3);
    fp12_exp(t3, t3, pi0.e); // A^e

    // e(h^tau,g)
    ep_null(t11);
    ep_new(t11);
    ep_mul(t11, ck.h, tau);
    pp_map_oatep_k12(t2, t11, ck.g2);

    // mathbb{U}
    fp12_mul(t3, t2, t3);
    
    bn_t ee;
    bn_new(ee);
    Hash(ee, t3, ck);
    b3 = (bn_cmp(ee, pi1.e) == RLC_EQ);

    if (b3)
    {
        if (outFlag)
        {
            printf("****************** Verification Passed ******************\n");
        }
        return 1;
    }
    else
    {
        if (outFlag)
        {
            printf("******************** Check 3 ERROR ********************\n");
        }
        return 0;
    }
    // End: Schnorr Proof Verification
}

void Hash(bn_t e, fp12_t U, CK ck)
{
    bn_null(e);
    bn_new(e);
    
    unsigned int seed = 31;
    unsigned int hash = 0;
    int len = 12 * RLC_FP_BYTES;

    uint8_t bin[len];
    fp12_write_bin(bin, sizeof(bin), U, 0);

    for (int i = 0; i < len; ++i)
    {
        hash = hash * seed + bin[i];
    }
    hash = hash & 0x7FFFFFFF;
    std::stringstream buffer;
    buffer << hash;
    string ss(buffer.str());
    const char *str = strdup(ss.c_str());
    bn_read_str(e, str, ss.length(), 10);
    bn_mod(e, e, ck.g1_order);
}