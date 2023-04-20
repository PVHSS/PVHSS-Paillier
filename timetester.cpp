#include "timetester.h"

void PVHSS_TIME_TEST(int msg_num, int degree_f, int cyctimes)
{
    Para param;
    EK ek0, ek1;
    PROOF pi0, pi1;
    param.msg_bits = 32;
    param.degree_f = degree_f;
    param.msg_num = msg_num;

    auto *Time = new double[cyctimes];
    double time, mean, stdev;

    // Start: Setup Test
    for (int i = 0; i < cyctimes; i++)
    {
        Para param00;
        EK ek000, ek100;
        time = GetTime();
        KeyGen(param00, ek000, ek100);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Setup algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
    // End: Setup Test

    // Start: True Setup
    KeyGen(param, ek0, ek1);
    // End: True Setup

    // Start: Random Input
    Vec<ZZ> X;
    X.SetLength(msg_num);
    for (int i = 0; i < msg_num; ++i)
    {
        RandomBits(X[i], param.msg_bits);
    }
    // End: Random Input

    // Start: Input Test
    for (int i = 0; i < cyctimes; i++)
    {
        Mat<ZZ> IIx;
        time = GetTime();
        ProbGen(IIx, param, X);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Input algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
    // End: Input Test

    // Start: True Input
    Mat<ZZ> Ix;
    ProbGen(Ix, param, X);
    // End: Tru Input

    // Start: Polynomial Generation
    vector<vector<int>> F_TEST;
    Random_Func(F_TEST, msg_num, degree_f);
    // F_TEST = generate_combinations(degree_f, msg_num);
    // End: Polynomial Generation

    bn_t y0, y1;
    // Start: Evaluation 0 Test
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        evaluate(y0, pi0, 0, param, ek0, Ix, F_TEST);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Evaluate0 algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
    // End: Evaluation 0 Test

    // Start: Evaluation 1 Test
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        evaluate(y1, pi1, 1, param, ek0, Ix, F_TEST);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Evaluate1 algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
    // End: Evaluation 1 Test

    // Start: Verification Test
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        verify(y0, y1, pi0, pi1, param, false);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Verify algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
    // End: Verfication Test
}

void NATIVE_TEST(int msg_num, int degree_f, int cyctimes)
{
    Para param;
    EK ek0, ek1;
    PROOF pi0, pi1;
    param.msg_bits = 32;
    param.degree_f = degree_f;
    param.msg_num = msg_num;
    KeyGen(param, ek0, ek1);

    ZZ y_native;
    Vec<ZZ> X;
    X.SetLength(msg_num);
    for (int i = 0; i < msg_num; ++i)
    {
        RandomBits(X[i], 32);
    }

    vector<vector<int>> F_TEST;
    Random_Func(F_TEST, msg_num, degree_f);

    auto *Time = new double[cyctimes];
    double time, mean, stdev;
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        NativeEval(y_native, degree_f, msg_num, X, param.ck.g1_order_ZZ, F_TEST);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Native algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
}

void PVHSS_ACC_TEST(int msg_num, int degree_f)
{
    Para param;
    EK ek0, ek1;
    PROOF pi0, pi1;
    param.msg_bits = 32;
    param.degree_f = degree_f;
    param.msg_num = msg_num;

    double time = GetTime();
    KeyGen(param, ek0, ek1);
    time = GetTime() - time;
    cout << "KeyGen algo time: " << time * 1000 << " ms\n";
    // cout << "PK size: ";
    // cout << sizeof(param.pk.N) + GetMemory_VecZZ(param.pk.D) + sizeof(param.ck.g1) + sizeof(param.ck.g2)
    // + sizeof(param.ck.gT_A) + sizeof(param.ck.h) + sizeof(param.ck.g1_order) << endl;

    // cout << "EK0 size: ";
    // cout << GetMemory_VecZZ(ek0.d_Bsk) + GetMemory_VecZZ(param.I_A) << endl;

    // cout << "EK1 size: ";
    // cout << GetMemory_VecZZ(ek1.d_Bsk) + GetMemory_VecZZ(param.I_A) << endl;

    Vec<ZZ> X;
    X.SetLength(msg_num);
    for (int i = 0; i < msg_num; ++i)
    {
        RandomBits(X[i], param.msg_bits);
    }

    Mat<ZZ> Ix;
    time = GetTime();
    ProbGen(Ix, param, X);
    time = GetTime() - time;
    cout << "ProbGen algo time: " << time * 1000 << " ms\n";

    vector<vector<int>> F_TEST;
    Random_Func(F_TEST, msg_num, degree_f);

    bn_t y0, y1;
    time = GetTime();
    evaluate(y0, pi0, 0, param, ek0, Ix, F_TEST);
    time = GetTime() - time;
    cout << "Eval0 algo time: " << time * 1000 << " ms\n";

    time = GetTime();
    evaluate(y1, pi1, 1, param, ek1, Ix, F_TEST);
    time = GetTime() - time;
    cout << "Eval1 algo time: " << time * 1000 << " ms\n";

    time = GetTime();
    verify(y0, y1, pi0, pi1, param);
    time = GetTime() - time;
    cout << "Veri algo time: " << time * 1000 << " ms\n";

    ZZ y_native, y_eval;
    bn_t y_eval_bn;
    bn_new(y_eval_bn);
    bn_sub(y_eval_bn, y1, y0);
    bn2ZZ(y_eval, y_eval_bn);
    y_eval = y_eval % param.ck.g1_order_ZZ;
    
    NativeEval(y_native, param.degree_f, msg_num, X, param.ck.g1_order_ZZ, F_TEST);
    cout << "True result: " << y_native << endl;
    cout << "Eval result: " << y_eval << endl;
    // cout << "g1_order_ZZ: " << param.ck.g1_order_ZZ << endl;
    core_clean();
}

void SNARK_TIME_TEST(int msg_num, int degree_f, int cyctimes)
{
    Para param;
    int msg_bits = 32;

    EK ek0, ek1;
    PROOF pi0, pi1;
    param.msg_bits = msg_bits;
    param.degree_f = degree_f;
    param.msg_num = msg_num;
    KeyGen(param, ek0, ek1);

    auto *Time = new double[cyctimes];
    double time, mean, stdev;

    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        BivPE_Gen(param.ck);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Gen algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    ZZ yb;
    RandomBits(yb, msg_bits);
    yb = yb % param.ck.g1_order_ZZ;
    int prf_key = 0;
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        BivPE_Prove(pi0, 0, yb, param.ck, prf_key);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Prove algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    BivPE_Prove(pi1, 1, yb, param.ck, prf_key);
    bn_t y0, y1;
    bn_new(y0);
    bn_new(y1);
    bn_rand_mod(y0, param.ck.g1_order);
    bn_rand_mod(y1, param.ck.g1_order);

    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        BivPE_Ver(y0, y1, pi0, pi1, param.ck, false);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Verify algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
}

void HSS_TIME_TEST(int msg_num, int degree_f, int cyctimes)
{
    PK pk;
    EK ek0, ek1;
    int msg_bits = 256;

    auto *Time = new double[cyctimes];
    double time, mean, stdev;
    for (int i = 0; i < cyctimes; i++)
    {
        HSS_Free(pk, ek0, ek1);
        time = GetTime();
        HSS_Gen(pk, ek0, ek1);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "HSS_Gen algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    for (int i = 0; i < cyctimes; i++)
    {
        ZZ X_single;
        RandomBits(X_single, msg_bits);
        vec_ZZ X_Input;
        time = GetTime();
        HSS_Input(X_Input, pk, X_single);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "HSS_Input algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    for (int i = 0; i < cyctimes; i++)
    {
        ZZ X_single;
        RandomBits(X_single, msg_bits);
        vec_ZZ X_Input, X_Memory;
        HSS_Input(X_Input, pk, X_single);
        int prf_key = 0;
        time = GetTime();
        HSS_ConvertInput(X_Memory, 0, pk, ek0, X_Input, prf_key);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "HSS_ConvertInput algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    for (int i = 0; i < cyctimes; i++)
    {
        ZZ X1, X2;
        RandomBits(X1, msg_bits);
        RandomBits(X2, msg_bits);
        vec_ZZ X_Input_1, X_Input_2;
        HSS_Input(X_Input_1, pk, X1);
        HSS_Input(X_Input_2, pk, X2);

        vec_ZZ z_I;
        time = GetTime();
        HSS_AddInput(z_I, pk, X_Input_1, X_Input_2);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "HSS_AddInput algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    for (int i = 0; i < cyctimes; i++)
    {
        ZZ X1, X2;
        RandomBits(X1, msg_bits);
        RandomBits(X2, msg_bits);
        vec_ZZ X_Input_1, X_Input_2;
        HSS_Input(X_Input_1, pk, X1);
        HSS_Input(X_Input_2, pk, X2);

        vec_ZZ z_Memory, X_Memory_1, X_Memory_2;
        int prf_key = 0;
        HSS_ConvertInput(X_Memory_1, 0, pk, ek0, X_Input_1, prf_key);
        HSS_ConvertInput(X_Memory_2, 0, pk, ek0, X_Input_2, prf_key);

        time = GetTime();
        HSS_AddMemory(z_Memory, pk, X_Memory_1, X_Memory_2);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "HSS_AddMemory algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    for (int i = 0; i < cyctimes; i++)
    {
        ZZ X;
        RandomBits(X, msg_bits);
        vec_ZZ X_Input;
        HSS_Input(X_Input, pk, X);

        vec_ZZ z_Memory, X_Memory;
        int prf_key = 0;
        HSS_ConvertInput(X_Memory, 0, pk, ek0, X_Input, prf_key);

        time = GetTime();
        HSS_Mul(z_Memory, pk, ek0, X_Input, X_Memory, prf_key);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "HSS_Mul algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    vector<vector<int>> F_TEST;
    Random_Func(F_TEST, msg_num, degree_f);

    Vec<ZZ> X;
    X.SetLength(msg_num);
    for (int i = 0; i < msg_num; ++i)
    {
        RandomBits(X[i], msg_bits);
    }

    Mat<ZZ> Ix;
    Ix.SetDims(msg_num, pk.l + 1);
    for (int i = 0; i < Ix.NumRows(); ++i)
    {
        HSS_Input(Ix[i], pk, X[i]);
    }

    Vec<ZZ> y0, y1;
    for (int i = 0; i < cyctimes; i++)
    {
        int prf_key = 0;
        time = GetTime();
        HSS_Evaluate(y0, 0, Ix, pk, ek0, prf_key, F_TEST);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Evaluate0 algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    // for (int i = 0; i < cyctimes; i++)
    // {
    //     int prf_key = 0;
    //     time = GetTime();
    //     HSS_Evaluate(y1, 1, Ix, pk, ek1, prf_key, F_TEST);
    //     Time[i] = GetTime() - time;
    // }
    // DataProcess(mean, stdev, Time, cyctimes);
    // cout << "Evaluate1 algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
}

void PAILLIER_TIME_TEST(int msg_num, int degree_f, int cyctimes)
{
    PK pk;
    ZZ d;

    auto *Time = new double[cyctimes];
    double time, mean, stdev;
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        Paillier_Gen(pk, d);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Paillier Gen algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    ZZ x, ct;
    RandomBits(x, 32);
    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        Pailler_Enc(ct, pk, x);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Paillier Enc algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";

    for (int i = 0; i < cyctimes; i++)
    {
        time = GetTime();
        Pailler_Dec(x, pk, d, ct);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean, stdev, Time, cyctimes);
    cout << "Paillier Dec algo time: " << mean * 1000 << " ms  RSD: " << stdev * 100 << "%\n";
}