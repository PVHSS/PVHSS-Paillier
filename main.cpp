#include "timetester.h"

int main(int, char **)
{
    int msg_num = 5;
    int degree_f = 15;
    int cyctimes = 5;

    for (int i = 1; i <= degree_f; ++i) {
        cout << "Degree: " << i << endl;
        PVHSS_TIME_TEST(msg_num, i, cyctimes);
    }
    // PVHSS_TIME_TEST(msg_num, degree_f, cyctimes);
    // PVHSS_ACC_TEST(msg_num, degree_f);
    // SNARK_TIME_TEST(msg_num, degree_f, cyctimes);
    // HSS_TIME_TEST(msg_num, degree_f, cyctimes);
    // PAILLIER_TIME_TEST(msg_num, degree_f, cyctimes);

    return 0;
}
