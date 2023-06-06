#include <iostream> 
#include "GnuP.h" 
#include <cmath> 

using namespace std;

const double eps = 0.00001;  // точность

double* Gauss(double** a, double* z, int n)
{
    double* r, ** I2;
    int k;
     
    r = new double[n];
    k = 0;
    I2 = new double* [5];
    for (int i = 0; i < 5; i++)
    {
        I2[i] = new double[6];
        for (int j = 0; j < 6; j++)
        {
            if (j < 5) I2[i][j] = a[i][j];
            if (j == 5) I2[i][j] = z[i];
        }
    }
    while (k < 5)
    {
        for (int i = k; i < 5; i++)
        {
            double temp = I2[i][k];
            if (abs(temp) < eps) continue;
            for (int j = 0; j < 6; j++) I2[i][j] = I2[i][j] / temp; //Приведене главной диагонали к 1
            if (i == k)  continue;
            for (int j = 0; j < 6; j++) I2[i][j] = I2[i][j] - I2[k][j]; //Приведение элементов под главной диагональю к 0
        }
        k++;
    }
    k = 4;
    while (k >= 0)
    {
        for (int i = k; i >= 0; i--)
        {
            double temp = I2[i][k];
            if (abs(temp) < eps) continue;
            if (i == k)  continue;
            for (int j = 5; j >= 0; j--) I2[i][j] = I2[i][j] - I2[k][j] * temp; //Приведение элементов над главной диагональю к 0
        }
        k--;
    }
    k = 5;
    for (int i = 0; i < n; i++)
    {
        if (i == 2) r[i] = 0;
        else r[i] = I2[i][k];
    }
    cout << "Gauss:" << endl;
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            if (j < 5) cout << I2[i][j] << " ";
            else if (j == 5) cout << " | " << r[i];
        }
        cout << endl;
    }
    return r;
}

int main()
{
    setlocale(LC_ALL, "rus");
    double* a, * a1, * z1, * z2, ** mat, * res, sum = 0, sumt = 0, sumz = 0, sumt2 = 0, k = 0, Mt = 0, Mz = 0, ch = 0, zn1 = 0, zn2 = 0, sr = 0;
    int n = 9, N = 5;
    GnuP gp;
    double t[9] = { 2, 2.13, 2.25, 2.38, 2.5, 2.63, 2.75, 2.88, 3 };
    double z[9] = { 12.57, 16.43, 19, 22.86, 26.71, 31.86, 37, 43.43, 49.86 };
    mat = new double* [5];
    res = new double[5];
    a = new double[5];
    a1 = new double[2];
    z1 = new double[9];
    z2 = new double[9];
    for (int i = 0; i < 5; i++) {
        mat[i] = new double[5];
        res[i] = 0;
        for (int j = 0; j < 5; j++) {
            mat[i][j] = 0;
        }
    }
    for (int i = 0; i < 9; i++) {
        mat[0][0] = 5; mat[0][1] += t[i]; mat[0][2] += 0; mat[0][3] += pow(t[i], 3); mat[0][4] += pow(t[i], 4); res[0] += z[i];
        mat[1][0] = pow(t[i], 1); mat[1][1] += pow(t[i], 2); mat[1][2] += 0; mat[1][3] += pow(t[i], 4); mat[1][4] += pow(t[i], 5); res[1] += z[i] * pow(t[i], 1);
        mat[2][0] = pow(t[i], 2); mat[2][1] += pow(t[i], 3); mat[2][2] += 0; mat[2][3] += pow(t[i], 5); mat[2][4] += pow(t[i], 6); res[2] += z[i] * pow(t[i], 2);
        mat[3][0] = pow(t[i], 3); mat[3][1] += pow(t[i], 4); mat[3][2] += 0; mat[3][3] += pow(t[i], 6); mat[3][4] += pow(t[i], 7); res[3] += z[i] * pow(t[i], 3);
        mat[4][0] = pow(t[i], 4); mat[4][1] += pow(t[i], 5); mat[4][2] += 0; mat[4][3] += pow(t[i], 7); mat[4][4] += pow(t[i], 8); res[4] += z[i] * pow(t[i], 4);
    }
    cout << "mat|res:" << endl;
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            if (j < 5) cout << mat[i][j] << " ";
            else if (j == 5) cout << " | " << res[i];
        }
        cout << endl;
    }
    a = Gauss(mat, res, N);
    cout << "a: ";
    for (int i = 0; i < N; i++)cout << a[i] << " ";
    cout << endl << "Линия регрессии: у= " << a[0] << " + " << a[1] <<  " * х" << " + " << a[2] << " * х^2" << " + " << a[3] << " * х^3" << " + " << a[4] << " * х^4" << endl;
    for (int i = 0; i < N; i++)
    {
        sum += t[i] * z[i];
        sumt += t[i];
        sumz += z[i];
        sumt2 += pow(t[i], 2);
    }
    a1[1] = (n * sum - sumt * sumz) / (n * sumt2 - pow(sumt, 2));
    a1[0] = sumz / n - a1[1] * sumt / n;
    for (int i = 0; i < 9; i++)
    {
        z1[i] = a1[0] + a1[1] * t[i];
        z2[i] = a[0] + a[1] * t[i] + a[2]*pow(t[i],2) + a[3] * pow(t[i], 3) + a[4] * pow(t[i], 4);
    }
    sum = 0;
    for (int i = 0; i < 9; i++) sum += pow(z[i] - z1[i], 2);
    sr = sum / n;
    cout << "Суммарная квадратичная ошибка: " << sum << endl << "Средняя квадратичная ошибка: " << sr << endl;
    for (int i = 0; i < n; i++) 
    {
        Mt += t[i];
        Mz += z[i];
    }
    Mt /= n;
    Mz /= n;
    for (int i = 0; i < n; i++) 
    {
        ch += pow((t[i] - Mt), 2);
        zn1 += pow((t[i] - Mt), 2);
        zn2 += pow((z[i] - Mz), 2);
    }
    k = ch / sqrt(zn1 * zn2);
    cout << "Коэффициент корреляции: r=" << k << endl;
    gp.plotArrayPar(9, t, z, 3, 2, 1, "z(t)");
    gp.plotArrayPar(9, t, z1, 3, 2, 6, "z1");
    gp.plotArrayPar(9, t, z2, 3, 2, 11, "z2");
    gp.plot();
    return 0;
}
