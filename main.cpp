// Степенной метод нахождения максимального и минимального собственного значения матрицы A
// A - СИММЕТРИЧНАЯ МАТРИЦА - запускать в режиме дебаг

#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace std;

void init_matrix(double** matrix, int N)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            cin >> matrix[i][j];
        }
    }
}

void is_matrix_symmetrical(int N, double** A)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if(A[i][j] != A[j][i])
            {
                cout << "Matrix is not symmetrical!!!" << endl;
                exit(0);
            }
        }
    }
}

void MatrVekt(int N, double **M, const double *V, double *R)
//N- размерность, M- матрица, V- вектор, R- результат
{
    for(int i=0; i<N; i++)
    {
        R[i]=0;
        for(int j=0; j<N; j++)
            R[i]+= M[i][j]*V[j];
    }
}

// вычисление первой нормы матрицы (корректно)
double matrix_norm(int n, double** a)
{
    double norma = 0;
    double* norms;
    norms = new double[n]();

    for (int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            norms[i] = norms[i] + abs(a[i][j]);
        }
    }

    for(int i = 0; i < n; i++)
    {
        if(norms[i] > norma) norma = norms[i];
    }

    delete[] norms;
    return norma;
}

// создать диагональную матрицу со значениями value на глав диагонали (корректно)
double** make_diag_matrix(int N, double value)
{
    double** matrix;
    matrix = new double*[N];
    for(int i = 0; i < N; i++)
    {
        matrix[i] = new double[N];
    }
    for(int i = 0; i < N; i++)
    {
        matrix[i][i] = value;
    }
    return matrix;
}

double find_eigenvalue(int N, double** A, double* x)
// Ax = l*x
{
    double eigenvalue;
    double* rez1;
    rez1 = new double[N]();
    MatrVekt(N, A, x, rez1);
    eigenvalue = rez1[0]/x[0];
    return eigenvalue;
}

// ВЫЧИСЛЕНИЕ ПЕРВОЙ НОРМЫ ВЕКТОРА (корректно)
double vector_norm(int n, const double* a)
{
    double norma = 0;
    for (int i = 0; i < n; i++)
        norma = norma + abs(a[i]);
    return norma;
}

double find_current_eps(int N, const double* x1, const double* x2)
{
    double* tempXX;
    tempXX = new double[N]();
    for(int i = 0; i < N; i++)
    {
        tempXX[i] = x1[i] - x2[i];
    }
    double norm = vector_norm(N, tempXX);
    return norm;
}

// вычитание матриц: E = E - A (КОРРЕКТНО)
void matrix_minus(int N, double** E, double** A)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            E[i][j] = E[i][j] - A[i][j];
        }
    }
}

int find_min_eigenvector(int N, double** A, double& min, double* x, double eps)
//N- размерность, A- матрица, min - минимальное собственное значение, eps- точность
{
    int count = 0;
    double temp_eps;
    double norm_A;
    double** E;
    double* rez;
    double* x_k;
    double norma;

    x_k = new double[N];
    for(int i = 0; i < N; i++) x_k[i] = 1; // задаём начальный вектор приближения единичным
    rez = new double[N]();
    norm_A = matrix_norm(N, A); // считаем ПЕРВУЮ норму матрицы A
    E = make_diag_matrix(N, norm_A); // создаем матрицу с нормой А на глав диагонали
    matrix_minus(N, E, A); // E = E(norm_A) - A

    // Основная итерационнная формула: x^(k+1) = (||A||_1 * E - A)*(x^k/||x^k||)
    do
    {
        MatrVekt(N, E, x_k, rez); // rez = A*x^k
        norma = vector_norm(N, x_k); // вычисляем первую норму вектора x^k
        for(int i = 0; i < N; i++)
        {
            x[i] = rez[i] * (1/norma);
        }

        temp_eps = find_current_eps(N, x, x_k);
        count++;

        for(int i = 0; i < N; i++) // переобозначение
        {
            x_k[i] = x[i];
        }
    } while (temp_eps >= eps);

    min = find_eigenvalue(N, A, x);

   return count;
}

int find_max_eigenvector(int N, double** A, double& max, double* x, double eps)
//N- размерность, A- матрица, max - минимальное собственное значение, x - собственный вектор, eps- точность
{
    int count = 0;
    double temp_eps;
    double norma;
    double* x_k;
    x_k = new double[N];
    for(int i = 0; i < N; i++) x_k[i] = 1;
    double* rez;
    rez = new double[N]();

    do
    {
        MatrVekt(N, A, x_k, rez); // A*x^k
        norma = vector_norm(N, x_k); // вычисляем первую норму вектора x^k
        for(int i = 0; i < N; i++) // основной итерационный процесс
        {
            x[i] = rez[i] * (1/norma);
        }

        count++;
        temp_eps = find_current_eps(N, x, x_k);

        for(int i = 0; i < N; i++) // переобозначение
        {
            x_k[i] = x[i];
        }
    } while (temp_eps >= eps);

    max = find_eigenvalue(N, A, x);

    delete[] x_k;
    delete[] rez;
    return count;
}

int main()
{
    double** A;
    double* x_min;
    double* x_max;
    int N;
    cin >> N; // matrix dimension
    x_min = new double[N]();
    x_max = new double[N]();
    A = new double*[N];
    for(int i = 0; i < N; i++)
    {
        A[i] = new double[N];
    }

    init_matrix(A,N);

    // is_matrix_symmetrical(N, A);

    double min;
    double max;
    int it_min; // amount of iterations
    int it_max; // amount of iterations
    double eps = 0.00000001;

    it_min = find_min_eigenvector(N, A, min, x_min, eps);
    it_max = find_max_eigenvector(N, A, max, x_max, eps);

    cout << "Iterations for min: " << it_min << endl;
    cout << "Iterations for max: " << it_max << endl;

    cout << "The min eigenvector: " << endl;
    for(int i = 0; i < N; i++ )
    {
        cout << fixed << setprecision(8) << x_min[i] << endl;
    }

    cout << "The max eigenvector: " << endl;
    for(int i = 0; i < N; i++ )
    {
        cout << fixed << setprecision(8) << x_max[i] << endl;
    }

    cout << endl << fixed << setprecision(8) << "Min eigenvalue = " << min << endl << "Max eigenvalue = " << max;

    for(int i = 0; i < N; i++)
    {
        delete A[i];
    }
    delete[] A;
    delete[] x_max;
    delete[] x_min;

    return 0;
}