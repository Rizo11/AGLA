#include <bits/stdc++.h>

using namespace std;

class Matrix {
private:
    // array stores 2D matrix
    vector<vector<double>> array = vector<vector<double>>();
    int rows = 0;
    int columns = 0;

public:

    /* returns the size of the row */
    int getRows() const {
        return rows;
    }

    /* returns the size of the columns */
    int getColumns() const {
        return columns;
    }

    /* constructor */
    Matrix(int rows, int columns) {
        this->rows = rows;
        this->columns = columns;

        // initialize matrix with zeros
        this->array = vector<vector<double>>(rows, vector<double> (columns, 0));
    }

    Matrix() {}

    /* returns current matrix with values of matrix B */
    Matrix* operator=(Matrix B) {
        for (int i = 0; i < this->rows; ++i) {
            for (int j = 0; j < this->columns; ++j) {
                this->array[i][j] = B(i, j);
            }
        }

        return this;
    }

    /* overloading () for putting value to matrix at given row, column */
    void operator()(int row, int column, double value) {
        array[row][column] = value;
    }

    /* indexing matrix at row, col*/
    double & operator()(int row, int column) {
        return array[row][column];
    }

    /* returns difference between current matrix and matrix B */
    Matrix operator-(Matrix B) {
        Matrix result = Matrix(this->rows, this->columns);
        for (int i = 0; i < this->rows; ++i) {
            for (int j = 0; j < this->columns; ++j) {
                result(i, j, this->array[i][j] - B.array[i][j] );
            }
        }

        return result;
    }

    /* returns sum of current matrix and matrix B */
    Matrix operator+(Matrix B) {
        Matrix result = Matrix(this->rows, this->columns);
        for (int i = 0; i < this->rows; ++i) {
            for (int j = 0; j < this->columns; ++j) {
                result(i, j, B(i, j) + this->array[i][j]);
            }
        }

        return result;
    }

    /* returns product of current matrix and matrix B */
    Matrix operator*(Matrix B) {
        Matrix result = Matrix(this->rows, B.columns);

        for (int i = 0; i < this->rows; ++i) {
            for (int j = 0; j < B.columns; ++j) {
                double product = 0;
                for (int k = 0; k < this->columns; ++k) {
                    double ths = this->array[i][k];
                    double b = B(k, j);
                    product += ths*b;
                }
                result(i, j, product);
            }
        }

        return result;
    }

    /* returns product of current matrix and matrix B */
    Matrix operator*(double c) {
        Matrix result = Matrix(this->getRows(), this->getColumns());

        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getColumns(); ++j) {
                result(i, j, this->array[i][j]*c);
            }
        }

        return result;
    }

    /* prints one line of Matrix m */
    void printLine(int line) {

        for (int j = 0; j < this->columns; ++j) {
            if(abs(this->array[line][j])<0.000001)
                cout << "0.0000 ";
            else {
                cout << fixed << setprecision(4) << this->array[line][j] << " ";
            }
        }
    }

    /* prints values of matrix A*/
    void printM() {
        for (int i = 0; i < this->getRows(); ++i) {
            this->printLine(i);
            cout << endl;
        }
    }

    bool diagonallyDominant() {
        int n = this->getRows();
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += abs(this->array[i][j]);
                }
            }
            if (abs(this->array[i][i]) <= sum) {
                return false;
            }
        }

        return true;
    }

    // matrix A should be square matrix
    // A will become A^(-1)
    void inverse(Matrix &a) {
        Matrix b = Matrix(a.rows, 1);

        int i, icol, irow, j, k, l, ll, n = a.getRows(), m = b.getColumns();

        double big, dum, pivinv;

        // to track pivoting
        vector<int> indxc(n), indxr(n), ipiv(n);

        for (j=0;j<n;j++) ipiv[j]=0;

        // main loop over the columns to be reduced
        for (i=0;i<n;i++) {
            big=0.0;
            // search for a pivot
            for (j=0;j<n;j++)
                if (ipiv[j] != 1)
                    for (k=0;k<n;k++) {
                        if (ipiv[k] == 0) {
                            double element = abs(a(j, k));
                            if (element >= big) {
                                big=element;
                                irow=j;
                                icol=k;
                            }
                        }
                    }

            ++(ipiv[icol]);

            // interchange rows, if needed, to put the pivot
            //element on the diagonal

            if (irow != icol) {
                for (l=0;l<n;l++) swap(a(irow, l), a(icol, l));
                swap(b(irow, 0),b(icol, 0));
            }
            indxr[i]=irow;
            indxc[i]=icol;
            // divide the pivot row by the
            //pivot element, located at irow and icol
            if (a(icol, icol) == 0.0) {
                if (100*abs(b(icol, 0)) != 0.0) {
                    cout << "NO\n";
                    exit (0);
                }
                cout << "INF\n";
                exit (0);
            }

            pivinv=1.0/a(icol, icol);
            a(icol, icol)=1.0;
            for (l=0;l<n;l++) {
                double d = a(icol, l);
                d *= pivinv;
                a(icol, l) *= pivinv;
            }
            b(icol, 0) *= pivinv;

            // reduce rows for pivots
            for (ll=0;ll<n;ll++) {
                if (ll != icol) {
                    dum=a(ll, icol);
                    a(ll, icol)=0.0;
                    for (l=0;l<n;l++) a(ll, l) -= a(icol, l) *dum;
                    b(ll, 0) -= b(icol, 0)*dum;
                }
            }

        }

        // reverse order that the permutation
        for (l=n-1;l>=0;l--) {
            if (indxr[l] != indxc[l])
                for (k=0;k<n;k++)
                    swap(a(k, indxr[l]),a(k, indxc[l]));
        }
    }

    // returns identity
    Matrix makeIdentity(int size) {
        Matrix I = Matrix(size, size);

        for (int i = 0; i < size; ++i) {
            I(i, i) = 1;
        }

        return I;
    }

    // returns upper triangular of a, excluding the diagonal elements
    Matrix getUpperTriangular(Matrix a) {
        Matrix upper = Matrix(a.getRows(), a.getRows());
        for (int i = 0; i < a.getRows(); ++i) {
            for (int j = i+1; j < a.getColumns(); ++j) {
                upper(i, j) = a(i, j);
            }
        }

        return upper;
    }

    // returns the lower triangular part of matrix alpha, including the diagonal elements
    Matrix getLowerTriangualar(Matrix a) {
        Matrix lower = Matrix(a.getRows(), a.getRows());
        for (int i = 0; i < a.getRows(); ++i) {
            for (int j = 0; j <= i; ++j) {
                lower(i, j) = a(i, j);
            }
        }

        return lower;
    }

};

double norm(Matrix a) {
    double sum = 0;
    for (int i = 0; i < a.getRows(); ++i) {
        sum += pow(a(i, 0), 2);
    }

    return ::sqrt(sum);
}

void gauss_seidel(Matrix A, Matrix b, double epsilon) {
    int n = A.getRows();
    Matrix alpha = Matrix(n, n);
    Matrix beta = Matrix(n, 1);

    Matrix x(n, 1), x_new(n, 1);

    // check if method is applicable
    if (!A.diagonallyDominant()) {
        cout << "The method is not applicable!" << endl;
        return;
    }

    // calculate alpha and beta
    for (int i = 0; i < n; i++) {
        beta(i, 0) = b(i, 0) / A(i, i);
        for (int j = 0; j < n; j++) {
            if (i == j) {
                alpha(i, j) = 0.0;
            } else {
                alpha(i, j) = -A(i, j) / A(i, i);
            }
        }
    }

    cout << "beta: " << endl;
    beta.printM();

    cout << "alpha:" << endl;
    alpha.printM();

    Matrix B = A.getLowerTriangualar(alpha);
    cout << "B:" << endl;
    B.printM();

    Matrix C = A.getUpperTriangular(alpha);
    cout << "C:" << endl;
    C.printM();

    Matrix I = A.makeIdentity(A.getRows());

    I = I - B;
    cout << "I-B:" << endl;
    I.printM();

    I.inverse(I);
    cout << "(I-B)_-1:" << endl;
    I.printM();

    // iterative Gauss-Seidel method
    int iter = 0;

    // output results
    x = beta;

    B = A.getLowerTriangualar(A);
    C = A.getUpperTriangular(A);
    B.inverse(B);
    Matrix T = (B*(-1))*C;
    Matrix K = B*b;

    x = beta;
    double error = 0.0;
    while (true) {
        cout << "x(" << iter << "):" << endl;
        x.printM();

        x_new = T*x + K;
        error = norm(x_new - x);
        cout << "e: " << error << endl;

        if (error < epsilon) {
            cout << "x(" << ++iter << "):" << endl;
            x_new.printM();
            break;
        }

        // update approximation
        x = x_new;
        iter++;
    }
}

int main() {
    int n, m;
    double accuracy;
    cin >> n;
    Matrix A = Matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A(i, j);
        }
    }


    cin >> m;
    Matrix b = Matrix(m, 1);
    for (int i = 0; i < m; i++) {
        cin >> b(i, 0);
    }

    cin >> accuracy;

    gauss_seidel(A, b, accuracy);
    return 0;
}
