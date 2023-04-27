#include <bits/stdc++.h>
#include <cstdio>

using namespace std;

#define GNUPLOT_NAME "gnuplot -persist"

class Matrix {
public:
    // array stores 2D matrix
    vector<vector<double>> array = vector<vector<double>>();
    int rows = 0;
    int columns = 0;

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

    // returns transpose of a
    Matrix transpose(Matrix a) {
        Matrix aT = Matrix(a.getColumns(), a.getRows());

        for (int i = 0; i < aT.getRows(); ++i) {
            for (int j = 0; j < aT.getColumns(); ++j) {
                aT(i, j) = a(j, i);
            }
        }

        return aT;
    }

};



int main() {

    int m = 0;

    // read matrix length
    cin >> m;

    vector<pair<int, int>> points = vector<pair<int, int>>();
    double x, y;
    //read points
    for (int i = 0; i < m; ++i) {
        cin >> x >> y;
        points.push_back(pair<double, double>(x, y));
    }

    int degree = 0;
    // read degree
    cin >> degree;

    Matrix A = Matrix(m, degree+1);
    Matrix b = Matrix(m, 1);

    // construct matrices A and b
    for (int i = 0; i < m; ++i) {

        x = points[i].first;
        y = points[i].second;
        for (int j = 0; j < A.getColumns(); ++j) {
            A(i, j) = pow(x, j);
        }

        b(i, 0) = y;
    }

    cout << "A: " << endl;
    A.printM();

    Matrix A_T = A.transpose(A);
    cout << "A_T:" << endl;

    A_T.printM();
    cout << "A_T*A:" << endl;
    Matrix ATA = (A_T*A);
    ATA.printM();


    cout << "(A_T*A)^-1:" << endl;
    ATA.inverse(ATA);
    ATA.printM();


    cout << "A_T*b:" << endl;
    Matrix A_Tb = A_T*b;
    (A_Tb).printM();

    cout << "x~:" << endl;
    Matrix xRes = (ATA*A_Tb);
    xRes.printM();

    FILE* pipe = popen(GNUPLOT_NAME, "w");
    const double npoints = 40;
    const double step = 0.25;
    double x1 = 0.0;
    double y1 = 0.0;

    if (pipe != NULL) {

        ::fprintf(pipe, "%s\n", "plot '-' using 1:2 title 'graph' with lines");

        for (int i = 0; i < npoints; ++i) {
            x1 = -5 + i*step;
            for (int j = 0; j < xRes.getRows(); ++j) {
                y1 += pow(x1, j) * xRes(j, 0);
            }
            cout << x1 << " " << y1 << endl;
            fprintf(pipe, "%f\t%f\n", x1, y1);
            y1 = 0.0;
        }
        fflush(pipe);


        pclose(pipe);
    } else {
        cout << "Could not open file" << endl;
    }

    return 0;
}
