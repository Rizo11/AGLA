#include <bits/stdc++.h>

using namespace std;

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

};


const double EPSILON = 1e-10;

double norm(Matrix a) {
    double sum = 0;
    for (int i = 0; i < a.rows; ++i) {
        sum += pow(a(i, 0), 2);
    }

    return ::sqrt(sum);
}

void jacobi(Matrix A, Matrix b, double epsilon) {
    int n = A.getRows();
    Matrix alpha = Matrix(n, n);
    Matrix beta = Matrix(n, 1);

    Matrix x(n, 1), x_new(n, 1);

    // check if method is applicable
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += abs(A(i, j));
            }
        }
        if (abs(A(i, i)) <= sum) {
            cout << "The method is not applicable!" << endl;
            return;
        }
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

    // iterative Jacobi method
    int iter = 0;


    // output results
    cout << "alpha:" << endl;
    alpha.printM();

    cout << "beta:" << endl;
    beta.printM();
    x = beta;
    while (true) {
        // calculate new approximation
        cout << "x(" << iter << "):" << endl;
        x.printM();
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += alpha(i, j) * x(j, 0);
                }
            }
            x_new(i, 0) = beta(i, 0) + sum;
        }

        // calculate accuracy
        double error = 0.0;
        error = norm(x_new-x);

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

    jacobi(A, b, accuracy);

    return 0;
}