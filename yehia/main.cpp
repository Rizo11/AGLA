
//Yehia Sobeh
//y.sobeh@innopolis.university
#include<bits/stdc++.h>

using namespace std;

#define go std::ios_base::sync_with_stdio(0); cin.tie(NULL); cout.tie(NULL)
#define ll long long
#define ull unsigned long long
#define db double
#define ld long double

class ColumnVector { // This class is defined to represent a column vector
protected:
    vector<double> v;
    int n;
public:
    ColumnVector() {}

    ColumnVector(int size) { // Constructor to initialize the vector with a given size
        v.resize(size);
        n=size;
    }
    int Size(){
        return v.size();}

    double &operator[](int index) { // Operator overloading to access the vector element by index
        return v[index];
    }

    double operator[](int index) const {
        return v[index];
    }

    ColumnVector operator+(const ColumnVector &other) const { // Operator overloading to add two column vectors
        ColumnVector result(v.size());
        for (int i = 0; i < v.size(); i++) {
            result[i] = v[i] + other[i];
        }
        return result;
    }
    ColumnVector operator-(const ColumnVector &other) const { // Operator overloading to subtract two column vectors
        ColumnVector result(v.size());
        for (int i = 0; i < v.size(); i++) {
            result[i] = v[i] - other[i];
        }
        return result;
    }

    ColumnVector operator*(double scalar) const { // Operator overloading to multiply a column vector with a scalar
        ColumnVector result(v.size());
        for (int i = 0; i < v.size(); i++) {
            result[i] = v[i] * scalar;
        }
        return result;
    }

    double norm() const { // Returns the norm of the column vector
        double sum = 0.0;
        for (int i = 0; i < v.size(); i++) {
            sum += v[i] * v[i];
        }
        return sqrt(sum);
    }
    friend ostream & operator << (ostream &out,  ColumnVector &ve) { // Operator overloading to print the matrix
        for (int i = 0; i < ve.Size(); i++) {
            cout << fixed << setprecision(4);
            out << ve[i] << " ";

            out << "\n";
        }
        return out;
    }
    friend istream & operator >> (istream &in,  ColumnVector  &ve)
    {
        for (int i = 0; i < ve.Size(); i++) {
            in >> ve[i] ;


        }
        return in;
    }
};


class Matrix { // This class is defined to represent a matrix
protected:
    vector<vector<double>> data;
public:
    Matrix() {}

    Matrix(int rows, int cols) { // Constructor to initialize the matrix with a given number of rows and columns
        data.resize(rows, vector<double>(cols));
    }
    Matrix(int rows, int cols,Matrix &m ,Matrix I ) { // Constructor to initialize the matrix with a given number of rows and columns
        data.resize(rows, vector<double>(2*cols));
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                data[i][j]=m[i][j];
                data[i][j+cols]=I[i][j];
            }

        }


    }

    vector<double> &operator[](int index) { // Operator overloading to access a matrix element by index
        return data[index];
    }

    const vector<double> &operator[](int index) const {
        return data[index];
    }

    int rows() const { // Returns the number of rows in the matrix
        return data.size();
    }

    int cols() const { // Returns the number of columns in the matrix
        if (rows() == 0)
            return 0;
        return data[0].size();
    }

    ColumnVector
    operator*(const ColumnVector &v) const { // Operator overloading to multiply a matrix with a column vector
        ColumnVector result(data.size());
        for (int i = 0; i < data.size(); i++) {
            double sum = 0.0;
            for (int j = 0; j < data[i].size(); j++) {
                sum += data[i][j] * v[j];
            }
            result[i] = sum;
        }
        return result;
    }
    Matrix operator/( ColumnVector &m)  {

        int row = rows(), col = cols();
        int m_row = m.Size(), m_col = 1;
        Matrix result = Matrix(row, m_col);
        if (col == m_row) {

            for (int i = 0; i < row; i++) {
                for (int j = 0; j < m_col; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < col; k++) {
                        sum += data[i][k] * m[k];
                    }
                    if(abs(sum)<1e-10)
                        result.data[i][j] = 0.00;
                    else
                        result.data[i][j] = sum;
                }
            }

            return result;
        } else
            cout << "Error: the dimensional problem occurred\n";
        return result;

    }

    Matrix operator*(const Matrix &m) const {
        int row = rows(), col = cols();
        int m_row = m.rows(), m_col = m.cols();
        Matrix result = Matrix(row, m_col);
        if (col == m_row) {

            for (int i = 0; i < row; i++) {
                for (int j = 0; j < m_col; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < col; k++) {
                        sum += data[i][k] * m.data[k][j];
                    }
                    if(abs(sum)<1e-10)
                        result.data[i][j] = 0.00;
                    else
                        result.data[i][j] = sum;
                }
            }

            return result;
        } else
            cout << "Error: the dimensional problem occurred\n";
        return result;

    }

    Matrix operator+(const Matrix &m) const {
        int row = rows(), col = cols();
        Matrix result = Matrix();
        if (row == m.rows() && col == m.cols()) {
            result = Matrix(rows(), cols());


            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) {
                    result.data[i][j] = data[i][j] + m.data[i][j];
                }
            }

            return result;
        } else
            cout << "Error: the dimensional problem occurred\n";
        //Matrix result2 = Matrix(0,0);
        return result;


    }

    Matrix operator-(const Matrix &m) const {
        int row = rows(), col = cols();
        Matrix result = Matrix();
        if (row == m.rows() && col == m.cols()) {
            result = Matrix(rows(), cols());

            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) {
                    result.data[i][j] = m.data[i][j] - data[i][j];
                }
            }
            return result;
        }
        cout << "Error: the dimensional problem occurred\n";


        return (result);


    }

    Matrix T() const {
        int row = rows(), col = cols();
        Matrix result(col, row);
        for (int i = 0; i < col; i++) {
            for (int j = 0; j < row; j++) {
                result.data[i][j] = data[j][i];
            }
        }
        return result;
    }

    Matrix operator=(const Matrix &m) {

        int row = m.rows(), col = m.cols();
        // make function resize
        data.resize(row, vector<double>(col));
        Matrix result = Matrix(row, col);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {

                data[i][j] = m.data[i][j];
                result.data[i][j] = m.data[i][j];
            }
        }


        return result;


    }

    Matrix operator=( ColumnVector &m) {

        int row = m.Size(), col = 1;
        // make function resize
        data.resize(row, vector<double>(col));
        Matrix result = Matrix(row, col);

        for (int i = 0; i < row; i++) {

            data[i][0] = m[i];
            result.data[i][0] = m[i];

        }


        return result;


    }


    void input() {
        int row = rows(), col = cols();
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                cin >> data[i][j];
            }
        }
    }

    friend ostream & operator << (ostream &out, const Matrix &m) { // Operator overloading to print the matrix
        for (int i = 0; i < m.rows(); i++) {
            for (int j = 0; j < m.cols(); j++) {
                cout << fixed << setprecision(4);
                out << m[i][j] << " ";
            }
            out << "\n";
        }
        return out;
    }
    friend istream & operator >> (istream &in,  Matrix  &m)
    {
        for (int i = 0; i < m.rows(); i++) {
            for (int j = 0; j < m.cols(); j++) {
                in >> m[i][j] ;
            }

        }
        return in;
    }
};
class Identity : public Matrix{

public:
    Identity() : Matrix(){};
    Identity(int n) : Matrix(n,n){
        for(int i = 0; i < n; i++){
            data[i][i] = 1;
        }
    }

    Identity(int n,int m) : Matrix(n,m){
        for(int i = 0; i < n; i++){
            for(int j=0;j<m;j++)
                data[i][j] = 1;
        }
    }
    Identity(int n,int m,int value) : Matrix(n,m){
        for(int i = 0; i < n; i++){
            data[i][i] = value;
        }
    }

};


class Permutation   : public Identity{

public:
    Permutation()  : Identity(){};
    Permutation (int n) : Identity(n){}

    void permutation(int line1 , int line2){

        line1--;
        line2--;
        data[line1][line1] = 0;
        data[line2][line2] = 0;
        data[line1][line2] = 1;
        data[line2][line1] = 1;
    }

};




class Elimination  : public Identity{
protected:
    int n;
    int m;
    Matrix A,id;
    ColumnVector x;
public:
    int sign =1;
    int step=1;
    Elimination() : Identity(){};
    Elimination(int n1) : Identity(n1){n = n1;}
    Elimination(Matrix &M) {A = M; n = M.rows(); Identity I = Identity(n); id =(Matrix)I;}
    //  Elimination(Matrix &M) {A = M; n = M.rows(); m =M.cols();Identity I = Identity(n,m); id =(Matrix)I;}
    Elimination(Matrix &M,ColumnVector& y) {x=y;A = M; n = M.rows(); Identity I = Identity(n); id =(Matrix)I;}
    void elimination (Matrix& m ,int line1 ,int line2){
        line1--,line2--;
        data[line2][line1] = -1.0*(m[line2][line1]/m[line1][line1]);///m[line1][line1]);

    }
    void pivoting(int f =0){


        for (int k = 0; k < n; k++) {

            double max_elem = abs(A[k][k]);
            int max_row = k;
            for (int i = k + 1; i < n; i++) {
                if (abs(A[i][k]) > max_elem) {
                    max_elem = abs(A[i][k]);
                    max_row = i;
                }
            }
            if (max_elem < 1e-10) {
                continue;

            }

            if(k!=max_row){
//cout<<"step #"<<step<<": permutation\n";
                step++;
                sign*=-1;
                Permutation P = Permutation(n);
                P.permutation(k+1,max_row+1);
                A = P * A;




                if(f==1){
                    id =P* id;
                    Matrix print =Matrix(n,n,A,id);
                    //     cout<<print;
                }
                else if(f==2){
                    x = P*x;
                    //cout<<A<<x;
                }

                else {
                    // cout<<A;
                }
            }
            for (int i = k + 1; i < n; i++) {

                if(A[i][k]!=0){
                    //  cout<<"step #"<<step<<": elimination\n";
                    step++;
                    Elimination E = Elimination(n);
                    E.elimination(A,k+1,i+1);

                    A = E * A;

                    if(f==1){
                        id = E * id;
                        Matrix print =Matrix(n,n,A,id);
                        //  cout<<print;
                    }
                    else if(f==2){
                        x = E*x;
                        //     cout<<A<<x;
                    }


                    /* else {
                      cout<<A;
                     }*/
                }
            }
        }


    }
    void upper_triangle(int f =3){
        for(int i =n-1 ;i>=0;i--){

            for(int j= i-1;j>=0;j--){
                if(A[j][i]!=0){
                    //   cout<<"step #"<<step<<": elimination\n";
                    step++;
                    Elimination E = Elimination(n);
                    E.elimination(A,i+1,j+1);

                    A = E * A;

                    if(f==1){
                        id = E *id;
                        Matrix print =Matrix(n,n,A,id);
                        //        cout<<print;
                    }
                    else if(f==2){
                        x = E*x;
                        //   cout<<A<<x;
                    }
                    /*else if(f==0){
                     cout<<A;
                    }*/
                }
            }
        }
    }
    void determinant(){double det =sign;
        for(int i=0;i<n;i++){
            det*=A[i][i];

        }
        cout<<"result:\n";
        cout<<det;}

    void diagonalNormalization(int f =3){
        Matrix Diagonal_normalization = Matrix(n,n);
        for(int i=0;i<n;i++){

            Diagonal_normalization[i][i] = 1/A[i][i];
        }
        A = Diagonal_normalization * A;

        //cout<<"Diagonal normalization:\n";
        if(f==1){
            id = Diagonal_normalization *id;
            Matrix print =Matrix(n,n,A,id);
            //    cout<<print;
        }
        else if(f==2){
            x = Diagonal_normalization*x;
            //  cout<<A<<x;
        }


        else if(f==0){
            // cout<<A;
        }






    }
    Matrix make_inverse (){

        pivoting(1);
        upper_triangle(1);
        diagonalNormalization(1);
        //cout<<id;
        return id;

    }
    void inverse(){
        cout<<"result:\n";
        cout<<id;
    }
    void solve(){
        cout<<"result:\n";
        cout<<x;
    }

};


void MatricesCalc() {
    int n1, n2, n3, m1, m2, m3;
    cin >> n1 >> m1;
    Matrix A = Matrix(n1, m1);

    cin>>A;

    cin >> n2 >> m2;
    Matrix B = Matrix(n2, m2);
    cin>>B;

    cin >> n3 >> m3;
    Matrix C = Matrix(n3, m3);
    cin>>C;
    Matrix D = Matrix();
    Matrix E = Matrix();
    Matrix F = Matrix();

    Matrix G = Matrix();
    D = A + B;
    cout << D;

    E = A - B;
    cout<<E;
    F = C * A;
    cout<<F;
    G = A.T();
    cout<<G;
}
void SquareMatricesCalc() {
    int n1, n2, n3, m1, m2, m3;
    cin >> n1 ;
    Matrix A = Matrix(n1, n1);

    cin>>A;

    cin >> n2 ;
    Matrix B = Matrix(n2, n2);
    cin>>B;

    cin >> n3 ;
    Matrix C = Matrix(n3, n3);
    cin>>C;
    Matrix D = Matrix();
    Matrix E = Matrix();
    Matrix F = Matrix();

    Matrix G = Matrix();
    D = A + B;
    cout << D;

    E = A - B;
    cout<<E;
    F = C * A;
    cout<<F;
    G = A.T();
    cout<<G;
}
void E_P_MatricesCalc() {
    int n1;
    cin >> n1 ;
    Matrix A = Matrix(n1, n1);
    cin>>A;
    Identity I =Identity(3);
    cout<<I;
    Matrix B = Matrix();
    Matrix C = Matrix();
    Elimination E21 = Elimination(n1);

    E21.elimination(A,1,2);
    cout<<E21;
    B = E21*A;
    cout<<B;
    Permutation P21 = Permutation(n1);
    P21.permutation(1,2);
    cout<<P21;
    C = P21*A;
    cout<<C;



}


void Determinant_MatricesCalc(){

    int n ;
    cin>>n;
    Matrix  A = Matrix(n,n);

    cin>>A;
    Elimination E = Elimination(A);
    E.pivoting(0);
    E.determinant();

}



void Inverse_MatricesCalc(){

    int n ;
    cin>>n;
    Matrix A = Matrix(n,n);
    cin>>A;

    Identity I = Identity(n);
    Matrix id = Matrix(n,n);
    id = (Matrix)I;
    Matrix augmentedMatrix = Matrix(n ,n, A,(Matrix)id);
    cout<<"step #0: Augmented Matrix\n";
    cout<<augmentedMatrix;
    cout<<"Direct way:\n";
    Elimination E = Elimination(A);
    E.pivoting(1);
    cout<<"Way back:\n";
    E.upper_triangle(1);
    E.diagonalNormalization(1);
    E.inverse();






}




void L_E_MatricesCalc(){

    int n ,m;
    cin>>n;
    Matrix A = Matrix(n,n);
    cin>>A;
    cin>>m;
    ColumnVector x = ColumnVector(m);
    cin>>x;
    cout<<"step #0:\n";
    cout<<A<<x;
    Elimination E = Elimination(A,x);
    E.pivoting(2);
    E.upper_triangle(2);
    E.diagonalNormalization(2);
    E.solve();


}
void least_square(){
    int m;cin>>m;
    ColumnVector t = ColumnVector(m);
    ColumnVector b = ColumnVector(m);
    Matrix temp = Matrix(m,2);
    cin>>temp;
    int degre;
    cin>>degre;
    degre++;
    Matrix A = Matrix(m,degre);
    for(int i=0;i<m;i++){
        A[i][0] = 1;
        A[i][1] = temp[i][0];
        int ac =A[i][1];
        for(int j=2;j<degre;j++){
            ac*=A[i][1];
            A[i][j]= ac;
        }
        t[i] = temp[i][0];
        b[i]= temp[i][1];

    }
    cout<<"A:\n"<<A;
    Matrix A_T = A.T();
//cout<<"A_T\n"<<A_T;

    Matrix A_T_A = A_T*A;
    cout<<"A_T*A:\n"<<A_T_A;
    Elimination E = Elimination(A_T_A);
/*E.pivoting(1);
E.upper_triangle(1);
E.diagonalNormalization(1);
E.inverse();*/

    Matrix A_T_A_inverse = E.make_inverse();
    cout<<"(A_T*A)^-1:\n"<<A_T_A_inverse;
    Matrix AT_b = A_T/b;
    cout<<"A_T*b:\n"<<AT_b;
    cout<<"x~:\n"<<A_T_A_inverse*AT_b;

}






void jacobiMethod(/*Matrix& A, ColumnVector& b, double epsilon*/) {
    int n ;//= A.Size();
    cin>>n;
    Matrix A = Matrix(n,n);
    cin>>A;
    int m;
    cin>>m;
    ColumnVector b = ColumnVector(n);
    cin>>b;
    double epsilon;
    cin>>epsilon;
    Matrix alpha(n, n);
    ColumnVector beta(n);
    ColumnVector x(n);
    ColumnVector prevX(n);
    double e = 1.0; // initial error
    int iteration = 0;

    // Initialize alpha and beta
    for (int i = 0; i < n; i++) {
        double diag = A[i][i];
        if (diag == 0) {
            cout << "The method is not applicable!\n" ;
            return;
        }
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += abs(A[i][j]);
            }
        }
        if (abs(diag) <= sum) {
            cout << "The method is not applicable!\n" ;
            return;
        }
    }


    for (int i = 0; i < n; i++) {
        double diag = A[i][i];
        if (diag == 0) {
            cout << "The method is not applicable!\n" ;
            return;
        }
        alpha[i][i] = 0;
        beta[i] = b[i] / diag;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                alpha[i][j] = -A[i][j] / diag;
            }
        }
    }

    cout<<"alpha:\n"<<alpha<<"beta:\n"<<beta;

    // Perform Jacobi iteration
    while (e > epsilon) {
        prevX = x; // Store previous x for error calculation

        for (int i = 0; i < n; i++) {
            x[i] = 0; // Initialize x to zero
            for (int j = 0; j < n; j++) {
                x[i] += alpha[i][j] * prevX[j];
            }
            x[i] += beta[i];
        }

        // Calculate error
        e = (x - prevX).norm();
        if (iteration > 0) {
            cout << "e: " << fixed << setprecision(4) << e << "\n";
        }

        // Print iteration results
        cout << "x(" << iteration << "):" << "\n";
        cout << x ;


        iteration++;
    }

}



void seidelMethod(/*Matrix& A, ColumnVector& b, double epsilon*/) {
    int n ;//= A.Size();
    cin>>n;
    Matrix A = Matrix(n,n);
    cin>>A;
    int m;
    cin>>m;
    Matrix b = Matrix(n,1);
    cin>>b;
    double epsilon;
    cin>>epsilon;
    Matrix alpha(n, n);
    Matrix beta(n,1);
    Matrix x(n,1);
    Matrix prevX(n,1);
    double e = 1.0; // initial error
    int iteration = 0;

    // Initialize alpha and beta
    for (int i = 0; i < n; i++) {
        double diag = A[i][i];
        if (diag == 0) {
            cout << "The method is not applicable!\n" ;
            return;
        }
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += abs(A[i][j]);
            }
        }
        if (abs(diag) <= sum) {
            cout << "The method is not applicable!\n" ;
            return;
        }
    }


    for (int i = 0; i < n; i++) {
        double diag = A[i][i];
        if (diag == 0) {
            cout << "The method is not applicable!\n" ;
            return;
        }
        alpha[i][i] = 0;
        beta[i][0] = b[i][0] / diag;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                alpha[i][j] = -A[i][j] / diag;
            }
        }
    }

    cout<<"alpha:\n"<<alpha<<"beta:\n"<<beta;
    Matrix B=Matrix(n,n);Matrix C=Matrix(n,n);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<i;j++){
            B[i][j]= alpha[i][j];
        }
        for(int j=i+1;j<n;j++){
            C[i][j]= alpha[i][j];
        }
    }
    cout<<"alpha\n"<<alpha<<"B\n"<<B<<"C\n"<<C;
    Identity id = Identity(n);
    Matrix I = (Matrix)id;
    cout<<"I\n"<<I;
    Matrix I_B= B-I;
    cout<<"I-B\n"<<I_B;

    Elimination E= Elimination(I_B);
    Matrix I_B_inverse =E.make_inverse();
    cout<<"(I-B)_-1:\n"<<I_B_inverse;
    // Perform Jacobi iteration
    Matrix LL=Matrix(n,n);Matrix KK=Matrix(n,n);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<=i;j++){
            B[i][j]= A[i][j];
        }
        for(int j=i+1;j<n;j++){
            C[i][j]= A[i][j];
        }
    }
    Identity negativ =Identity(n,n,-1);
    Elimination E1 = Elimination(B);
    B =E1.make_inverse();
    B= negativ*B;
    Matrix L =B*C;
    B= negativ*B;
    Matrix K =B*b;
    x=beta;
    Matrix X_new(n,1);
    while (e > epsilon) {
        /*prevX = x; // Store previous x for error calculation

        for (int i = 0; i < n; i++) {
            x[i] = 0; // Initialize x to zero
            for (int j = 0; j < n; j++) {
                x[i] += alpha[i][j] * prevX[j];
            }
            x[i] += beta[i];
            prevX[i]= x[i];
        }*/

        //    Matrix tmp =L/x;
        X_new = L*x+K;
        // Calculate error
        //     e = (x - prevX).norm();
        ColumnVector a1(n),a2(n);
        for(int i=0;i<n;i++){
            a1[i] =X_new[i][0];
            a2[i] =x[i][0];

        }
        if (iteration > 0) {
            e = (a1-a2).norm();
            cout << "e: " << fixed << setprecision(4) << e << "\n";
        }

        // Print iteration results
        cout << "x(" << iteration << "):" << "\n";
        cout << x ;
        x= X_new;
        iteration++;
    }

}
int main() {
    go;
    seidelMethod();
    //   jacobiMethod();
    /*  int t = 6;
      //cin>>t;
      if (t == 1) {
          MatricesCalc();
          return 0;
      }
      else if (t == 2) {
          SquareMatricesCalc();
          return 0;
      }
      else if (t == 3) {
          E_P_MatricesCalc();
          return 0;
      }
      else if (t == 4) {
          Determinant_MatricesCalc();
          return 0;
      }
      else if(t == 5){
           Inverse_MatricesCalc();
          return 0;
      }
      else if(t == 6){
           L_E_MatricesCalc();
          return 0;
      }
  */



}