import java.io.IOException;

public class Matrix {
    private int m; // number of rows
    private int n; // number of columns
    private ComplexNumber[][] data;

    // create a matrix with zeros (ComplexNumber)
    public Matrix(int m, int n) {
        this.m = m;
        this.n = n;
        this.data = new ComplexNumber[m][n];
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                this.data[a][b] = new ComplexNumber();
            }
        }
    }

    public Matrix(ComplexNumber[][] thatData) {
        this.m = thatData.length;
        this.n = thatData[0].length;
        this.data = new ComplexNumber[m][n];
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                this.data[a][b] = new ComplexNumber(thatData[a][b]);
            }
        }
    }

    public Matrix(double[][] thatData) {
        this.m = thatData.length;
        this.n = thatData[0].length;
        this.data = new ComplexNumber[m][n];
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                this.data[a][b] = new ComplexNumber(thatData[a][b]);
            }
        }
    }

    public Matrix(Matrix that) {
        this(that.data);
    }

    // Note: static can be called without an instance
    public static Matrix identity(int d) {
        Matrix I = new Matrix(d, d);
        for(int a = 0; a < d; a++) {
            I.data[a][a] = new ComplexNumber(1,0);
        }
        return I;
    }

    public static Matrix sigz() {
        double[][] sigz_data = {{-1, 0, 0}, {0, 0, 0}, {0, 0, 1}};
        return new Matrix(sigz_data);
    }

    // return the traspose of a matrix
    public Matrix transpose() {
        Matrix t = new Matrix(n, m);
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                t.data[b][a] = this.data[a][b];
            }
        }
        return t;
    }

    public Matrix complexConjugate() {
        Matrix ct = new Matrix(n, m);
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                ct.data[b][a] = this.data[a][b].complexConjugate();
            }
        }
        return ct;
    }

    public ComplexNumber trace() {
        ComplexNumber t = new ComplexNumber(0, 0);
        if(m != n) {
            throw new RuntimeException("Not a square Matrix.");
        }else {
            for(int a = 0; a < m; a++) {
                t = t.add(this.data[a][a]);
            }
        }
        return t;
    }

    public Matrix add(Matrix that) {
        Matrix A = this;
        Matrix B = that;
        Matrix C = new Matrix(m, n);
        if ((A.m != B.m) || (A.n != B.m)) {
            throw new RuntimeException("Illegal matrix dimension for plus operation.");
        }else {
            for(int a = 0; a < m; a++) {
                for(int b = 0; b < n; b++) {
                    C.data[a][b] = A.data[a][b].add(B.data[a][b]);
                }
            }
        }
        return C;
    }

    public Matrix add(double that) {
        Matrix A = this;
        double B = that;
        Matrix C = new Matrix(m, n);
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                C.data[a][b] = A.data[a][b].add(B);
            }
        }
        return C;
    }

    public Matrix subtract(Matrix that) {
        Matrix A = this;
        Matrix B = that;
        Matrix C = new Matrix(m, n);
        if((A.m != B.m) || (A.n != B.n)) {
            throw new RuntimeException("Illegal matrix dimension for minus operation.");
        }else {
            for(int a = 0; a < m; a++) {
                for(int b = 0; b < n; b++) {
                    C.data[a][b] = A.data[a][b].subtract(B.data[a][b]);
                }
            }
        }
        return C;
    }

    public Matrix subtract(double that) {
        Matrix A = this;
        double B = that;
        Matrix C = new Matrix(m, n);
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                C.data[a][b] = A.data[a][b].subtract(B);
            }
        }
        return C;
    }

    public Matrix multiplyNotUsing(Matrix that) {
        Matrix C = new Matrix(this.m, that.n);
        if(this.n != that.m) {
            throw new RuntimeException("Illegal matrix dimensions for multiply operation.");
        }else {
            Runnable[][] runnables = new Runnable[this.m][that.n];
            Thread[][] threads = new Thread[this.m][that.n];
            int l = this.n;
            for(int a = 0; a < C.m; a++) {
                for(int b = 0; b < C.n; b++) {
                    threads[a][b] = new Thread(new CommonMethods.MatrixMultiThread(a, b, l, C, this, that));
                    threads[a][b].run();
                }
            }
            for(int a = 0; a < C.m; a++) {
                for (int b = 0; b < C.n; b++) {
                    try {
                        threads[a][b].join();
                    }catch(Exception ex) {
                        System.out.println(ex);
                    }
                }
            }
        }
        return C;
    }

    public Matrix multiply(Matrix that) {
        Matrix C = new Matrix(this.m, that.n);
        if(this.n != that.m) {
            throw new RuntimeException("Illegal matrix dimensions for multiply operation.");
        }else {
            ComplexNumber[][] data = new ComplexNumber[C.m][C.n];
            int l = this.n;
            for(int a = 0; a < C.m; a++) {
                for(int b = 0; b < C.n; b++) {
                    for(int id = 0; id < l; id++) {
                        C.getData()[a][b] = C.getData()[a][b].add(this.getData()[a][id].multiply(that.getData()[id][b]));
                    }
                }
            }

        }
        return C;
    }

    public Matrix multiply(ComplexNumber that) {
        Matrix A = this;
        ComplexNumber B = that;
        Matrix C = new Matrix(m, n);
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                C.data[a][b] = A.data[a][b].multiply(B);
            }
        }
        return C;
    }

    public Matrix multiply(double that) {
        Matrix A = this;
        double B = that;
        Matrix C = new Matrix(m, n);
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                C.data[a][b] = A.data[a][b].multiply(B);
            }
        }
        return C;
    }

    public int getM() {
        return this.m;
    }

    public int getN() {
        return this.n;
    }

    public boolean equals(Matrix that) {
        if(this.m != that.m || this.n != that.n) {
            throw new RuntimeException("Matrices Dimensions don't match.");
//            return false;
        }else {
            for(int a = 0; a < m; a++) {
                for(int b = 0; b < n; b++) {
                    if (!(this.data[a][b].equals(that.data[a][b]))) {
                        return false;
                    }
                }
            }
            return true;
        }
    }

    public Matrix kron(Matrix that) {
        Matrix A = this;
        Matrix B = that;
        Matrix C = new Matrix(A.m*B.m, A.n*B.n);
        for(int a = 0; a < A.m; a++) {
            for(int b = 0; b < A.n; b++) {
                for(int c = 0; c < B.m; c++) {
                    for(int d = 0; d < B.n; d++) {
                        C.data[a*B.m+c][b*B.n+d] = A.data[a][b].multiply(B.data[c][d]);
                    }
                }
            }
        }
        return C;
    }

    public double norm2() {
        double sum = 0;
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                double re = data[a][b].reOutput();
                double im = data[a][b].imOutput();
                sum += re * re + im * im;
            }
        }
        return Math.sqrt(sum);
    }

    public ComplexNumber[][] getData() {
        return this.data;
    }

    public void setData(int r, int c, ComplexNumber C) {
        this.data[r][c] = C;
    }

    public void show() {
        System.out.println();
        System.out.println("" + m + " by " + n + " matrix");
        for(int a = 0; a < m; a++) {
            for(int b = 0; b < n; b++) {
                System.out.print("      " + this.data[a][b].toString());
            }
            System.out.println();
        }
    }
}
