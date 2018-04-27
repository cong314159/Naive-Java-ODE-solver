public class CommonMethods {
    // brute force method fail
    public static Matrix expmBruteForce(Matrix m) {
        double currentDenominator = 1;
        Matrix currentNumerator = Matrix.identity(m.getM());
        Matrix etothem = Matrix.identity(m.getM());
        for(int k = 1; k < 100; k++) {
            currentDenominator *= k;
            currentNumerator = currentNumerator.multiply(m);
            etothem = etothem.add(currentNumerator.multiply(1 / currentDenominator));
        }
//        etothem.show();
        return etothem;
    }

    public static Matrix commutator(Matrix a, Matrix b) {
        return a.multiply(b).subtract(b.multiply(a));
    }

    public static Matrix bipartiteTraceoutElectron(Matrix rho_ev, int N_v, int N_e) {
        ComplexNumber[][] rho_e_data = new ComplexNumber[N_e][N_e];
        for(int idx = 0; idx < N_e; idx++) {
            for(int idxx = 0; idxx < N_e; idxx++) {
                rho_e_data[idx][idxx] = new ComplexNumber(0, 0);
            }
        }
        for(int id = 0; id < N_e; id++) {
            for(int idx = 0; idx < N_e; idx++) {
                for(int idxx = 0; idxx < N_v; idxx++) {
                    rho_e_data[id][idx] = rho_e_data[id][idx].add(rho_ev.getData()[id*N_v + idxx][idx*N_v + idxx]);
                }
            }
        }
        return new Matrix(rho_e_data);
    }

    public static Matrix bipartiteTraceoutVibron(Matrix rho_ev, int N_v, int N_e) {
        ComplexNumber[][] rho_v_data = new ComplexNumber[N_v][N_v];
        for(int idx = 0; idx < N_v; idx++) {
            for(int idxx = 0; idx < N_v; idxx++) {
                rho_v_data[idx][idxx] = new ComplexNumber(0, 0);
            }
        }
        for(int id = 0; id < N_v; id++) {
            for(int idx = 0; id < N_v; idx++) {
                for(int idxx = 0; idxx < N_e; idxx++) {
                    rho_v_data[id][idx] = rho_v_data[id][idx].add(rho_ev.getData()[id*N_e + idxx][idx*N_e + idxx]);
                }
            }
        }
        return new Matrix(rho_v_data);
    }

    public static double errorNorm2(Matrix a, Matrix b) {
        Matrix e = a.subtract(b);
        return e.norm2();
    }

    public static class MatrixMultiThread implements Runnable {

        private static int r;
        private static int c;
        private static int l;
        private static Matrix t;
        private static Matrix a;
        private static Matrix b;

        public MatrixMultiThread(int r, int c, int l, Matrix t, Matrix a, Matrix b) {
            this.r = r;
            this.c = c;
            this.l = l;

            this.t = t;
            this.a = a;
            this.b = b;
        }

        public void run() {
            ComplexNumber ij = new ComplexNumber(0, 0);
            for(int idx = 0; idx < l; idx ++) {
                ij = ij.add(a.getData()[r][idx].multiply(b.getData()[idx][c]));
            }
            t.setData(r, c, ij);
        }
    }
}





























