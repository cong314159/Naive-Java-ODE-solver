
public class RHS {
    private static HamiltonianElectronic hamiltonianElectronic;
    private static Matrix H_v;
    private static Matrix H_ev;
    private static Matrix[] L;
    private static double hbar_eVs = 6.582E-16;
    private static double t;
    private static double gamma;
    private static int N_v;
    private static int N_e;

    public RHS(HamiltonianElectronic hamiltonianElectronic, Matrix H_v, Matrix H_ev, Matrix[] L, double gamma, int N_v, int N_e) {
        this.H_v = H_v;
        this.H_ev = H_ev;
        this.L = L;
        this.gamma = gamma;
        this.N_v = N_v;
        this.N_e = N_e;
        this.hamiltonianElectronic = hamiltonianElectronic;
    }

    public Matrix getDrhoDt(double t, Matrix rho) {
        Matrix H_e = hamiltonianElectronic.getHamiltonianElectronic(t);
        Matrix H_total = H_e.kron(Matrix.identity(N_v)).add(Matrix.identity(N_e).kron(H_v)).add(H_ev);
        int L_size = L.length;
        int d = N_e * N_v;
        double[][] Lindbladian_data = new double[d][d];
        for(int id = 0; id < d; id++) {
            for(int idx = 0; idx < d; idx++) {
                Lindbladian_data[id][idx] = 0;
            }
        }
        Matrix Lindbladian = new Matrix(Lindbladian_data);
        for(int id = 0; id < L_size; id++) {
            Matrix LdL_half = L[id].complexConjugate().multiply(L[id]).multiply((double) - 1 / 2);
            Matrix adder1 = L[id].multiply(rho).multiply(L[id].complexConjugate());
            Matrix adder2 = LdL_half.multiply(rho);
            Matrix adder3 = rho.multiply(LdL_half);
            Lindbladian = Lindbladian.add(adder1).add(adder2).add(adder3);
        }
        ComplexNumber i = new ComplexNumber(0, 1);
        Matrix drhodt = CommonMethods.commutator(H_total, rho);
        drhodt = drhodt.multiply(i.multiply((double)- 1 / hbar_eVs));
        drhodt = drhodt.add(Lindbladian);
//        drhodt.show();
        return drhodt;
    }
}
