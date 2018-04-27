public class HamiltonianTotal {
    private static Matrix H_ev;
    private static Matrix H_e;
    private static Matrix H_v;
    private static int N_e;
    private static int N_v;
    private static Matrix H;

    public HamiltonianTotal(Matrix H_ev, Matrix H_e, Matrix H_v, int N_e, int N_v) {
        this.H_ev = H_ev;
        this.H_e = H_e;
        this.H_v = H_v;
        this.N_e = N_e;
        this.N_v = N_v;
        this.calculation();
    }

    public void calculation() {
        Matrix H_ev_component = H_ev;
        Matrix H_e_component = H_e.kron(Matrix.identity(N_v));
        Matrix H_v_component = Matrix.identity(N_e).kron(H_v);
        this.H = H_e_component.add(H_v_component).add(H_ev_component);
    }

    public Matrix getH() {
        return this.H;
    }
}
