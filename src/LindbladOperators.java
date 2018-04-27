public class LindbladOperators {
    // constants
    private static double qe = 1.602e-19;
    private static double amu2eV = 931.5e6;
    private static double epsilon0 = 8.854e-12;
    private static double c_cms = 2.997E10;
    private static double c_nms = 2.997E17;
    private static double hbar_eVs = 6.582E-16;
    private static double kb = 8.6173303e-5;
    // parameters
    private static double temp;
    private static double f_inv_cms;
    private static double m_amu;
    private static double lambda; // [meV]
    private static double gamma; // [eV]
    private static double a; //[nm]
    private static int N_v;
    private static int N_e;
    // result
    private Matrix[] L;

    public LindbladOperators(double a, double temp, double f_inv_cms, double m_amu, double lambda, double gamma, int N_v, int N_e) {
        this.a = a;
        this.temp = temp;
        this.f_inv_cms = f_inv_cms;
        this.m_amu = m_amu;
        this.lambda = lambda;
        this.gamma = gamma;
        this.N_v = N_v;
        this.N_e = N_e;
        this.calculation();
    }

    public void calculation() {
        double f = f_inv_cms * c_cms;
        double tau = 1 / f;
        double T_Lindblad = tau / 2;
        double omega = f * 2 * Math.PI;
        double hbw = hbar_eVs * omega;
        Matrix[] Lindblad = new Matrix[2];
        ComplexNumber[][] anni_data = new ComplexNumber[N_v][N_v];
        for(int idx = 0; idx < N_v; idx++) {
            for(int idxx = 0; idxx < N_v; idxx++) {
                if(idx+1 == idxx) {
                    anni_data[idx][idxx] = new ComplexNumber(Math.sqrt(idx+1));
                }else {
                    anni_data[idx][idxx] = new ComplexNumber(0);
                }
            }
        }
        Matrix anni = new Matrix(anni_data);
        Matrix anni_d = anni.complexConjugate();

        Lindblad[0] = anni.multiply(Math.sqrt(1 / T_Lindblad));
        Lindblad[0] = Matrix.identity(N_e).kron(Lindblad[0]);
        Lindblad[1] = anni_d.multiply(Math.exp(-hbw / (2 * kb * temp)) * Math.sqrt(1 / T_Lindblad));
        Lindblad[1] = Matrix.identity(N_e).kron(Lindblad[1]);

        this.L = Lindblad;
    }

    public Matrix[] getL() {
        return L;
    }
}
