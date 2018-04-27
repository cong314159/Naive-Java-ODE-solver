public class HamiltonianVibron {
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
    private Matrix H_v;

    public HamiltonianVibron(double a, double temp, double f_inv_cms, double m_amu, double lambda, double gamma, int N_v, int N_e) {
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
        double m = (amu2eV * m_amu) / (c_nms * c_nms);
        double f = f_inv_cms * c_cms;
        double omega = Math.PI * 2 * f;
        // data for annihilation operator;
        ComplexNumber[][] anni_data = new ComplexNumber[N_v][N_v];
        for(int a = 0; a < N_v; a++) {
            for(int b = 0; b < N_v; b++) {
                if(a+1 == b) {
                    anni_data[a][b] = new ComplexNumber(Math.sqrt(a+1));
                }else {
                    anni_data[a][b] = new ComplexNumber(0);
                }
            }
        }
        Matrix anni = new Matrix(anni_data);
        Matrix anni_d = anni.complexConjugate();
        double hbw = hbar_eVs * omega;
        Matrix ada = anni_d.multiply(anni);
        Matrix id = Matrix.identity(N_v);
        H_v = ada.add(id.multiply(0.5)).multiply(hbw);
    }

    public Matrix getH_v() {
        return H_v;
    }
}
