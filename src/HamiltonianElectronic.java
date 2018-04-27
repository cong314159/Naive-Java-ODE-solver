public class HamiltonianElectronic {
    // constants
    private static double qe = 1.602e-19;
    private static double amu2eV = 931.5e6;
    private static double epsilon0 = 8.854e-12;
    private static double c_cms = 2.997E10;
    private static double c_nms = 2.997E17;
    private static double hbar_eVs = 6.582E-16;
    private static double kb = 8.6173303e-5;
    // parameters
    private static double a;
    private static double T_clk;
    private static double T_drv;
    private static double gamma;
    private static double amp_clk;
    private static int N_v;
    private static int N_e;
    private static ParameterClockingField PCF;
    private static ParameterDriverActivation PDA;
    private static ParameterDriverPolarizationSinusoid PDPS;

    public HamiltonianElectronic(ParameterClockingField PCF, ParameterDriverActivation PDA, ParameterDriverPolarizationSinusoid PDPS, double a, double T_clk, double T_drv, double gamma, double amp_clk, int N_v, int N_e) {
        this.a = a;
        this.T_clk = T_clk;
        this.T_drv = T_drv;
        this.gamma = gamma;
        this.amp_clk = amp_clk;
        this.N_v = N_v;
        this.N_e = N_e;
        this.PCF = PCF;
        this.PDA = PDA;
        this.PDPS = PDPS;
    }

    public Matrix getHamiltonianElectronic(double t) {
        // get clocking field at time t
        double Ez_t = PCF.getCurrentField(t);
//        System.out.println(Ez_t);

        double A_drv_t = PDA.getCurrentActivation(t);
//        System.out.println(A_drv_t);

        double P_drv_t = PDPS.getCurrentPolarization(t);
//        System.out.println(P_drv_t);

        double[][] H_clk_data = {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}};
        Matrix H_clk = new Matrix(H_clk_data);
        H_clk = H_clk.multiply(-Ez_t * a / 2);
//        H_clk.show();

        double[][] H_kin_data = {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}};
        Matrix H_kin = new Matrix(H_kin_data);
        H_kin = H_kin.multiply(-gamma);
//        H_kin.show();

        double q0_drv_t = A_drv_t * qe * (1 - P_drv_t) / 2;
        double q1_drv_t = A_drv_t * qe * (1 + P_drv_t) / 2;
        double qn_drv_t = - A_drv_t * qe;

        double E1 = 0;
        E1 += qe * q0_drv_t / (4 * Math.PI * Math.sqrt(2) * a * 1e-9 * epsilon0);
        E1 += qe * q1_drv_t / (4 * Math.PI * a * 1e-9 * epsilon0);
        E1 += qe * qn_drv_t / (4 * Math.PI * Math.sqrt((double)3/2) * a * 1e-9 * epsilon0);
        E1 += - qe * q0_drv_t / (4 * Math.PI * Math.sqrt((double)3/2) * a * 1e-9 * epsilon0);
        E1 += - qe * q1_drv_t / (4 * Math.PI * Math.sqrt((double)3/2) * a * 1e-9 * epsilon0);
        E1 += - qe * qn_drv_t / (4 * Math.PI * a * 1e-9 * epsilon0);
        E1 = E1 / Math.abs(qe);

        double E0 = 0;
        E0 += qe * q1_drv_t / (4 * Math.PI * Math.sqrt(2) * a * 1e-9 * epsilon0);
        E0 += qe * q0_drv_t / (4 * Math.PI * a * 1e-9 * epsilon0);
        E0 += qe * qn_drv_t / (4 * Math.PI * Math.sqrt((double)3/2) * a * 1e-9 * epsilon0);
        E0 += - qe * q1_drv_t / (4 * Math.PI * Math.sqrt((double)3/2) * a * 1e-9 * epsilon0);
        E0 += - qe * q0_drv_t / (4 * Math.PI * Math.sqrt((double)3/2) * a * 1e-9 * epsilon0);
        E0 += - qe * qn_drv_t / (4 * Math.PI * a * 1e-9 * epsilon0);
        E0 = E0 / Math.abs(qe);
//        System.out.println(E0);

        double En = 0;

        double[][] H_nei_data = {{E0, 0, 0}, {0, En, 0}, {0, 0, E1}};
        Matrix H_nei = new Matrix(H_nei_data);

        Matrix H_e = H_clk.add(H_kin).add(H_nei);
        return H_e;
    }
}
