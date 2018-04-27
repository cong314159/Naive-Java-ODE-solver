import java.io.*;

public class ODEsolver {
    private static double qe = 1.602e-19;
    private static double amu2eV = 931.5e6;
    private static double epsilon0 = 8.854e-12;
    private static double c_cms = 2.997E10;
    private static double c_nms = 2.997E17;
    private static double hbar_eVs = 6.582E-16;
    private static double kb = 8.6173303e-5;

    public static void main(String[] args) {
        System.out.println("ODEsolver starting up...");
        System.out.println("Doing preliminary calculation...");

        double a = 1;
        double temp = 300;
        double f_inv_cms = 298;
        double m_amu = 5.58;
        double lambda = 440;
        double gamma = 0.1;
        double amp_clk = 1;
        double a_drv = 0.9;
        int N_v = 21;
        int N_e = 3;
        double f = f_inv_cms * c_cms;
        double tau = 1 / f;
        double T_clk = 10 * tau;
        double T_drv = 20 * tau;
        double n = 0.5; // [] number of clock cycles
        double T_total = n * T_clk;
        double errorTolerance = 0.01;

        System.out.println("Setting up parameter objects...");
        ParameterClockingField PCF = new ParameterClockingField(T_clk, amp_clk);
        ParameterDriverActivation PDA = new ParameterDriverActivation(a_drv);
        ParameterDriverPolarizationSinusoid PDPS = new ParameterDriverPolarizationSinusoid(T_drv);

        HamiltonianVibron hamiltonianVibron = new HamiltonianVibron(a, temp, f_inv_cms, m_amu, lambda, gamma, N_v, N_e);
        Matrix H_v = hamiltonianVibron.getH_v();
//        H_v.show();

        HamiltonianCoupling hamiltonianCoupling = new HamiltonianCoupling(a, temp, f_inv_cms, m_amu, lambda, gamma, N_v, N_e);
        Matrix H_ev = hamiltonianCoupling.getH_ev();
//        H_ev.show();

        HamiltonianElectronic hamiltonianElectronic = new HamiltonianElectronic(PCF, PDA, PDPS, a, T_clk, T_drv, gamma, amp_clk, N_v, N_e);
        Matrix H_e_initial = hamiltonianElectronic.getHamiltonianElectronic(0);
//        H_e.show();
//        HamiltonianElectronic hamiltonianElectronic = new HamiltonianElectronic();

        HamiltonianTotal hamiltonianTotal = new HamiltonianTotal(H_ev, H_e_initial, H_v, N_e, N_v);
        Matrix H_initial = hamiltonianTotal.getH();
//        H.show();

        // Calculate thermal equilibrium density matrix


        // preparation for Lindblad solver
        LindbladOperators lindbladOperators = new LindbladOperators(a, temp, f_inv_cms, m_amu, lambda, gamma, N_v, N_e);
        Matrix[] L = lindbladOperators.getL();
//        L[0].show();
//        L[1].show();
        double tstart = System.currentTimeMillis();
        // read thermal equilibrium density matrix from Matlab calculation
        System.out.println("Reading thermal equilibrium density matrix from file...");
        double[][] rho_eq_data = new double[N_e * N_v][N_e * N_v];
        try {
            String path = System.getProperty("user.dir");
            File file = new File(path + "/src/rho_eq.txt");
            FileReader fileReader = new FileReader(file);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            String line;
            int idx = 0;
            while ((line = bufferedReader.readLine()) != null) {
                String[] lineArray = line.split("  ");
                for(int idxx = 1; idxx < lineArray.length; idxx++) {
                    rho_eq_data[idx][idxx-1] = Double.parseDouble(lineArray[idxx]);
                }
                idx++;
            }
            fileReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        Matrix rho_eq = new Matrix(rho_eq_data);
//        rho_eq.show();

        System.out.println("Setting up RightHandSide function object...");
        RHS rhs = new RHS(hamiltonianElectronic, H_v, H_ev, L, gamma, N_v, N_e);

        double[] timeSpan = {0, T_total};
        System.out.println("Target time is: " + T_total);

//        ODE45
        try {
            PrintWriter writer = new PrintWriter("p_vs_t45.txt", "UTF-8");

            ODE45 ode45 = new ODE45(rhs, rho_eq, timeSpan, T_clk, errorTolerance);
            double timeOutput = ode45.getCurrentTime();
            double ptgtOutput = CommonMethods.bipartiteTraceoutElectron(ode45.getCurrentRho(), N_v, N_e).multiply(Matrix.sigz()).trace().reOutput();
            writer.println(timeOutput + " " + ptgtOutput);
            System.out.println("Percentage of finish: " + timeOutput / T_total);
            while(ode45.getCurrentTime() < T_total) {
                ode45.ODE45SolveStep();
                timeOutput = ode45.getCurrentTime();
                ptgtOutput = CommonMethods.bipartiteTraceoutElectron(ode45.getCurrentRho(), N_v, N_e).multiply(Matrix.sigz()).trace().reOutput();
                writer.println(timeOutput + " " + ptgtOutput);
                System.out.println("Percentage of finish: " + timeOutput / T_total);
            }
            writer.close();
            System.out.println(ode45.getCurrentTime());
        }catch (IOException e) {
            e.printStackTrace();
        }

        double tend = System.currentTimeMillis();

        System.out.println("Total time elapsed: " + (tend - tstart) + "s");

//        ODE23
//        try {
//            PrintWriter writer = new PrintWriter("p_vs_t23.txt", "UTF-8");
//
//            ODE23 ode23 = new ODE23(rhs, rho_eq, timeSpan, T_clk, errorTolerance);
//            double timeOutput = ode23.getCurrentTime();
//            double ptgtOutput = CommonMethods.bipartiteTraceoutElectron(ode23.getCurrentRho(), N_v, N_e).multiply(Matrix.sigz()).trace().reOutput();
//            writer.println(timeOutput + " " + ptgtOutput);
//            System.out.println("Percentage of finish: " + timeOutput / T_total);
//            while(ode23.getCurrentTime() < T_total) {
//                ode23.ODE45SolveStep();
//                timeOutput = ode23.getCurrentTime();
//                ptgtOutput = CommonMethods.bipartiteTraceoutElectron(ode23.getCurrentRho(), N_v, N_e).multiply(Matrix.sigz()).trace().reOutput();
//                writer.println(timeOutput + " " + ptgtOutput);
//                System.out.println("Percentage of finish: " + timeOutput / T_total);
//            }
//            writer.close();
//            System.out.println(ode23.getCurrentTime());
//        }catch (IOException e) {
//            e.printStackTrace();
//        }

//        ODERK4
//        try {
//            PrintWriter writer = new PrintWriter("p_vs_tRK4.txt", "UTF-8");
//
//            ODERK4 oderk4 = new ODERK4(rhs, rho_eq, timeSpan, T_clk, errorTolerance);
//            double timeOutput = oderk4.getCurrentTime();
//            double ptgtOutput = CommonMethods.bipartiteTraceoutElectron(oderk4.getCurrentRho(), N_v, N_e).multiply(Matrix.sigz()).trace().reOutput();
//            writer.println(timeOutput + " " + ptgtOutput);
//            System.out.println("Percentage of finish: " + timeOutput / T_total);
//            while(oderk4.getCurrentTime() < T_total) {
//                oderk4.ODERK4SolveStep();
//                timeOutput = oderk4.getCurrentTime();
//                ptgtOutput = CommonMethods.bipartiteTraceoutElectron(oderk4.getCurrentRho(), N_v, N_e).multiply(Matrix.sigz()).trace().reOutput();
//                writer.println(timeOutput + " " + ptgtOutput);
//                System.out.println("Percentage of finish: " + timeOutput / T_total);
//            }
//            writer.close();
//            System.out.println(oderk4.getCurrentTime());
//        }catch (IOException e) {
//            e.printStackTrace();
//        }
    }
}
