public class ODE45 {
    private static double c1 = (double)0;
    private static double c2 = (double)1 / 5;
    private static double c3 = (double)3 / 10;
    private static double c4 = (double)4 / 5;
    private static double c5 = (double)8 / 9;
    private static double c6 = (double)1;
    private static double c7 = (double)1;

    private static double a21 = (double)1 / 5;

    private static double a31 = (double)3 / 40;
    private static double a32 = (double)9 / 40;

    private static double a41 = (double)44 / 45;
    private static double a42 = (double)- 56 / 15;
    private static double a43 = (double)32 / 9;

    private static double a51 = (double)19372 / 6561;
    private static double a52 = (double)-25360 / 2187;
    private static double a53 = (double)64448 / 6561;
    private static double a54 = (double)- 212 / 729;

    private static double a61 = (double)9017 / 3168;
    private static double a62 = (double)-355 / 33;
    private static double a63 = (double)46732 / 5247;
    private static double a64 = (double)49 / 176;
    private static double a65 = (double)- 5103 / 18656;

    private static double a71 = (double)35 / 384;
    private static double a72 = (double)0;
    private static double a73 = (double)500 / 1113;
    private static double a74 = (double)125 / 192;
    private static double a75 = (double)- 2187 / 6784;
    private static double a76 = (double)11 / 84;

    private static double b1 = (double)35 / 384;
    private static double b2 = (double)0;
    private static double b3 = (double)500 / 1113;
    private static double b4 = (double)125 / 192;
    private static double b5 = (double)- 2187 / 6784;
    private static double b6 = (double)11 / 84;
    private static double b7 = (double)0;

    private static double be1 = (double)5179 / 57600;
    private static double be2 = (double)0;
    private static double be3 = (double)7571 / 16695;
    private static double be4 = (double)393 / 640;
    private static double be5 = (double)- 92097 / 339200;
    private static double be6 = (double)187 / 2100;
    private static double be7 = (double)1 / 40;

    private static RHS rhs;
    private static double timeSpanStart;
    private static double timeSpanEnd;
    private static double errorTolerance;
    private double h;
    private static double stepIncrease = 1.05;
    private static double stepDecrese = 0.5;

    private double currentTime;
    private Matrix currentRho;
    private double currentPtgt;
    private Matrix FSAL;

    public ODE45(RHS rhs, Matrix rho, double[] timeSpan, double T_clk, double errorTolerance) {
        System.out.println("ODE45 starting up...");
        System.out.println("Initializing ODE45 parameter...");
        this.rhs = rhs;
        this.timeSpanStart = timeSpan[0];
        this.timeSpanEnd = timeSpan[1];
        this.currentTime = timeSpanStart;
        this.currentRho = rho;
        this.FSAL = rhs.getDrhoDt(currentTime, rho);
        this.errorTolerance = errorTolerance;
        this.h = T_clk / 1000;
        System.out.println("Solving with Dormand-Prince Method 5(4)...");
        System.out.println("...\n...\n...\n...\n");
    }

    public void ODE45SolveStep(){
        Matrix some = new Matrix(currentRho.getM(), currentRho.getN());
        Matrix k1 = FSAL;
//        k1.show();
        some = k1.multiply(a21);
        some = some.multiply(h);

        Matrix k2 = rhs.getDrhoDt(currentTime + c2 * h, currentRho.add(some));
        some = k1.multiply(a31);
        some = some.add(k2.multiply(a32));
        some = some.multiply(h);

        Matrix k3 = rhs.getDrhoDt(currentTime + c3 * h, currentRho.add(some));
        some = k1.multiply(a41);
        some = some.add(k2.multiply(a42));
        some = some.add(k3.multiply(a43));
        some = some.multiply(h);

        Matrix k4 = rhs.getDrhoDt(currentTime + c4 * h, currentRho.add(some));
        some = k1.multiply(a51);
        some = some.add(k2.multiply(a52));
        some = some.add(k3.multiply(a53));
        some = some.add(k4.multiply(a54));
        some = some.multiply(h);

        Matrix k5 = rhs.getDrhoDt(currentTime + c5 * h, currentRho.add(some));
        some = k1.multiply(a61);
        some = some.add(k2.multiply(a62));
        some = some.add(k3.multiply(a63));
        some = some.add(k4.multiply(a64));
        some = some.add(k5.multiply(a65));
        some = some.multiply(h);

        Matrix k6 = rhs.getDrhoDt(currentTime + c6 * h, currentRho.add(some));
        some = k1.multiply(a71);
        some = some.add(k2.multiply(a72));
        some = some.add(k3.multiply(a73));
        some = some.add(k4.multiply(a74));
        some = some.add(k5.multiply(a75));
        some = some.add(k6.multiply(a76));
        some = some.multiply(h);

        Matrix k7 = rhs.getDrhoDt(currentTime + c7 * h, currentRho.add(some));

        some = k1.multiply(b1);
        some = some.add(k2.multiply(b2));
        some = some.add(k3.multiply(b3));
        some = some.add(k4.multiply(b4));
        some = some.add(k5.multiply(b5));
        some = some.add(k6.multiply(b6));
        some = some.add(k7.multiply(b7));
        some = some.multiply(h);
        Matrix estimatedRho = currentRho.add(some);

        some = k1.multiply(be1);
        some = some.add(k2.multiply(be2));
        some = some.add(k3.multiply(be3));
        some = some.add(k4.multiply(be4));
        some = some.add(k5.multiply(be5));
        some = some.add(k6.multiply(be6));
        some = some.add(k7.multiply(be7));
        some = some.multiply(h);
        Matrix estimatedRho_error = currentRho.add(some);

        double error = CommonMethods.errorNorm2(estimatedRho, estimatedRho_error);
//        System.out.println("Estimated Error: " + error);
//        System.out.println("Error Tolerance is currently: " + errorTolerance);
        if(error > errorTolerance) {
            this.h = this.h * stepDecrese;
            this.ODE45SolveStep();
//            System.out.println("Step size decreased for error control...");
        }else {
            this.FSAL = k7;
            this.currentRho = estimatedRho;
            this.currentTime = currentTime + this.h;
            this.currentPtgt = CommonMethods.bipartiteTraceoutElectron(estimatedRho, 21, 3).multiply(Matrix.sigz()).trace().reOutput();
            this.h = this.h * stepIncrease;
//            System.out.println(this.h);
//            System.out.println("Step size increased for efficiency...");
//            System.out.println("...\n...\n...\n...\n");
//            CommonMethods.bipartiteTraceoutElectron(estimatedRho, 21, 3).show();
        }
    }

    public double getCurrentTime() {
        return currentTime;
    }

    public Matrix getCurrentRho() {
        return currentRho;
    }

    public double getCurrentPtgt() {
        return currentPtgt;
    }
}























