public class ODE23 {
    private static double c2 = (double) 1 / 2;
    private static double c3 = (double) 3 / 4;
    private static double c4 = (double) 1;

    private static double a21 = (double) 1 / 2;

    private static double a32 = (double) 3 / 4;

    private static double a41 = (double) 2 / 9;
    private static double a42 = (double) 1 / 3;
    private static double a43 = (double) 4 / 9;

    private static double b1 = (double) 2 / 9;
    private static double b2 = (double) 1 / 3;
    private static double b3 = (double) 4 / 9;
    private static double b4 = (double) 0;

    private static double be1 = (double) 7 / 24;
    private static double be2 = (double) 1 / 4;
    private static double be3 = (double) 1 / 3;
    private static double be4 = (double) 1 / 8;

    private static RHS rhs;
    private static double timeSpanStart;
    private static double timeSpanEnd;
    private static double errorTolerance;
    private double h;
    private static double stepIncrease = 1.05;
    private static double stepDecrese = 0.6;

    private double currentTime;
    private Matrix currentRho;
    private double currentPtgt;
    private Matrix FSAL;

    public ODE23(RHS rhs, Matrix rho, double[] timeSpan, double T_clk, double errorTolerance) {
        System.out.println("ODE23 starting up...");
        System.out.println("Initializing ODE23 parameter...");
        this.rhs = rhs;
        this.timeSpanStart = timeSpan[0];
        this.timeSpanEnd = timeSpan[1];
        this.currentTime = timeSpanStart;
        this.currentRho = rho;
        this.FSAL = rhs.getDrhoDt(currentTime, rho);
        this.errorTolerance = errorTolerance;
        this.h = T_clk / 1000;
        System.out.println("Solving with Bogacki-Shampine Method 3(2)...");
        System.out.println("...\n...\n...\n...\n");
    }

    public void ODE45SolveStep(){
        Matrix some = new Matrix(currentRho.getM(), currentRho.getN());
        Matrix k1 = FSAL;
        some = k1.multiply(a21);
        some = some.multiply(h);

        Matrix k2 = rhs.getDrhoDt(currentTime + c2 * h, currentRho.add(some));
        some = k2.multiply(a32);
        some = some.multiply(h);

        Matrix k3 = rhs.getDrhoDt(currentTime + c3 * h, currentRho.add(some));
        some = k1.multiply(a41);
        some = some.add(k2.multiply(a42));
        some = some.add(k3.multiply(a43));
        some = some.multiply(h);

        Matrix k4 = rhs.getDrhoDt(currentTime + c4 * h, currentRho.add(some));
        some = k1.multiply(b1).add(k2.multiply(b2)).add(k3.multiply(b3)).add(k4.multiply(b4));
        some = some.multiply(h);
        Matrix estimatedRho = currentRho.add(some);

        some = k1.multiply(be1).add(k2.multiply(be2)).add(k3.multiply(be3)).add(k4.multiply(be4));
        some = some.multiply(h);
        Matrix estimatedRho_error = currentRho.add(some);

        double error = CommonMethods.errorNorm2(estimatedRho, estimatedRho_error);
//        System.out.println("Estimated Error: " + error);
//        System.out.println("Error Tolerance is currently: " + errorTolerance);
        if(error > errorTolerance) {
            System.out.println("Step size decreased for error control...");
            this.h = this.h * stepDecrese;
            this.ODE45SolveStep();
        }else {
            this.FSAL = k4;
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























