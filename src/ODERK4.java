// Not viable

public class ODERK4 {
    private static double c2 = (double) 1 / 2;
    private static double c3 = (double) 1 / 2;
    private static double c4 = (double) 1;

    private static double a21 = (double) 1 / 2;

    private static double a32 = (double) 1 / 2;

    private static double a43 = (double) 1;

    private static double b1 = (double) 1 / 6;
    private static double b2 = (double) 2 / 6;
    private static double b3 = (double) 2 / 6;
    private static double b4 = (double) 1 / 6;

    private static RHS rhs;
    private static double timeSpanStart;
    private static double timeSpanEnd;
    private static double errorTolerance;
    private double h;

    private double currentTime;
    private Matrix currentRho;
    private double currentPtgt;

    public ODERK4(RHS rhs, Matrix rho, double[] timeSpan, double T_clk, double errorTolerance) {
        System.out.println("ODERK4 uses Runge Kutta Classical 4th order method...");
        System.out.println("Not recommended...");
        System.out.println("ODERK4 starting up...");
        System.out.println("Initializing ODERK4 parameter...");
        this.rhs = rhs;
        this.timeSpanStart = timeSpan[0];
        this.timeSpanEnd = timeSpan[1];
        this.currentTime = timeSpanStart;
        this.currentRho = rho;
        this.errorTolerance = errorTolerance;
        this.h = T_clk / 1000;
        System.out.println("Solving with Runge Kutta 4th order Method 3(2)...");
        System.out.println("...\n...\n...\n...\n");
    }

    public void ODERK4SolveStep(){
        Matrix some = new Matrix(currentRho.getM(), currentRho.getN());
        Matrix k1 = rhs.getDrhoDt(currentTime, currentRho);
        some = k1.multiply(a21);
        some = some.multiply(h);

        Matrix k2 = rhs.getDrhoDt(currentTime + c2 * h, currentRho.add(some));
        some = k2.multiply(a32);
        some = some.multiply(h);

        Matrix k3 = rhs.getDrhoDt(currentTime + c3 * h, currentRho.add(some));
        some = k3.multiply(a43);
        some = some.multiply(h);

        Matrix k4 = rhs.getDrhoDt(currentTime + c4 * h, currentRho.add(some));
        some = k1.multiply(b1).add(k2.multiply(b2)).add(k3.multiply(b3)).add(k4.multiply(b4));
        some = some.multiply(h);
        Matrix estimatedRho = currentRho.add(some);

        this.currentRho = estimatedRho;
        this.currentTime += h;
        this.currentPtgt = CommonMethods.bipartiteTraceoutElectron(estimatedRho, 21, 3).multiply(Matrix.sigz()).trace().reOutput();
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























