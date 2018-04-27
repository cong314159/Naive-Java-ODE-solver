public class ParameterDriverPolarizationSinusoid {
    private double T;

    public ParameterDriverPolarizationSinusoid(double T) {
        this.T = T;
    }

    public double getCurrentPolarization(double t) {
        return Math.sin(2 * Math.PI * t / T);
    }
}
