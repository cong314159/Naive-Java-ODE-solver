public class ParameterClockingField {
    private double T; // [s]
    private double amp = 1; // [V/nm]

    public ParameterClockingField(double T, double amp) {
        this.T = T;
        this.amp = amp;
    }

    public ParameterClockingField(double T) {
        this.T = T;
    }

    public double getCurrentField(double t) {
        return amp * Math.cos((2 * Math.PI * t) / T);
    }
}
