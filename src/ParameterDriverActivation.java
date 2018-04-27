public class ParameterDriverActivation {
    private double activation;

    // set a constant activation
    public ParameterDriverActivation(double activation) {
        this.activation = activation;
    }

    public double getCurrentActivation(double t) {
        return activation;
    }

    public double getCurrentActivation() {
        return activation;
    }
}
