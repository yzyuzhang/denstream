public class ExponentialDecayFunc {
    double base;
    double lambda;

    public ExponentialDecayFunc(double base, double lambda) {
        this.base = base;
        this.lambda = lambda;
    }

    public double decayValue(double t_0, double t) {
        return Math.pow(this.base, -1 * this.lambda * (t - t_0));
    }
}
