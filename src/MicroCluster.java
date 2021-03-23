import java.util.ArrayList;
import java.util.List;

public class MicroCluster extends CFCluster {

    private long lastEditT = -1;
    private long creationTimestamp = -1;
    private double lambda;
    private long currentTimestamp;

    public MicroCluster(double[] center, long creationTimestamp, double lambda, long currentTimestamp) {
        super(center);
        this.creationTimestamp = creationTimestamp;
        this.lastEditT = creationTimestamp;
        this.lambda = lambda;
        this.currentTimestamp = currentTimestamp;
    }

    public void insert(CFCluster cfc, long timestamp) {
        N += 1;
        super.setWeight(super.getWeight() + 1);
        this.lastEditT = timestamp;

        for (int i = 0; i < cfc.LS.length; i++) {
            LS[i] += cfc.LS[i];
            SS[i] += cfc.LS[i] * cfc.LS[i];
        }
    }

    public void insert(Point p, long timestamp) {
        N += p.initWeight;
        super.setWeight(super.getWeight() + p.initWeight);
        this.lastEditT = timestamp;

        for (int i = 0; i < p.pos.length; i++) {
            LS[i] += p.pos[i];
            SS[i] += p.pos[i] * p.pos[i];
        }
    }

    public long getLastEditTimestamp() {
        return lastEditT;
    }

    public double[] calcCF2(long dt) {
        double[] cf2 = new double[SS.length];
        for (int i = 0; i < SS.length; i++) {
            cf2[i] = Math.pow(2, -lambda * dt) * SS[i];
        }
        return cf2;
    }

    public double[] calcCF1(long dt) {
        double[] cf1 = new double[LS.length];
        for (int i = 0; i < LS.length; i++) {
            cf1[i] = Math.pow(2, -lambda * dt) * LS[i];
        }
        return cf1;
    }

    @Override
    public double getWeight() {
        return getWeight(currentTimestamp);
    }

    public double getWeight(long timestamp) {
        long dt = timestamp - lastEditT;
        return (N * Math.pow(2, -lambda * dt));
    }

    public long getCreationTime() {
        return creationTimestamp;
    }

    @Override
    public double[] getCenter() {
        return getCenter(currentTimestamp);
    }

    public double[] getCenter(long timestamp) {
        long dt = timestamp - lastEditT;
        double w = getWeight(timestamp);
        double[] res = new double[LS.length];
        for (int i = 0; i < LS.length; i++) {
            res[i] = LS[i];
            res[i] *= Math.pow(2, -lambda * dt);
            res[i] /= w;
        }
        return res;
    }

    public double getRadius() {
        return getRadius(currentTimestamp);
    }

    public double getRadius(long timestamp) {
        long dt = timestamp - lastEditT;
        double[] cf1 = calcCF1(dt);
        double[] cf2 = calcCF2(dt);
        double w = getWeight(timestamp);

        double res = Math.sqrt(getVectorAbs(cf2) / w - Math.pow(getVectorAbs(cf1) / w, 2));
        return res;
    }

    public MicroCluster copy() {
        MicroCluster copy = new MicroCluster(this.LS.clone(), this.getCreationTime(), this.lambda, this.currentTimestamp);
        copy.setWeight(this.N + 1);
        copy.N = this.N;
        copy.SS = this.SS.clone();
        copy.LS = this.LS.clone();
        copy.lastEditT = this.lastEditT;
        return copy;
    }

    public CFCluster getCF() {
        CFCluster cf = copy();
        double w = getWeight();
        cf.setWeight(w);
        return cf;
    }
}

//public class MicroCluster {
//    int d; // dimension
//    int n; // number of points in micro-cluster
//    double t0;  // latest update time
//    double lambda;
//    double[] cf1;
//    double[] cf2;
//    double w;  // total weight of micro-cluster
//
//    double euclidDist(double[] p1, double[] p2) {
//        double sum = 0;
//        for (int i=0; i<p1.length; i++) {
//            sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
//        }
//        return Math.sqrt(sum);
//    }
//
//    double absolute(double[] vector) {
//        int d = vector.length;
//        double[] zero = new double[d];
//        for (int i=0; i<d; i++) {
//            zero[i] = 0;
//        }
//        return euclidDist(vector, zero);
//    }
//
//    double fadingFunction(double t_0, double t) {
//        return Math.pow(2.0, -1 * lambda * (t - t_0));
//    }
//
//    double radius() {
//        return Math.sqrt((absolute(cf2) / w) - (absolute(cf1) / w) * (absolute(cf1) / w));
//    }
//
//    // p-micro-cluster c_p = (CF^2, CF^1, w).
//    CMC cfVector(double t) {
//        double fadeFactor = fadingFunction(this.t0, t);
//        for (int i=0; i<d; i++) {
//            cf1[i] *= fadeFactor;
//            cf2[i] *= fadeFactor;
//        }
//        w = fadeFactor * w;
//
//        this.t0 = t;
//        return new CMC(cf1, cf2, w);
//    }
//
//    // Merge with point p with weight p_w.
//    void merge(double[] p, double p_w) {
//        for (int i=0; i<d; i++) {
//            cf2[i] += p[i] * p[i];
//            cf1[i] += p[i];
//        }
//        w += p_w;
//        n++;
//    }
//}
