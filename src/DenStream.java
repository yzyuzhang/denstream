import java.util.ArrayList;
import java.util.List;

public class DenStream {
    private double weightThreshold = 0.01;
    double lambda;
    double epsilon;
    double minWeight;
    double mu;
    double beta;

    List<Point> initBuffer;
    List<MicroCluster> p_micro_clusters;
    List<MicroCluster> o_micro_clusters;

    boolean initialized;
    private long timestamp = -1;
    long currentTimestamp = 0;
    // Tp: time period to check clustering.
    long tp;


    public DenStream(double lambda, double epsilon, double minWeight, double mu, double beta) {
        currentTimestamp = 0;

        this.lambda = lambda;
        this.epsilon = epsilon;
        this.minWeight = minWeight;
        this.mu = mu;
        this.beta = beta;

        initialized = false;
        p_micro_clusters = new ArrayList<>();
        o_micro_clusters = new ArrayList<>();
        initBuffer = new ArrayList<>();
        tp = Math.round(1 / lambda * Math.log((beta * mu) / (beta * mu - 1))) + 1;
    }

//    private List<Point> getNeighbors(Point point, List<Point> points, double eps) {
//        ArrayList<Integer> neighbourIDs = new ArrayList<Integer>();
//        for (int p = 0; p < points.size(); p++) {
//            Point npoint = points.get(p);
//            if (!npoint.covered) {
//                double dist = distance(point.toDoubleArray(), points.get(p).toDoubleArray());
//                if (dist < eps) {
//                    neighbourIDs.add(p);
//                }
//            }
//        }
//        return neighbourIDs;
//    }
//    public void initialDBScan() {
//        for (int p = 0; p < initBuffer.size(); p++) {
//            Point point = initBuffer.get(p);
//            List<Point> neighbors = getNeighbors(point, initBuffer, epsilon);
//            double neighborWeight = 0;
//            for (Point nei : neighbors) {
//                neighborWeight += nei.initWeight;
//            }
//            if (neighborWeight > beta*mu) {
//                MicroCluster mc = new MicroCluster(point, point.numAttributes(), timestamp, lambda, currentTimestamp);
//                expandCluster(mc, initBuffer, neighbourhood);
//                p_micro_clusters.add(mc);
//            }
//        }
//    }

    // Updates on receiving new point (pos, weight, true_label) at timestamp.
    public void trainOnNewPoint(double[] pos, double weight, int trueLabel, long time) {
        currentTimestamp = time;
        Point point = new Point(pos, weight, trueLabel, time);
        // Algorithm 1 Merging(p).
        // Try to merge p into its nearest p-micro-cluster c_p;
        boolean merged = false;
        if (!p_micro_clusters.isEmpty()) {
            MicroCluster x = nearestCluster(point, p_micro_clusters);
            MicroCluster xCopy = x.copy();
            xCopy.insert(point, timestamp);
            if (xCopy.getRadius(timestamp) <= epsilon) {
                x.insert(point, timestamp);
                merged = true;
            }
        }
        if (!merged && !o_micro_clusters.isEmpty()) {
            MicroCluster x = nearestCluster(point, o_micro_clusters);
            MicroCluster xCopy = x.copy();
            xCopy.insert(point, timestamp);

            if (xCopy.getRadius(timestamp) <= epsilon) {
                x.insert(point, timestamp);
                merged = true;
                // Remove x from outlier-buffer and create a new p-micro-cluster.
                if (x.getWeight() > beta * mu) {
                    o_micro_clusters.remove(x);
                    p_micro_clusters.add(x);
                }
            }
        }
        // Create a new o-micro-cluster by p and insert it into the outlier-buffer.
        if (!merged) {
            o_micro_clusters.add(new MicroCluster(point.pos, timestamp, lambda, timestamp));
        }

        // Periodic cluster removal.
        if (timestamp % tp == 0) {
            List<MicroCluster> removalList = new ArrayList<>();
            for (MicroCluster c : p_micro_clusters) {
                if (c.getWeight() < beta * mu) {
                    removalList.add(c);
                }
            }
            for (MicroCluster c : removalList) {
                p_micro_clusters.remove(c);
            }

            removalList.clear();
            for (MicroCluster c : o_micro_clusters) {
                long t0 = c.getCreationTime();
                double xsi1 = Math.pow(2, (-lambda * (timestamp - t0 + tp))) - 1;
                double xsi2 = Math.pow(2, -lambda * tp) - 1;
                double xsi = xsi1 / xsi2;
                if (c.getWeight() < xsi) {
                    removalList.add(c);
                }
            }
            for (MicroCluster c : removalList) {
                o_micro_clusters.remove(c);
            }
        }
    }

    // Generate clusters.
    public List<MicroPoint> generateClusters(long currentTimestamp) {
        List<MicroPoint> pointsToCluster = new ArrayList<>();
        for (MicroCluster mc : p_micro_clusters) {
            MicroPoint mPoint = new MicroPoint(mc.getCenter(currentTimestamp), mc.getRadius(currentTimestamp),
                    mc.getWeight(currentTimestamp));
            pointsToCluster.add(mPoint);
        }
        GenerateClusters dbscan = new GenerateClusters(epsilon, mu);
        return dbscan.cluster(pointsToCluster);
    }

    // Finds the nearest p-micro-cluster.
    public MicroCluster nearestCluster(Point p, List<MicroCluster> cl) {
        double minDist = 0;
        MicroCluster minCluster = null;
        for (MicroCluster x : cl) {
            if (minCluster == null) {
                minCluster = x;
            }
            double dist = distance(p.pos, x.getCenter(timestamp));
//            dist -= x.getRadius(timestamp);
            if (dist < minDist) {
                minDist = dist;
                minCluster = x;
            }
        }
        return minCluster;
    }

    private double distance(double[] pointA, double[] pointB) {
        double distance = 0.0;
        for (int i = 0; i < pointA.length; i++) {
            double d = pointA[i] - pointB[i];
            distance += d * d;
        }
        return Math.sqrt(distance);
    }
}
