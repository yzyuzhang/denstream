
// MicroPoint abstracts each micro-cluster as a point.
public class MicroPoint {
    double[] pos;  // center of micro-cluster
    double radius;
    double weight;
    int clusterId;
    boolean visited;

    public MicroPoint(double[] pos, double radius, double weight) {
        this.pos = pos;
        this.radius = radius;
        this.weight = weight;
        this.visited = false;
        this.clusterId = -1;  // initially a noise point
    }

    public MicroPoint(MicroPoint p) {
        this(p.pos, p.radius, p.weight);
    }

    public double euclidDist(MicroPoint p1) {
        double sumDistSq = 0.0;
        for (int i = 0; i < pos.length; i++) {
            sumDistSq += Math.pow((pos[i] - p1.pos[i]), 2);
        }
        return Math.sqrt(sumDistSq);
    }
}
