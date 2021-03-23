public class Point {
    int trueLabel;
    int testLabel;
    double[] pos;
    double initWeight;
    boolean visited;
    long initTimestamp;

    public Point(double[] pos, double weight, int trueLabel, long initTimestamp) {
        this.pos = pos;
        this.initWeight = weight;
        this.trueLabel = trueLabel;
        this.visited = false;
        this.initTimestamp = initTimestamp;
    }

    public Point(Point p1) {
        this.pos = p1.pos.clone();
        this.initWeight = p1.initWeight;
        this.trueLabel = p1.trueLabel;
        this.visited = p1.visited;
        this.initTimestamp = p1.initTimestamp;
    }

    public double euclidDist(Point p1) {
        double sumDistSq = 0.0;
        int d = pos.length;
        for (int i = 0; i < d; i++) {
            sumDistSq += Math.pow(pos[i] - p1.pos[i], 2);
        }
        return Math.sqrt(sumDistSq);
    }

    public double getWeight(long timestamp, double lambda) {
        long dt = timestamp - this.initTimestamp;
        return (initWeight * Math.pow(2, -lambda * dt));
    }
}
