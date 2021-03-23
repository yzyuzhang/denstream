
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Batch DBSCAN clusterer.
 */
public class DBSCAN {
    double eps;
    double mu;
    double lambda;
    int clusterGlobalID;  // cluster unique ID, start from 0

    public DBSCAN(final double eps, final double mu, final double lambda) {
        this.eps = eps;
        this.mu = mu;
        this.lambda = lambda;
        clusterGlobalID = 0;
    }

    /**
     * Batch DBSCAN clustering algorithm.
     * The clustering result is each point is labelled
     * with either a cluster index or noise.
     *
     * @param points points to cluster
     */
    public List<Point> cluster(final List<Point> points, long currentTimestamp) {
        resetCluster();

        // copy points to a new list to cluster (deep copy)
        List<Point> pointsToCluster = new ArrayList<>();
        for (Point p : points) {
            pointsToCluster.add(new Point(p));
        }

        for (Point point : pointsToCluster) {
            if (point.visited) {
                continue;
            }
            point.visited = true;
            final List<Point> neighbors = getNeighbors(point, pointsToCluster);

            double totalWeight = 0;
            for (Point nei : neighbors) {
                totalWeight += nei.getWeight(currentTimestamp, lambda);
            }

            if (totalWeight >= mu) {
                point.trueLabel = clusterGlobalID;
                expandCluster(point, neighbors, pointsToCluster, clusterGlobalID, currentTimestamp);
                clusterGlobalID++;
            } else {
                // noise point temporarily, may become border point later
                point.trueLabel = -1;
            }
        }

        return pointsToCluster;
    }

    private void resetCluster() {
        clusterGlobalID = 0;
    }


    // Expands the cluster to include all density-reachable points.
    private void expandCluster(final Point point, final List<Point> neighbors,
                               final List<Point> points, int clusterId, long currentTimestamp) {
        List<Point> seeds = new ArrayList<>(neighbors);
        int index = 0;
        while (index < seeds.size()) {
            final Point current = seeds.get(index);
            // only check non-visited points
            if (!current.visited) {
                current.visited = true;
                current.trueLabel = clusterId;
                final List<Point> currentNeighbors = getNeighbors(current, points);

                double totalWeight = 0;
                for (Point nei : neighbors) {
                    totalWeight += nei.getWeight(currentTimestamp, lambda);
                }

                // current point is a density-connected core point
                if (totalWeight >= mu) {
                    for (Point currentNbr : currentNeighbors) {
                        seeds.add(currentNbr);
                    }
                }
            }

            // assign cluster ID to boarder point
            if (current.trueLabel == -1) {
                current.visited = true;
                current.trueLabel = clusterId;
            }

            index++;
        }
    }

    // Return a list of density-reachable neighbors of a point.
    private List<Point> getNeighbors(final Point point, final List<Point> points) {
        final List<Point> neighbors = new ArrayList<>();
        for (final Point p : points) {
            // include point itself
            if (point.euclidDist(p) <= eps) {
                neighbors.add(p);
            }
        }
        return neighbors;
    }
}
