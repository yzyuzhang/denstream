import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Batch DBSCAN clusterer for DenStream.
 */
public class GenerateClusters {
    private final double eps;
    private final double mu;
    private int clusterGlobalId;  // increasing cluster unique ID, start from 0.
    List<MicroPoint> pointsToCluster;

    public GenerateClusters(final double eps, final double mu) {
        this.eps = eps;
        this.mu = mu;
        clusterGlobalId = 0;
        pointsToCluster = new ArrayList<>();
    }

    /**
     * Batch DBSCAN clustering algorithm.
     * The clustering result is each point is labelled
     * with either a cluster index or noise.
     *
     * @param points points to cluster
     */
    public List<MicroPoint> cluster(final List<MicroPoint> points) {
        resetCluster();

        // copy points to a new list to cluster (deep copy)
        for (MicroPoint p : points) {
            pointsToCluster.add(new MicroPoint(p));
        }

        for (MicroPoint point : pointsToCluster) {
            if (point.visited) {
                continue;
            }
            point.visited = true;
            List<MicroPoint> neighbors = getNeighbors(point, pointsToCluster);

            double totalWeight = 0;
            for (MicroPoint p : neighbors) {
                totalWeight += p.weight;
            }
            if (totalWeight >= mu) {
                point.clusterId = clusterGlobalId;
                expandCluster(point, neighbors, pointsToCluster, clusterGlobalId);
                clusterGlobalId++;
            } else {
                // noise point temporarily, may become border point later
                point.clusterId = -1;
            }
        }

        return pointsToCluster;
    }

    /**
     * Clear cluster state.
     */
    private void resetCluster() {
        clusterGlobalId = 0;
        pointsToCluster.clear();
    }

    /**
     * Expands the cluster to include all density-reachable points.
     *
     * @param point     starting core point
     * @param neighbors point's neighbors
     * @param points    all point set
     * @param clusterId new cluster id
     */
    private void expandCluster(final MicroPoint point, final List<MicroPoint> neighbors,
                               final List<MicroPoint> points, int clusterId) {
        List<MicroPoint> seeds = new ArrayList<>(neighbors);
        int index = 0;
        while (index < seeds.size()) {
            final MicroPoint current = seeds.get(index);
            // Only check non-visited points.
            if (!current.visited) {
                current.visited = true;
                current.clusterId = clusterId;
                final List<MicroPoint> currentNeighbors = getNeighbors(current, points);

                // Current point is a density-connected core point.
                double totalWeight = 0;
                for (MicroPoint p : currentNeighbors) {
                    totalWeight += p.weight;
                }
                if (totalWeight >= mu) {
                    for (MicroPoint currentNbr : currentNeighbors) {
                        seeds.add(currentNbr);
                    }
                }
            }

            // Assign cluster ID to boarder point.
            if (current.clusterId == -1) {
                current.visited = true;
                current.clusterId = clusterId;
            }

            index++;
        }
    }

    /**
     * Return a list of density-reachable neighbors of a {@code point}
     *
     * @param p  the point to look for
     * @param points all points
     * @return neighbors (including point itself)
     */
    private List<MicroPoint> getNeighbors(final MicroPoint p,
                                          final List<MicroPoint> points) {
        final List<MicroPoint> neighbors = new ArrayList<>();
        for (final MicroPoint p1 : points) {
            double dist = p.euclidDist(p1);
            // Include point itself.
            if (dist <= 2 * eps && dist <= p.radius + p1.radius) {
                neighbors.add(p);
            }
        }
        return neighbors;
    }
}
