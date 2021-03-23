import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

public class Main {

    double eps = 1.5;
    double mu = 3;
    double lambda = 2;

    public static void main(String[] args) {
        String inFn = "./input/sample_data.txt";
        String denStreamOut = "./output/denstream-accuracy.txt";
        int numPointsToRead = 2000;
        List<Point> stream = readData(inFn, numPointsToRead);

        // Number of points received.
        int count = 0;
        List<Point> points = new ArrayList<>();
        int interval = 10000;
        for (Point p : stream) {
            count++;
            points.add(p);
            long currentTime = p.initTimestamp;
            if (count % interval == 0) {
                // Run DenStream.
                List<Point> denStreamResult = new ArrayList<>();

                // Run batch DBSCAN.
                List<Point> batchResult = new Main().runBatch(points, currentTime);

                // Compute the accuracy.
                double accuracy = computeAccuracy(batchResult, denStreamResult);
            }
        }
    }

    private static double computePurity(List<Integer> list) {
        HashMap<Integer, Integer> map = new HashMap<>();
        for (int e : list) {
            map.put(e, map.getOrDefault(e, 0) + 1);
        }

        int maxCount = 0;
        for (int key : map.keySet()) {
            if (map.get(key) > maxCount) {
                maxCount = map.get(key);
            }
        }
        double pur = ((double) maxCount) / list.size();
        return pur;
    }

    // Compute purity of the clusters.
    private static double computeAccuracy(List<Point> batch, List<Point> denStream) {
        for (int i=0; i<batch.size(); i++) {
            batch.get(i).testLabel = denStream.get(i).testLabel;
        }
        HashMap<Integer, List<Integer>> map = new HashMap<>();
        for (Point p : batch) {
            if (!map.containsKey(p.testLabel)) {
                map.put(p.testLabel, new ArrayList<Integer>());
            }
            map.get(p.testLabel).add(p.trueLabel);
        }

        double pur = 0;
        for (int key : map.keySet()) {
            List<Integer> values = map.get(key);
            pur += computePurity(values);
        }
        return pur/map.size();
    }

    private List<Point> runBatch(List<Point> input, long currentTime) {
        // Run batch DBSCAN to get the true label.
        DBSCAN dbscan = new DBSCAN(eps, mu, lambda);
        List<Point> res = dbscan.cluster(input, currentTime);
        return res;
    }

    private static List<Point> readData(String inFileName, int numPointsToRead) {
        List<Point> points = new ArrayList<>();
        try {
            Scanner sc = new Scanner(new File(inFileName));
            int numPoints = 0;
            while (sc.hasNextLine() && numPoints < numPointsToRead) {
                String line = sc.nextLine();
                String[] strs = line.split(",");
                int d = strs.length;
                double[] pos = new double[d];
                for (int i = 0; i < d; i++) {
                    pos[i] = Double.parseDouble(strs[i]);
                }
                double initWeight = 1.0;
                long timestamp = 0;
                Point p = new Point(pos, initWeight, -1, timestamp);
                points.add(p);
                numPoints++;
            }
            sc.close();
        } catch (Exception e) {
            System.out.println("read input exception");
        }
        System.out.println("Reading data complete, number of points: " + points.size() + "\n");
        return points;
    }
}
