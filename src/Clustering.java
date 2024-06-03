import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class Clustering {

    public static void main(String[] args) {
        try{
            Scanner scanner = new Scanner(System.in);
            boolean tryAgain = true;

            while(tryAgain) {
                System.out.print("Enter the filename to use (Content root path): ");
                String filename = scanner.nextLine();

                System.out.println("Select one of the following algorithms to run (1 to 4):");
                System.out.println("1. Algorithm 1: Gonzalez + LLoyd's Algorithm");
                System.out.println("2. Algorithm 2: Range-based Initialization + Lloyd's Algorithm");
                System.out.println("3. Algorithm 3: LOF + a heuristic to solve k-median");
                System.out.println("4. Algorithm 4: Binary Search + Greedy Algorithm");
                System.out.println("Note that the output files will be saved as alg.sol. The printing of the results will also be displayed on the console.");

                int choice = scanner.nextInt();
                scanner.nextLine();

                ArrayList<Points> inputs = readInputs(filename);
                int k = inputs.get(0).getX();
                int m = inputs.get(0).getY();
                inputs.remove(0);
                System.out.println("Instance file: " + filename + " has been read.");

                switch (choice) {
                    case 1:
                        // Method 1
                        solveKMedianProblem1_instance(inputs, k, m, filename);
                        break;
                    case 2:
                        // Method 2
                        solveKMedianProblem2_instance(inputs, k, m, filename);
                        break;
                    case 3:
                        //Method 3
                        System.out.println("You have chosen LOF + a heuristic to solve k-median with outliers. Choose the heuristic (1 or 2):");
                        System.out.println("1. Algorithm 1: Gonzalez + LLoyd's Algorithm");
                        System.out.println("2. Algorithm 2: Range-based Initialization + Lloyd's Algorithm");
                        ArrayList<Points> cleanedData = removeOutliersLOF(inputs, k, m);
                        int heuristicChoice = scanner.nextInt();
                        scanner.nextLine();
                        if(heuristicChoice == 1)
                            solveKMedianWithOutliers1_instance(cleanedData, k, m, filename);
                        else if(heuristicChoice == 2)
                            solveKMedianProblem2_instance(cleanedData, k, m, filename);
                        else
                            System.out.println("Invalid choice.");
                        break;
                    case 4:
                        // Method 4
                        solveKMedianWithOutliers2_instance(inputs, k, m, filename);
                        break;
                    default:
                        System.out.println("Invalid choice.");
                        break;
                }
                System.out.println();
                System.out.print("Do you want to try a different file or different algorithm? (yes/no): ");
                String tryAgainInput = scanner.nextLine().toLowerCase();
                tryAgain = tryAgainInput.equals("yes");

            }
            scanner.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    //Algorithm 1 to solve k-median problem (initialize k centers farthest from each other)
    public static void solveKMedianProblem1_instance(ArrayList<Points> inputs, int k, int m, String filename) {
        //Initialization of k centers
        //ArrayList to store centers of clusters
        ArrayList<Points> centers = new ArrayList<>();
        boolean[] isCenter = new boolean[inputs.size()];   //isCenter default = false

        //Step 1 Initialize: First center by randomizing
        int randomInt = new Random().nextInt(inputs.size());
        centers.add(inputs.get(randomInt));
        isCenter[randomInt] = true;
        //Step 2 Initialize: Remaining centers by calculating distance from each point to nearest center
        //Calculate minimum distance from each point to already chosen center
        while (centers.size() < k) {
            double[] minDistance = new double[inputs.size()];
            for (int i = 0; i < inputs.size(); i++) {
                if (!isCenter[i]) {
                    minDistance[i] = Double.MAX_VALUE;
                    for (Points center : centers) {
                        double distance = inputs.get(i).getDistance(center);
                        if (distance < minDistance[i]) {
                            minDistance[i] = distance;
                        }
                    }
                }
            }

            //Find the point furthest from the nearest center
            double maxDistance = minDistance[0];
            int maxIndex = 0;  //the index of the point furthest from the nearest center
            for (int i = 1; i < inputs.size(); i++) {
                if (minDistance[i] > maxDistance) {
                    maxDistance = minDistance[i];
                    maxIndex = i;
                }
            }
            //Add the point to the list of centers
            centers.add(inputs.get(maxIndex));
            isCenter[maxIndex] = true;
        }

        //Assign remaining points to the nearest center to form a cluster
        boolean converged = false;
        int iterations = 0;

        //Initialize clusters
        ArrayList<ArrayList<Points>> clusters = new ArrayList<>();
        for (int i = 0; i < centers.size(); i++) {
            clusters.add(new ArrayList<>());
        }

        while (!converged) {
            iterations++;
            //Clear the clusters for new iteration
            for (ArrayList<Points> cluster : clusters) {
                cluster.clear();
            }


            //Assign each point to the nearest center
            for (int i = 0; i < inputs.size(); i++) {
                if (!isCenter[i]) {
                    double minDist = Double.MAX_VALUE;
                    int clusterIndex = 0;
                    for (int j = 0; j < centers.size(); j++) {
                        double dist = inputs.get(i).getDistance(centers.get(j));
                        if (dist < minDist) {
                            minDist = dist;
                            clusterIndex = j;
                        }
                    }
                    clusters.get(clusterIndex).add(inputs.get(i));
                }
            }
            converged = true;

            //Update new centers
            for (int i = 0; i < clusters.size(); i++) {
                Points newCenter = calculateNewCenter(clusters.get(i), inputs);
                if (!newCenter.equals(centers.get(i))) {
                    centers.set(i, newCenter);
                    converged = false;
                }
            }
        }

        printCenters(centers);
        System.out.println("Number of iterations: " + iterations);
        double totalDistance = calculateTotalClusterDistances(clusters, centers);
        System.out.println("Total distance of points to their centers: " + totalDistance);

        //Writing the solution to a file
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(filename + ".alg1.sol"));
            writer.write("k " + k);
            writer.newLine();
            writer.write("m " + m);
            writer.newLine();
            for (Points center : centers) {
                writer.write(center.getX() + " " + center.getY());
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //Algorithm 2 to solve k-median problem (initialize k centers by selecting random points in each interval)
    public static void solveKMedianProblem2_instance(ArrayList<Points> inputs, int k, int m, String filename) {
        // Find the range of the data
        int xmax = Integer.MIN_VALUE;
        int xmin = Integer.MAX_VALUE;
        int ymax = Integer.MIN_VALUE;
        int ymin = Integer.MAX_VALUE;

        //Find the extreme values of x and y
        for (Points point : inputs) {
            if (point.getX() > xmax) xmax = point.getX();
            if (point.getX() < xmin) xmin = point.getX();
            if (point.getY() > ymax) ymax = point.getY();
            if (point.getY() < ymin) ymin = point.getY();
        }

        // Determine the larger range to set intervals (range of x or range of y)
        int range = Math.max(xmax - xmin, ymax - ymin);
        double intervalSize = range / (double) k;    //calculate the interval size given k

        ArrayList<Points> centers = new ArrayList<>();
        boolean[] isCenter = new boolean[inputs.size()];
        Random random = new Random();

        // Initialize centers
        boolean useX = (xmax - xmin) > (ymax - ymin);   //boolean to use x if x has larger range
        // Initialize k centers by selecting a random point in each interval
        for (int i = 0; i < k; i++) {
            //Calculate the start and end of the interval
            double start = (useX ? xmin : ymin) + i * intervalSize;  //ternary: if useX is true, start = xmin, else start = ymin
            double end = start + intervalSize;
            ArrayList<Integer> eligibleIndices = new ArrayList<>();

            for (int j = 0; j < inputs.size(); j++) {
                //Check if the point is within the interval and has not been selected as a center
                double value = useX ? inputs.get(j).getX() : inputs.get(j).getY();   //ternary: if useX is true, value = x, else value = y
                if (value >= start && value <= end && !isCenter[j]) {
                    eligibleIndices.add(j);   //add the index of the point to the list of eligible indices
                }
            }

            if (!eligibleIndices.isEmpty()) {
                //Select a random point from the list of eligible indices as center
                int selectedIdx = eligibleIndices.get(random.nextInt(eligibleIndices.size()));   //random point within interval
                centers.add(inputs.get(selectedIdx));
                isCenter[selectedIdx] = true;  //update boolean isCenter
            }
        }

        // Use the initialized centers to proceed with clustering
        boolean converged = false;
        int iterations = 0;

        //Initialize clusters
        ArrayList<ArrayList<Points>> clusters = new ArrayList<>();
        for (int i = 0; i < centers.size(); i++) {
            clusters.add(new ArrayList<>());
        }

        while (!converged) {
            iterations++;
            //Clear the clusters for new iteration
            for (ArrayList<Points> cluster : clusters) {
                cluster.clear();
            }

            //Assign each point to the nearest center
            for (int i = 0; i < inputs.size(); i++) {
                if (!isCenter[i]) {
                    double minDist = Double.MAX_VALUE;
                    int clusterIndex = 0;
                    for (int j = 0; j < centers.size(); j++) {
                        double dist = inputs.get(i).getDistance(centers.get(j));
                        if (dist < minDist) {
                            minDist = dist;
                            clusterIndex = j;
                        }
                    }
                    clusters.get(clusterIndex).add(inputs.get(i));
                }
            }
            converged = true;

            //Update new centers
            for (int i = 0; i < clusters.size(); i++) {
                Points newCenter = calculateNewCenter(clusters.get(i), inputs);
                if (!newCenter.equals(centers.get(i))) {
                    centers.set(i, newCenter);
                    converged = false;
                }
            }
        }

        printCenters(centers);
        System.out.println("Number of iterations: " + iterations);
        double totalDistance = calculateTotalClusterDistances(clusters, centers);
        System.out.println("Total distance of points to their centers: " + totalDistance);

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(filename + ".alg2.sol"));
            writer.write("k " + k);
            writer.newLine();
            writer.write("m " + m);
            writer.newLine();
            for (Points center : centers) {
                writer.write(center.getX() + " " + center.getY());
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //Algorith 3 to solve k-median problem with outliers (remove outliers using LOF)
    public static void solveKMedianWithOutliers1_instance(ArrayList<Points> inputs, int k, int m, String filename) {
        //Repeat same clustering method as in solveKMedianProblem1_instance to find clusters to cleaned data
        ArrayList<Points> centers = new ArrayList<>();
        boolean[] isCenter = new boolean[inputs.size()];
        int randomInt = new Random().nextInt(inputs.size());
        centers.add(inputs.get(randomInt));
        isCenter[randomInt] = true;

        while (centers.size() < k) {
            double[] minDistance = new double[inputs.size()];
            for (int i = 0; i < inputs.size(); i++) {
                if (!isCenter[i]) {
                    minDistance[i] = Double.MAX_VALUE;
                    for (Points center : centers) {
                        double distance = inputs.get(i).getDistance(center);
                        if (distance < minDistance[i]) {
                            minDistance[i] = distance;
                        }
                    }
                }
            }
            double maxDistance = minDistance[0];
            int maxIndex = 0;
            for (int i = 1; i < inputs.size(); i++) {
                if (minDistance[i] > maxDistance) {
                    maxDistance = minDistance[i];
                    maxIndex = i;
                }
            }
            centers.add(inputs.get(maxIndex));
            isCenter[maxIndex] = true;
        }
        boolean converged = false;
        int iterations = 0;
        ArrayList<ArrayList<Points>> clusters = new ArrayList<>();
        for (int i = 0; i < centers.size(); i++) {
            clusters.add(new ArrayList<>());
        }

        while (!converged) {
            iterations++;
            for (ArrayList<Points> cluster : clusters) {
                cluster.clear();
            }
            for (int i = 0; i < inputs.size(); i++) {
                if (!isCenter[i]) {
                    double minDist = Double.MAX_VALUE;
                    int clusterIndex = 0;
                    for (int j = 0; j < centers.size(); j++) {
                        double dist = inputs.get(i).getDistance(centers.get(j));
                        if (dist < minDist) {
                            minDist = dist;
                            clusterIndex = j;
                        }
                    }
                    clusters.get(clusterIndex).add(inputs.get(i));
                }
            }
            converged = true;

            //Update new centers
            for (int i = 0; i < clusters.size(); i++) {
                Points newCenter = calculateNewCenter(clusters.get(i), inputs);
                if (!newCenter.equals(centers.get(i))) {
                    centers.set(i, newCenter);
                    converged = false;
                }
            }
        }

        printCenters(centers);
        System.out.println("Number of iterations: " + iterations);
        double totalDistance = calculateTotalClusterDistances(clusters, centers);
        System.out.println("Total distance of points to their centers: " + totalDistance);
        // Calculate the largest minimum distance of all points to the centers
        double largestMinDistance = calculateLargestMinDistance(inputs, centers);
        System.out.println("The radius is:  " + largestMinDistance);

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(filename + ".alg3.sol"));
            writer.write("k " + k);
            writer.newLine();
            writer.write("m " + m);
            writer.newLine();
            writer.write("r " + largestMinDistance);
            writer.newLine();
            for (Points center : centers) {
                writer.write(center.getX() + " " + center.getY());
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    //Algorithm 4 to solve k-median problem with outliers (binary search to find optimal r and greedy algorithm to for k clusters)
    public static void solveKMedianWithOutliers2_instance(ArrayList<Points> inputs, int k, int m, String filename){
        double low = 0;
        double high = 0;
        ArrayList<Points> bestCenters = new ArrayList<>();

        // Find the largest distance between any two points
        for (int i = 0; i < inputs.size(); i++) {
            for (int j = i + 1; j < inputs.size(); j++) {
                double distance = inputs.get(i).getDistance(inputs.get(j));
                if (distance > high) {
                    high = distance;
                }
            }
        }

        //Binary search to find the optimal radius
        //Initialize r to the largest distance between any two points
        double bestRadius = high;
        while (high - low > 1e-5) {  // Set precision limit
            double mid = (low + high) / 2;
            ArrayList<Points> currentCenters = new ArrayList<>();
            //Check if new radius can cover m points with k centers
            if (canCoverMPoints(inputs, currentCenters, k, m, mid)) {
                bestRadius = mid;
                //If can cover m points with k centers, reduce the radius by updating upperbound to be mid
                high = mid;
                bestCenters = new ArrayList<>(currentCenters); // Update best centers
            } else {
                //If cannot cover m points with k centers, increase the radius by updating lower bound to be mid
                low = mid;
            }
        }
        printResults(bestRadius, bestCenters);

        // Calculate the covered count
        int coveredCount = 0;
        for (Points point : inputs) {
            double minDist = Double.MAX_VALUE;
            for (Points center : bestCenters) {
                double dist = point.getDistance(center);
                if (dist <= bestRadius) {
                    minDist = dist;
                    coveredCount++;
                    break; // Break if the point is covered by a center
                }
            }
        }

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(filename + ".alg4.sol"));
            writer.write("k " + k);
            writer.newLine();
            writer.write("m " + coveredCount);
            writer.newLine();
            writer.write("r " + bestRadius);
            writer.newLine();
            for (Points center : bestCenters){
                writer.write(center.getX() + " " + center.getY());
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // Check if m points can be covered by k centers with radius r using greedy algorithm
    public static boolean canCoverMPoints(ArrayList<Points> inputs, ArrayList<Points> centers, int k, int m, double r) {
        boolean[] covered = new boolean[inputs.size()];
        int coveredCount = 0; // Number of points covered by the centers
        centers.clear();

        for (int i = 0; i < k && coveredCount < m; i++) {
            Points bestCenter = null;
            int bestCoverage = -1;
            // Find the best new center that covers the most points
            for (Points candidate : inputs) {
                int currentCoverage = 0;
                for (int j = 0; j < inputs.size(); j++) {
                    if (!covered[j] && candidate.getDistance(inputs.get(j)) <= r) {
                        currentCoverage++;
                    }
                }
                // Update the best center if the current candidate covers more points
                if (currentCoverage > bestCoverage) {
                    bestCenter = candidate;
                    bestCoverage = currentCoverage;
                }
            }
            if (bestCenter == null) {
                return false;
            }
            // Mark points as covered
            for (int j = 0; j < inputs.size(); j++) {
                if (bestCenter.getDistance(inputs.get(j)) <= r) {  // Check if the point is within the radius
                    if (!covered[j]) {
                        covered[j] = true;
                        coveredCount++;  // Increment the number of points covered and add the point to the cluster
                    }
                }
            }
            centers.add(bestCenter);
        }
        return coveredCount >= m;
    }

    //Using Local Outlier Factor to remove outliers
    public static ArrayList<Points> removeOutliersLOF(ArrayList<Points> data, int k, int m) {
        // Calculate LOF for each data point
        double[] lof = calculateLOF(data, k);

        // Sort the LOF values
        ArrayList<Double> lofList = new ArrayList<>();
        for (double value : lof) {
            lofList.add(value);
        }
        Collections.sort(lofList, Collections.reverseOrder());

        // Remove (numPoints - m) outliers
        int numOutliersToRemove = data.size() - m;
        for (int i = 0; i < numOutliersToRemove; i++) {
            double outlierLof = lofList.get(i);
            for (int j = 0; j < data.size(); j++) {
                if (lof[j] == outlierLof) {
                    data.remove(j);
                    lofList.remove(i);
                    break;
                }
            }
        }

        // return cleaned data
        return data;
    }

    //Read inputs from a file
    public static ArrayList<Points> readInputs(String filename)throws FileNotFoundException {
        Scanner scanner = new Scanner(new File(filename));

        ArrayList<Points> inputs = new ArrayList<>();
        if(scanner.hasNext()) {
            scanner.next(); //Consume string "k"
            int k = scanner.nextInt();
            scanner.next(); //Consume string "m"
            int m = scanner.nextInt();
            inputs.add(new Points (k, m)); //The first point in inputs is just to store k and m
        }

        while(scanner.hasNextInt()) {
            int x = scanner.nextInt();
            int y = scanner.nextInt();
            inputs.add(new Points(x, y));
        }
        return inputs;
    }

    public static double calculateLargestMinDistance(ArrayList<Points> inputs, ArrayList<Points> centers) {
        double largestMinDistance = Double.MIN_VALUE;
        for (Points point : inputs) {
            double minDistance = Double.MAX_VALUE;
            for (Points center : centers) {
                double distance = point.getDistance(center);
                if (distance < minDistance) {
                    minDistance = distance;
                }
            }
            if (minDistance > largestMinDistance) {
                largestMinDistance = minDistance;
            }
        }
        return largestMinDistance;
    }

    public static double calculateLocalReachabilityDensity(ArrayList<Points> data, Points point, int k) {
        int dataSize = data.size();
        ArrayList<Double> distances = new ArrayList<>();

        // Calculate distances to all other points
        for (int j = 0; j < dataSize; j++) {
            Points otherPoint = data.get(j);
            double distance = point.getDistance(otherPoint);
            distances.add(distance);
        }

        // Sort distances and find k-nearest neighbors
        Collections.sort(distances);
        double kDistance = distances.get(k - 1);

        // Calculate reachability distance
        double sumReachabilityDistance = 0.0;
        for (int j = 0; j < k; j++) {
            double reachabilityDistance = Math.max(distances.get(j), kDistance);
            sumReachabilityDistance += reachabilityDistance;
        }

        // Calculate local reachability density
        double avgReachabilityDistance = sumReachabilityDistance / k;
        double localReachabilityDensity = 1.0 / avgReachabilityDistance;

        return localReachabilityDensity;
    }

    // Calculate the Local Outlier Factor (LOF) for each data point
    public static double[] calculateLOF(ArrayList<Points> data, int k) {
        int dataSize = data.size();
        double[] lof = new double[dataSize];

        for (int i = 0; i < dataSize; i++) {
            Points point = data.get(i);
            ArrayList<Double> distances = new ArrayList<>();

            // Calculate distances to all other points
            for (int j = 0; j < dataSize; j++) {
                if (i != j) {
                    Points otherPoint = data.get(j);
                    double distance = calculateDistance(point, otherPoint);
                    distances.add(distance);
                }
            }

            // Sort distances and find k-nearest neighbors
            Collections.sort(distances);
            double kDistance = distances.get(k - 1);

            // Calculate reachability distance and local reachability density
            double sumReachabilityDistance = 0.0;
            for (int j = 0; j < k; j++) {
                double reachabilityDistance = Math.max(distances.get(j), kDistance);
                sumReachabilityDistance += reachabilityDistance;
            }
            double avgReachabilityDistance = sumReachabilityDistance / k;
            double localReachabilityDensity = 1.0 / avgReachabilityDistance;

            // Calculate LOF
            double sumLof = 0.0;
            for (int j = 0; j < k; j++) {
                double neighborLrd = calculateLocalReachabilityDensity(data, data.get(j), k);
                sumLof += neighborLrd / localReachabilityDensity;
            }
            double avgLof = sumLof / k;
            lof[i] = avgLof;
        }

        return lof;
    }

    public static double calculateDistance(Points p1, Points p2) {
        // Implement distance calculation between two points (e.g., Euclidean distance)
        return Math.sqrt(Math.pow(p1.getX() - p2.getX(), 2) + Math.pow(p1.getY() - p2.getY(), 2));
    }

    // Print the optimal radius and centers for algorithm 2 k-median with outliers
    private static void printResults(double radius, ArrayList<Points> centers) {
        System.out.println("Algorithm 4: Optimal Radius: " + radius);
        System.out.println("Centers:");
        for (Points center : centers) {
            System.out.println("(" + center.getX() + ", " + center.getY() + ")");
        }
    }

    // Calculate the total distance of all points to their respective centers
    public static double calculateTotalClusterDistances(ArrayList<ArrayList<Points>> clusters, ArrayList<Points> centers) {
        double totalDistance = 0;

        // Iterate over each cluster and its corresponding center
        for (int i = 0; i < clusters.size(); i++) {
            ArrayList<Points> cluster = clusters.get(i);
            Points center = centers.get(i);

            // Calculate the distance of each point in the cluster to the center
            for (Points point : cluster) {
                totalDistance += point.getDistance(center);
            }
        }

        return totalDistance;
    }

    // Calculate the new center of a cluster
    public static Points calculateNewCenter(ArrayList<Points> cluster, ArrayList<Points> inputs) {
        if (cluster.isEmpty()) return null; // Return null if the cluster is empty

        double sumX = 0, sumY = 0;
        // Calculate the center by averaging the coordinates of all points in the cluster
        for (Points p : cluster) {
            sumX += p.getX();
            sumY += p.getY();
        }
        int x = (int) (sumX / cluster.size());
        int y = (int) (sumY / cluster.size());
        Points center = new Points(x, y);

        // Find the actual point in inputs that is closest to the calculated center
        Points closestPoint = null;
        double minDistance = Double.MAX_VALUE;
        for (Points p : inputs) {
            double distance = center.getDistance(p);
            if (distance < minDistance) {
                minDistance = distance;
                closestPoint = p;
            }
        }

        return closestPoint; // Return the actual closest point from the original inputs
    }

    //Print the centers
    public static void printCenters(ArrayList<Points> centers) {
        int count = 1;
        for (Points center : centers) {
            System.out.println("Center " +count+  " :(" + center.getX() + ", " + center.getY() + ")");
            count++;
        }
    }

    //Print the clusters
    public static void printClusters(ArrayList<ArrayList<Points>> clusters) {
        for (int i = 0; i < clusters.size(); i++) {
            System.out.println("Cluster " + (i + 1) + ":");
            for (Points point : clusters.get(i)) {
                System.out.println("\tPoint: (" + point.getX() + ", " + point.getY() + ")");
            }
        }
    }
}
