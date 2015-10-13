package pb.g2;

import pb.sim.Asteroid;
import pb.sim.Point;

import java.util.ArrayList;
import java.util.Hashtable;

public class CollisionChecker {
    public static class CollisionPair {
        public int asteroid1_index, asteroid2_index;
    }

    public static long checkCollision(Asteroid a1, Asteroid a2, long max_time, long real_time, long time_limit) {
        // avoid allocating a new Point object for every position
        Point p1 = new Point(), p2 = new Point();
        // search for collision with other asteroids
        double r = a1.radius() + a2.radius();
        for (long ft = 0 ; ft != max_time ; ++ft) {
            long t = real_time + ft;
            if (t >= time_limit) break;
            a1.orbit.positionAt(t - a1.epoch, p1);
            a2.orbit.positionAt(t - a2.epoch, p2);
            // if collision, return push to the simulator
            if (Point.distance(p1, p2) < r) {
                return t;
            }
        }
        return -1; // No collision within deadline
    }

    public static Hashtable<Long, ArrayList<CollisionPair>> checkCollision(
            Asteroid[] asteroids, long max_time, long real_time, long time_limit) {
        // Avoid allocating a new Point object for every position
        int n = asteroids.length;
        Point[] asteroid_future_positions = new Point[n];

        // Search for collision with other asteroids
        double[][] radius_matrix = new double[n][n];

        // Not we're only computing total radius for all i<j to save computation time
        for (int i = 0; i < n; i++) {
            for (int j = i+1; j < n; j++) {
                radius_matrix[i][j] = asteroids[i].radius() + asteroids[j].radius();
            }
        }

        Hashtable<Long, ArrayList<CollisionPair>> collisions = new Hashtable<Long, ArrayList<CollisionPair>>();

        for (long dt = 0; dt <= max_time; dt++) {
            long t = real_time + dt;
            if (t >= time_limit) {
                break;
            }

            ArrayList<CollisionPair> collisions_at_t = new ArrayList<CollisionPair>();

            for (int i = 0; i < n; i++) {
                long epoch = asteroids[i].epoch;
                asteroids[i].orbit.positionAt(t - epoch, asteroid_future_positions[i]);
            }

            // Find lowest energy collision and return push to the simulator
            // Not we're only checking distance for all i<j to save computation time
            for (int i = 0; i < n; i++) {
                for (int j = i+1; j < n; j++) {
                    if (Point.distance(asteroid_future_positions[i], asteroid_future_positions[j])
                            < radius_matrix[i][j]) {
                        CollisionPair collision = new CollisionPair();
                        collision.asteroid1_index = i;
                        collision.asteroid2_index = j;
                        collisions_at_t.add(collision);
                    }
                }
            }

            // If there was a collision at this time, add it to the table
            if (!collisions_at_t.isEmpty()) {
                // TODO: Sort collisions by some metric
                collisions.put(t, collisions_at_t);
            }
        }

        return collisions;
    }
}
