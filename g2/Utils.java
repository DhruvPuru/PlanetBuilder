package pb.g2;

import pb.sim.Asteroid;

public class Utils {
    /**
     * Returns index of asteroid of largest radius.
     */
    public static int largestOrbitRadius(Asteroid[] asteroids) {
        int largest = 0;
        for (int i = 1; i < asteroids.length; i++) {
            if (asteroids[i].orbit.a > asteroids[largest].orbit.a) {
                largest = i;
            }
        }
        return largest;
    }

    public static int largestMass(Asteroid[] asteroids) {
        int largest = 0;
        for (int i = 1; i < asteroids.length; i++) {
            if (asteroids[i].mass > asteroids[largest].mass) {
                largest = i;
            }
        }
        return largest;
    }
}
