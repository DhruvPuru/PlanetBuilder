package pb.g6;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Random;

/**
 * Created by llxtt on 10/10/15.
 */


public class Utils {
    public Utils() {
    }

    // Euclidean distance between two points
    public static double getPerpendicularAngle(double oriAngle)
    {
        // used to pick push angle perpendicular to v.direction() inward or outward
        Random random = new Random();
        double newAngle = oriAngle + (random.nextDouble() - 0.5) * Math.PI * 0.5;
        return newAngle;
    }
}
