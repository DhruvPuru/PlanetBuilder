package pb.g4;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Random;

public class Player implements pb.sim.Player {

    // used to PIck asteroid and velocity boost randomly
    private Random random = new Random();

    // current time, time limit
    private long time = -1;
    private long time_limit = -1;

    // time until next push
    private long time_of_push = 0;

    // number of retries
    private int retries_per_turn = 1;
    private int turns_per_retry = 100;

    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit)
    {
        if (Orbit.dt() != 24 * 60 * 60)
            throw new IllegalStateException("Time quantum is not a day");
        this.time_limit = time_limit;
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids,
            double[] energy, double[] direction)
    {
        if (++time == time_of_push) {
            double min = Double.MAX_VALUE;
            int minIndex1 = -1;
            int minIndex2 = -1;
            // if not yet time to push do nothing
            for (int i = 0; i < asteroids.length; i++) {
                Asteroid a = asteroids[i];
                Point aPoint = a.orbit.positionAt(++time);
                for (int j = i+1; j < asteroids.length; j++) {
                    Asteroid b = asteroids[j];
                    Point bPoint = b.orbit.positionAt(time);
                    double dist = Point.distance(aPoint, bPoint);
                    if (dist < min) {
                        minIndex1 = i;
                        minIndex2 = j;
                        min = dist;
                    }
                }
            }

            Asteroid furthestFromSun = asteroids[minIndex1];
            Asteroid closerToSun = asteroids[minIndex2];
            int indexToPush = minIndex1;

            Point origin = new Point(0, 0);
            Point furthest = furthestFromSun.orbit.positionAt(time);
            Point closest = closerToSun.orbit.positionAt(time);

            if (Point.distance(furthest, origin) < Point.distance(closest, origin)) {
                Asteroid temp = furthestFromSun;
                furthestFromSun = closerToSun;
                closerToSun = temp;
                indexToPush = minIndex2;
            }

            furthest = furthestFromSun.orbit.positionAt(time);
            closest = closerToSun.orbit.positionAt(time);

            double mass = furthestFromSun.mass;
            double furthestAngle = Math.atan2(closest.y - furthest.y, closest.x - furthest.x);

            Point v = furthestFromSun.orbit.velocityAt(time);
            // add 5-50% of current velocity in magnitude
            double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
            double v2 = v1 * 0.2 + 0.05;

            for (double angle = furthestAngle - Math.PI/9; angle < furthestAngle + Math.PI/9; angle += Math.PI/36) {
                for (double velocity = v2; velocity < v1 * 0.5; velocity += v2 * 0.10) {
                    double pushEnergy = 0.05 * mass * velocity * velocity * 0.5;
                    if (prediction(furthestFromSun, closerToSun, pushEnergy, angle)) {
                        System.out.println("Collision!");
                        energy[indexToPush] = pushEnergy;
                        direction[indexToPush] = angle;
                        return;
                    }
                }
            }
            time_of_push = time + turns_per_retry;
            System.out.println("No collision");
        }

        // if (++time <= time_of_push) 
        // 	return;
        // for (int retry = 1 ; retry <= retries_per_turn ; ++retry) {
        // 	// PIck a random asteroid and get its velocity
        // 	int i = random.nextInt(asteroids.length);
        // 	Point v = asteroids[i].orbit.velocityAt(time);
        // 	// add 5-50% of current velocity in magnitude
        // 	double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
        // 	double v2 = v1 * (random.nextDouble() * 0.45 + 0.05);
        // 	// apply push at -π/8 to π/8 of current angle
        // 	double d1 = Math.atan2(v.y, v.x);
        // 	double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
        // 	// compute energy
        // 	double E = 0.5 * asteroids[i].mass * v2 * v2;
        // 	// try to push asteroid
        // 	Asteroid a1 = null;
        // 	try {
        // 		a1 = Asteroid.push(asteroids[i], time, E, d2);
        // 	} catch (InvalidOrbitException e) {
        // 		System.out.println("  Invalid orbit: " + e.getMessage());
        // 		continue;
        // 	}
        // 	// avoid allocating a new Point object for every position
        // 	Point p1 = v, p2 = new Point();
        // 	// search for collision with other asteroids
        // 	for (int j = 0 ; j != asteroids.length ; ++j) {
        // 		if (i == j) continue;
        // 		Asteroid a2 = asteroids[j];
        // 		double r = a1.radius() + a2.radius();
        // 		// look 10 years in the future for collision
        // 		for (long ft = 0 ; ft != 3650 ; ++ft) {
        // 			long t = time + ft;
        // 			if (t >= time_limit) break;
        // 			a1.orbit.positionAt(t - a1.epoch, p1);
        // 			a2.orbit.positionAt(t - a2.epoch, p2);
        // 			// if collision, return push to the simulator
        // 			if (Point.distance(p1, p2) < r) {
        // 				energy[i] = E;
        // 				direction[i] = d2;
        // 				// do not push again until collision happens
        // 				time_of_push = t + 1;
        // 				return;
        // 			}
        // 		}
        // 	}
        // 	System.out.println("  No collision ...");
        // }
        // time_of_push = time + turns_per_retry;
    }

    public boolean prediction(Asteroid source, Asteroid target, double energy, 
            double direction) {
        try {
            source = Asteroid.push(source, time, energy, direction);
        } catch (InvalidOrbitException e) {
            System.out.println("  Invalid orbit: " + e.getMessage());
        }
        // avoid allocating a new Point object for every position
        // search for collision with other asteroids
        Point p1 = source.orbit.velocityAt(time);
        Point p2 = new Point();
        double r = source.radius() + target.radius();
        // look 10 years in the future for collision
        for (long ft = 0 ; ft != 3650; ++ft) {
            long t = time + ft;
            if (t >= time_limit) 
                break;
            source.orbit.positionAt(t - source.epoch, p1);
            target.orbit.positionAt(t - target.epoch, p2);
            // if collision, return push to the simulator
            if (Point.distance(p1, p2) < r) {
                System.out.println("Collision predicted at time " + t);
                time_of_push = t + 1;
                return true;
            }
        }
        return false;
    }
}
