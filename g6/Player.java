package pb.g6;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.ArrayList;
import java.util.Random;

public class Player implements pb.sim.Player {

	// used to pick asteroid and velocity boost randomly
	private Random random = new Random();

	// current time, time limit
	private long time = -1;
	private long time_limit = -1;

	// time until next push
	private long time_of_push = 0;

	// number of retries
	private int retries_per_turn = 1;
	private int turns_per_retry = 10;

    private int magnet_index = -1;

    private int get_biggest_asteroid_index(Asteroid[] asteroids){
        int max_mass_index = -1;
        double max_mass = 0;
        for(int i = 0; i < asteroids.length; i++){
            if(asteroids[i].mass > max_mass){
                max_mass_index = i;
                max_mass = asteroids[i].mass;
            }
        }
        return max_mass_index;
    }

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit) {
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
	}

	// try to push asteroid
	public void play(Asteroid[] asteroids,
					 double[] energy, double[] direction) {
		Collision collision = new Collision(asteroids);
		// if not yet time to push do nothing
		if (++time <= time_of_push) return;
		System.out.println("Year: " + (1 + time / 365));
		System.out.println("Day: " + (1 + time % 365));

        int a2_index = get_biggest_asteroid_index(asteroids);
        Asteroid a2 = asteroids[a2_index];

		for (int retry = 1; retry <= retries_per_turn; ++retry) { 
            
            if (a2_index ==-1) {
                System.out.println("NOOOOOOOOOOO INDEX ERROR");
                System.exit(0);
            }
            for (int i = 0; i != asteroids.length; ++i) {
                // TODO: to check magnet here!!!!
                if (i == a2_index) continue;
    			// pick a random asteroid and get its velocity
    			// int i = random.nextInt(asteroids.length);
    			Point v = asteroids[i].orbit.velocityAt(time);
    			// add 2-10% of current velocity in magnitude
    			System.out.println("Try: " + retry + " / " + retries_per_turn);
    			double v1 = v.magnitude();
    			double v2 = v1 * (random.nextDouble() * 0.08 + 0.02);
    			System.out.println("  Speed: " + v1 + " +/- " + v2);
    			double d1 = v.direction();
    			// double d2 = Utils.getPerpendicularAngle(d1);
                double d2 = d1;
    			System.out.println("  Angle: " + d1 + " -> " + d2);
    			// compute energy
    			double E = 0.5 * asteroids[i].mass * v2 * v2;
    			// try to push asteroid
    			Asteroid a1 = null;
    			try {
    				a1 = Asteroid.push(asteroids[i], time, E, d2);
    			} catch (InvalidOrbitException e) {
    				System.out.println("  Invalid orbit: " + e.getMessage());
    				continue;
    			}
    			// avoid allocating a new Point object for every position
    			Point p1 = v, p2 = new Point();

    			// search for collision with other asteroids
			
				// Asteroid a2 = asteroids[j];

                double r = a1.radius() + a2.radius();
				for (long ft = 0; ft != 3650; ++ft) {
					long t = time + ft;
					if (t >= time_limit) break;
					a1.orbit.positionAt(t - a1.epoch, p1);
					a2.orbit.positionAt(t - a2.epoch, p2);
					// if collision, return push to the simulator
					if (Point.distance(p1, p2) < r) {
						energy[i] = E;
						direction[i] = d2;
						// do not push again until collision happens
						time_of_push = t + 1;
						System.out.println("  Collision prediction !");
						System.out.println("  Year: " + (1 + t / 365));
						System.out.println("  Day: " + (1 + t % 365));
						return;
					}
				}

			}
			System.out.println("  No collision ...");
		}
		time_of_push = time + turns_per_retry;
	}
}