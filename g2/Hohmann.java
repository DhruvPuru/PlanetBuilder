package pb.g2;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;

public class Hohmann {

	/**
	 * Generate a Hohmann Transfer push of asteroid a1 into a2.
	 * Assumes both are currently in circular orbits.
	 */
	public static Push generatePush(Asteroid a1, int asteroid_to_push_index, Asteroid a2, long time) {
		double r1 = a1.orbit.a;
		double r2 = a2.orbit.a;
		double dv = Math.sqrt(pb.sim.Orbit.GM / r1) * (Math.sqrt(2 * r2 / (r1 + r2)) - 1);
		long expected_collision_time = (long) Math.ceil(Math.PI * Math.sqrt(Math.pow(r1 + r2, 3) / (8 * Orbit.GM)) / Orbit.dt());
		double energy = a1.mass * Math.pow(dv, 2) / 2;
		double direction = a1.orbit.velocityAt(time - a1.epoch).direction();
		if (dv < 0)
			direction += Math.PI;
        //System.out.println("expected collision time: " + expected_collision_time);
		return new Push(a1, asteroid_to_push_index, energy, direction, expected_collision_time);
	}

	/**
	 * Generate a Hohmann Transfer push to correct an elliptical orbit
	 * to a circular one with radius = semi-major axis length of ellipse.
	 */
	public static Push generateCorrection(Asteroid asteroid, int asteroid_to_push_index, long time) {
		Point p = asteroid.orbit.positionAt(time - asteroid.epoch);

		Point v1 = new Orbit(p).velocityAt(0); // Velocity for round

		Point v = asteroid.orbit.velocityAt(time - asteroid.epoch);
		Point dv = new Point(v1.x - v.x, v1.y - v.y);

		double energy = asteroid.mass * Math.pow(dv.magnitude(), 2) / 2;
		double direction = dv.direction();
		return new Push(asteroid, asteroid_to_push_index, energy, direction, -1);
	}

	/**
	 * Return the time after current time when asteroid a can be pushed to b with Hohmann Transfer.
 	 */
	public double timeToPush(long time, Asteroid a, Asteroid b) {
		double ra = a.orbit.a;
		double rb = b.orbit.a;

		double angle = Math.PI * (1 - Math.pow(1 + ra / rb, 1.5) / Math.sqrt(8));

		double aa = a.orbit.positionAt(time - a.epoch).direction();
		double ba = b.orbit.positionAt(time - b.epoch).direction();

		double angle_now = ba - aa;
		if (angle_now < 0) angle_now += Math.PI * 2;

		double alphaa = a.orbit.velocityAt(time - a.epoch).magnitude() / ra;
		double alphab = a.orbit.velocityAt(time - b.epoch).magnitude() / rb;

		if (angle_now >= angle) {
			if (alphaa >= alphab) {
				return (angle_now - angle) / (alphaa - alphab);
			} else {
				return (Math.PI * 2 - angle_now + angle) / (alphab - alphaa);
			}
		} else {
			if (alphab >= alphaa) {
				return (angle - angle_now) / (alphaa - alphab);
			} else {
				return (Math.PI * 2 - angle + angle_now) / (alphab - alphaa);
			}
		}

	}
}