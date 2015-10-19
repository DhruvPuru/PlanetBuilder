package pb.sim;

public class Asteroid {

	// the density of all asteroids (default is Sun density)
	static double density = 1410.0;

	// unique ids for asteroids
	private static long serial_id = -1;

	// mass of the asteroid
	public final double mass;

	// creation time (used to find position)
	public final long epoch;

	// unique ID of asteroid
	public final long id;

	// orbit of the asteroid
	public final Orbit orbit;

	// internal constructor that sets all parameters including id
	private Asteroid(Orbit orbit, double mass, long epoch, long id)
	{
		if (orbit == null || mass <= 0 || epoch < 0)
			throw new IllegalArgumentException();
		this.orbit = orbit;
		this.mass  = mass;
		this.epoch = epoch;
		this.id = id;
	}

	// create a asteroid by setting all parameters
	public Asteroid(Orbit orbit, double mass, long epoch)
	{
		this(orbit, mass, epoch, ++Asteroid.serial_id);
	}

	// radius of the asteroid
	public double radius()
	{
		double volume = mass / density;
		return Math.cbrt(0.75 * volume / Math.PI);
	}

	// find which (sets of) asteroids collide at given time (ignore nulls)
	public static int[][] test_collision(Asteroid[] asteroids, long time)
	{
		if (time < 0) throw new IllegalArgumentException();
		Point[]  center = new Point  [asteroids.length];
		double[] radius = new double [asteroids.length];
		for (int i = 0 ; i != asteroids.length ; ++i) {
			if (asteroids[i] == null) continue;
			// ignore asteroids not yet created
			long t = time - asteroids[i].epoch;
			if (t < 0) continue;
			// store position and radius
			center[i] = new Point();
			asteroids[i].orbit.positionAt(t, center[i]);
			radius[i] = asteroids[i].radius();
		}
		return Point.overlaps(center, radius);
	}

	// force collision of asteroids without distance checks
	public static Asteroid force_collision(Asteroid[] asteroids, long time)
	{
		double m = 0, mr_x = 0, mr_y = 0, mv_x = 0, mv_y = 0;
		// compute weighted average of position and velocity
		Point r = new Point();
		Point v = new Point();
		for (int i = 0 ; i != asteroids.length ; ++i) {
			long t = time - asteroids[i].epoch;
			asteroids[i].orbit.positionAt(t, r);
			asteroids[i].orbit.velocityAt(t, v);
			// weighted average of position is the barycenter
			mr_x += r.x * asteroids[i].mass;
			mr_y += r.y * asteroids[i].mass;
			// weighter average of velocity is the momentum
			mv_x += v.x * asteroids[i].mass;
			mv_y += v.y * asteroids[i].mass;
			// sum the mass
			m += asteroids[i].mass;
		}
		r.x = mr_x / m;  r.y = mr_y / m;
		v.x = mv_x / m;  v.y = mv_y / m;
		return new Asteroid(new Orbit(r, v), m, time);
	}

	// push asteroid to specific direction
	public static Asteroid push(Asteroid asteroid, long time,
	                            double energy, double direction)
	{
		if (Double.isNaN(energy) || Double.isInfinite(energy)
		                         || energy < 0.0)
			throw new IllegalArgumentException("Invalid energy");
		if (Double.isNaN(direction) || Double.isInfinite(direction))
			throw new IllegalArgumentException("Invalid direction");
		// find current position and velocity of asteroid
		long t = time - asteroid.epoch;
		Point r = asteroid.orbit.positionAt(t);
		Point v = asteroid.orbit.velocityAt(t);
		// translate push energy to velocity and combine
		double magnitude = Math.sqrt(2.0 * energy / asteroid.mass);
		v.x += magnitude * Math.cos(direction);
		v.y += magnitude * Math.sin(direction);
		// return (new object of) the "same" asteroid with new orbit
		Orbit orbit = new Orbit(r, v);
		return new Asteroid(orbit, asteroid.mass, time, asteroid.id);
	}
}
