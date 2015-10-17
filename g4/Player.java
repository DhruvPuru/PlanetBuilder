package pb.g4;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Random;
import java.util.*;

public class Player implements pb.sim.Player{

	public static double dt = 24 * 60 * 60;
	public static Point sun = new Point(0, 0);
	public static final long NUMBER_OF_YEARS_AHEAD = 100;
	public static final long TIME_THRESHOLD = 5000;

	// used to PIck asteroid and velocity boost randomly
	private Random random = new Random();

	// current time, time limit
	private long time = -1;
	private long time_limit = -1; 

	// time until next push
	private long time_of_push = 0;

	// number of retries
	private int retries_per_turn = 1;
	private int turns_per_retry = 3;

	private int timeSincePush = 0;
	private double collisionTime = -1;

	private Asteroid furthestFromSun;
	private Asteroid closerToSun;
	private Asteroid largest_asteroid;
	private int indexToPush = -1;
	private int indexToHit = -1;
	private Set<Asteroid> asteroidOrder;
	private double fifty_percent_mass;

	private int num_closest_asteroids = 3;
	private int initial_number_of_asteroids;

	//stores asteroid masses
	private HashMap<Asteroid, Double> cached_asteroid_masses = new HashMap<Asteroid, Double>();

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit)
	{
		initial_number_of_asteroids = asteroids.length;
		storeMass(asteroids);
		System.out.println("Init");

		//Printing out info in case it tells me (Dhruv) something 
		for (Asteroid a : asteroids) {
			// System.out.println("Orbital radius to body radius ratio:"+ 
			// 	Point.distance(a.orbit.positionAt(0), sun) / a.radius());

			double dt = 24 * 60 * 60;
			double t = (Math.PI + Math.PI) * Math.sqrt(a.orbit.a / Orbit.GM) * a.orbit.a;
			double t_dt = (long) Math.ceil(t / dt);
			// System.out.println("Time period: " + t_dt);
		}

		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
	}

	public void storeMass(Asteroid[] asteroids) {
		double mass_sum = 0;
		for(Asteroid asteroid: asteroids)
		{
			cached_asteroid_masses.put(asteroid, asteroid.mass);
			mass_sum += asteroid.mass;
			System.out.println(asteroid.mass);
		}
		fifty_percent_mass = 0.5*mass_sum;
		System.out.println("50% mass: " + fifty_percent_mass);
	} 

	public void updateMass(Asteroid asteroid1, Asteroid asteroid2, Asteroid[] asteroids) {
		cached_asteroid_masses.remove(asteroid1);
		cached_asteroid_masses.remove(asteroid2);
		for(Asteroid asteroid: asteroids)
		{
			if(!cached_asteroid_masses.containsKey(asteroid))
			{
				cached_asteroid_masses.put(asteroid, asteroid.mass);
			}
		}
	} 

	private void printMassVelocity(Asteroid[] asteroids) {
		for (Asteroid asteroid: asteroids) {
			System.out.println("mass, velocity:" + asteroid.mass + ", " + asteroid.orbit.velocityAt(time));
		}
	}

	// try to push asteroid
	public void play(Asteroid[] asteroids,
		double[] energy, double[] direction) {
		num_closest_asteroids = Math.min(num_closest_asteroids, asteroids.length - 1);
		if (++time%100 == 0 && time > collisionTime) {
			push_closest_to_largest(asteroids, energy, direction);
		}
	}

	public int findLargestAsteroidIndex(Asteroid[] asteroids)
	{
		double max_mass = 0;
		int idx = 0;
		for(int i = 0; i < asteroids.length; i++)
		{
			if(asteroids[i].mass > max_mass)
			{
				max_mass = asteroids[i].mass;
				idx = i;
			}
		}
		return idx;
	}

	public void push_closest_to_largest(Asteroid[] asteroids, double[] energy, double[] direction)
	{
		if (time % 365 == 0) {
			System.out.println("Year: " + time / 365);
		}
		PriorityQueue<Asteroid> heap = new PriorityQueue<Asteroid>(asteroids.length, new AsteroidComparator());
		int largest_asteroid_idx = findLargestAsteroidIndex(asteroids);
		double min = Double.MAX_VALUE;
		int minIndex = -1;
		// if not yet time to push do nothing
		largest_asteroid = asteroids[largest_asteroid_idx];
		Point largestAsteroidPosition = largest_asteroid.orbit.positionAt(time);
		double largestAsteroidDistFromSun = Point.distance(largestAsteroidPosition, sun);

		for (int i = 0; i < asteroids.length; i++) {
			if (i != largest_asteroid_idx) 
				heap.add(asteroids[i]);
		}

		Asteroid lowestEnergyAsteroid = asteroids[0];
		double leastEnergy = Double.MAX_VALUE;
		double bestAngle = 0;

		for (int i = 0; i < num_closest_asteroids; i++) {
			Asteroid other_asteroid = heap.poll();
			Point closest = other_asteroid.orbit.positionAt(time);
			Point v = other_asteroid.orbit.velocityAt(time);
			double distBetweenPointAndSun = Point.distance(closest, sun);
			double pushAngle;

			double mass = other_asteroid.mass;
			double arc = Math.atan2(closest.x - largestAsteroidPosition.x, closest.y - largestAsteroidPosition.y);

			// add 5-50% of current velocity in magnitude
			double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
			double v2 = v1 * 0.20 + 0.05;
			collisionTime = time;

			for (pushAngle = arc - Math.PI/18; pushAngle < arc + Math.PI/18; pushAngle += Math.PI/36) {
				for (double velocity = v2; velocity < 0.4 * v1; velocity += v2 * 0.05) {
					double pushEnergy = 0.05 * mass * velocity * velocity * 0.5;
					long predictedTimeOfCollission = improvedPrediction(other_asteroid, largest_asteroid, 
						time, pushEnergy, pushAngle);
					if (predictedTimeOfCollission > 0) {
						System.out.println("collision predicted" + " at energy: "
							+ pushEnergy + " and direction: " + pushAngle + " at year: " + time / 365);
						timeSincePush = 0;

						if (pushEnergy < leastEnergy) {
							collisionTime = predictedTimeOfCollission;
							lowestEnergyAsteroid = other_asteroid;
							leastEnergy = pushEnergy;
							bestAngle = pushAngle;
						}
					}
				}
			}
		}

		int indexToPush = 0;
		for (int j = 0; j < asteroids.length; j++) {
			if (lowestEnergyAsteroid == asteroids[j]) {
				indexToPush = j;
			}
		}

		if (leastEnergy == Double.MAX_VALUE) {
			leastEnergy = 0;
		}
		energy[indexToPush] = leastEnergy;
		direction[indexToPush] = bestAngle;
	}

	public long prediction(Asteroid source, Asteroid target, long time, double energy, 
		double direction) {
		try {
			source = Asteroid.push(source, time, energy, direction);
		} catch (InvalidOrbitException e) {
			e.printStackTrace();
			return -1;
		}

		// avoid allocating a new Point object for every position
		// search for collision with other asteroids
		Point p1 = source.orbit.velocityAt(time);
		Point p2 = new Point();
		double r = source.radius() + target.radius();
			// look 10 years in the future for collision
		for (long ft = 0 ; ft != 3650 / 2; ++ft) {
			long t = time + ft;
			if (t >= time_limit) 
				break;
			source.orbit.positionAt(t - source.epoch, p1);
			target.orbit.positionAt(t - target.epoch, p2);
				// if collision, return push to the simulator
			if (Point.distance(p1, p2) < r) {
				collisionTime = t;
				return t;
			}
		}
		return -1;
	}

	public long improvedPrediction(Asteroid source, Asteroid target, long time, double energy, 
		double direction) {
		try {
			source = Asteroid.push(source, time, energy, direction);
		} catch (InvalidOrbitException e) {
			e.printStackTrace();
			return -1;
		}

		Point p1 = new Point();
		double r = source.radius() + target.radius();

		double T = (Math.PI + Math.PI) * Math.sqrt(source.orbit.a / Orbit.GM) * source.orbit.a;
		long T_dt = (long) Math.ceil(T / dt);
		// System.out.println("Time period for asteroid: " + source.id + " is: " + T_dt);

		//Find points at which the source will be on the target's orbit
		HashMap<Long, Point> timesOnOrbit = new HashMap<Long, Point>();
		for (long ft = 0; ft < Math.min(T_dt, TIME_THRESHOLD); ft++) {
			long t = time + ft;
			if (t >= time_limit || timesOnOrbit.keySet().size() == 4) 
				break;

			source.orbit.positionAt(t - source.epoch, p1);
			if (isOnOrbit(target.orbit, p1, 0.01)) {
				timesOnOrbit.put(t, p1);
			}
		}

		Point p2 = new Point();
		for (Map.Entry<Long, Point> entry : timesOnOrbit.entrySet()) {
			long timeOfIntersection = entry.getKey();
			Point pointOfIntersection = entry.getValue();

			for (int year = 1; year < NUMBER_OF_YEARS_AHEAD; year++) {
				long timeToCheck = timeOfIntersection + year * T_dt;
				target.orbit.positionAt(timeToCheck - target.epoch, p2);
				// if collision, return push to the simulator
				if (Point.distance(p1, p2) < r) {
					collisionTime = timeToCheck;
					return timeToCheck;
				}
			}
		}
		return -1;
	}

	public boolean isOnOrbit(Orbit o, Point p, double threshold) {
		double e = Math.sqrt(1 - (o.b*o.b)/(o.a*o.a));
		double c1 = o.a * e;
		double cD = o.A + Math.PI;

		double xC = c1 * Math.cos(cD);
		double yC = c1 * Math.sin(cD);

		double xVal = Math.pow((p.x - xC) * Math.cos(o.A) + (p.y-yC) * Math.sin(o.A), 2);
		double yVal = Math.pow((p.x - xC) * Math.sin(o.A) - (p.y-yC) * Math.cos(o.A), 2);

		double finalValToCheck = xVal/(o.a*o.a) + yVal/(o.b*o.b);

		if (Math.abs(1 - finalValToCheck) < threshold) {
			System.out.println("Found orbital intersection: " + finalValToCheck);
			return true;
		}
		return false;
	}

	public class AsteroidComparator implements Comparator<Asteroid> {
		public int compare(Asteroid a1, Asteroid a2) {
			double a1Distance = Point.distance(a1.orbit.positionAt(time), largest_asteroid.orbit.positionAt(time));
			double a2Distance = Point.distance(a2.orbit.positionAt(time), largest_asteroid.orbit.positionAt(time));

			if (a1Distance < a2Distance ) {
				return -1;
			} else if (a2Distance < a1Distance ) {
				return 1;
			} else {
				return 0;
			}
		}
	}
}
