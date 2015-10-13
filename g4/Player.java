package pb.g4;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Random;
import java.util.*;

public class Player implements pb.sim.Player{
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

	private int num_closest_asteroids = 5;
	private int initial_number_of_asteroids;

	//stores asteroid masses
	private HashMap<Asteroid, Double> cached_asteroid_masses = new HashMap<Asteroid, Double>();

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit)
	{
		initial_number_of_asteroids = asteroids.length;
		asteroidOrder = new HashSet<Asteroid>();
		storeMass(asteroids);
		System.out.println("Init");
		dynamicProgramming(asteroids);
		System.out.println(asteroidOrder.size());
		for(Asteroid a: asteroidOrder)
		{
			System.out.println(a.mass);
		}
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
	}
	public void storeMass(Asteroid[] asteroids)
	{
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

	public void updateMass(Asteroid asteroid1, Asteroid asteroid2, Asteroid[] asteroids)
	{
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

	private void printMassVelocity(Asteroid[] asteroids)
	{
		for(Asteroid asteroid: asteroids)
		{
			System.out.println("mass, velocity:" + asteroid.mass + ", " + asteroid.orbit.velocityAt(time));
		}

	}

	public void dynamicProgramming(Asteroid[] asteroids)
	{
		System.out.println("Starting Dynamic Programming");
		ArrayList<Double> masses = new ArrayList<Double>(); //exp
		ArrayList<Double> energies = new ArrayList<Double>(); //stam
		masses.add(0.0);
		energies.add(0.0);
		for(double i = Math.pow(10, 35); i <= Math.pow(10, 39); i += Math.pow(10,35))
		{
			boolean write = false;
			for(int j = 0; j < asteroids.length; j++)
			{
				Point v = asteroids[j].orbit.velocityAt(time);
				Double velocity = Math.sqrt(v.x * v.x + v.y * v.y);
				Double mass = asteroids[j].mass;
				Double energy = 0.5* mass * velocity * velocity;
				// System.out.println("energy: " + energy);
				if(energy >= i && energy < (i + Math.pow(10, 35)))
				{
					masses.add(asteroids[j].mass);
					energies.add(energy);
					write = true;
				}
			}
			if (write == false)
			{
				masses.add(0.0);
				energies.add(0.0);
			}
		}
		optimize(masses, asteroids);
	}
	public void optimize(ArrayList<Double> mass, Asteroid[] asteroids)
	{
		System.out.println("Optimizing");
		ArrayList<Double> r = new ArrayList<Double>();
		String results = null;
		ArrayList<String[]> r1 = new ArrayList<String[]>();
		Set<Asteroid> r2 = new HashSet<Asteroid>();
		ArrayList<String> r3 = new ArrayList<String>();
		for(double i: mass)
		{
			r3.add("0");
			r.add(0.0);
		}
		for (int i = 0; i < mass.size(); i++)
		{
			double q = 0.0;
			String q1 = "0 0 0";
			for (int j = 0; j <= i; j++ )
			{
				if(q < (mass.get(j) + r.get(i-j)))
				{
					q =  mass.get(j) + r.get(i-j);
					q1 = mass.get(j)+ " " + r3.get(i-j);
					if( mass.get(j)!= 0 && r.get(i-j)!=0)
					{
						results = mass.get(j) + " " + r.get(i-j) + " " + q + "\n";
						r1.add(results.split(" "));
					}
				}
			}
			r.set(i,q);
			r3.set(i,q1);
		}
		String[] last = r3.get(r3.size()-1).split(" ");
		for(int k = 0; k < asteroids.length; k++)
		{
			for(int j = 0; j < last.length; j++)
			{
				if(Double.parseDouble(last[j]) == asteroids[k].mass)
				{
					r2.add(asteroids[k]);
				}
			}
		}
		asteroidOrder = r2;
	}

	// try to push asteroid
	public void play(Asteroid[] asteroids,
		double[] energy, double[] direction)
	{
		num_closest_asteroids = Math.min(num_closest_asteroids, asteroids.length - 1);
		if (++time%10 == 0 && time > collisionTime) {
			// if(cached_asteroid_masses.size() < asteroids.length && indexToPush != -1 && indexToHit != -1)
			// {
			// 	updateMass(furthestFromSun, closerToSun, asteroids);
			// }
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

	public void perturb(Asteroid[] asteroids, double[] energy, double[] direction)
	{
		Asteroid a1;
		//It's been a long time since a push, so perturb the system
		if (timeSincePush > 7300) {
			boolean validOrbitNotFound = true;
			while (validOrbitNotFound) {
				int i = random.nextInt(asteroids.length);
				Point v = asteroids[i].orbit.velocityAt(time);
				// add 5-50% of current velocity in magnitude
				double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
				double v2 = v1 * (random.nextDouble() * 0.45 + 0.05);
				// apply push at -π/8 to π/8 of current angle
				double d1 = Math.atan2(v.y, v.x);
				double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
				// compute energy
				double E = 0.5 * asteroids[i].mass * v2 * v2;
				try {
					a1 = Asteroid.push(asteroids[i], time, E, d2);
					validOrbitNotFound = false;
				} catch (InvalidOrbitException e) {
					System.out.println("  Invalid orbit: " + e.getMessage());
					validOrbitNotFound = false;
					continue;
				}
				energy[i] = E;
				direction[i] = d2;
			}
		}
	}

	public void push_closest_to_largest(Asteroid[] asteroids, double[] energy, double[] direction)
	{
		PriorityQueue<Asteroid> heap = new PriorityQueue<Asteroid>(asteroids.length, new AsteroidComparator());
		int largest_asteroid_idx = findLargestAsteroidIndex(asteroids);
		double min = Double.MAX_VALUE;
		int minIndex = -1;
		// if not yet time to push do nothing
		largest_asteroid = asteroids[largest_asteroid_idx];
		Point aPoint = largest_asteroid.orbit.positionAt(time);

		for (int i = 0; i < asteroids.length; i++) {
			if (i != largest_asteroid_idx) 
				heap.add(asteroids[i]);
		}

		Asteroid lowestEnergyAsteroid = asteroids[0];
		double leastEnergy = Double.MAX_VALUE;
		double bestAngle = 0;

		for (int i = 0; i < num_closest_asteroids; i++) {
			Asteroid other_asteroid = heap.poll();

			Point origin = new Point(0, 0);
			Point largest = largest_asteroid.orbit.positionAt(time);
			Point closest = other_asteroid.orbit.positionAt(time);

			double mass = other_asteroid.mass;
			double arc = Math.atan2(closest.x - largest.x, closest.y - largest.y);

			Point v = other_asteroid.orbit.velocityAt(time);
			// add 5-50% of current velocity in magnitude
			double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
			double v2 = v1 * 0.2 + 0.05;

			int loopCount = 0;
			for (double angle = arc - Math.PI/9; angle < arc + Math.PI/9; angle += Math.PI/36) {
				for (double velocity = v2; velocity < 0.5 * v1; velocity += v2 * 0.10) {
					double pushEnergy = 0.05 * mass * velocity * velocity * 0.5;
					loopCount++;
					if (prediction(other_asteroid, largest_asteroid, time, pushEnergy, angle)) {
						System.out.println("collision " + " at energy: "
							+ pushEnergy + " and direction: " + angle + " at year: " + time / 365);
						timeSincePush = 0;

						if (pushEnergy < leastEnergy) {
							lowestEnergyAsteroid = other_asteroid;
							leastEnergy = pushEnergy;
							bestAngle = angle;
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

	public boolean prediction(Asteroid source, Asteroid target, long time, double energy, 
		double direction) {
		try {
			source = Asteroid.push(source, time, energy, direction);
		} catch (InvalidOrbitException e) {
			e.printStackTrace();
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
				collisionTime = t;
				return true;
			}
		}
		return false;
	}

	public class AsteroidComparator implements Comparator<Asteroid> {
		public int compare(Asteroid a1, Asteroid a2) {
			double a1Distance = Point.distance(a1.orbit.positionAt(time), largest_asteroid.orbit.positionAt(time));
			double a2Distance = Point.distance(a2.orbit.positionAt(time), largest_asteroid.orbit.positionAt(time));

			if ( a1Distance < a2Distance ) {
				return -1;
			} else if ( a2Distance < a1Distance ) {
				return 1;
			} else {
				return 0;
			}
		}
	}
}
