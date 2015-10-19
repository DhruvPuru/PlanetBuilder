package pb.g9;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.*;

public class Player implements pb.sim.Player {

	// iteration number
	int iteration = 0;

	// used to pick asteroid and velocity boost randomly
	private Random random = new Random();

	// current time, time limit
	private long time = -1;
	private long time_limit = -1;

	// time until next push
	private long time_of_push = 0;

	// number of retries
	private int retries_per_turn = 1;
	private int turns_per_retry = 0;

	private Point origin = new Point(0, 0);

	private onePush bestpush;
	private long wait_time = 100;
	private double average_energy = Double.MAX_VALUE/20.0;
	private int push_times = 0;

	private ArrayList<Point> find_intersection(Asteroid a, Asteroid b, HashSet<Long> timelist){
		ArrayList<Point> intersection_list = new ArrayList<Point>();
		long period = a.orbit.period();
		if (period / 365.0 > 50) {
			return intersection_list;
		}
		if (period > 3650)
			period = 3650;
		double r = a.radius() + b.radius();
		for (long ft = 0 ; ft <= period ; ++ft) {
			long t = time + wait_time + ft;
			Point p1 = new Point();
			Point c = new Point();
			a.orbit.positionAt(t - a.epoch, p1);
			b.orbit.center(c);
			Point foci = new Point(c.x*2,c.y*2);
			double dist = Point.distance(p1,foci)+Point.distance(p1,new Point(0,0));
			if (Math.abs(dist - b.orbit.a*2) < r){
				intersection_list.add(p1);
				timelist.add(t);
			}
				
		}
		return intersection_list;
	}
	
	private ArrayList<Point> fast_find_intersection(Asteroid a, Asteroid b, HashSet<Long> timelist){
		ArrayList<Point> intersection_list = new ArrayList<Point>();
		long period = a.orbit.period();
		long origin_period = a.orbit.period();
		if (period / 365.0 > 20) {
			return intersection_list;
		}
		//System.out.println("  fast_find_intersection begin !");
		if (period > 3650)
			period = 3650;
		double r = a.radius() + b.radius();
		Point c = new Point();
		b.orbit.center(c);
		Point foci = new Point(c.x*2,c.y*2);
		long step_length = Math.round(Math.sqrt(period));
		double[] diff = new double[3];
		
		long t = time + wait_time;
		Point p1 = new Point();
		a.orbit.positionAt(t - a.epoch, p1);

		double dist = Point.distance(p1,foci)+Point.distance(p1,new Point(0,0));
		diff[1] = dist - b.orbit.a * 2;
		long ft;
		for (ft = step_length ; ft <= period ; ft+=step_length) {
			t = time + wait_time + ft;
			p1 = new Point();
			a.orbit.positionAt(t - a.epoch, p1);
			dist = Point.distance(p1,foci)+Point.distance(p1,new Point(0,0));
			diff[2] = dist - b.orbit.a * 2;
			if (diff[1] >= 0 && diff[2] <0 || diff[1] <= 0 && diff[2] >0){
				for(long ft1 = ft-step_length;ft1 < ft;ft1++){
					//System.out.print(ft1+" ("+ft+") ");
					t = time + wait_time + ft1;
					p1 = new Point();
					a.orbit.positionAt(t - a.epoch, p1);
					dist = Point.distance(p1,foci)+Point.distance(p1,new Point(0,0));
					if (Math.abs(dist - b.orbit.a*2) < r){
						intersection_list.add(p1);
						timelist.add(t%origin_period);
						break;
					}
				}
			}
			diff[0] = diff[1];
			diff[1] = diff[2];
		}

		//System.out.println("  fast_find_intersection mid !");
		t = time + wait_time + period;
		p1 = new Point();
		a.orbit.positionAt(t - a.epoch, p1);		
		dist = Point.distance(p1,foci)+Point.distance(p1,new Point(0,0));
		diff[2] = dist - b.orbit.a * 2;
		if (diff[1] >= 0 && diff[2] <0 || diff[1] <= 0 && diff[2] >0){
			for(long ft1 = ft-step_length;ft1 < period;ft1++){
				t = time + wait_time + ft1;
				p1 = new Point();
				a.orbit.positionAt(t - a.epoch, p1);
				dist = Point.distance(p1,foci)+Point.distance(p1,new Point(0,0));
				if (Math.abs(dist - b.orbit.a*2) < r){
					intersection_list.add(p1);
					timelist.add(t%origin_period);
				}
			}
		}
		//System.out.println("  fast_find_intersection end !");
		return intersection_list;
	}

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit) {

		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
		bestpush = new onePush(-1, Double.MAX_VALUE, 0, 0, 0, 0, 0);
	}

	// try to push asteroid
	public void play(Asteroid[] asteroids, double[] energy, double[] direction) {

		time++;
		if (time == bestpush.time && bestpush.energy< 5 * average_energy){
			push_times++;
			average_energy = average_energy*((push_times-1.0)/push_times)+bestpush.energy/push_times;
		 	int i = bestpush.i;
		 	energy[i] = bestpush.energy;
		 	direction[i] = bestpush.direction;

		 	 // do not push again until collision happens
		 	time_of_push = bestpush.collision_time + 1;
		 	bestpush = new onePush(-1, Double.MAX_VALUE,0, 0, 0, 0, 0);
		 	iteration++;
		 	System.out.println("average energy: "+ average_energy);
		 	return;
		}
		// if not yet time to push do nothing
		if (time <= time_of_push)
			return;
		//System.out.println("Year: " + (1 + time / 365));
		//System.out.println("Day: " + (1 + time % 365));
		for (int retry = 1; retry <= retries_per_turn; ++retry) {
			// pick a random asteroid and get its velocity
			// int i = random.nextInt(asteroids.length);

			//System.out.println("Try: " + retry + " / " + retries_per_turn);

			int j;
			if (iteration == 1) {
				List<Integer> desiredOrbits = findMiddleOrbits(asteroids);
				j = getHeaviestAsteroidAmong(asteroids, desiredOrbits);
			}
			else
				j = getHeaviestAsteroid(asteroids);

			// for (int i = 0; i != asteroids.length; ++i) 

			// 10% of the total number of asteroids
			int asteroidsToConsider = asteroids.length/10;
			if (asteroidsToConsider < 1)
				asteroidsToConsider = 1;

			List<Integer> favorableAsteroids = getKHighestWeight(asteroids, asteroidsToConsider);
			for (int i : favorableAsteroids) {
				if (i == j)
					continue;

				Point v = asteroids[i].orbit.velocityAt(time + wait_time
						- asteroids[i].epoch);
				// add 5-50% of current velocity in magnitude
				double v1 = Math.sqrt(v.x * v.x + v.y * v.y);

//				double E_i= Orbit.GM/(2*asteroids[i].orbit.a);
//				double E_j= Orbit.GM/(asteroids[j].orbit.a+asteroids[i].orbit.a);
//				double E_diff = asteroids[i].mass*Math.abs(E_i-E_j);
				for (double k = 0.05; k < .8; k = k + 0.005) {
					//for (double k2 = 0; k2 < 5; k2 = k2 + 0.1) {
						double v2 = v1 * k;
						//System.out.println("  Speed: " + v1 + " +/- " + v2);

						// apply push at -/8 to /8 of current angle
						double d1 = Math.atan2(v.y, v.x);
						double d2;
						if (asteroids[j].orbit.a < asteroids[i].orbit.a)
							//d2 = Math.PI + d1 + (k2 - 0.5) * Math.PI * 0.25;
							d2 = Math.PI + d1;
						else
							//d2 = d1 + (k2 - 0.5) * Math.PI * 0.25;
							d2 = d1;
						//System.out.println("  Angle: " + d1 + " -> " + d2);

						// compute energy
						double E = 0.5 * asteroids[i].mass * v2 * v2;
						//double E = k*E_diff;

						// try to push asteroid
						Asteroid a1 = null;
						try {
							a1 = Asteroid.push(asteroids[i], time + wait_time
									- asteroids[i].epoch, E, d2);
						} catch (InvalidOrbitException e) {
//							System.out.println("  Invalid orbit: "
//							+ e.getMessage());
							continue;
						}

						// avoid allocating a new Point object for every
						// position
						Point p1 = v, p2 = new Point();

						// search for collision with any other asteroids
						Asteroid a2 = asteroids[j];
						HashSet<Long> timelist = new HashSet<Long>();

						ArrayList<Point> intersections = fast_find_intersection(a1,
								a2, timelist);
						if (intersections.size() == 0)
							continue;
						double r = a1.radius() + a2.radius();

						// look 10 years in the future for collision
						for (long ft = 0; ft < 3650; ++ft) {

							long t = time + wait_time + ft;
							if (t >= time_limit)
								break;
							if (!timelist.contains(t%a1.orbit.period()))
								continue;
							a1.orbit.positionAt(t - a1.epoch, p1);
							a2.orbit.positionAt(t - a2.epoch, p2);

							// if collision, return push to the simulator
							if (Point.distance(p1, p2) < r) {

								onePush push = new onePush(time + wait_time, E,
										d2, i, a1.id, t, a2.id);
									
								if (E < bestpush.energy) {
								//	if (CheckIncidentalCollisions(push, asteroids)){
										bestpush = push;
								//	}
								}
//								energy[i] = E;
//								direction[i] = d2;
//
//								// do not push again until collision happens
//								time_of_push = t + 1;
//								iteration++;
//								System.out.println("  Collision prediction !");
//								System.out.println("  Year: " + (1 + t / 365));
//								System.out.println("  Day: " + (1 + t % 365));

//								return;
							}
						}
					//}
				}
				// System.out.println("  No collision ...");
			}
		}
		time_of_push = time + turns_per_retry;
	}

	private List<Integer> findMiddleOrbits(Asteroid[] asteroids) {
		Map<Double, Integer> radiusToIndexUnsorted = new HashMap<Double, Integer>();

		for (int i = 0; i < asteroids.length; i++) {
			radiusToIndexUnsorted.put(asteroids[i].radius(), i);
		}

		Map<Double, Integer> radiusToIndex = new TreeMap<Double, Integer>(radiusToIndexUnsorted);

		int i = 0;
		int n = asteroids.length;
		Iterator<Map.Entry<Double, Integer>> iter = radiusToIndex.entrySet().iterator();
		List<Integer> desiredOrbits = new ArrayList<Integer>();
		while (i < (n/2) - (n/10)) {
			if (iter.hasNext())
				iter.next();
			i++;
		}
		while (i <= (n/2) + (n/10)) {
			if (iter.hasNext())
				desiredOrbits.add(iter.next().getValue());
			i++;
		}
		return desiredOrbits;
	}

	private int getHeaviestAsteroidAmong(Asteroid asteroids[], List<Integer> desiredOrbits) {
		int max = 0;
		double mass = 0.0;
		double maxMass = Double.MIN_VALUE;

		for (int i: desiredOrbits) {
			mass = asteroids[i].mass;
			if (mass > maxMass) {
				max = i;
				maxMass = mass;
			}
		}
		return max;
	}

	private int getHeaviestAsteroid(Asteroid asteroids[]) {
		int max = 0;
		double mass = 0.0;
		double maxMass = Double.MIN_VALUE;

		for (int i = 0; i < asteroids.length; i++) {
			mass = asteroids[i].mass;
			if (mass > maxMass) {
				max = i;
				maxMass = mass;
			}
		}
		return max;
	}

	private int getHighestWeightAsteroid(Asteroid asteroids[]) {
		/** 
		* Returns the index of the asteroid with the maximum distance to mass ratio.
		* (Higher the distance mass ratio, more "pushable" the asteroid - costs less energy)
		**/

		int max = 0;
		double distance = 0.0;
		double mass = 0.0;
		double weight = 0.0;
		double maxWeight = Double.MIN_VALUE;

		for (int i = 0; i < asteroids.length; i++) {

			Asteroid a = asteroids[i];
			Point a_center = new Point();
			a.orbit.positionAt(time - a.epoch, a_center);
			distance = Double
					.valueOf(Point.distance(a_center, new Point(0, 0)));
			mass = asteroids[i].mass;
			weight = distance / mass;
			if (weight > maxWeight) {
				max = i;
				maxWeight = weight;
			}
		}
		return max;
	}

	private List<Integer> getKLowestWeight(Asteroid[] asteroids, int k) {
		double[] weights = new double[asteroids.length];
		double distance = 0.0;
		double mass = 0.0;

		for (int i = 0; i < asteroids.length; i++) {
			Asteroid a = asteroids[i];
			Point a_center = new Point();
			a.orbit.positionAt(time - a.epoch, a_center);
			distance = Double
					.valueOf(Point.distance(a_center, new Point(0, 0)));
			mass = asteroids[i].mass;
			weights[i] = distance / mass;
		}

		List<Integer> lowestKIndices;
		double[] temp = new double[k];
		double maxWeight = Integer.MIN_VALUE;
		int maxIndex = -1;
		Integer[] tempIndexInWeights = new Integer[k];

		// put k weights in temp, and keep track of max seen so far
		for (int i = 0; i < k; i++) {
			temp[i] = weights[i];
			tempIndexInWeights[i] = i;
			if (temp[i] > maxWeight) {
				maxWeight = temp[i];
				maxIndex = i;
			}
		}

		// check remaining weights and maintain lowest k in temp
		for (int i = k; i < weights.length; i++) {
			if (weights[i] < maxWeight) {
				// replace the max so far with this value
				temp[maxIndex] = weights[i];
				tempIndexInWeights[maxIndex] = i;

				//find new max
				maxWeight = Double.MIN_VALUE;
				maxIndex = -1;

				for (int j = 0; j< k; j++) {
					if (temp[j] > maxWeight) {
						maxWeight = temp[j];
						maxIndex = j;
					}
				}
			}
		}
		lowestKIndices = new ArrayList<Integer>(Arrays.asList(tempIndexInWeights));
		return lowestKIndices;
	}
	
	private List<Integer> getKHighestWeight(Asteroid[] asteroids, int k) {
		double[] weights = new double[asteroids.length];
		double distance = 0.0;
		double mass = 0.0;

		for (int i = 0; i < asteroids.length; i++) {
			Asteroid a = asteroids[i];
			Point a_center = new Point();
			a.orbit.positionAt(time - a.epoch, a_center);
			distance = Double
					.valueOf(Point.distance(a_center, new Point(0, 0)));
			mass = asteroids[i].mass;
			weights[i] = distance / mass;
		}

		List<Integer> highestKIndices;
		double[] temp = new double[k];
		double minWeight = Integer.MAX_VALUE;
		int minIndex = -1;
		Integer[] tempIndexInWeights = new Integer[k];

		// put k weights in temp, and keep track of max seen so far
		for (int i = 0; i < k; i++) {
			temp[i] = weights[i];
			tempIndexInWeights[i] = i;
			if (temp[i] < minWeight) {
				minWeight = temp[i];
				minIndex = i;
			}
		}

		// check remaining weights and maintain lowest k in temp
		for (int i = k; i < weights.length; i++) {
			if (weights[i] > minWeight) {
				// replace the max so far with this value
				temp[minIndex] = weights[i];
				tempIndexInWeights[minIndex] = i;

				//find new max
				minWeight = Double.MAX_VALUE;
				minIndex = -1;

				for (int j = 0; j< k; j++) {
					if (temp[j] < minWeight) {
						minWeight = temp[j];
						minIndex = j;
					}
				}
			}
		}
		highestKIndices = new ArrayList<Integer>(Arrays.asList(tempIndexInWeights));
		return highestKIndices;
	}

	private int getLowestWeightAsteroid(Asteroid[] asteroids) {
		/** 
		* Returns the index of the asteroid with the minimum distance to mass ratio.
		* (Higher the distance mass ratio, more "pushable" the asteroid - costs less energy)
		* 
		* Hence, it is not good to push the one with the lowest distance to mass ratio.
		* because it is too close to the sun, and too heavy.
		**/

		int min = 0;
		double distance = 0.0;
		double mass = 0.0;
		double weight = 0.0;
		double minWeight = Double.MAX_VALUE;

		for (int i = 0; i < asteroids.length; i++) {

			Asteroid a = asteroids[i];
			Point a_center = new Point();
			a.orbit.positionAt(time - a.epoch, a_center);
			distance = Double
					.valueOf(Point.distance(a_center, new Point(0, 0)));
			mass = asteroids[i].mass;
			weight = distance / mass;
			if (weight < minWeight) {
				min = i;
				minWeight = weight;
			}
		}
		return min;
	}

	private int getMinDistanceTime(Asteroid a, Asteroid  b) {
		/**
		* Returns the time at which the Euclidean distance between asteroids is minimized.
		**/

		double minDistance = -1;
		int minDistTime = Integer.MAX_VALUE;
		for(int t = 0; t < 20*365; t++){
			double distance = Point.distance(a.orbit.positionAt((long)t - a.epoch),b.orbit.positionAt((long)t - b.epoch));
			if(distance < minDistance){
				minDistance = distance;
				minDistTime = t;
			}			
		}
		return minDistTime;
	} 
	
	private void attemptCollisions(Asteroid a, Asteroid b, int minDistTime){
		// Calculate which needs to be pushed
		Point posA = a.orbit.positionAt((long)minDistTime - a.epoch - 365);
		Point posB = b.orbit.positionAt((long)minDistTime - b.epoch);		
		// For now, just push "a"
		boolean didPush = false;
		// Create vector that draws a line from location of push to where I 
		// want the intersection to take place
		Point dirVector = new Point(posB.x - posA.x, posB.y - posA.y);
		double dir = dirVector.direction();
		Point v = a.orbit.velocityAt((long)minDistTime - a.epoch - 365);
		double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
		double r = a.radius() + b.radius();
		Point newPos;
		for(int i = 0; i < 500; i++){
			double v2 = v1 * ((1+i/100f));
			double E = 0.5 * a.mass * v2 * v2;
			Asteroid testAst = pushTest(a, minDistTime - 365, E, dir);
			newPos = testAst.orbit.positionAt((long)minDistTime - testAst.epoch);
			if (Point.distance(posB, newPos) < r){
				System.out.println("Good push!");
				Asteroid.push(a, minDistTime - 365, E, dir);
				didPush = true;
				break;
			} 
			
		}
		if(didPush == false){
			System.out.println("Could not find good push!");			
		}
	}
	
	private boolean CheckIncidentalCollisions(onePush push, Asteroid[] asteroids){
		Asteroid pushAst = null;
		for(int i = 0; i < asteroids.length; i++){
			if(asteroids[i].id == push.asteroidPushId){
				pushAst = asteroids[i];
				break;
			}
		}
		if(pushAst == null){
			System.out.println("ERRROR");
		}
		Asteroid newAst = pushTest(pushAst, push.time, push.energy,
				push.direction);
		
		Point newAstPos;
		Point otherAstPos;
		double r;
		for(long t = push.time; t< push.collision_time; t++){
			newAstPos = newAst.orbit.positionAt(t);
			for(int j = 0; j<asteroids.length; j++){
				if(asteroids[j].id != push.asteroidPushId && asteroids[j].id != push.asteroidCollidedId){					
					otherAstPos = asteroids[j].orbit.positionAt(t);
					r = newAst.radius() + asteroids[j].radius();
					if(Point.distance(newAstPos, otherAstPos) < r){
						System.out.println("Bad Push: Incidental Collision Detected with Asteroid "+j);
						return false;
					}
					
				}
				
			}
		}
		return true;
		
	}
	public static Asteroid pushTest(Asteroid asteroid, long time,
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
		return new Asteroid(new Orbit(r, v), asteroid.mass, time);
	}

	
}

class onePush {

	// the time the push is made
	long time;

	// energy used to make the push
	double energy;

	// the direction of the push
	double direction;

	// the time at which the push will cause a collision
	long collision_time;
	
	// index of the asteroid to push in the array asteroids[]
	int i;
	
	// id of the asteroid to push
	long asteroidPushId;
	
	// id of the asteroid to push towards
    long asteroidCollidedId;
	
	public onePush(long time, double energy, double direction, int i, long pushId,
			long collision_time, long collidedId) {
		this.time = time;
		this.energy = energy;
		this.direction = direction;
		this.i = i;
		this.asteroidPushId = pushId;
		this.collision_time = collision_time;
		this.asteroidCollidedId = collidedId;
		
	}
}
