package pb.g1;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Arrays;
import java.util.Comparator;
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
	private int turns_per_retry = 3;
	
	private double velocity_factor = 0.00;
	
	private boolean directionFlag = true;
	
	

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
		
		// if not yet time to push do nothing
		if (++time <= time_of_push) return;
		
		
		if(directionFlag){
		velocity_factor = velocity_factor +0.1;
		}
		
		if(!directionFlag){
			velocity_factor = velocity_factor - 0.1;
		}
		
		if(velocity_factor > 3.5){
			velocity_factor = 3.5;
			directionFlag = !directionFlag;
		}
		else if(velocity_factor < 0.1){
			velocity_factor = 0.1;
			directionFlag = !directionFlag;
		}
		
		/*velocity_factor = velocity_factor +0.1;
		if(velocity_factor > 2.5){
			velocity_factor = 2.5;
			
		}*/
		
		

		AsteroidIndex[] asteroidsWrapper = new AsteroidIndex[asteroids.length];
		double min_energy = Double.MAX_VALUE;
		double min_direction = 0;
		long min_energy_time=0;
		
		AsteroidIndex asteroidIndex = null;
		for (int i=0; i<asteroids.length; i++)
		{
			asteroidIndex = new AsteroidIndex(i, asteroids[i].mass);
			asteroidsWrapper[i] = asteroidIndex;
		}
		
		Arrays.sort(asteroidsWrapper, new AsteroidComparator());
		
		for (int i=0; i<asteroidsWrapper.length; i++)
		{
			System.out.println(asteroidsWrapper[i].index + " : mass : " + asteroidsWrapper[i].mass );			
		}
		
		System.out.println("---------------------------------------");
		
		System.out.println("Year: " + (1 + time / 365));
		System.out.println("Day: "  + (1 + time % 365));
		for (int retry = 1 ; retry <= retries_per_turn ; ++retry) {
			// pick a random asteroid and get its velocity
			//int i = random.nextInt(asteroids.length);
			
			int i= getLightestAsteroid(asteroids);
			
			System.out.println("Index " + i);
			
			
			Point v = asteroids[i].orbit.velocityAt(time);
			// add 5-50% of current velocity in magnitude
			System.out.println("Try: " + retry + " / " + retries_per_turn);
			double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
			double v2 = v1 * (random.nextDouble() * 0.45 + 0.05);
			//double v2 = v1 * velocity_factor;
			System.out.println("  Speed: " + v1 + " +/- " + v2);
			// apply push at -π/8 to π/8 of current angle
			double d1 = Math.atan2(v.y, v.x);
			double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
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
			for (int j = 0 ; j != asteroids.length ; ++j) {
				if (i == j) continue;
				Asteroid a2 = asteroids[j];
				//Asteroid a2 = asteroids[getHeaviestAsteroid(asteroids)];
				double r = a1.radius() + a2.radius();
				// look 20 years in the future for collision
				for (long ft = 0 ; ft != 3650 ; ++ft) {
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
						/*if(E < min_energy){
							min_energy = E;
							min_direction = d2;
							min_energy_time = t + 1;
						}*/
						System.out.println("  Collision prediction !");
						System.out.println("  Year: " + (1 + t / 365));
						System.out.println("  Day: "  + (1 + t % 365));
						return;
					}
				}
				
			}
			
		/*	if(min_energy != Double.MAX_VALUE){
			energy[i] = min_energy;
			direction[i] = min_direction;
			// do not push again until collision happens
			time_of_push = min_energy_time;
			return;
			}*/
			
			System.out.println("  No collision ...");
		}
		
		
		time_of_push = time + turns_per_retry;
		
	}
	
	public static int getLightestAsteroid(Asteroid[] asteroids){
		double min_weight = Double.MAX_VALUE;
		int min_weight_index = 0;
		for(int i=0;i<asteroids.length;i++){
			if (asteroids[i].mass< min_weight){
				min_weight = asteroids[i].mass;
				min_weight_index = i;
			}
		}
		
		return min_weight_index;
		
	}
	
	public static int getHeaviestAsteroid(Asteroid[] asteroids){
		double max_weight = 0;
		int max_weight_index = 0;
		for(int i=0;i<asteroids.length;i++){
			if (asteroids[i].mass> max_weight){
				max_weight = asteroids[i].mass;
				max_weight_index = i;
			}
		}
		
		return max_weight_index;
		
	}
	
	
	public static int getFarthestAsteroid(Asteroid[] asteroids){
		System.out.println("Method Called");
		double max_major_radius = 0;
		int max_major_radius_index = 0;
		for(int i=0;i<asteroids.length;i++){
			if (asteroids[i].orbit.a > max_major_radius){
				System.out.println("Asteroid "+i+" asteroids[i].orbit.a ");
				max_major_radius = asteroids[i].orbit.a ;
				max_major_radius_index = i;
			}
		}
		
		return max_major_radius_index;
		
	}
	
	public boolean doesElipseIntersect(Asteroid a1, Asteroid a2){
		
		while(a1.orbit.iterator().hasNext()){
			
			Point p = a1.orbit.iterator().next();
			
			
			
			double foci = a2.orbit.a * 2;
			
			
		}
			
		
		return true;
	}
	
	public double getEccentricity(Asteroid a){
		
		Point v = a.orbit.velocityAt(time);
		
		double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
		
		Point centre = a.orbit.center();
		
		
		
		return 0.0;
	}
	
	class AsteroidComparator implements Comparator<AsteroidIndex>
	{
	    public int compare(AsteroidIndex h1, AsteroidIndex h2)
	    {
	    	int cmp = 0;
	    	
	    	if (h1.mass < h2.mass)
	    	{
	    		return 1;
	    	}
	    	if (h1.mass >= h2.mass)
	    	{
	    		return -1;
	    	}
	    	return cmp;
	    }
	}
}
