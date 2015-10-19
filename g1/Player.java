package pb.g1;


import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.ArrayList;
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
	
	
	private int prevAsteroidsLength = 0;
	private int initAsteroidsLength = 0;
	
	
	private boolean correctionInProgress;
	
	
	// Select the asteroid that is closest to the given mass
	private double asteroidCombinedMass;
	
	int playCount = 0;
	
	
	

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit)
	{
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
		
		this.prevAsteroidsLength = asteroids.length;
		this.initAsteroidsLength = asteroids.length;
		
	}
	
	
	private double norm (Point point) 
	{
		return Math.sqrt(point.x*point.x + point.y*point.y);
	}
	
	
	private void orbitCorrection(int i, Asteroid[] asteroids, double[] energy, double[] direction)
	{
		Point location = new Point();
		asteroids[i].orbit.positionAt(time - asteroids[i].epoch, location);
		Point velocity = new Point();
		asteroids[i].orbit.velocityAt(time - asteroids[i].epoch, velocity);
		double orbitRadius = norm(location);
	    double circularSpeed = Math.sqrt(Orbit.GM / orbitRadius);
	    double targetAngle = Math.PI/2 + Math.atan2(location.y, location.x);
		double targetX = circularSpeed * Math.cos(targetAngle);
		double targetY = circularSpeed * Math.sin(targetAngle);
		double pushX = targetX - velocity.x;
		double pushY = targetY - velocity.y;
		double pushAngle = Math.atan2(pushY, pushX);
		double pushSpeed = norm(new Point(pushX, pushY));
		System.out.println("Correction push ");
	    energy[i] = 0.5 * asteroids[i].mass * pushSpeed * pushSpeed;
	    direction[i] = pushAngle;
	}
	
	
	

	
	// try to push asteroid
	public void play(Asteroid[] asteroids,
	                 double[] energy, double[] direction)
	{
		
		for (int i=0; i<asteroids.length; i++)
		{
			System.out.println(asteroids[i].orbit.a);
		}
		
		//Try to solve the problem thinking that there are just 2 asteroids available ..
		time++;
		boolean allCorrect = true;
		
		
		int currCount = asteroids.length;
		if (currCount != prevAsteroidsLength)
		{
			
		}
		

		for (int i=0; i<asteroids.length; i++)
		{
			// What should be the push ?
			if (Math.abs(asteroids[i].orbit.a - asteroids[i].orbit.b) > 2)
			{
				System.out.println("Something is wrong !!! " + i);
				allCorrect = false;
			}
		}
		
		if (allCorrect)
		{
			//System.out.println("Yippieee");
		}
		
		
		int currAsteroidsLength = asteroids.length;
		
		if (currAsteroidsLength != prevAsteroidsLength)
		{
			
			//Correct the orbit of the asteroid 
			  for (int i=0; i<asteroids.length; i++)
			  {
				// What should be the push ?
				//if (asteroids[i].orbit.a != asteroids[i].orbit.b)
				if (Math.abs(asteroids[i].orbit.a - asteroids[i].orbit.b) > 100)
				{
					orbitCorrection(i, asteroids, energy, direction);
				    return; 
			
				}
		
		}
		}
		
		
		//prevAsteroidsLength = asteroids.length;
		
		if (allCorrect)
		{
			/*int innerIndex = 0;
			int outerIndex = 0;
			//double r1 = Math.min(asteroids[0].orbit.a, asteroids[1].orbit.a)
			if (asteroids[0].orbit.a < asteroids[1].orbit.a )
			{
				innerIndex = 0;
				outerIndex = 1;
			}
			else
			{
				innerIndex = 1;
				outerIndex = 0;
			}
			
			double r1 = asteroids[innerIndex].orbit.a;
			double r2 = asteroids[outerIndex].orbit.a;
			
			double omegaOuter = Math.sqrt(Orbit.GM/ Math.pow(r2, 3));
			Point velocityInner = asteroids[innerIndex].orbit.velocityAt(time);
			double angleInner = Math.atan2(velocityInner.y, velocityInner.x);
			
			Point velocityOuter = asteroids[outerIndex].orbit.velocityAt(time);
			double angleOuter = Math.atan2(velocityOuter.y, velocityOuter.x);
			
			double timeTransfer = Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM));
			
			double sumRadii = asteroids[innerIndex].radius() + asteroids[outerIndex].radius();
			
			double alignThreshold = sumRadii/r2;
			
			double alignment = Math.abs(Math.PI - timeTransfer*omegaOuter - (angleOuter - angleInner)); 
			
			if (alignment < alignThreshold)
			{
				double velocityNew = Math.sqrt(Orbit.GM / r1) * (Math.sqrt( 2*r2 / (r1+r2)) - 1);
				energy[innerIndex] = 0.5*asteroids[innerIndex].mass * velocityNew * velocityNew;
				direction[innerIndex] = angleInner;
			}
			
			*/
			
			
			//Try to solve the problem thinking that there are just 2 asteroids available ..
			//time++;
			
			
			int innerIndex = 0;
			int outerIndex = 0;
			if (asteroids[0].orbit.a < asteroids[1].orbit.a )
			{
				innerIndex = 0;
				outerIndex = 1;
			}
			else
			{
				innerIndex = 1;
				outerIndex = 0;
			}
			
			double r1 = asteroids[innerIndex].orbit.a;
			double r2 = asteroids[outerIndex].orbit.a;
			
			
			double omegaOuter = Math.sqrt(Orbit.GM/ Math.pow(r1, 3));
			Point velocityInner = asteroids[innerIndex].orbit.velocityAt(time - asteroids[innerIndex].epoch);
			double angleInner = Math.atan2(velocityInner.y, velocityInner.x);
			
			Point velocityOuter = asteroids[outerIndex].orbit.velocityAt(time - asteroids[outerIndex].epoch);
			double angleOuter = Math.atan2(velocityOuter.y, velocityOuter.x);
			
			double timeTransfer = Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM));
			
			double sumRadii = asteroids[innerIndex].radius() + asteroids[outerIndex].radius();
			
			double alignThreshold = sumRadii/r1;
			
			double alignment = Math.abs(Math.PI - timeTransfer*omegaOuter - (angleInner - angleOuter)); 
			
			if (alignment < alignThreshold)
			{
				double anglePush = 0.0;
				if (angleOuter > 0)
				{
					anglePush = angleOuter - Math.PI;
				}
				else
				{
					anglePush = Math.PI + angleOuter;
				}
				//double velocityNew = Math.sqrt(Orbit.GM / r1) * (Math.sqrt( 2*r2 / (r1+r2)) - 1);
				double velocityNew = Math.sqrt(Orbit.GM / r2) * (1 - Math.sqrt( 2*r1 / (r1+r2)) );
				
				energy[outerIndex] = 0.5*asteroids[outerIndex].mass * velocityNew * velocityNew;
				direction[outerIndex] = anglePush;
				
				
				
				//store mass
				//asteroidCombinedMass = (asteroids[outerIndex].mass + asteroids[innerIndex].mass);
			} 
			
			
			//double r1 = Math.min(asteroids[0].orbit.a, asteroids[1].orbit.a)
			
			
			
			
		}
		
		prevAsteroidsLength = asteroids.length;
		//playCount++;
		

			
	} 
	
	
	

	// try to push asteroid
	public void play2(Asteroid[] asteroids,
	                 double[] energy, double[] direction)
	{
		
		time++;
		
		if ( time !=0 && time%100 == 0 )
		{
			// Get the nearest asteroid
			int nearestAsteroidIndex = getNearestAsteroid(asteroids);
			// Push the nearest asteroid
			
			Point v = asteroids[nearestAsteroidIndex].orbit.velocityAt(time);
			// add 5-50% of current velocity in magnitude
			double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
			
			double v2 = 0.05*v1;
			
			// apply push at -π/8 to π/8 of current angle
			double angle1 = Math.atan2(v.y, v.x);
			
			double angle2 = angle1 - Math.PI/2;
			
			if (angle2 < -Math.PI)
			{
				angle2 = 2*Math.PI + angle2;
			}
			
			double E = 0.5 * asteroids[nearestAsteroidIndex].mass * v2 * v2;
			
			energy[nearestAsteroidIndex] = E;
			direction[nearestAsteroidIndex] = angle2;
			
			
			
			
			
			
		} 
		
		
		/*for (int i=0; i<asteroids.length; i++)
		{
			Point velocity = asteroids[i].orbit.velocityAt(time);
			double v1 = Math.sqrt(velocity.x * velocity.x + velocity.y * velocity.y);
			double timePeriod = ((2*Math.PI*asteroids[i].orbit.a)/v1)/(24*60*60); 
			System.out.println("Time period of asteroid " + i + " : " + timePeriod);
			double d1 = Math.atan2(velocity.y, velocity.x);
			
			double degrees = Math.toDegrees(d1);
			//System.out.println("Degrees of first asteroid " + degrees);
			
			System.out.println("Degrees of asteroid " + i + " : " + degrees);
			
		}
		
		
		AsteroidIndex[] asteroidsWrapper = new AsteroidIndex[asteroids.length];
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
		
		
		// if not yet time to push do nothing
		if (++time <= time_of_push) return;
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
			System.out.println("  Speed: " + v1 + " +/- " + v2);
			// apply push at -π/8 to π/8 of current angle
			double d1 = Math.atan2(v.y, v.x);
			
			//double degrees = Math.toDegrees(d1);
			
			//System.out.println("Degrees of first asteroid " + degrees);
			
			
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
				double r = a1.radius() + a2.radius();
				// look 10 years in the future for collision
				for (long ft = 0 ; ft != 7320 ; ++ft) {
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
						System.out.println("  Day: "  + (1 + t % 365));
						return;
					}
				}
			}
			System.out.println("  No collision ...");
		}
		time_of_push = time + turns_per_retry; */
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
	
	
	
	public static int getNearestAsteroid(Asteroid[] asteroids){
		System.out.println("Method Called");
		double min_major_radius = Double.MAX_VALUE;
		int min_major_index = 0;
		for(int i=0;i<asteroids.length;i++){
			if (asteroids[i].orbit.a < min_major_radius){
				System.out.println("Asteroid "+i+" asteroids[i].orbit.a ");
				min_major_radius = asteroids[i].orbit.a ;
				min_major_index = i;
			}
		}
		return min_major_index;
	}
	
	
	
	class AsteroidComparator implements Comparator<AsteroidIndex>
	{
	    public int compare(AsteroidIndex h1, AsteroidIndex h2)
	    {
	    	int cmp = 0;
	    	
	    	/*if (h1.mass < h2.mass)
	    	{
	    		return 1;
	    	}
	    	if (h1.mass >= h2.mass)
	    	{
	    		return -1;
	    	}*/
	    	if (h1.radius < h2.radius)
	    	{
	    		return -1;
	    	}
	    	if (h1.radius >= h2.radius)
	    	{
	    		return 1;
	    	}
	    	return cmp;
	    }
	}
	
	
	
	
}
