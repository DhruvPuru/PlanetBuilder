package pb.g7;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Random;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Comparator;
import java.util.Iterator;

public class Player implements pb.sim.Player {

    // used to pick asteroid and velocity boost randomly
    private Random random = new Random();

    // current time, time limit
    private long time = -1;
    private long time_limit = -1;

    // time until next push
    //private long time_of_collision = 0;
    //private boolean collisionStarted = false;

    private int prevNumAsteroids = -1;

    private static final double G = Orbit.GM / Orbit.M;//6.67382967579392e-11;
    private static final double M = Orbit.M;//1.98855e30;

    private long accumulatorID;

    private HashMap<Long,Long> collidingAsteroids;

    HashSet<push_move> asteroids_in_phase;

    private int init_num_asteroids;
    // number of retries



    HashMap<Long, push_move> push_queue = new HashMap<Long, push_move>();

    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit)
    {
	if (Orbit.dt() != 24 * 60 * 60)
	    throw new IllegalStateException("Time quantum is not a day");
	this.time_limit = time_limit;
	prevNumAsteroids = asteroids.length;
	init_num_asteroids = asteroids.length;
	asteroids_in_phase = getAsteroidsInPhase(asteroids, this.time_limit - 100*365, asteroids.length-1);
	this.collidingAsteroids = calc_asteroids(asteroids, asteroids_in_phase, asteroids.length-1);

    }

    private HashSet<push_move> getAsteroidsInPhase(Asteroid[] asteroids, long time_limit, int accumulator_index) {
	int n = asteroids.length;
	ArrayList<asteroid_index> asteroids_ordered = new ArrayList <asteroid_index> ();
	Point astr_location = new Point();
	for (int i=0; i<n; i++) {
	    asteroids[i].orbit.positionAt(time - asteroids[i].epoch, astr_location);
	    asteroids_ordered.add(new asteroid_index(i, l2norm(astr_location)));
	}
	Collections.sort(asteroids_ordered);	
	Asteroid accumulator = asteroids[asteroids_ordered.get(accumulator_index).index];
	this.accumulatorID = accumulator.id;
	Point v1 = new Point();
	Point v2 = new Point();
	double r2 = 0.5*(accumulator.orbit.a + accumulator.orbit.b);
	double omega2 = Math.sqrt(Orbit.GM / Math.pow(r2,3));

	HashSet<push_move> asteroids_in_phase = new HashSet<push_move>();
	for (int i=0; i<n; i++) {
	    if (i == accumulator_index) {continue;}
	    Asteroid push_asteroid = asteroids[asteroids_ordered.get(i).index];
	    double r1 = 0.5*(push_asteroid.orbit.a+push_asteroid.orbit.b);
	    double thresh = push_asteroid.radius() + accumulator.radius();
	    double tH = Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM));

	    for (long t = 0 ; t < .8*time_limit ; ++t) {
		push_asteroid.orbit.velocityAt(t, v1);
		double theta1 = Math.atan2(v1.y,v1.x);
		accumulator.orbit.velocityAt(t, v2);
		double theta2 = Math.atan2(v2.y,v2.x);		    
		if ( Math.abs(theta1 + Math.PI - theta2 - tH*omega2) < thresh / r2) {
		    double deltav = Math.sqrt(Orbit.GM / r1) * (Math.sqrt( 2*r2 / (r1+r2)) - 1);
		    double E = 0.5*asteroids[asteroids_ordered.get(i).index].mass * (deltav * deltav + 2 * l2norm(v1) * deltav);
		    //double E = 0.5*asteroids[i].mass*(deltav*deltav);
		    asteroids_in_phase.add(new push_move(i,push_asteroid.id, E, theta1, time+t, E / push_asteroid.mass));
		    break;
		}
	    }
	}	      
	return asteroids_in_phase;
    }

    	public double get_total_mass(Asteroid[] asteroids)
    {
    	double total_mass = 0.0;
    	for(Asteroid each_asteroid : asteroids)
    		total_mass += each_asteroid.mass;
    	return total_mass;
    }

    public HashMap<Long, Long> calc_asteroids(Asteroid[] asteroids,
					      HashSet<push_move> colliding_asteroids, int accumulator)
    {
    	double total_mass = get_total_mass(asteroids);
    	HashMap<Long,Long> toPush = new HashMap<Long,Long>();
    	double accumulated_mass = asteroids[accumulator].mass;

    	Comparator<push_move> density_comparator = new Comparator<push_move>(){
		public int compare(push_move pm1, push_move pm2)
		{
			double d1 = pm1.density;
			double d2 = pm2.density;
			return Double.compare(d1, d2);
		} };
    	

	Comparator<push_move> energy_comparator = new Comparator<push_move>(){
		public int compare(push_move pm1, push_move pm2)
		{
			double e1 = pm1.energy;
			double e2 = pm2.energy;
			return Double.compare(e1, e2);

		    //if(Math.abs(e1-e2)<.00000000000001)
		    //	return 0;
		    //else if(e1 > e2)
		    //	return 1;
		    //else
		    //	return -1;
		} };


	PriorityQueue<push_move> min_density_q = new PriorityQueue<push_move>(colliding_asteroids.size(), density_comparator);
	PriorityQueue<push_move> min_energy_q = new PriorityQueue<push_move>(colliding_asteroids.size(), energy_comparator);

	for(push_move pm : colliding_asteroids)
	    {
		min_density_q.add(pm);
		min_energy_q.add(pm);
	    }

	    /*PriorityQueue<push_move> min_E_q_copy = min_energy_q.clone();
	    while(!min_E_q_copy.isEmpty())
	    {
	    	System.out.println(min_E_q_copy.poll().energy);
	    }*/

	while(accumulated_mass < total_mass*.55)
	    {
		if(min_energy_q.isEmpty() || min_density_q.isEmpty())
		    {
			System.out.println("Asteroids don't have enough mass");
			break;
		    }
		push_move minE = min_energy_q.peek();
		push_move minD = min_density_q.peek();
		double massE = asteroids[minE.index].mass; //density = Energy/Mass 
		double massD = asteroids[minD.index].mass;
		if((accumulated_mass + massE)>= total_mass/2.0)
		    {

			accumulated_mass+=massE;
			toPush.put(minE.id, minE.time);
			min_energy_q.remove(minE);
			min_density_q.remove(minE);
			//System.out.println("Density= " + minE.density);
			//System.out.println("Energy = " + minE.energy);

		    }
		else
		    {
			accumulated_mass+=massD;
			toPush.put(minD.id, minD.time);
			min_energy_q.remove(minD);
			min_density_q.remove(minD);
			//System.out.println("Density= " + minD.density);
			//System.out.println("Energy = " + minD.energy);
		    }

	    }
	//System.out.println("accumulated mass = " + accumulated_mass);
	//	System.out.println("total mass = " + total_mass);
	//	System.out.println("number of asteroids to push = " + toPush.size());
	return toPush;

    }
    
    private double l2norm(Point p) {return Math.sqrt(p.x*p.x+p.y*p.y);}
    private double l2norm(double x, double y) {return Math.sqrt(x*x+y*y);}

    private ArrayList<Double> calculatePush(double speed, double angle, double targetSpeed, double targetAngle) {
	double vx = speed * Math.cos(angle);
	double vy = speed * Math.sin(angle);
	double target_vx = targetSpeed * Math.cos(targetAngle);
	double target_vy = targetSpeed * Math.sin(targetAngle);
	double push_vx = target_vx - vx;
	double push_vy = target_vy - vy;
	double push_angle = Math.atan2(push_vy, push_vx);
	double push_speed = l2norm(push_vx, push_vy);

	ArrayList<Double> parameters = new ArrayList<Double>();
	parameters.add(new Double(push_speed));
	parameters.add(new Double(push_angle));
	return parameters;
    }

    public boolean pathFree(Asteroid[] asteroids, int push_index, int target_index, double push_energy, double push_dir, long push_time) {
	Asteroid newAsteroid = null;
	try {
	    newAsteroid = Asteroid.push(asteroids[push_index], push_time, push_energy, push_dir);
	}
	catch (pb.sim.InvalidOrbitException e){
	    System.out.println(" Invalid orbit: " + e.getMessage());
	    return false;
	}
	long t = push_time;
	int n = asteroids.length;
	Point p1 = new Point();
	Point p2 = new Point();
	while (t < time_limit) {
	    newAsteroid.orbit.positionAt(t - newAsteroid.epoch, p1);
	    for (int i=0; i<n; i++) {
		if (i == push_index) {continue;}
		asteroids[i].orbit.positionAt(t - asteroids[i].epoch, p2);
		if (l2norm(p1.x-p2.x,p1.y-p2.y) < asteroids[i].radius() + newAsteroid.radius()) {
		    if (i == target_index) {
			return true;
		    }
		    else {
			return false;
		    }
		}
	    }
	    t++;
	}
	System.out.println("No collision...");
	return false;
    } 

    // try to push asteroid
    public void play(Asteroid[] asteroids,
		     double[] energy, double[] direction)
    {
	time++;
	int n = asteroids.length;
	ArrayList<asteroid_index> asteroids_ordered = new ArrayList <asteroid_index> ();
	Point astr_location = new Point();
	for (int i=0; i<n; i++) {
	    asteroids[i].orbit.positionAt(time - asteroids[i].epoch, astr_location);
	    asteroids_ordered.add(new asteroid_index(i, l2norm(astr_location)));
	}
	Collections.sort(asteroids_ordered);

	if (n < prevNumAsteroids) { //collision occurred, correct orbit of outermost asteroid

	    boolean collided = true;
	    for (int i=0; i<asteroids.length; i++) {
		if (asteroids[i].id == accumulatorID) {
		    collided = false;
		    break;
		}
	    }

	    if (collided) {
		accumulatorID = -1;
		int accumulatorIndex = -1;
		for (int i=0; i<asteroids.length; i++) {
		    if (asteroids[i].id > accumulatorID) {
			accumulatorID = asteroids[i].id;
			accumulatorIndex = i;
		    }
		}
		
		Point location = new Point();
		Asteroid accumulator = asteroids[accumulatorIndex];
		accumulator.orbit.positionAt(time - accumulator.epoch, location);
		double r2 = l2norm(location); 
		double orbit_speed = Math.sqrt(Orbit.GM / r2);
		double tangent_theta = Math.PI/2 + Math.atan2(location.y, location.x);
		double omega2 = Math.sqrt(Orbit.GM / Math.pow(r2,3));
		Point v2 = accumulator.orbit.velocityAt(time - accumulator.epoch);
		double normv2 = l2norm(v2);
		double theta2 = Math.atan2(v2.y,v2.x);
		ArrayList<Double> parameters = calculatePush(normv2, theta2, orbit_speed, tangent_theta);
		energy[accumulatorIndex] = 0.5 * accumulator.mass * Math.pow(parameters.get(0),2);
		direction[accumulatorIndex] = parameters.get(1);
		prevNumAsteroids = n;
		//System.out.println(accumulatorID);
		//collisionStarted = false;
		return;
	    }
	}

	int outerIndex = -1;
	for (int i=0; i<asteroids.length; i++) {
	    if (asteroids[i].id == accumulatorID) {
		outerIndex = i;
	    }
	}

	if (outerIndex == -1) {
	    System.out.println("No valid accumulator asteroid");
	    System.exit(0);
	}
	
	Asteroid outerAsteroid = asteroids[outerIndex];
	double r2 = asteroids_ordered.get(n-1).getRadius();// - asteroids[asteroids_ordered.get(n-1).index].radius();
	double omega2 = Math.sqrt(Orbit.GM / Math.pow(r2,3));
	Point v2 = outerAsteroid.orbit.velocityAt(time - outerAsteroid.epoch);
	double theta2 = Math.atan2(v2.y,v2.x);
	double normv2 = l2norm(v2);	    
	
	int min_index = -1;
	long min_time = 0;
	for (int i=n-1; i>=0; i--) {
	    int innerIndex = asteroids_ordered.get(i).index;
	    if (asteroids[innerIndex].id == accumulatorID) {continue;}
	    Asteroid innerAsteroid = asteroids[innerIndex];
	    if (innerAsteroid.orbit.a != innerAsteroid.orbit.b) {continue;} //if its not a circle, continue
	    double r1 = asteroids_ordered.get(i).getRadius();

	    long t = time;
	    Point v1 = innerAsteroid.orbit.velocityAt(t - innerAsteroid.epoch);
	    double normv1 = l2norm(v1);
	    double theta1 = Math.atan2(v1.y,v1.x);
	    double tH = Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM));
	    double thresh = innerAsteroid.radius() + outerAsteroid.radius();
	    v2 = outerAsteroid.orbit.velocityAt(t - outerAsteroid.epoch);
	    theta2 = Math.atan2(v2.y,v2.x);

	    if(collidingAsteroids.containsKey(innerAsteroid.id))
	    {
	    	if(time == collidingAsteroids.get(innerAsteroid.id))
	    	{
		    //System.out.println("Expecting to see asteroid " + innerAsteroid.id + " pushed now.");

	    	}
	    }
	    if ( Math.abs(theta1 + Math.PI - theta2 - tH*omega2) < thresh / r2) {
		if (collidingAsteroids.containsKey(innerAsteroid.id) || time_limit - time < 50*365){// && Math.abs(collidingAsteroids.get(innerAsteroid.id) - time) < 2*365 ) {
		    collidingAsteroids.remove(innerAsteroid.id);
		    double deltav = Math.sqrt(Orbit.GM / r1) * (Math.sqrt( 2*r2 / (r1+r2)) - 1);
		    double E = 0.5*asteroids[innerIndex].mass * deltav * deltav;
		    energy[innerIndex] = E;
		    direction[innerIndex] = theta1;

		    if (time_limit - time < 50*365) {System.out.println("Running out of time...");}
		    //System.out.println("Asteroid " + innerAsteroid.id + " was just pushed.");
		    //collisionStarted = true;
		    //System.out.println(collidingAsteroids.size() + " remaining asteroids to collide.");
		    /*if(collidingAsteroids.size() == 1)
		    {
		    	for(Long id : collidingAsteroids.keySet())
		    	{
			    System.out.println("time of predicted collision = DAY " + (1+collidingAsteroids.get(id)%365) + "YEAR " + (1+collidingAsteroids.get(id)/365));
		    	}
			}*/
		}
	    }
	}
    }
}
