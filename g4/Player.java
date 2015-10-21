package pb.g4;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.*;

public class Player implements pb.sim.Player{
  public static final int FUTURE_DAYS_TO_SEARCH = 5;

  private Point sun = new Point(0, 0);
  private long time = -1;
  private long time_limit = -1;
  private long collisionTime = -1;

  private Asteroid furthestFromSun;
  private Asteroid closerToSun;
  private Asteroid largestAsteroid;
  private int indexToPush = -1;
  private int indexToHit = -1;
  private Set<Asteroid> asteroidOrder;
  private double fifty_percent_mass;

  private int num_otherAsteroidLocation_asteroids = 4;
  private int numAsteroids;
  private boolean firstCollision = true;
  private boolean[] colliding;

  private double bestEnergy;
  private int bestAsteroidIndex;
  private double bestAngle;
  private long timeOfNextPush;

  //stores asteroid masses
  private HashMap<Asteroid, Double> cached_asteroid_masses = new HashMap<Asteroid, Double>();

  // print orbital information
  public void init(Asteroid[] asteroids, long time_limit) {
    firstCollision = false;
    colliding = new boolean[asteroids.length];
    numAsteroids = asteroids.length;
    asteroidOrder = new HashSet<Asteroid>();
    System.out.println("Init");
    if (Orbit.dt() != 24 * 60 * 60)
      throw new IllegalStateException("Time quantum is not a day");
    this.time_limit = time_limit;
  }

  // try to push asteroid
  public void play(Asteroid[] asteroids,
      double[] energy, double[] direction) {
    if (time % 365 == 0) {
      System.out.println("Year: " + time / 365);
    }
    correctCollidedOrbits(asteroids, energy, direction);
    if (asteroids.length != numAsteroids) {
      colliding = new boolean[asteroids.length];
      firstCollision = false;
      numAsteroids = asteroids.length;
    }
    else if (++time%10 == 0 && time > collisionTime) {
      pushAllToLargest(asteroids, energy, direction);
    }
    if (timeOfNextPush == time && bestEnergy < Double.MAX_VALUE) {
      System.out.println("Push at time: " + time + " bestEnergy " + bestEnergy + " bestAngle: " + bestAngle);
      colliding[bestAsteroidIndex] = true;
      energy[bestAsteroidIndex] = bestEnergy;
      direction[bestAsteroidIndex] = bestAngle;
    }
  }

  public void correctCollidedOrbits(Asteroid[] asteroids, double[] energy, double[] direction) {

    for (int i = 0; i < asteroids.length; i++ ) {
      if (Math.abs(asteroids[i].orbit.b - asteroids[i].orbit.a) > asteroids[i].radius() && !colliding[i]) {
        Point position = asteroids[i].orbit.positionAt(time - largestAsteroid.epoch);
        //Velocity for a hypothetical circular orbit at this position
        Point circularVelocity = new Orbit(position).velocityAt(0);
        Point currentVelocity = asteroids[i].orbit.velocityAt(time - largestAsteroid.epoch);
        Point dv = new Point(circularVelocity.x - currentVelocity.x, circularVelocity.y - currentVelocity.y);

        double pushEnergy = asteroids[i].mass * Math.pow(dv.magnitude(), 2) / 2;
        double pushAngle = dv.direction();

        energy[i] = pushEnergy;
        direction[i] = pushAngle;
      }
    }
  }

  public int findLargestAsteroidIndex(Asteroid[] asteroids) {
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

  public void pushAllToLargest(Asteroid[] asteroids, double[] energy, double[] direction) {
    int largestAsteroid_idx;
    if (!firstCollision) {
      largestAsteroid_idx = findLargestAsteroidIndex(asteroids);
    }
    else {
      largestAsteroid_idx = getBestAsteroid(asteroids, 0.5);
    }

    largestAsteroid = asteroids[largestAsteroid_idx];
    double r2 = largestAsteroid.orbit.a;

    //Reset best energy
    bestEnergy = Double.MAX_VALUE;

    for (int i = 0; i < asteroids.length; i++) {
      if (i != largestAsteroid_idx) {
        Asteroid otherAsteroid = asteroids[i];
        double mass = otherAsteroid.mass;
        double pushAngle = otherAsteroid.orbit.velocityAt(time - otherAsteroid.epoch).direction();

        double r1 = otherAsteroid.orbit.a;
        double dv = Math.sqrt(Orbit.GM / r1)
          * (Math.sqrt((2 * r2)/(r1 + r2)) - 1);

        if (dv < 0) pushAngle += Math.PI;

        double pushEnergy = mass * dv * dv * 0.5;

        for (int days = 0; days < FUTURE_DAYS_TO_SEARCH; days++) {
          long predictedTimeOfCollision = prediction(otherAsteroid, largestAsteroid, time + days, pushEnergy, pushAngle);

          //If collision predicted and it's better than previous best collision
          if (predictedTimeOfCollision > 0 && pushEnergy < bestEnergy) {
            System.out.println("collision predicted" + " at energy: "
                + pushEnergy + " and direction: " + pushAngle + " at year: " + time + days / 365);
            bestEnergy = pushEnergy;
            bestAngle = pushAngle;
            bestAsteroidIndex = i;
            timeOfNextPush = time + days;
            collisionTime = predictedTimeOfCollision;
          }
        }
      }
    }

    System.out.println("Time diff: " + (timeOfNextPush - time));
    // if (bestEnergy < Double.MAX_VALUE) {	
    // }
  }

  public long prediction(Asteroid source, Asteroid target, long time, double energy, 
      double direction) {
    try {
      source = Asteroid.push(source, time, energy, direction);
    } catch (InvalidOrbitException e) {
      System.out.println("Invalid orbit predicted with energy " + energy + " and angle " + direction);
    }

    Point p1 = source.orbit.velocityAt(time - source.epoch);
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
        System.out.println("Collision predicted at year " + t / 365 + " day: " + t % 365 + " with energy: " + energy);
        return t;
      }
    }
    return -1;
  }

  public class AsteroidComparator implements Comparator<Asteroid> {
    public int compare(Asteroid a1, Asteroid a2) {
      double a1Distance = Point.distance(a1.orbit.positionAt(time - a1.epoch), largestAsteroid.orbit.positionAt(time - largestAsteroid.epoch));
      double a2Distance = Point.distance(a2.orbit.positionAt(time - a2.epoch), largestAsteroid.orbit.positionAt(time - largestAsteroid.epoch));

      if ( a1Distance < a2Distance ) {
        return -1;
      } else if ( a2Distance < a1Distance ) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  public class OrbitComparator implements Comparator<Asteroid>
  {
    public int compare(Asteroid a1, Asteroid a2)
    {
      return (int)(getDistanceFromSun(a1) - getDistanceFromSun(a2));
    }


    public double getDistanceFromSun(Asteroid a)
    {
      return Point.distance(a.orbit.positionAt(time), sun);
    }
  }

  public int getBestAsteroid(Asteroid[] asteroids, double percent)
  {
    ArrayList<Asteroid> center_asteroids = new ArrayList<Asteroid>();
    ArrayList<Asteroid> sorted_asteroids = new ArrayList<Asteroid>(Arrays.asList(asteroids));
    Collections.sort(sorted_asteroids, new OrbitComparator());

    int num_asteroids = sorted_asteroids.size();
    int select_asteroids = (int)(percent*num_asteroids);
    double max_mass = 0;
    Asteroid max_asteroid;
    for (int i = (int)(0.5*num_asteroids - 0.5*select_asteroids); i < (int)(0.5*num_asteroids + 0.5*select_asteroids); i++)
    {
      if(sorted_asteroids.get(i).mass > max_mass)
      {
        max_mass = asteroids[i].mass;
        max_asteroid = asteroids[i];
      }

    }
    for(int i = 0; i < asteroids.length; i++)
    {
      if(max_mass == asteroids[i].mass)
      {
        return i;
      }
    }
    return -1;
  }

}
