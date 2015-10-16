package pb.g8;


import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;
import pb.sim.Orbit;
import pb.sim.Point;

public class Player implements pb.sim.Player {

    //iteration number
    int iteration = 1;

    // used to pick asteroid and velocity boost randomly
    private Random random = new Random();

    // current time, time limit
    private long time = -1;
    private long time_limit = -1;

    private Point origin = new Point(0,0);

    // time until next push
    private long time_of_push = 0;

    // number of retries
    private int retries_per_turn = 1;
    private int turns_per_retry = 3;

    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit) 
    {
        if (Orbit.dt() != 24 * 60 * 60)
            throw new IllegalStateException("Time quantum is not a day");
        this.time_limit = time_limit;
    }
    // try to push asteroid
    public void play(Asteroid[] asteroids, double[] energy, double[] direction) 
    {
        // if not yet time to push do nothing
        if (++time <= time_of_push) return;
        long time_left_per_asteroid = (time_limit - time)/asteroids.length;
        time_left_per_asteroid = Math.max(time_left_per_asteroid, 3650);
        //System.out.println("Year: " + (1 + time / 365));
        //System.out.println("Day: "  + (1 + time % 365));
        List<Push> pushes = new ArrayList<Push>();
        for (int retry = 1 ; retry <= retries_per_turn ; ++retry) 
        {
        	if(iteration == 1) 
            {
        		int i = getFarthestAsteroid(asteroids);
        		Point v = asteroids[i].orbit.velocityAt(time);
        		double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
                double v2 = 0.0;
                for (double k=0; k< 1; k=k+0.005) 
                {
                	v2 = v1 * (k * 0.45 + 0.05); //Can't change this
                    //System.out.println("  Speed: " + v1 + " +/- " + v2);
                    // apply push at -pi/8 to pi/8 of current angle
                    double d1 = Math.atan2(v.y, v.x);
                    double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
                    //System.out.println("  Angle: " + d1 + " -> " + d2);
                    // compute energy
                    double E = 0.5 * asteroids[i].mass * v2 * v2;
                    // try to push asteroid
                    Asteroid a1 = null;
                    try {
                        a1 = Asteroid.push(asteroids[i], time, E, d2);
                    } catch (InvalidOrbitException e) {
                        //System.out.println("  Invalid orbit: " + e.getMessage());
                        continue;
                    }
                    // avoid allocating a new Point object for every position
                    Point p1 = v, p2 = new Point();
                 // search for collision with other asteroids
                    for (int j = 0; j != asteroids.length; ++j) 
                     {
                         if (i == j) continue;
                         Asteroid a2 = asteroids[j];
                         
                         // look in the future for collision         
                         boolean willCollide = willCollide(a1, a2, time_left_per_asteroid, p1, p2);
                         
                         if (willCollide) 
                         {
                         	debug("COllide !!" + willCollide);
                            pushes.add(new Push(i, E, d2, time_of_push));
//                            energy[i] = E;
//                            direction[i] = d2;
                            // do not push again until collision happens
                            long t = time_of_push - 1;
                            System.out.println("  Collision prediction !");
                            System.out.println("  Year: " + (1 + t / 365));
                            System.out.println("  Day: " + (1 + t % 365));
                            debug("current time: " + time);
                            debug("origin time: " + t + " || ");
                             //iteration number breaks here
                            iteration++;
                            //System.out.println("will overlap: " + willOverlap(p1, a1.radius(), p2, a2.radius()));                         
                            break;
                }
            }
                }
                if(pushes.size() > 0) 
                {
                    Push min_push = pushes.get(0);
                    for(Push push: pushes) 
                    {
                    	System.out.println("Energy: " + push);
                        if(push.energy < min_push.energy) 
                        {
                            min_push = push;
                        }
                    }
                    System.out.println("Min Energy: " + min_push);
                    energy[min_push.asteroid_id] = min_push.energy;
                    direction[min_push.asteroid_id] = min_push.direction;
                    time_of_push = min_push.time_of_push;
                    return;
                }
                time_of_push = time + turns_per_retry;

            }
            if(iteration > 1)
            {
            	 for(int i = 0; i<asteroids.length; i++)
            	 {
            		 int j = getHeaviestAsteroid(asteroids);
            		 Asteroid a2 = asteroids[j];
	                 if (i == j) continue;
	                 Point v = asteroids[i].orbit.velocityAt(time);
	         		 double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
	                 double v2 = 0.0;
	                for (double k=0; k< 1; k=k+0.005) 
	                {
	                	v2 = v1 * (k * 0.45 + 0.05); //Can't change this
	                    //System.out.println("  Speed: " + v1 + " +/- " + v2);
	                    // apply push at -pi/8 to pi/8 of current angle
	                    double d1 = Math.atan2(v.y, v.x);
	                    double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
	                    //System.out.println("  Angle: " + d1 + " -> " + d2);
	                    // compute energy
	                    double E = 0.5 * asteroids[i].mass * v2 * v2;
	                    // try to push asteroid
	                    Asteroid a1 = null;
	                    try {
	                        a1 = Asteroid.push(asteroids[i], time, E, d2);
	                    } catch (InvalidOrbitException e) {
	                        //System.out.println("  Invalid orbit: " + e.getMessage());
	                        continue;
	                    }
	                    // avoid allocating a new Point object for every position
	                    Point p1 = v, p2= new Point();
	                 // search for collision with biggest asteroid
	                     boolean willcollide = willCollide(a1, a2, time_left_per_asteroid, p1, p2);
	                    if(willcollide)
	                     {
	                    	debug("COllide !!" + willcollide);
                            pushes.add(new Push(i, E, d2, time_of_push));
//                            energy[i] = E;
//                            direction[i] = d2;
                            // do not push again until collision happens
                            long t = time_of_push - 1;
                            System.out.println("  Collision prediction !");
                            System.out.println("  Year: " + (1 + t / 365));
                            System.out.println("  Day: " + (1 + t % 365));
                            debug("current time: " + time);
                            debug("origin time: " + t + " || ");
                             //iteration number breaks here
                            iteration++;
                            //System.out.println("will overlap: " + willOverlap(p1, a1.radius(), p2, a2.radius()));                         
                            break;	                     
	                     }
            }
            	 }
	                if(pushes.size() > 0) 
	                {
	                    Push min_push = pushes.get(0);
	                    for(Push push: pushes) 
	                    {
	                    	System.out.println("Energy: " + push);
	                        if(push.energy < min_push.energy) 
	                        {
	                            min_push = push;
	                        }
	                    }
	                    System.out.println("Min Energy: " + min_push);
	                    energy[min_push.asteroid_id] = min_push.energy;
	                    direction[min_push.asteroid_id] = min_push.direction;
	                    time_of_push = min_push.time_of_push;
	                    return;
	                }
	                time_of_push = time + turns_per_retry;

        }
        }
        }
        public int getFarthestAsteroid(Asteroid[] asteroids) 
        {
            int index = 0;
            double maxDistance = 0;
            Point point = new Point();
            for(int aa = 0; aa < asteroids.length; aa++) {
                Asteroid a1 = asteroids[aa];
                a1.orbit.positionAt(time - a1.epoch, point);
                double distance = Point.distance(point, origin);
                if(distance > maxDistance) {
                    maxDistance = distance;
                    index = aa;
                }
            }
            return index;
        }

        public boolean willCollide(Asteroid a11, Asteroid a22, long timeInterval,Point p1, Point p2)
        {    	
        	Asteroid a1, a2;
            //Make sure a1 is always the one has bigger period to short the loop time
            if(a11.orbit.period() > a22.orbit.period())
            {
                a1 = a22;
                a2 = a11;
            }
            else
            {
            	a1 = a11;
            	a2 = a22;
            }
            long threshold = 20;
            ArrayList<Long> cts = getCollisonPoints(a1, a2);
           
            a1.orbit.positionAt(time - a1.epoch, p1);
            a2.orbit.positionAt(time - a2.epoch, p2);
            if(willOverlap(p1, a1.radius(), p2, a2.radius()))
            {
                debug("overlap ***********");
                time_of_push = time + 2;
                return true;
            }

            for(int i = 0; i < cts.size(); i++)
            {
                long startTime = cts.get(i) - threshold;
                long endTime = cts.get(i) + threshold;

                while(startTime < timeInterval + time){
                    //debug("checked Time: " + startTime + " to " + endTime);
                    endTime = endTime > timeInterval + time? timeInterval + time: endTime;
                    for(int t = (int)startTime; t <= endTime; t++){
                        a1.orbit.positionAt(t - a1.epoch, p1);
                        a2.orbit.positionAt(t - a2.epoch, p2);
                        if(willOverlap(p1, a1.radius(), p2, a2.radius()))
                        {
                        	debug("will overlap");
    	                    	time_of_push = t + 1;
    	                        return true;
                        }
                    }
                    startTime = startTime + a1.orbit.period();
                    endTime = endTime + a1.orbit.period();
                }
            }
            return false;

        }

        public boolean willOverlap(Point p1, double r1, Point p2, double r2)
        {
            double distance = Point.distance(p1, p2);
            if(distance < (r1 + r2))
                return true;
            //debug("will not overlap");
            return false;
        }

        public ArrayList<Long> getCollisonPoints(Asteroid a1, Asteroid a2)
        {
            double period = a1.orbit.period();
            Point tmp = new Point();
            Point center = a2.orbit.center();
            double a = a2.orbit.a;
            double threshold = 2 * a1.radius();
            Point foci1 = new Point(0, 0);
            Point foci2 = new Point(2 * center.x, 2 * center.y);
            ArrayList<Long> cts = new ArrayList<Long>();
            for(int i = 0; i < period; i++){
                a1.orbit.positionAt(time + i - a1.epoch, tmp);
                double distance =Point.distance(foci1, tmp) + Point.distance(foci2, tmp);
                if(Math.abs(2 * a - distance) <= threshold)
                {
                    Long ct = time + i;
                    cts.add(ct);
                }
            }
            ArrayList<Long> ret = new ArrayList<Long>();
            if(cts.size() > 0){
                ret.add(cts.get(0));
                for(int i = 1; i < cts.size(); i++)
                {
                    if(cts.get(i) - cts.get(i - 1) > 1)
                    {
                        ret.add(cts.get(i));
                    }
                }
            }

            return ret;
        }

        public void debug(String str)
        {
            //System.out.print("debug: " + str + "\n");
        }

        private int getLightestAsteroid(Asteroid asteroids[]) 
        {
            int min = 1;
            for (int i = 0; i < asteroids.length; i++)
            {
                if (asteroids[i].mass < asteroids[min].mass)
                    min = i;
            }
            return min;
        }

        private int getHeaviestAsteroid(Asteroid asteroids[])
        {
            int max = 1;
            for (int i = 0; i < asteroids.length; i++)
            {
                if (asteroids[i].mass > asteroids[max].mass)
                    max = i;
            }
            return max;
        }

    }
