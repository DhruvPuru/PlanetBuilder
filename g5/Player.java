package pb.g5;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Arrays;
import java.util.Random;

public class Player implements pb.sim.Player {

	// used to pick asteroid and velocity boost randomly

	// current time, time limit
	private long time = -1;
	private long time_limit;
	
	private boolean collisionImminent;
	private long nextCollision;
	private int n_asteroids;
	
	private Integer[] closestPair;
	
	int pushed = -1;

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit)
	{
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
		collisionImminent = false;
		nextCollision = Long.MIN_VALUE;
		n_asteroids = asteroids.length;
		closestPair = null;
	}
	
	private Integer[] getClosestPairAtTime(Asteroid[] asteroids, long time) {
		double closestDist = Double.MAX_VALUE;
		Integer[] pair = new Integer[2];
		Asteroid a1, a2;
		a1 = a2 = null;
		for(int i = 0; i < asteroids.length; ++i) {
			a1 = asteroids[i];
			Point p1 = a1.orbit.positionAt(time - a1.epoch);
			for(int j = i+1; j < asteroids.length; ++j) {
				a2 = asteroids[j];
				Point p2 = a2.orbit.positionAt(time - a2.epoch);
				double dist = Point.distance(p1,p2);
				if(dist < closestDist) {
					closestDist = dist;
					pair[0] = i;
					pair[1] = j;
				}
			}
		}
		return pair;
	}
	
	private Integer[] getClosestApproachToTargetWithinTime(Asteroid[] asteroids, int target, long time) {
		double closestDist = Double.MAX_VALUE;
		Integer[] pair = new Integer[2];
		Asteroid a1, a2;
		a1 = asteroids[target];
		a2 = null;
		Point p1, p2;
		p1 = a1.orbit.positionAt(time - a1.epoch);
		p2 = null;
		for(long t = 0; t <= 1000; ++t) {
			for(int i = 0; i < asteroids.length; ++i) {
				if(i == target)
					continue;
				a2 = asteroids[i];
				p2 = a2.orbit.positionAt(time+t - a2.epoch);
				double dist = Point.distance(p1,p2);
				if(dist < closestDist) {
					closestDist = dist;
					pair[0] = i;
					pair[1] = target;
				}
			}
		}
		return pair;
	}
	
	private int getHeaviestAsteroid(Asteroid[] asteroids) {
		double heaviest = Double.MIN_VALUE;
		int index = -1;
		for(int i = 0; i < asteroids.length; ++i) {
			if(heaviest < asteroids[i].mass) {
				heaviest = asteroids[i].mass;
				index = i;
			}
		}
		return index;
	}

	// try to push asteroid
	public void play(Asteroid[] asteroids,
	                 double[] energy, double[] direction)
	{
		++time;
				
		if(time % 365 == 0) {
			System.out.println("Year: " + (1 + time / 365));
			System.out.println("Day: "  + (1 + time % 365));
		}

		if(asteroids.length < n_asteroids) {
			collisionImminent = false;
			--n_asteroids;
		}
				
		if(collisionImminent && time == nextCollision) {
			Asteroid a1 = asteroids[closestPair[0]];
			Point p1 = a1.orbit.positionAt(time-a1.epoch);
			Asteroid a2 = asteroids[closestPair[1]];
			Point p2 = a2.orbit.positionAt(time-a2.epoch);
			if(Point.distance(p1, p2) > a1.radius()+a2.radius()) {
				System.out.println(Point.distance(p1,p2));
				System.exit(0);
			}
		}
			
		if(collisionImminent) {
			return;
		}
			
		int heaviest = getHeaviestAsteroid(asteroids);
		closestPair = getClosestApproachToTargetWithinTime(asteroids, heaviest, time+1000);
		Asteroid a1 = asteroids[closestPair[0]];
		Asteroid a2 = asteroids[closestPair[1]];
		int pushedIndex = closestPair[0];
		if(a2.mass < a1.mass) {
			Asteroid temp = a1;
			a1 = a2;
			a2 = temp;
			pushedIndex = closestPair[1];
		}
		GradientDescent gd = new GradientDescent(a1, a2, time);
		Push p = gd.tune();
		if(p != null) {
			int i = pushedIndex;
			energy[i] = p.energy;
			direction[i] = p.direction;
			collisionImminent = true;
			nextCollision = time + gd.predictedTime;
			return;
		} else if (1 > 0) {
			return;
		}
	}
}
