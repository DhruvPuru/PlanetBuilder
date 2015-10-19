package pb.g5;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

public class Player implements pb.sim.Player {

	// used to pick asteroid and velocity boost randomly

	// current time, time limit
	private long time = -1;
	private long time_limit;
	
	private boolean haveQueuedPush;
	private PushWithTime queuedPush;
	private boolean collisionImminent;
	private long nextCollision;
	private int n_asteroids;
	private int starting_n_asteroids;
	
	private Integer[] closestPair;
	
	private HashMap<Long, Long> highEnergyCollisionAsteroids;
	
	long pushed = -2;

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit)
	{
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
		collisionImminent = false;
		nextCollision = Long.MIN_VALUE;
		n_asteroids = asteroids.length;
		starting_n_asteroids = asteroids.length;
		closestPair = null;
		
		highEnergyCollisionAsteroids = new HashMap<Long, Long>();
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
				if(highEnergyCollisionAsteroids.containsKey(a2.id)) {
					continue;
				}
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
	
	private int findAsteroidToPush(Asteroid[] asteroids, int targetIndex, long time) {
		Asteroid target = asteroids[targetIndex];
		Point targetLocation = target.orbit.positionAt(time - target.epoch);
		double targetTheta = Math.atan2(targetLocation.y, targetLocation.x);
		Integer[] sortedByPhase = new Integer[asteroids.length];
		for(int i = 0; i < sortedByPhase.length; ++i) {
			sortedByPhase[i] = i;
		}
		
		Arrays.sort(sortedByPhase, new Comparator<Integer>() {
			public int compare(final Integer o1, final Integer o2) {
				Asteroid a1 = asteroids[o1];
				Point p1 = a1.orbit.positionAt(time - a1.epoch);
				double theta1 = targetTheta - Math.atan2(p1.y,p1.x);
				Asteroid a2 = asteroids[o2];
				Point p2 = a2.orbit.positionAt(time - a2.epoch);
				double theta2 = targetTheta - Math.atan2(p2.y,p2.x);
				while(theta1 < 0) theta1 += Math.PI*2;
				while(theta2 < 0) theta2 += Math.PI*2;
				return Double.compare(theta1, theta2);
			}
		});
		
		
		double minTheta = 10*Math.PI/180;
		double maxTheta = 30*Math.PI/180;
		
		double minDist = Double.MAX_VALUE;
		int minID = -1;
		for(int i = 1; i < sortedByPhase.length; ++i) {
			int idx = sortedByPhase[i];
			Asteroid a1 = asteroids[idx];
			Point p1 = a1.orbit.positionAt(time - a1.epoch);
			double theta1 = targetTheta - Math.atan2(p1.y, p1.x);
			while(theta1 < 0) theta1 += Math.PI*2;
			if(theta1 >= minTheta && theta1 <= maxTheta) {
				double dist = Point.distance(p1, targetLocation);
				if(dist < minDist) {
					minDist = dist;
					minID = idx;
				}
			}
		}
		return minID;
	}
	
	private int getKthAsteroid(Asteroid[] asteroids, int k) {
		Integer[] sortedByOrbitRadius = new Integer[asteroids.length];
		for(int i = 0; i < sortedByOrbitRadius.length; ++i) {
			sortedByOrbitRadius[i] = i;
		}
		Arrays.sort(sortedByOrbitRadius, new Comparator<Integer>() {
			public int compare(final Integer o1, final Integer o2) {
				return Double.compare(asteroids[o1].orbit.a, asteroids[o2].orbit.a);
			}
		});
		for(int i = 0; i < sortedByOrbitRadius.length && time == 0; ++i) {
//			System.out.println(asteroids[sortedByOrbitRadius[i]].orbit.a);
		}
//		if(time == 0) System.out.println(asteroids[sortedByOrbitRadius[k]].orbit.positionAt(0).magnitude());
		return sortedByOrbitRadius[k];
	}
	
	private int getOutermostAsteroid(Asteroid[] asteroids, int targetIndex) {
		double maxA = Double.MIN_VALUE;
		int idx = -1;
		for(int i = 0; i < asteroids.length; ++i) {
			if(i == targetIndex)
				continue;
			double a = asteroids[i].orbit.a;
			if(a > maxA) {
				maxA = a;
				idx = i;
			}
		}
		return idx;
	}
	
	private boolean withinTheta(Asteroid target, Asteroid toPush, long time, double minDeg, double maxDeg) {
		Point t = target.orbit.positionAt(time-target.epoch);
		Point p = toPush.orbit.positionAt(time-toPush.epoch);
		double tTheta = Math.atan2(t.y, t.x);
		double pTheta = Math.atan2(p.y, p.x);
		double dt = (tTheta - pTheta) * 180 / Math.PI;
		while(dt < 0) dt += 360;
		if(minDeg <= dt && dt <= maxDeg)
			return true;
		return false;
	}
	
	private int getEllipticalAsteroid(Asteroid[] asteroids) {
		for(int i = 0; i < asteroids.length; ++i) {
			if(Math.abs(asteroids[i].orbit.a - asteroids[i].orbit.b) > 1e-9)
				return i;
		}
		return -1;
	}

	// try to push asteroid
	public void play(Asteroid[] asteroids,
	                 double[] energy, double[] direction)
	{
		++time;
				
		if(time % 365 == 0 || time/365==70) {
			System.out.println("Year: " + (1 + time / 365));
			System.out.println("Day: "  + (1 + time % 365));
		}
		
		if(time/365==70 && time%365==158) {
			int a = 0;
			int b = a+1;
		}
		
		if(haveQueuedPush) {
			if(queuedPush.time == time) {
				energy[queuedPush.index] = queuedPush.push.energy;
				direction[queuedPush.index] = queuedPush.push.direction;
				collisionImminent = true;
				nextCollision = queuedPush.time + queuedPush.predictedTime;
				pushed = queuedPush.time;
				haveQueuedPush = false;
			} else {
				return;
			}
		}
		
		if(asteroids.length < n_asteroids) {
			collisionImminent = false;
			--n_asteroids;
		}
		
		/*if(time == pushed + 1) {
			Scanner s = new Scanner(System.in);
			System.out.print("Type g: ");
			while(!s.next().equals("g")) {
				System.out.print("Type g: ");
			}
		}*/

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
		
		ArrayList<Long> toRemove = new ArrayList<Long>();
		for(Long id : highEnergyCollisionAsteroids.keySet()) {
			long lastCheckTime = highEnergyCollisionAsteroids.get(id);
			if(time - lastCheckTime > 365) {
				toRemove.add(id);
			}
		}
		for(Long id : toRemove) {
//			System.out.println("freeing "+id);
			highEnergyCollisionAsteroids.remove(id);
		}
		
		int targetIndex = getHeaviestAsteroid(asteroids);
		if(starting_n_asteroids == asteroids.length)
			targetIndex = getKthAsteroid(asteroids, (int) (0.8 * asteroids.length));
//		if(time == 0) {
//			targetIndex = getKthAsteroid(asteroids, Math.round(asteroids.length-1));
//		}
		//int targetIndex = getOutermostAsteroid(asteroids, -1);
		
		closestPair = getClosestApproachToTargetWithinTime(asteroids, targetIndex, time+1000);
		Asteroid a1 = asteroids[closestPair[0]];
		Asteroid a2 = asteroids[closestPair[1]];
		int pushedIndex = closestPair[0];
		if(a2.mass < a1.mass) {
			Asteroid temp = a1;
			a1 = a2;
			a2 = temp;
			pushedIndex = closestPair[1];
		}
/*		closestPair = new Integer[2];
		closestPair[1] = targetIndex;
		closestPair[0] = getOutermostAsteroid(asteroids, targetIndex);//findAsteroidToPush(asteroids, targetIndex, time);
		Asteroid a1 = asteroids[closestPair[0]];
		Asteroid a2 = asteroids[closestPair[1]];
//		if(!withinTheta(a2, a1, time, 10, 45))
//			return;*/
//		int pushedIndex = closestPair[0];
		long pushTime = time+3000;
		GradientDescent gd = new GradientDescent(a1, a2, pushTime);
		Push p = gd.tune();
		if(p != null) {
//			System.out.println("Found push looking for better");
			PushWithTime pt = tunePushTime(a1, a2, new PushWithTime(p, pushedIndex, pushTime, gd.predictedTime), time);
			Point p1 = a1.orbit.positionAt(pt.time-a1.epoch);
			Point p2 = a2.orbit.positionAt(pt.time-a2.epoch);
			double t1 = Math.atan2(p1.y, p1.x)*180/Math.PI;
			double t2 = Math.atan2(p2.y, p2.x)*180/Math.PI;
//			System.out.println("Angle:"+(t2-t1));
			if(pt.push.energy > 5*Math.pow(10, 36)) {
				//System.out.println("energy too high skipping");
				highEnergyCollisionAsteroids.put(a1.id, time);
				//System.out.println("locking "+a1.id);
				return;
			}
			if(pt != null) {
				haveQueuedPush = true;
				queuedPush = pt;
//				System.out.println("found better setting for future: "+pt.time);
				return;
			}
//			System.out.println("Did not find better");
			energy[pushedIndex] = p.energy;
			direction[pushedIndex] = p.direction;
			collisionImminent = true;
			nextCollision = time + gd.predictedTime;
			pushed = time;
			return;
		}
	}
	
	private PushWithTime tunePushTime(Asteroid a1, Asteroid a2, PushWithTime resultantPush, long timeNow) {
		long bestTime = resultantPush.time;
		Point p1 = a1.orbit.positionAt(time-a1.epoch);
		Point p2 = a2.orbit.positionAt(time-a2.epoch);
		long delta = 300;
		Push bestPush = resultantPush.push;
		long predictedTime = resultantPush.predictedTime;
		int numInDirection = 0;
		while(delta > 1) {
			if(Math.abs(numInDirection) > 4) {
				delta += delta;
				numInDirection = 0;
			}
			
			GradientDescent l = new GradientDescent(a1, a2, bestTime+delta);
			GradientDescent e = new GradientDescent(a1, a2, bestTime-delta);
			
			Push later = l.tune();
			Push earlier = e.tune();
			
			if(later != null && later.energy < bestPush.energy) {
				bestTime += delta;
				bestPush = later;
				predictedTime = l.predictedTime;
				if(numInDirection < 0) numInDirection = 1;
				else ++numInDirection;
			} else if(earlier != null && earlier.energy < bestPush.energy) {
				bestTime -= delta;
				bestPush = earlier;
				predictedTime = e.predictedTime;
				if(numInDirection < 0) numInDirection = -1;
				else --numInDirection;
			} else {
				delta /= 2;
			}
		}
		
		if(bestTime < timeNow) {
//			System.out.println("tried to push too soon: "+bestTime+"\tnow: "+timeNow);
			System.exit(0);
		}
//		System.out.println("StartingE: "+resultantPush.push.energy+"\tFinalE: "+bestPush.energy);
		return new PushWithTime(bestPush, resultantPush.index, bestTime, predictedTime);
	}
}
