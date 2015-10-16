package pb.g5;

import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;
import pb.sim.Point;

public class GradientDescent {
	Asteroid a;
	Asteroid target;
	long timeOfStart, predictedTime;

	public GradientDescent(Asteroid a, Asteroid target, long time) {
		super();
		this.a = a;
		this.target = target;
		timeOfStart = time;
		this.predictedTime = 0;
	}
	
	public Push tune() {
		double direction = 0;
		double energy = varyEnergy(direction);
		for(int i = 0; i < 3; ++i) {
			//System.out.println(direction + "," + energy);
			direction = varyDirection(energy);
			energy = varyEnergy(direction);
		}
//		System.out.println(direction + "," + energy);
		Asteroid pushed = Asteroid.push(a, timeOfStart, energy, direction);
		predictedTime = getTimeOfClosestDistance(pushed, timeOfStart);
		double dist = getClosestDistanceVaryTime(pushed, timeOfStart);
//		System.out.println("DIST : "+dist);
//		System.out.println("RAD  : "+(a.radius()+target.radius()));
//		System.out.println("DELTA: "+getDistance(a, timeOfStart));
//		System.out.println("TIME : "+predictedTime);
		if(dist < a.radius() + target.radius()) {
			return new Push(energy, direction);
		}
		return null;
	}
	
	
	private double getClosestDistanceVaryTime(Asteroid a, long timeOfStart) {
		int delta = 200;
		long timeDel = 1000;
		double current = getDistance(a, timeOfStart + timeDel);
		while(delta >= 1) {
			double later = getDistance(a,timeOfStart + timeDel + delta);
			double earlier = getDistance(a,timeOfStart + timeDel - delta);
			if(later < current) {
				timeDel += delta;
				current = later;
			} else if(earlier < current) {
				timeDel -= delta;
				current = earlier;
			} else {
				if(delta >= 10)
					delta /= 10;
				else if(delta > 1)
					delta = 1;
				else
					delta = 0;
			}
			while(timeDel <= delta) {
				delta /= 2;
			}
		}
		return current;
	}
	
	private long getTimeOfClosestDistance(Asteroid a, long timeOfStart) {
		int delta = 200;
		long timeDel = 1000;
		double current = getDistance(a, timeOfStart + timeDel);
		while(delta >= 1) {
			double later = getDistance(a,timeOfStart + timeDel + delta);
			double earlier = getDistance(a,timeOfStart + timeDel - delta);
			if(later < current) {
				timeDel += delta;
				current = later;
			} else if(earlier < current) {
				timeDel -= delta;
				current = earlier;
			} else {
				if(delta >= 10)
					delta /= 10;
				else if(delta > 1)
					delta = 1;
				else
					delta = 0;
			}
			while(timeDel <= delta) {
				delta /= 2;
			}
		}
		return timeDel;
	}

	
	private double varyDirection(double pushEnergy) {
		double currentDirection = 0;
		double delta = Math.PI/6;
		double current = Double.MAX_VALUE;
		try {
			Asteroid pushed = Asteroid.push(a, timeOfStart, pushEnergy, currentDirection);
			current = getClosestDistanceVaryTime(pushed, timeOfStart);
		} catch (InvalidOrbitException e) {
			// do nothing;
		}
		double angleDel = 0;
		while(delta >= Math.PI/288) {
			double dh = Double.MAX_VALUE;
			double dl = Double.MAX_VALUE;
			
			try {
				Asteroid higher = Asteroid.push(a, timeOfStart, pushEnergy, currentDirection + delta);
				dh = getClosestDistanceVaryTime(higher, timeOfStart);
			} catch(InvalidOrbitException e) {
				// do nothing;
			}
			try {
				Asteroid lower = Asteroid.push(a, timeOfStart, pushEnergy, currentDirection - delta);
				dl = getClosestDistanceVaryTime(lower, timeOfStart);
			} catch (InvalidOrbitException e) {
				// do nothing;
			}
			
			if(dh < current) {
				angleDel += delta;
				current = dh;
			} else if(dl < current) {
				angleDel -= delta;
				current = dl;
			} else {
				if(delta >= Math.PI/36)
					delta /= 2;
				else if(delta > Math.PI/288)
					delta = Math.PI/288;
				else
					delta = 0;
			}
		}
//		System.out.println("DirectionCount: "+count);
		return angleDel;
	}
	
	private double varyEnergy(double pushDirection) {
		Point vel = a.orbit.velocityAt(timeOfStart - a.epoch);
		double asteroidEnergy = Math.sqrt(vel.x*vel.x + vel.y*vel.y) * a.mass * 0.5;
		double currentEnergy = asteroidEnergy * 0.65;
		double delta = asteroidEnergy * 0.25;
		Asteroid pushed = null;
		try{
			pushed = Asteroid.push(a, timeOfStart, currentEnergy, pushDirection);
		} catch (Exception e) {
			System.out.println(timeOfStart+","+currentEnergy+","+pushDirection);
		}
		double current = getClosestDistanceVaryTime(pushed, timeOfStart);
		while(delta >= asteroidEnergy * 0.01) {
			double dh = Double.MAX_VALUE;
			double dl = Double.MAX_VALUE;
			
			try {
				Asteroid higher = Asteroid.push(a, timeOfStart, currentEnergy + delta, pushDirection);
				dh = getClosestDistanceVaryTime(higher, timeOfStart);
			} catch(InvalidOrbitException e) {
				// do nothing;
			}
			try {
				Asteroid lower = Asteroid.push(a, timeOfStart, currentEnergy - delta, pushDirection);
				dl = getClosestDistanceVaryTime(lower, timeOfStart);
			} catch (InvalidOrbitException e) {
				// do nothing;
			}

			if(dl <= current) {
				currentEnergy -= delta;
				current = dl;
			} else if(dh < current) {
				currentEnergy += delta;
				current = dh;
				delta += delta;
			} else {
				if(delta >= asteroidEnergy * 0.02)
					delta /= 2;
				else if(delta > asteroidEnergy * 0.01)
					delta = asteroidEnergy * 0.01;
				else
					delta = 0;
			}
			while(delta >= currentEnergy) {
				delta = currentEnergy/2;
			}
		}
//		System.out.println("EnergyCount: "+count);
		return currentEnergy;
	}
	
	private double getDistance(Asteroid a, long time) {
		Point p1 = target.orbit.positionAt(time - target.epoch);
		Point p2 = a.orbit.positionAt(time - a.epoch);
		double sumRadii = target.radius() + a.radius();
		double dist = Point.distance(p1, p2);
		if(dist <= sumRadii)
			return 0;
		return dist;
	}
}
