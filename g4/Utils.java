// public class Utils {
// 	public static void optimize(ArrayList<Double> mass, Asteroid[] asteroids) {
// 		System.out.println("Optimizing");
// 		ArrayList<Double> r = new ArrayList<Double>();
// 		String results = null;
// 		ArrayList<String[]> r1 = new ArrayList<String[]>();
// 		Set<Asteroid> r2 = new HashSet<Asteroid>();
// 		ArrayList<String> r3 = new ArrayList<String>();
// 		for(double i: mass)
// 		{
// 			r3.add("0");
// 			r.add(0.0);
// 		}
// 		for (int i = 0; i < mass.size(); i++)
// 		{
// 			double q = 0.0;
// 			String q1 = "0 0 0";
// 			for (int j = 0; j <= i; j++ )
// 			{
// 				if(q < (mass.get(j) + r.get(i-j)))
// 				{
// 					q =  mass.get(j) + r.get(i-j);
// 					q1 = mass.get(j)+ " " + r3.get(i-j);
// 					if( mass.get(j)!= 0 && r.get(i-j)!=0)
// 					{
// 						results = mass.get(j) + " " + r.get(i-j) + " " + q + "\n";
// 						r1.add(results.split(" "));
// 					}
// 				}
// 			}
// 			r.set(i,q);
// 			r3.set(i,q1);
// 		}

// 		String[] last = r3.get(r3.size()-1).split(" ");
// 		for(int k = 0; k < asteroids.length; k++)
// 		{
// 			for(int j = 0; j < last.length; j++)
// 			{
// 				if(Double.parseDouble(last[j]) == asteroids[k].mass)
// 				{
// 					r2.add(asteroids[k]);
// 				}
// 			}
// 		}
// 		asteroidOrder = r2;
// 	}

// 	public static void dynamicProgramming(Asteroid[] asteroids)
// 	{
// 		System.out.println("Starting Dynamic Programming");
// 		ArrayList<Double> masses = new ArrayList<Double>(); //exp
// 		ArrayList<Double> energies = new ArrayList<Double>(); //stam
// 		masses.add(0.0);
// 		energies.add(0.0);
// 		for(double i = Math.pow(10, 35); i <= Math.pow(10, 39); i += Math.pow(10,35))
// 		{
// 			boolean write = false;
// 			for(int j = 0; j < asteroids.length; j++)
// 			{
// 				Point v = asteroids[j].orbit.velocityAt(time);
// 				Double velocity = Math.sqrt(v.x * v.x + v.y * v.y);
// 				Double mass = asteroids[j].mass;
// 				Double energy = 0.5* mass * velocity * velocity;
// 				// System.out.println("energy: " + energy);
// 				if(energy >= i && energy < (i + Math.pow(10, 35)))
// 				{
// 					masses.add(asteroids[j].mass);
// 					energies.add(energy);
// 					write = true;
// 				}
// 			}
// 			if (write == false)
// 			{
// 				masses.add(0.0);
// 				energies.add(0.0);
// 			}
// 		}
// 		optimize(masses, asteroids);
// 	}

// 	public void perturb(Asteroid[] asteroids, double[] energy, double[] direction)
// 	{
// 		Asteroid a1;
// 		//It's been a long time since a push, so perturb the system
// 		if (timeSincePush > 7300) {
// 			boolean validOrbitNotFound = true;
// 			while (validOrbitNotFound) {
// 				int i = random.nextInt(asteroids.length);
// 				Point v = asteroids[i].orbit.velocityAt(time);
// 				// add 5-50% of current velocity in magnitude
// 				double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
// 				double v2 = v1 * (random.nextDouble() * 0.45 + 0.05);
// 				// apply push at -π/8 to π/8 of current angle
// 				double d1 = Math.atan2(v.y, v.x);
// 				double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
// 				// compute energy
// 				double E = 0.5 * asteroids[i].mass * v2 * v2;
// 				try {
// 					a1 = Asteroid.push(asteroids[i], time, E, d2);
// 					validOrbitNotFound = false;
// 				} catch (InvalidOrbitException e) {
// 					System.out.println("  Invalid orbit: " + e.getMessage());
// 					validOrbitNotFound = false;
// 					continue;
// 				}
// 				energy[i] = E;
// 				direction[i] = d2;
// 			}
// 		}
// 	}
// }