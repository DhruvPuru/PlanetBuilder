/*public class Utils {
	public static void optimize(ArrayList<Double> mass, Asteroid[] asteroids) {
		System.out.println("Optimizing");
		ArrayList<Double> r = new ArrayList<Double>();
		String results = null;
		ArrayList<String[]> r1 = new ArrayList<String[]>();
		Set<Asteroid> r2 = new HashSet<Asteroid>();
		ArrayList<String> r3 = new ArrayList<String>();
		for(double i: mass)
		{
			r3.add("0");
			r.add(0.0);
		}
		for (int i = 0; i < mass.size(); i++)
		{
			double q = 0.0;
			String q1 = "0 0 0";
			for (int j = 0; j <= i; j++ )
			{
				if(q < (mass.get(j) + r.get(i-j)))
				{
					q =  mass.get(j) + r.get(i-j);
					q1 = mass.get(j)+ " " + r3.get(i-j);
					if( mass.get(j)!= 0 && r.get(i-j)!=0)
					{
						results = mass.get(j) + " " + r.get(i-j) + " " + q + "\n";
						r1.add(results.split(" "));
					}
				}
			}
			r.set(i,q);
			r3.set(i,q1);
		}

		String[] last = r3.get(r3.size()-1).split(" ");
		for(int k = 0; k < asteroids.length; k++)
		{
			for(int j = 0; j < last.length; j++)
			{
				if(Double.parseDouble(last[j]) == asteroids[k].mass)
				{
					r2.add(asteroids[k]);
				}
			}
		}
		asteroidOrder = r2;
	}

	public static void dynamicProgramming(Asteroid[] asteroids)
	{
		System.out.println("Starting Dynamic Programming");
		ArrayList<Double> masses = new ArrayList<Double>(); //exp
		ArrayList<Double> energies = new ArrayList<Double>(); //stam
		masses.add(0.0);
		energies.add(0.0);
		for(double i = Math.pow(10, 35); i <= Math.pow(10, 39); i += Math.pow(10,35))
		{
			boolean write = false;
			for(int j = 0; j < asteroids.length; j++)
			{
				Point v = asteroids[j].orbit.velocityAt(time);
				Double velocity = Math.sqrt(v.x * v.x + v.y * v.y);
				Double mass = asteroids[j].mass;
				Double energy = 0.5* mass * velocity * velocity;
				// System.out.println("energy: " + energy);
				if(energy >= i && energy < (i + Math.pow(10, 35)))
				{
					masses.add(asteroids[j].mass);
					energies.add(energy);
					write = true;
				}
			}
			if (write == false)
			{
				masses.add(0.0);
				energies.add(0.0);
			}
		}
		optimize(masses, asteroids);
	}

	  public void storeMass(Asteroid[] asteroids) {
    double mass_sum = 0;
    for(Asteroid asteroid: asteroids)
    {
      cached_asteroid_masses.put(asteroid, asteroid.mass);
      mass_sum += asteroid.mass;
      System.out.println(asteroid.mass);
    }
    fifty_percent_mass = 0.5*mass_sum;
    System.out.println("50% mass: " + fifty_percent_mass);
  } 

  public void updateMass(Asteroid asteroid1, Asteroid asteroid2, Asteroid[] asteroids) {
    cached_asteroid_masses.remove(asteroid1);
    cached_asteroid_masses.remove(asteroid2);
    for(Asteroid asteroid: asteroids)
    {
      if(!cached_asteroid_masses.containsKey(asteroid))
      {
        cached_asteroid_masses.put(asteroid, asteroid.mass);
      }
    }
  } 

    private void printMassVelocity(Asteroid[] asteroids) {
    for (Asteroid asteroid: asteroids) {
      System.out.println("mass, velocity:" + asteroid.mass + ", " + asteroid.orbit.velocityAt(time));
    }
  }
}*/
