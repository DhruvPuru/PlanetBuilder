package pb.sim;

public class InvalidOrbitException extends RuntimeException {

	public InvalidOrbitException(String message)
	{
		super(message);
	}

	public static final long serialVersionUID = 0xdeadbeefcafebabeL;
}
