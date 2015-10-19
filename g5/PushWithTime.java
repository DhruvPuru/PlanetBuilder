package pb.g5;

public class PushWithTime {
	Push push;
	int index;
	long time;
	long predictedTime;
	
	public PushWithTime(Push push, int index, long time, long predictedTime) {
		super();
		this.push = push;
		this.index = index;
		this.time = time;
		this.predictedTime = predictedTime;
	}
}
