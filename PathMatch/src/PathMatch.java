
public class PathMatch {
	static int graphSize;
	static int[][] graph;
	static int[][] distance;
	static int[][] path;
	
	
	public static void readGraph() {}
	public static void readQuery() {}
	public static void readCorrespondence() {}
	public static void generateQueryCorrespondences() {}
	public static void calculateQueryOrders() {}
	public static void calculateQueryWeights() {}
	public static void queryTopPaths() {}
	public static void floyd() {
		// initialize
		for(int i = 0; i < graphSize; i++) {
			distance[i][i] = 0;
			for(int j = 0; j < graphSize; j++) {
				path[i][j] = (i == j) ? i : -1;
				distance[i][j] = graph[i][j];
			}
		}
			
		// k - iterations
		for(int k = 0; k < graphSize; k++) {
			for(int i = 0; i < graphSize; i++) {
				for(int j = 0; j < graphSize; j++) {
					if(distance[i][k] != Short.MAX_VALUE &&
						distance[k][j] != Short.MAX_VALUE &&
						distance[i][k] + distance[k][j] < distance[i][j]) {
						
						distance[i][j] = distance[i][k] + distance[k][j];
						path[i][j] = path[k][j];
					}
				}
			}
		}
	}
	
}
