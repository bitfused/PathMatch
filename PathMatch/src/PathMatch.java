import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;


public class PathMatch {
	static String INPUT_NAME = "input";
	static String QUERY_NAME = "query";
    static String SCORE_NAME = "corr";
	static int graphSize;
	static int[][] graph;
	static int[][] distance;
	static int[][] path;
	static List<String> names;
	static List<String> query;


    private static void readCorr() throws IOException {
        


    }



    public static void readGraph() throws IOException {
        // Input file contains the input graph
        // Each line defines all the edges from the first vertex to all the other vertices.
        // The edges/lines represents the vertices that are related (proteins or metabolic processes that are related).

        // v1 v2 v3 v4 ==> 2 directed edges:  (v1, v3) and (v1, v4)
        //                 1 undirected edge: (v1, v2)
        // v2 v1 v3 v4 ==> 2 directed edges:  (v2, v3) and (v2, v4)
        //                 1 undirected edge: (v2, v1)

		names = new LinkedList<>();
		FileReader input = new FileReader(INPUT_NAME);
		BufferedReader bufRead = new BufferedReader(input);
		String myLine = null;

		while ((myLine = bufRead.readLine()) != null) {

            // Add the first vertex of each line into the linked list
			names.add(myLine.substring(0, myLine.indexOf(" ")));
		}
		
		// initialize graph adj. matrix once we know the number of nodes
		graphSize = names.size();
		graph = new int[graphSize][graphSize];
		
		for(int i = 0; i < graphSize; i++) {
			for(int j = 0; j < graphSize; j++) {
				graph[i][j] = Short.MAX_VALUE;
			}
		}
		
		bufRead.close();
		input.close();

        // Reopen the input file
		input = new FileReader(INPUT_NAME); // for some reason reset is not supported...
		bufRead = new BufferedReader(input);
		while ((myLine = bufRead.readLine()) != null) {    
			String[] line = myLine.split(" ");

            // For each line (first vertex), add the vertices (2nd vertex onwards) that it connects to.
			int i = names.indexOf(line[0]);
			for (int n = 1; n < line.length; n++) {
				int j = names.indexOf(line[n]);
				if (j != -1) {
					graph[i][j] = 1;
				}
			}
		}
		bufRead.close();
		input.close();
	}
	
	public static void readQuery() throws IOException {
        // Query file contains contains the query path in sorted order (line by line).
        // Each line contains the name of the vertex.

		query = new LinkedList<>();
		FileReader fr = new FileReader(QUERY_NAME);
		BufferedReader bufRead = new BufferedReader(fr);
		String myLine = null;

		while ((myLine = bufRead.readLine()) != null) {   
			query.add(myLine);
		}
		bufRead.close();
		fr.close();
	}
	
	public static void readCorrespondence() {
		
		
	}
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
	
	public static void main(String [] args) {
		try {
			readGraph();
			readQuery();
            readCorr();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//System.out.println(Arrays.deepToString(graph));
		System.out.println(query.toString());		
	}

}
