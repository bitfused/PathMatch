import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.Buffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.jgrapht.alg.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;


public class PathMatch {

    /* CONSTANTS */
    static int NUM_TOP_PATHS = 10;                  // Number of pathways with top scores that are required to be in output
    static int MAX_NUM_MISMATCHES_AND_GAPS = 1;    // Maximum number of mismatches and gaps allowed
    static double INDEL_PENALTY = 23.0;             // Penalty for mismatches and gaps
    static double CORR_E_VALUE_CUTOFF= 1e200;      // Cutoff for determining similarity

    static String INPUT_NAME = "input";
	static String QUERY_NAME = "query";
    static String CORR_NAME = "corr";
	static int graphSize;
    static double[][] correspondence;
	static int[][] graph;
	static int[][] distance;
    static List<List<Vertex>> dag;
    static int[][] path;
	static List<String> names;           // List of nodes in input graph - G
	static List<String> query;           // List of nodes in query path  - Q

    public static class Vertex {
        String fromNode;
        double weight;
        List<Edge> edges;

        public Vertex(String name, double wt) {
            fromNode = name;
            weight = wt;
            edges = new ArrayList<>();
        }
        public void addEdge(Edge e) {
            edges.add(e);
        }
        public String getNode() {
            return fromNode;
        }
    }

    public static class Edge {
        String toNode;
        double weight;

        public Edge(String name, double wt) {
            toNode = name;
            weight = wt;
        }
    }


    private static void readCorrespondence() throws IOException {
        // Correspondence file contains the correspondences (or substitution scores)
        // Each line contains the correspondences between vertices in the query path and the input graph
        // Log the correspondences if they are below the cutoff value and are from blast2seq pathway.

        // Initialize the correspondence matrix
        correspondence = new double[query.size()][graphSize];


        // Parse the file and fill in the graph
        FileReader corr = new FileReader(CORR_NAME);
        BufferedReader br = new BufferedReader(corr);
        String myLine = null;

        String v1;
        String v2;
        double value;

        while ((myLine = br.readLine()) != null) {

            String[] line = myLine.split(" +");

            v1 = line[0];
            if (line[1].indexOf("|") == -1) {
                v2 = line[1];
            } else {
                v2 = line[1].substring(0, line[1].indexOf("|"));
            }
            value = Double.parseDouble(line[2]);

            int i = query.indexOf(v1);
            int j = names.indexOf(v2);

            if (value <= CORR_E_VALUE_CUTOFF) {
                if(value <= 1e-299){
                    correspondence[i][j]=Math.log(1e-299);
                } else {
                    correspondence[i][j]=Math.log(value);
                }
            }
        }
            br.close();
            corr.close();
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
            String full = myLine.substring(0, myLine.indexOf(" "));

            if (full.indexOf("|") == -1) {
                names.add(myLine.substring(0, myLine.indexOf(" ")));
            } else {
                names.add(full.substring(0, full.indexOf("|")));
            }
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
            int i;
            if (line[0].indexOf("|") == -1) {
                i = names.indexOf(line[0]);
            } else {
                i = names.indexOf(line[0].substring(0, line[0].indexOf("|")));
            }
			for (int n = 1; n < line.length; n++) {
                int j;
			    if (line[n].indexOf("|") == -1) {
                    j = names.indexOf(line[n]);
                } else {
                    j = names.indexOf(line[n].substring(0, line[n].indexOf("|")));
                }
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

    public static void createDAG() {
        // Create and initialize DAG G' from the query path Q and the input graph G
        // Add the vertices in the correct "levels"

        int totalNodesInDAG = 0;

        // Initialize the number of levels in the G'
        int numLevels = query.size() + 2;
        dag = new ArrayList<>(numLevels);

        // Add the source level to G'
        List<Vertex> sLevel = new ArrayList<Vertex>(1);
        Vertex s = new Vertex("Source", 0);
        sLevel.add(s);
        dag.add(sLevel);

        for (int i=0; i<query.size(); i++) {
            // Create a new list for the current level
            List<Vertex> level = new ArrayList<Vertex>();

            // Find all the associated vertices to the current query path "level"
            // Add to the "level" list.
            for (int j=0; j<names.size(); j++) {
                if (correspondence[i][j] != 0) {
                    int lvl = i + 1;
                    String fromName = names.get(j) + "_Lvl:" + lvl;
                    Vertex v = new Vertex(fromName, correspondence[i][j]);

                    // Add the edge "from Vertex to Sink"
                    double wt = INDEL_PENALTY * ((query.size()-1) - i);
                    Edge e1 = new Edge("Sink", wt);
                    v.addEdge(e1);

                    // Add the edge "from Source to Vertex"
                    wt = correspondence[i][j] + INDEL_PENALTY * i;
                    Edge e2 = new Edge(fromName, wt);
                    s.addEdge(e2);

                    level.add(v);
                    totalNodesInDAG++;
                }
            }

            // Add the completed "level" list to G'
            dag.add(level);
        }

        // Add the sink level to G'
        List<Vertex> tLevel = new ArrayList<Vertex>(1);
        Vertex t = new Vertex("Sink", 0);
        tLevel.add(t);
        dag.add(tLevel);

        // Add the edge weights
        calculateDAGWeights();
    }

	public static void calculateDAGWeights() {
        // Calculate the edge weights for all vertices in the DAG G'

        int levels = dag.size();

        // Calculate the edge weights for non-source and non-sink vertices
        for (int i=1; i<levels-1; i++) {

            // For each vertex v in the current level of G':
            for (Vertex v : dag.get(i)) {
                String fromNode = v.getNode().split("_Lvl:")[0];

                int levelsDown = 1;
                // Repeat process for however many levels below as allowed by the max number of gaps
                while ( levelsDown <= MAX_NUM_MISMATCHES_AND_GAPS+1  && i+levelsDown < levels-1 ) {
                    // Try to connect an edge to the level below:
                    List<Vertex> nextLevel = dag.get(i+levelsDown);

                    // Get each vertex nextV from the level below
                    for (Vertex nextV : nextLevel) {
                        String toNode = nextV.getNode().split("_Lvl:")[0];

                        // Find all the reachable neighbours for v
                        // Adjacent neighbour, hence no gap penalty
                        // Edge weight will be the weight(v)
                        int dist = distance[names.indexOf(fromNode)][names.indexOf(toNode)];
                        if (dist <= MAX_NUM_MISMATCHES_AND_GAPS + 1) {
                            if (dist > levelsDown) {
                                int level = i + levelsDown;
                                Edge e = new Edge(toNode+"_Lvl:"+level, nextV.weight + INDEL_PENALTY*(dist-1));
                                if (e.weight < 0) {
                                    System.out.println("Negative edge weight.");
                                }
                                v.addEdge(e);
                            } else {
                                int level = i + levelsDown;
                                Edge e = new Edge(toNode+"_Lvl:"+level, nextV.weight + INDEL_PENALTY*(levelsDown-1));
                                if (e.weight < 0) {
                                    System.out.println("Negative edge weight.");
                                }
                                v.addEdge(e);
                            }


                        }
                        // "toNode" is not reachable.
                        // Edge weight will be MAX_FLT
                        else {
//                            Edge e = new Edge(toNode+"_Lvl:"+level, MAX_FLT);
//                            v.addEdge(e);
                        }
                    }
                    levelsDown++;
                }
            }
        }
        System.out.println(dag.size());
    }

    public static void floyd() {
    	distance = new int[graphSize][graphSize];
    	path = new int[graphSize][graphSize];
		// initialize

        distance = new int[graphSize][graphSize];
        path = new int[graphSize][graphSize];

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
    
    public static void kShortestPaths() {
        // Create the template for storing the formatted DAG
    	SimpleDirectedWeightedGraph<String, DefaultWeightedEdge> fDAG = new SimpleDirectedWeightedGraph<String, DefaultWeightedEdge> (DefaultWeightedEdge.class);

        // All the vertices in the DAG can be obtained from the "source" level
        Vertex source = dag.get(0).get(0);
        fDAG.addVertex(source.fromNode);                                        // Add the source vertex

        for (Edge e : source.edges) {
            fDAG.addVertex(e.toNode);                                           // Add the rest of the vertices in DAG
            DefaultWeightedEdge e1 = fDAG.addEdge(source.fromNode,e.toNode);
            fDAG.setEdgeWeight(e1, e.weight);                                   // Add the corresponding edge weight "source -> vertex"
        }

        Vertex sink = dag.get(dag.size()-1).get(0);
        fDAG.addVertex(sink.fromNode);                                          // Add the sink vertex


        // Add the remaining edges to all non-sink/non-source vertices
        for (int level=1; level<dag.size()-1; level++) {
            for (Vertex v : dag.get(level)) {
                for (Edge e : v.edges) {
                    DefaultWeightedEdge edge = fDAG.addEdge(v.fromNode,e.toNode);
                    fDAG.setEdgeWeight(edge, e.weight);
                }
            }
        }

//        List shortest_path = DijkstraShortestPath.findPathBetween(fDAG, source.fromNode, sink.fromNode);
//        System.out.println(shortest_path);

    	KShortestPaths<String, DefaultWeightedEdge> kpath = new KShortestPaths<String, DefaultWeightedEdge>(fDAG, source.fromNode, NUM_TOP_PATHS);

        List<GraphPath<String,DefaultWeightedEdge>> output = kpath.getPaths(sink.fromNode);

        System.out.println("K-Shortest Paths:");
        for(GraphPath<String, DefaultWeightedEdge> lst : output) {
        	System.out.println(lst.getEdgeList());

            double score = 0;
            for (DefaultWeightedEdge e : lst.getEdgeList()) {
                String start = e.toString().split(" : ")[0].trim();
                String end = e.toString().split(" : ")[1].trim();
                System.out.println(start + " --> " + end + " -- score = " + fDAG.getEdgeWeight(e));
                score += fDAG.getEdgeWeight(e);
            }
            System.out.println("Score: " + score + "\n");

        }
    }
    

	
	public static void main(String [] args) {
		try {
            readGraph();
			readQuery();
            readCorrespondence();
            floyd();
            createDAG();
            kShortestPaths();
            
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
