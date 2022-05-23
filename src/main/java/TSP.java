import org.apache.commons.io.IOUtils;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.alg.connectivity.ConnectivityInspector;
import org.jgrapht.alg.cycle.PatonCycleBase;
import org.jgrapht.graph.*;
import org.json.*;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Comparator;

public class TSP {

    public static void main(String[] args) {
        TSP tsp = new TSP("src/main/resources/gr120.json");
        System.out.println("Lower bound using Minimal Spanning Tree: " + tsp.getLowerBoundMinSpanTree());
        System.out.println("Lower bound using 1-Tree: " + tsp.getLowerBound1Tree());
        System.out.println("Lower bound using the doubling edges technique: " + tsp.getLowerBoundDoubleEdges());

    }

    private final Graph<Integer, DefaultWeightedEdge> graph;

    public TSP(String path) {
        JSONObject obj = null;
        try (InputStream inputStream = new FileInputStream(path)) {
            String jsonText = IOUtils.toString(inputStream, StandardCharsets.UTF_8);
            obj = new JSONObject(jsonText);
        } catch (IOException e) {
            e.printStackTrace();
        }
        graph = new DefaultUndirectedWeightedGraph<>(DefaultWeightedEdge.class);
        if (obj != null) {
            JSONArray edgeWeights = obj.getJSONArray("edge_weights");
            //Add vertices
            for (int i = 0; i < edgeWeights.length(); i++) {
                graph.addVertex(i);
            }
            //Add Edges
            for (int i = 0; i < edgeWeights.length(); i++) {
                JSONArray edgeWeightsForVertex = edgeWeights.getJSONArray(i);
                for (int j = 0; j < edgeWeightsForVertex.length(); j++) {
                    if (edgeWeightsForVertex.getInt(j) != 0) {
                        DefaultWeightedEdge e = new DefaultWeightedEdge();
                        graph.addEdge(i, j, e);
                        graph.setEdgeWeight(e, edgeWeightsForVertex.getInt(j));
                    }
                }
            }
        }
    }


    public int getLowerBoundMinSpanTree() {
        //Perform Kruskal algorithm
        Graph<Integer, DefaultWeightedEdge> minSpanTree = Kruskal(graph);
        //Get total weight and return it
        int weight = 0;
        for (DefaultWeightedEdge e : minSpanTree.edgeSet()) {
            weight += graph.getEdgeWeight(e);
        }
        return weight;
    }

    public int getLowerBound1Tree() {
        //Get vertex with the highest combined weight of two edges
        int vertexToRemove = Integer.MAX_VALUE;
        double maxEdgesWeight = 0;
        DefaultWeightedEdge[] edgesToAddLater = new DefaultWeightedEdge[2];
        for (Integer v : graph.vertexSet()) {
            //Get weight of the two minimal edges
            DefaultWeightedEdge[] minEdges = new DefaultWeightedEdge[2];
            for (DefaultWeightedEdge e : graph.edgesOf(v)) {
                if (minEdges[0] == null) {
                    minEdges[0] = e;
                    continue;
                }
                if (minEdges[1] == null) {
                    if (graph.getEdgeWeight(e) < graph.getEdgeWeight(minEdges[0])) {
                        minEdges[1] = minEdges[0];
                        minEdges[0] = e;
                    } else {
                        minEdges[1] = e;
                    }
                    continue;
                }
                if (graph.getEdgeWeight(e) < graph.getEdgeWeight(minEdges[0])) {
                    minEdges[0] = minEdges[1];
                    minEdges[0] = e;
                }
            }
            double weight = graph.getEdgeWeight(minEdges[0]) + graph.getEdgeWeight(minEdges[1]);

            if (weight > maxEdgesWeight) {
                maxEdgesWeight = weight;
                vertexToRemove = v;
                edgesToAddLater = minEdges;
            }
        }
        //Create temporary copy without the max vertex
        Graph<Integer, DefaultWeightedEdge> tempGraph = new DefaultUndirectedWeightedGraph<>(DefaultWeightedEdge.class);
        Graphs.addGraph(tempGraph, graph);
        tempGraph.removeVertex(vertexToRemove);
        //Compute minimum spanning tree of altered graph
        Graph<Integer, DefaultWeightedEdge> oneTree = Kruskal(tempGraph);
        //Add the removed vertex with the two edges with the highest weight
        oneTree.addVertex(vertexToRemove);
        oneTree.addEdge(graph.getEdgeSource(edgesToAddLater[0]), graph.getEdgeTarget(edgesToAddLater[0]), edgesToAddLater[0]);
        oneTree.addEdge(graph.getEdgeSource(edgesToAddLater[1]), graph.getEdgeTarget(edgesToAddLater[1]), edgesToAddLater[1]);
        //Get total weight and return it
        int weight = 0;
        for (DefaultWeightedEdge e : oneTree.edgeSet()) {
            weight += graph.getEdgeWeight(e);
        }
        return weight;
    }

    public int getLowerBoundDoubleEdges() {
        //Create distance matrix
        int[][] distances = new int[graph.vertexSet().size()][graph.vertexSet().size()];
        for (int i = 0; i < graph.vertexSet().size(); i++) {
            for (int j = 0; j < graph.vertexSet().size(); j++) {
                if (i == j) {
                    distances[i][j] = Integer.MAX_VALUE;
                } else {
                    distances[i][j] = (int) graph.getEdgeWeight(graph.getEdge(i, j));
                }
            }
        }

        //Go through rows and deduct the minimum
        int u = 0;
        for (int i = 0; i < distances.length; i++) {
            //Find row minimum
            int min = Integer.MAX_VALUE;
            for (int j = 0; j < distances[i].length; j++) {
                min = Math.min(distances[i][j], min);
            }
            //Deduct minimum from all values in the same row
            for (int j = 0; j < distances[i].length; j++) {
                distances[i][j] -= min;
            }
            u += min;
        }
        //Do the same again for each column
        int v = 0;
        for (int j = 0; j < distances[0].length; j++) {
            //Find column minimum
            int min = Integer.MAX_VALUE;
            for (int i = 0; i < distances.length; i++) {
                min = Math.min(distances[i][j], min);
            }
            //Deduct minimum from all values in the same column
            for (int i = 0; i < distances.length; i++) {
                distances[i][j] -= min;
            }
            v += min;
        }
        return u + v;
    }

    private Graph<Integer, DefaultWeightedEdge> Kruskal(Graph<Integer, DefaultWeightedEdge> inputGraph) {
        ArrayList<DefaultWeightedEdge> edges = new ArrayList<>(inputGraph.edgeSet());
        //Sort edges ascending by weight
        edges.sort(Comparator.comparingInt(e -> (int) inputGraph.getEdgeWeight(e)));

        Graph<Integer, DefaultWeightedEdge> minSpanTree = new DefaultUndirectedWeightedGraph<>(DefaultWeightedEdge.class);
        //Add all vertices
        for (int v : inputGraph.vertexSet()) {
            minSpanTree.addVertex(v);
        }

        for (DefaultWeightedEdge e : edges) {
            //Add the edge to the tree
            minSpanTree.addEdge(inputGraph.getEdgeSource(e), inputGraph.getEdgeTarget(e), e);
            //If this created a cycle, remove the edge again
            PatonCycleBase<Integer, DefaultWeightedEdge> cycleBase = new PatonCycleBase<>(minSpanTree);
            if (!cycleBase.getCycleBasis().getCycles().isEmpty())
                minSpanTree.removeEdge(e);

            //If the graph is connected, the algorithm is done
            ConnectivityInspector<Integer, DefaultWeightedEdge> conn = new ConnectivityInspector<>(minSpanTree);
            if (conn.isConnected())
                break;
        }
        return minSpanTree;
    }
}
