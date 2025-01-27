#include "../../Src/graph_class.hpp"


int main() 
{
    // Create an obstacle
    Point p1(0, 0);
    Point p2(1, 0);
    Point p3(1, 1);
    Point p4(0, 1);
    Segment s1(p1, p2);
    Segment s2(p2, p3);
    Segment s3(p3, p4);
    Segment s4(p4, p1);
    vector<Arc*> obstacle_arcs_1 = {&s1, &s2, &s3, &s4};
    Obstacle obstacle1(obstacle_arcs_1);


    // Create a new obstacle
    Point p5(2, 0);
    Point p6(3, 0);
    Point p7(3, 1);
    Point p8(2, 1);
    Segment s5(p5, p6);
    Segment s6(p6, p7);
    Segment s7(p7, p8);
    Segment s8(p8, p5);
    vector<Arc*> obstacle_arcs_2 = {&s5, &s6, &s7, &s8};
    Obstacle obstacle2(obstacle_arcs_2);


    // Create a graph with the obstacle
    Point start1(2, 2);
    Point goal1(2, 3);
    vector<Obstacle> obstacles = {obstacle1};
    Graph graph1(obstacles, start1, goal1);
    string graph_filename = "graph_output.txt";
    graph1.exporte(graph_filename);
    cout << "Graph exported to " << graph_filename << endl;


    // Add a new obstacle to the graph
    graph1.add_obstacle(obstacle2);
    // Export the updated graph
    string updated_graph_filename = "updated_graph_output.txt";
    graph1.exporte(updated_graph_filename);
    cout << "Updated graph exported to " << updated_graph_filename << endl;


    // Try the padding function to see if there is an extension of the obstacle
    // Create an obstacle that will be the moving piece
    Point p9(-3, 0);
    Point p10(-2, 0);
    Point p11(-3, 1);
    Point p12(-2, 1);
    Segment s9(p9, p10);
    Segment s10(p10, p11);
    Segment s11(p11, p12);
    Segment s12(p12, p9);
    vector<Arc*> moving_piece_arcs = {&s9, &s10, &s11, &s12};
    Obstacle moving_piece(moving_piece_arcs);
    // Create a graph with the moving piece
    Point goal = Point(3, 0.5);
    Graph graph2(obstacles, goal, moving_piece);
    graph2.apply_padding_to_obstacles();
    string padded_graph_filename = "padded_graph_output.txt";
    graph2.exporte(padded_graph_filename);
    cout << "Padded graph exported to " << padded_graph_filename << endl;


    return 0;
}
