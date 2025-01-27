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
    vector<Arc*> obstacle_arcs = {&s1, &s2, &s3, &s4};
    Obstacle obstacle(obstacle_arcs);
    // Create the start and goal points
    Point start(-3, 1);
    Point goal(2, 0);
    vector<Obstacle> obstacles = {obstacle};
    Graph graph(obstacles, start, goal);
    // Export the graph
    string graph_filename = "graph.txt";
    graph.exporte(graph_filename);
    cout<<"Graph exported to "<<graph_filename<<endl;
    vector<vector<double>> adj_matrix = graph.build_adjacency_matrix();
    // show the adjacency matrix
    cout<<"Adjacency matrix: "<<endl;
    for (int i = 0; i < adj_matrix.size(); ++i){
        for (int j = 0; j < adj_matrix[i].size(); ++j){
            cout<<adj_matrix[i][j]<<" ";
        }
        cout<<endl;
    }
    // find the shortest path
    graph.find_shortest_path_and_exporte("shortest_path.txt");


    // Second test with two obstacles
    // Create an obstacle
    Point p5(0, -1);
    Point p6(1, -1);
    Point p7(0, -2);
    Point p8(1, -2);
    Segment s5(p5, p6);
    Segment s6(p6, p7);
    Segment s7(p7, p8);
    Segment s8(p8, p5);
    vector<Arc*> obstacle_arcs_2 = {&s5, &s6, &s7, &s8};
    Obstacle obstacle1(obstacle_arcs);
    Obstacle obstacle2(obstacle_arcs_2);
    vector<Obstacle> obstacles2 = {obstacle1, obstacle2};
    Point start2(-1, 0.5);
    Point goal2(2, -0.5);
    Graph graph2(obstacles2, start2, goal2);
    // Export the graph
    string graph_filename2 = "graph2.txt";
    graph2.exporte(graph_filename2);
    cout<<"Graph exported to "<<graph_filename2<<endl;
    // find the shortest path
    graph2.find_shortest_path_and_exporte("shortest_path2.txt");


    // third one with 3 obstacles
    // Create an obstacle
    Point p9(2, 0);
    Point p10(3, 0);
    Point p11(2, -1);
    Point p12(3, -1);
    Segment s9(p9, p10);
    Segment s10(p10, p11);
    Segment s11(p11, p12);
    Segment s12(p12, p9);
    vector<Arc*> obstacle_arcs_3 = {&s9, &s10, &s11, &s12};
    Obstacle obstacle3(obstacle_arcs_3);
    vector<Obstacle> obstacles3 = {obstacle1, obstacle2, obstacle3};
    Point start3(-1, 0.5);
    Point goal3(4, -0.5);
    Graph graph3(obstacles3, start3, goal3);
    // Export the graph
    string graph_filename3 = "graph3.txt";
    graph3.exporte(graph_filename3);
    cout<<"Graph exported to "<<graph_filename3<<endl;
    // find the shortest path
    graph3.find_shortest_path_and_exporte("shortest_path3.txt");


    // testing while having two obstacles with side by side
    // Create an obstacle
    Point p13(0, 0);
    Point p14(1, 0);
    Point p15(1, -1);
    Point p16(0, -1);
    Segment s13(p13, p14);
    Segment s14(p14, p15);
    Segment s15(p15, p16);
    Segment s16(p16, p13);
    vector<Arc*> obstacle_arcs_4 = {&s13, &s14, &s15, &s16};
    Obstacle obstacle4(obstacle_arcs_4);
    vector<Obstacle> obstacles4 = {obstacle1, obstacle4};
    Point start4(-1, 0.5);
    Point goal4(2, -0.5);
    Graph graph4(obstacles4, start4, goal4);
    // Export the graph
    string graph_filename4 = "graph4.txt";
    graph4.exporte(graph_filename4);
    cout<<"Graph exported to "<<graph_filename4<<endl;
    // find the shortest path
    graph4.find_shortest_path_and_exporte("shortest_path4.txt");


    return 0;
}
