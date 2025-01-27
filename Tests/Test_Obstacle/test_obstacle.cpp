#include "../../Src/obstacle_class.hpp"


int main()
{
    // Test polygon with segments
    vector<Arc*> sides1;
    // Create some segments (as an example)
    Segment s1(Point{0, 0}, Point{4, 0});
    Segment s2(Point{4, 0}, Point{4, 6});
    Segment s3(Point{4, 6}, Point{0, 6});
    Segment s4(Point{0, 6}, Point{0, 0});
    // Add segments to the obstacle
    sides1.push_back(&s1);
    sides1.push_back(&s2);
    sides1.push_back(&s3);
    sides1.push_back(&s4);
    // Create obstacle
    Obstacle obstacle1(sides1);
    // Test the area and barycenter
    Point centroid = barycenter(obstacle1);
    double obstacle_area = area(obstacle1);
    cout<<"Obstacle1 Arcs:"<<endl;
    for (Arc* arc : obstacle1.arcs){
        cout<<*arc<<endl;
    }
    cout<<"Obstacle1 Barycenter: ("<<centroid.x<<", "<<centroid.y<<")"<<endl;
    cout<<"Obstacle1 Area: "<<obstacle_area<<endl;


    // Test polygon with segments
    vector<Arc*> sides2;
    // Create some segments (as an example)
    Segment s5(Point{-1, -1}, Point{4, 0});
    Segment s6(Point{4, 0}, Point{5, 5});
    Segment s7(Point{5, 5}, Point{0, 4});
    Segment s8(Point{0, 4}, Point{-1, -1});
    // Add segments to the obstacle
    sides2.push_back(&s5);
    sides2.push_back(&s6);
    sides2.push_back(&s7);
    sides2.push_back(&s8);
    // Create obstacle
    Obstacle obstacle2(sides2);
    // Test the area and barycenter
    Point centroid2 = barycenter(obstacle2);
    double obstacle_area2 = area(obstacle2);
    fstream file = fstream("obstacle2.txt", ios::out);
    obstacle2.exporte(file);
    cout<<"Obstacle2 Arcs:"<<endl;
    for (Arc* arc : obstacle2.arcs){
        cout<<*arc<<endl;
    }
    cout<<"Obstacle2 Barycenter: ("<<centroid2.x<<", "<<centroid2.y<<")"<<endl;
    cout<<"Obstacle2 Area: "<<obstacle_area2<<endl;
    

    // Test to see if the intersection of a segment with an obstacle is detected
    Segment s9(Point{-1, 0}, Point{5, 6});
    bool intersect1 = does_intersect(s9, obstacle1);
    cout<<"Does segment9 intersect with obstacle1: "<<intersect1<<endl; // intersection
    Segment s10(Point{5, -1}, Point{5, 6});
    bool intersect2 = does_intersect(s10, obstacle1);
    cout<<"Does segment10 intersect with obstacle1: "<<intersect2<<endl; // no intersection
    Segment s11(Point{4, -1}, Point{4, 7});
    bool intersect3 = does_intersect(s11, obstacle1);
    cout<<"Does segment11 intersect with obstacle1: "<<intersect3<<endl; // there is a point of intersection on the border 
    Segment s12(Point{4, -1}, Point{4, 0});
    bool intersect4 = does_intersect(s12, obstacle1);
    cout<<"Does segment12 intersect with obstacle1: "<<intersect4<<endl; // there is no intersection if there only is one point of intersection and it is a corner of the obstacle


    return 0;
}