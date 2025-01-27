#include "../../Src/arc_class.hpp"


const double pi=4*atan(1.);
Param_Fun fsin([](double t){return Point(-sin(0.5*pi*(1+t)),t);},"x(t) = -sin(0.5*pi*(1+t)), y(t) = t");


int main()
{
    Point A={0.,1.}; 
    Point B(1.,0.);
    // Segment class validation
    Segment S={A,B};
    cout<<S<<" length="<<S.length()<<endl;
    cout<<"S.reverse() -> "<<S.reverse()<<endl;
    // ArcCircle class validation
    ArcCircle AC(Point(0.,0.),1.,0,-pi/2);
    cout<<AC<<" length="<<AC.length()<<endl;
    cout<<"AC.reverse() -> "<<AC.reverse()<<endl;
    // ArcParameter class validation
    ArcParameter AP(fsin,-1,1, 10);
    cout<<AP<<" length="<<AP.length()<<endl;
    cout<<"AP.reverse() -> "<<AP.reverse()<<endl;


    // intersection tests
    // creating points
    Point p1(0, 0);
    Point p2(1, 1);
    Point p3(2, 0);
    Point p4(1, -1);
    Point p5(1, 0);
    Point center(0, 0);
    
    // create Segment: from p1 to p2
    Segment seg1(p1, p2);
    // create Segment: from p1 to p5
    Segment seg2(p1, p5);
    // create ArcCircle: center at (0, 0), radius 1, angle range from 0 to pi/2
    ArcCircle arc_circle1(center, 1, 0, M_PI / 2);
    // create parametric arc with a simple parametrization function
    Param_Fun param_fun1([](double t){return Point(t, -1+t*t);}, "x(t) = t, y(t) = -1+t*t"); // Example function x = t, y = -1+t^2
    ArcParameter param_arc1(param_fun1, 0, 1, 10);
    // create another parametric arc
    Param_Fun param_fun2([](double t){return Point(t, -t*t);}, "x(t) = t, y(t) = -t*t"); // Example function x = t, y = -t^2
    ArcParameter param_arc2(param_fun2, 0, 1, 10);
    // create another ArcCircle: from center (0, 0), radius 1, angle range from pi/2 to pi
    ArcCircle arc_circle2(center, 1, M_PI / 2, M_PI);

    // check if two arc circles intersect
    cout<<"\nArcCircle1 intersects with ArcCircle2: "<<(do_arc_circle_intersect(arc_circle1, arc_circle2) ? "Yes" : "No")<<endl;
    // check if a segment intersects with an arc circle
    cout<<"\nSegment1 intersects with ArcCircle1: "<<(do_segment_intersect_arc_circle(seg1, arc_circle1) ? "Yes" : "No")<<endl;
    // check if an arc circle intersects with a parametric arc
    cout<<"\nArcCircle1 intersects with Parametric Arc1: "<<(do_arc_circle_intersect_parametric(arc_circle1, param_arc1) ? "Yes" : "No")<<endl;
    // check if two parametric arcs intersect
    cout<<"\nParametric Arc1 intersects with Parametric Arc2: "<<(do_parametric_intersect(param_arc1, param_arc2) ? "Yes" : "No")<<endl;
    // check if a segment intersects with an arc circle
    cout<<"\nSegment2 intersects with ArcCircle1: "<<(do_segment_intersect_arc_circle(seg2, arc_circle1) ? "Yes" : "No")<<endl;
    

    // additional tests on segment intersections
    Point P1(0, 0);
    Point P2(1, 1);
    Point P3(2, 2);
    Segment S1(P1, P3);
    Segment S2(P2, P3);
    Segment S3(P1, P2);
    Point P4(1, 2);
    Segment S4(P4, P2);
    cout<<"\nSegment S1: "<<S1<<endl;
    cout<<"Segment S2: "<<S2<<endl;
    cout<<"Segment S3: "<<S3<<endl;
    cout<<"Segment S4: "<<S4<<endl;
    cout<<"\nSegment S1 overlaps with Segment S2: "<<(are_overlapping(S1, S2) ? "Yes" : "No")<<endl;
    cout<<"\nSegment S2 intersects with Segment S3: "<<(do_segments_intersect(S2, S3) ? "Yes" : "No")<<endl;
    cout<<"\nSegment S2 overlaps with Segment S3: "<<(are_overlapping(S2, S3) ? "Yes" : "No")<<endl;
    cout<<"\nSegment S1 intersects with Segment S4: "<<(do_segments_intersect(S1, S4) ? "Yes" : "No")<<endl;
    cout<<"\nSegment S3 overlaps with Segment S4: "<<(are_overlapping(S3, S4) ? "Yes" : "No")<<endl;


    return 0;
}