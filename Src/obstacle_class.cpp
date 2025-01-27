#include "obstacle_class.hpp"



//==========================================================================
// class Obstacle : polygonal obstacle
//==========================================================================

// member functions
//==========================================================================
// default constructor
Obstacle::Obstacle()
{

}
// obstacle with the list of sides and their properties
Obstacle::Obstacle(const vector<Arc*>& sides) : arcs(sides)
{

}
// apply padding to the obstacle
void Obstacle::apply_padding(double padding_distance)
{
    Point center = barycenter(*this);
    for (auto& arc : arcs) {
        for (auto& point : arc->points) {
            // expand each point away from the center of the obstacle
            Point direction = point - center;
            double magnitude = sqrt(direction.x * direction.x + direction.y * direction.y);
            if (magnitude > 0) {
                direction.x /= magnitude;
                direction.y /= magnitude;
                point.x += padding_distance * direction.x;
                point.y += padding_distance * direction.y;
            }
        }
    }
}
// export the obstacle to a file
void Obstacle::exporte(ostream& out) 
{
    for (Arc* arc : arcs) {
        arc->exporte(out);
    }
}

// external functions
//==========================================================================
// check if a point is inside the obstacle
bool is_inside(const Point& P, const Obstacle& obstacle)
{
    int nb_intersections = 0;

    for (Arc* arc : obstacle.arcs){
        if (const ArcCircle* arc_circle = dynamic_cast<const ArcCircle*>(arc)){
            // Check if P is inside the sector of the circle
            double dx = P.x - arc_circle->center.x;
            double dy = P.y - arc_circle->center.y;
            double distance_squared = dx * dx + dy * dy;
            double radius_squared = arc_circle->radius * arc_circle->radius;

            // P is within the circle (i.e., distance from center is less than the radius)
            if (distance_squared <= radius_squared){
                // Angle between point P and center
                double angle = atan2(dy, dx);
                if (angle < 0){
                    angle += 2 * M_PI; // Ensure angle is positive
                }
                // Check if the angle is within the arc's angular range
                if (arc_circle->angle1 <= angle && angle <= arc_circle->angle2){
                    nb_intersections++;  // Point P is within the arc's sector
                }
            }
        }
        else if (const ArcParameter* arc_param = dynamic_cast<const ArcParameter*>(arc)){
            // Approximate the point using the parametric function
            // You can sample the parameter space and check if the point matches
            double step = (arc_param->tf - arc_param->t0) / arc_param->nbt;
            for (double t = arc_param->t0; t <= arc_param->tf; t += step){
                Point sampled_point = arc_param->parametrisation(t);
                // Check if the sampled point is close to P (within a tolerance)
                if (fabs(sampled_point.x - P.x) < 1e-6 && fabs(sampled_point.y - P.y) < 1e-6){
                    nb_intersections++;  // Point P is on the parametric arc
                    break;
                }
            }
        }
        else if (const Segment* seg = dynamic_cast<const Segment*>(arc)){
            // for a segment
            Point p1 = seg->points[0];
            Point p2 = seg->points[1];
            // Check if the ray from P intersects the segment (ray-casting)
            if ((p1.y > P.y) != (p2.y > P.y) && P.x < (p2.x - p1.x) * (P.y - p1.y) / (p2.y - p1.y) + p1.x){
                nb_intersections++;
            }
        }
    }

    // Point is inside if the number of intersections is odd
    return nb_intersections % 2 != 0;
}
// Check if the given arc intersects with any arc in the obstacle
bool does_intersect(const Arc& arc, const Obstacle& obstacle)
{
    for (Arc* obstacle_arc : obstacle.arcs){
        if (const Segment* seg1 = dynamic_cast<const Segment*>(&arc)){
            if (const Segment* seg2 = dynamic_cast<const Segment*>(obstacle_arc)){
                // Check if two segments intersect
                if (do_segments_intersect(*seg1, *seg2)){
                    return true;
                }
            }
            else if (const ArcCircle* arc_circle = dynamic_cast<const ArcCircle*>(obstacle_arc)){
                // Check if the segment intersects with the arc circle
                if (do_segment_intersect_arc_circle(*seg1, *arc_circle)){
                    return true;
                }
            }
            else if (const ArcParameter* param_arc = dynamic_cast<const ArcParameter*>(obstacle_arc)){
                // Check if the segment intersects with the parametric arc
                if (do_segment_intersect_parametric(*seg1, *param_arc)){
                    return true;
                }
            }
        }
        else if (const ArcCircle* arc_circle = dynamic_cast<const ArcCircle*>(&arc)){
            if (const Segment* seg2 = dynamic_cast<const Segment*>(obstacle_arc)){
                // Check if the arc circle intersects with the segment
                if (do_segment_intersect_arc_circle(*seg2, *arc_circle)){
                    return true;
                }
            }
            else if (const ArcCircle* arc_circle2 = dynamic_cast<const ArcCircle*>(obstacle_arc)){
                // Check if two arc circles intersect
                if (do_arc_circle_intersect(*arc_circle, *arc_circle2)){
                    return true;
                }
            }
            else if (const ArcParameter* param_arc = dynamic_cast<const ArcParameter*>(obstacle_arc)){
                // Check if the arc circle intersects with the parametric arc
                if (do_arc_circle_intersect_parametric(*arc_circle, *param_arc)){
                    return true;
                }
            }
        }
        else if (const ArcParameter* param_arc = dynamic_cast<const ArcParameter*>(&arc)){
            if (const Segment* seg2 = dynamic_cast<const Segment*>(obstacle_arc)){
                // Check if the parametric arc intersects with the segment
                if (do_segment_intersect_parametric(*seg2, *param_arc)){
                    return true;
                }
            }
            else if (const ArcCircle* arc_circle = dynamic_cast<const ArcCircle*>(obstacle_arc)){
                // Check if the parametric arc intersects with the arc circle
                if (do_arc_circle_intersect_parametric(*arc_circle, *param_arc)){
                    return true;
                }
            }
            else if (const ArcParameter* param_arc2 = dynamic_cast<const ArcParameter*>(obstacle_arc)){
                // Check if two parametric arcs intersect
                if (do_parametric_intersect(*param_arc, *param_arc2)) {
                    return true;
                }
            }
        }
    }
    return false;
}
// barycentre of the obstacle (only for polygonal obstacles for now)
Point barycenter(const Obstacle& obstacle)
{
    for (Arc* arc : obstacle.arcs){
        const Segment* seg = dynamic_cast<const Segment*>(arc);
        if (!seg){
            cout << "The obstacle is not a polygon" << endl;
            exit(-1);
        }
    }
    // getting all the obstacle points (vertices)
    vector<Point> vertices;
    for (Arc* arc : obstacle.arcs){
        const Segment* seg = dynamic_cast<const Segment*>(arc);
        vertices.push_back(seg->points[0]);
    }

    Point centroid(0,0);
    for (int i = 0; i < vertices.size(); i++){
        centroid.x += vertices[i].x;
        centroid.y += vertices[i].y;
    }
    centroid.x /= vertices.size();
    centroid.y /= vertices.size();
    return centroid;
}
// obstacle area  (only for polygonal obstacles for now)
double area(const Obstacle& obstacle)
{
    if (obstacle.arcs.size() < 3){
        cout << "The obstacle is not a polygon" << endl;
        exit(-1);
    }

    // getting all the obstacle points (vertices)
    vector<Point> vertices;
    for (Arc* arc : obstacle.arcs){
        const Segment* seg = dynamic_cast<const Segment*>(arc);
        vertices.push_back(seg->points[0]);
    }

    const Point& G = barycenter(obstacle); // handle if the obstacle is not a polygon
    double total_area = 0.0;
    for (int i=0; i<vertices.size(); i++){
        const Point& A = G;
        const Point& B = vertices[i];
        const Point& C = vertices[(i+1)%vertices.size()];
        total_area += surface(A, B, C);
    }
    return total_area;
}
