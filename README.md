# Trajectory planification
> I carried out a little trajectory planification project. The aim was to optimize the trajectory of an object in a complex environment full of obstacles. 
I needed this kind of code for a school project and it also a good way of learning about modelling and finding the shortest path. Doing in C++ was also a good training to improve my skills in this language.
This is just a first version where a lot of things have been simplified  (I'll make a list of possible improvements and points to work on below) for an autonomous vehicule that we don't yet know right now what it will look like.


#### Tables of contents
* [Path tree](#path-tree)
* [Direct links to folders](#direct-links-to-folders)  
* [Simplifications](#simplification)  
* [Improvements to be made](#improvements_to_be_made)


## Path tree
```
Trajectory_Planification/
├── Src/
│   ├── arc_class       // define what an arc is             
│   ├── graph_class     // define what an graph is 
│   ├── obstacle_class  // define what an obstacle is  
│   └── point_class     // define what a point is 
│
└── Tests/
    ├── Test_Arc/           
    ├── Test_Djikstra/         
    ├── Test_Graph/         
    ├── Test_Obstacle/         
    ├── Test_Param/         
    └── Test_Point/         
```


## Direct links to folders 
* [Src](./Src/) : contains the codes and differents class and functions associated
* [Tests](./Tests/) : contains the differents tests for each part of the code 
    * [Test_Arc](./Tests/Test_Arc/) : contains the test of the arcs functionalities
    * [Test_Djikstra](./Tests/Test_Djikstra/) : contains the test of the search for the shortest path
    * [Test_Arc](./Tests/Test_Graph/) : contains the test of the graph caracteristics
    * [Test_Arc](./Tests/Test_Obstacle/) : contains the test of the obstacles class
    * [Test_Arc](./Tests/Test_Param/) : contains the test of the association of the function object and its string formula
    * [Test_Arc](./Tests/Test_Point/) : contains the test of the point class


## Simplifications
The object that need to be moved can only be a point or a rectangle travelling through segments with obstacles that are polygonal


## Improvements to be made
About the code restriction : 
- The arcs have already been set up to account for circular arcs and parameterized arcs.
- Overlapping occurs when using non-decimal numbers, leading to minor overlapping issues on nodes and arcs in the graph.
- Obstacles are assumed not to intersect, though they can be in contact.
- Currently, there must be enough space between obstacles to allow for padding. Future improvements could handle intersecting obstacles by merging them into larger obstacles.
- Interpolating a shape using a polygon? → This could be considered if no solution is found for computing the barycenter and area with parameterized arcs and circular arcs.
- Adjust margins and sample sizes for testing on a real track with parameterized arcs for arcs intersections, if necessary (requirements will vary depending on the scale of the problem).
What can be added : 
- Obstacles could be circular, movable, or even cease to be obstacles, with varying alert levels depending on time dynamics, making the graph of the problem constantly different.
- Adapt arcs and padding for non-polygonal shapes (requires research on how to properly achieve this). Improve padding optimization beyond simple rectangles, considering whether an object can rotate and, later, incorporating the physical properties of the environment to make the padding more realistic.
- Model elevation changes by extending distances and adding associated cost penalties. Consider modeling with adjacency matrices based on elevation to determine cost and integrate track datas.
- Incorporate the physical characteristics of the vehicle to establish physics-based constraints on obstacle traversal.
- Modify the search algorithm if arcs are no longer simple segments; introduce functions to account for arc costs, among other necessary adaptations.
