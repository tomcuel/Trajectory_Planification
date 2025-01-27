#include "../../Src/point_class.hpp"


int main()
{
    // Test de la classe Point
    Point P1;
    cout<<"\nP1="<<P1<<endl;

    Point P2(6,10);
    cout<<"P2="<<P2<<endl;
    cout<<"norme P2="<<norme(P2)<<endl;
    Point P3 = P2 + 2 + P1;
    cout<<"P3="<<P3<<endl;
    cout<<"\n";


    return 0;
}