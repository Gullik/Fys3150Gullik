#include <circle.h>
#include <cmath>


Circle::Circle(double radius)
{
    double m_radius = radius;

}

double Circle::Area()
{
    return m_radius*m_radius*3.14;
}
