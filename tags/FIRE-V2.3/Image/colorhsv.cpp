/**
 * \file colorhsv.cc
 * \author David Rybach <rybach@i6.informatik.rwth-aachen.de>
 * \date 2004-03-13
 */
#include <vector>
#include <cmath>
#include <iostream>
#include "colorhsv.hpp"
#include <math.h>

ColorHSV& ColorHSV::operator=(const ColorHSV &rhs)
{
    h = rhs.Hue();
    s = rhs.Saturation();
    v = rhs.Value();
    return *this;
}

/**
 * Convert RGB-Value to HSV
 * Code from Thomas Deselaers, who got it from http://www.cs.rit.edu/~ncs/color/t_convert.html
 */
void ColorHSV::convertRGBtoHSV(double red, double green, double blue)
{
    double min, max, delta;

    min=red;
    if(green<min) min=green;
    if(blue<min) min=blue;

    max=red;
    if(green>max) max=green;
    if(blue>max) max=blue;

    v = max;
    delta = max - min;

    if( max != 0 )
        s = delta / max;        // s
    else
    {
        s = 0;
        h = 0;
        return;
    }

    if(delta==0) h=0;
    else
    {
        if( red == max )
        {
            if(delta==0) std::cerr << "FEHLER" << std::endl;
            h = ( green - blue ) / delta;        // between yellow & magenta
        }
        else if( green == max )
        {
            if(delta==0) std::cerr << "FEHLER" << std::endl;
            h = 2 + ( blue - red ) / delta;    // between cyan & yellow
        }
        else
        {
            if(delta==0) std::cerr << "FEHLER" << std::endl;
            h = 4 + ( red - green ) / delta;    // between magenta & cyan
        }
        h *= 60;                // degrees
    }
    if( h < 0 )
        h += 360;
}

void ColorHSV::complexPixel(std::vector<double> &pixel) const
{
    pixel[0] = s * cos(h);
    pixel[1] = s * sin(h);
    pixel[2] = v;
    // std::cout << "complex: " << h << "," << s << "," << v << " -> " << ret(0) << "," << ret(1) << "," << ret(2) << std::endl;
}
