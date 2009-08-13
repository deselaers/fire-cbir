/**
 * \file colorhsv.h
 * \author David Rybach <rybach@i6.informatik.rwth-aachen.de>
 * \date 2004-03-13
 */
#ifndef COLORHSV_H
#define COLORHSV_H

#include <vector>

/**
 * HSV-Pixel
 */
class ColorHSV
{
public:
    /**
     * Create HSV-Pixel from RGB-Pixel
     */
    ColorHSV(double red, double green, double blue) { convertRGBtoHSV(red, green, blue); }
    /**
     * Create empty HSV-Pixel
     */
    ColorHSV() : h(0), s(0), v(0) {}
    
    /**
     * Convert RGB-Pixel to HSV-Pixel
     */
    ColorHSV& operator=(const ColorHSV &rhs);
    
    /**
     * Calculate complex values for HSV-Values
     */
    void complexPixel(std::vector<double> &pixel) const;
    
    /**
     * Hue
     */
    double Hue() const { return h; }
    /**
     * Saturation
     */
    double Saturation() const { return s; }
    /**
     * Value
     */
    double Value() const { return v; }
    
private:
    /**
     * Calcualte HSV-Values from RGB-Values
     */
    void convertRGBtoHSV(double red, double green, double blue);
    double h, s, v;
};

#endif
