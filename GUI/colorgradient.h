//color gradient class derived from "Code - heatmaps and color gradients" by Andrew Noske
#ifndef COLORGRADIENT_H
#define COLORGRADIENT_H
#include <QVector>

/**
* \class ColorGradient
* \brief Simple class to generate heatmaps and color gradients for sound pressure and phase visualization.
*
* This class is derived from "Code - heatmaps and color gradients" by Andrew Noske (https://www.andrewnoske.com/wiki/Code_-_heatmaps_and_color_gradients)
*/
class ColorGradient
{
private:
    /**
    * \struct ColorPoint
    * \brief Defines a color point in the color gradient.
    */
    struct ColorPoint
    {
        float r,g,b;
        float val;
        ColorPoint(){}

        /**
        * \brief Constructs a new color point with the given color and value.
        * \param red The red component of the color.
        * \param green The green component of the color.
        * \param blue The blue component of the color.
        * \param value The value of the color point in the gradient that has to be visualized.
        */
        explicit ColorPoint(float red, float green, float blue, float value)
        {
            this->r = red;
            this->g = green;
            this->b = blue;
            this->val = value;
        }
    };
    QVector<ColorPoint> color;

public:
    ColorGradient(){ defaultGradient();}

  //-- Inserts a new color point into its correct position:
    void addColorPoint(float red, float green, float blue, float value)
        {
        for(int i=0; i<color.size(); i++)
        {
            if(value < color[i].val)
            {
                color.insert(color.begin()+i, ColorPoint(red,green,blue, value));
                return;
            }
        }
        color.push_back(ColorPoint(red,green,blue, value));
    }
    void clear(){color.clear();}

    QVector<ColorPoint> getColor(){return color;}

    /**
    * \brief Sets the default gradient for sound pressure visualization.
    */
    void defaultGradient()
    {
        color.clear();
        color.push_back(ColorPoint(0, 0, 1,   0.0f));      // blau.
        color.push_back(ColorPoint(0, 1, 1,   0.25f));     // cyan.
        color.push_back(ColorPoint(0, 1, 0,   0.5f));      // gruen.
        color.push_back(ColorPoint(1, 1, 0,   0.75f));     // gelb.
        color.push_back(ColorPoint(1, 0, 0,   1.0f));      // rot.
    }

    /**
    * \brief Sets the gradient for phase visualization.
    */
    void phaseGradient()
    {
        color.clear();
        color.push_back(ColorPoint(0, 0, 1,   0.0));      // blau.
        color.push_back(ColorPoint(0, 1, 1,   0.25f * global::PI));     // cyan.
        color.push_back(ColorPoint(0, 1, 0,   0.5f * global::PI));      // gruen.
        color.push_back(ColorPoint(1, 1, 0,   0.75f * global::PI));     // gelb.
        color.push_back(ColorPoint(1, 0, 0,   1.0f * global::PI));      // rot.
        color.push_back(ColorPoint(1, 1, 0,   1.25f * global::PI));     // gelb.
        color.push_back(ColorPoint(0, 1, 0,   1.5f * global::PI));      // gruen.
        color.push_back(ColorPoint(0, 1, 1,   1.75f * global::PI));     // cyan.
        color.push_back(ColorPoint(0, 0, 1,   2.0f * global::PI));      // blau.
    }

  void getColorAtValue(const float value, float &red, float &green, float &blue)
  {
    if(color.size()==0)
      return;

    for(int i=0; i<color.size(); i++)
    {
        ColorPoint &currC = color[i];
        if(value < currC.val)
        {
            ColorPoint &prevC  = color[ std::max(0,i-1) ];
            float valueDiff    = (prevC.val - currC.val);
            float fractBetween = (valueDiff==0) ? 0 : (value - currC.val) / valueDiff;
            red   = (prevC.r - currC.r)*fractBetween + currC.r;
            green = (prevC.g - currC.g)*fractBetween + currC.g;
            blue  = (prevC.b - currC.b)*fractBetween + currC.b;
            return;
        }
    }
    red   = color.back().r;
    green = color.back().g;
    blue  = color.back().b;
    return;
  }
};
#endif // COLORGRADIENT_H

