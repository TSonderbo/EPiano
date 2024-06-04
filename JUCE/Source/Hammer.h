/*
  ==============================================================================

    Hammer.h
    Created: 8 Mar 2024 10:34:19am
    Author:  Sonderbo

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "Configuration.h"

class Hammer
{
public:

    Hammer();
    //==============================================================================
    void prepareToPlay(double k, const std::vector<float> contact, int N);
    void beginHammer(float velocity);
    void stopHammer();
    float calculateForce(std::vector<float>::iterator tine_u);
    float getMass();
private:

    float k; //Sampling period
    int N;

    float K = 1.5f * pow(10, 12); //Hammer stiffness/Spring constant
    float alpha = 2.8f; //Local geometry of impact coefficient - typically in range[1.5, 3.5]
    float mu = 0.6f;
    float m = 1.1f * pow(10, -2); //Hammer mass

    const float uIn = -1.0f * pow(10, -4); //Initial hammer position
    float u;
    float uPrev;

    float x;
    float xPrev;

    bool active;
    bool hasMadeContact;

    std::vector<float> contact;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Hammer)
};