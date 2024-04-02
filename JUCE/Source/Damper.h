/*
  ==============================================================================

    Damper.h
    Created: 19 Mar 2024 3:00:26pm
    Author:  Sonderbo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

class Damper
{
public:

    void prepareToPlay(double sampleRate, int N);
    void calculateDamping(const float* tine_u);

private:

    float K = 1.5 * pow(10, 11); //Damper stiffness/Spring constant
    float alpha = 2.5; //Local geometry of impact coefficient - typically in range[1.5, 3.5]
    float M = 1.1 * pow(10, -2); //Damper mass
    float b = -0.00001f; //Damper location

    //States
    std::vector<std::vector<float>> psiStates; //3 by N+1 state matrix
    std::vector<float*> psi; //uStates pointer
};