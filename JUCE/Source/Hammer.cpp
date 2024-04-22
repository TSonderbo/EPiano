/*
  ==============================================================================

    Hammer.cpp
    Created: 8 Mar 2024 10:34:19am
    Author:  Sonderbo

  ==============================================================================
*/

#include "Hammer.h"

Hammer::Hammer()
{

}

void Hammer::prepareToPlay(double sampleRate, std::vector<float> contact, int N)
{
    k = 1.0f / (sampleRate * config::oversampling);

    this->contact = contact;
    this->N = N;
    //TODO set stiffness based on hammer/tine number
}

void Hammer::beginHammer(float velocity)
{
    active = true;

    u = uIn + k * velocity;
    uPrev = uIn ;

    xPrev = 0;
    x = 0;

}

void Hammer::stopHammer()
{
    active = false;
}

float Hammer::calculateForce(std::vector<float>::iterator tine_u)
{
    if (active == false)
    {
        return 0.0f;
    }

    //Calculate compression
    x = 0.0f;
    for (int i = 0; i <= N; i++)
    {
        x += contact[i] * (u - tine_u[i]);
    }

    //Calculate Force
    float force = 0.0f;
    if (x <= 0.0f)
    {
        force = 0.0f;
        x = 0;
    }
    else
    {
        float v = (1 / k) * (x - xPrev);
        force = K * powf(x, alpha) * (1 + mu * v);
    }
    //Update states
    float uNext = 2 * u - uPrev - (force * k * k);
    uPrev = u;
    u = uNext;

    xPrev = x;

    if (u <= uIn)
    {
        stopHammer();
    }

    return force;
}

float Hammer::getMass()
{
    return m;
}