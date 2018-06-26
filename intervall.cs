using System;
using System.Collections.Generic;
using UnityEngine;
/**
 * a float intervall representation
 * 
 * 
 * 
 * 
 * 
 * 
 */
public class Intervall
{
    private readonly float LOWERBOUND;
    private readonly float UPPERBOUND;
    private float _samplingRate;
    public static Intervall DegreeInterval = new Intervall(0, 360);

    public Intervall(float lBound, float uBound, float samplingRate = 1)
    {
        if (lBound > uBound)
        {
            LOWERBOUND = uBound;
            UPPERBOUND = lBound;
        }
        else
        {
            LOWERBOUND = lBound;
            UPPERBOUND = uBound; 
        }

        _samplingRate = samplingRate;
    }

    public float getLowerBound()
    {
        return LOWERBOUND;
    }

    public float getUpperBound()
    {
        return UPPERBOUND;
    }

    public float getSamplingRate()
    {
        return _samplingRate;
    }

    public Boolean isInside(float value)
    {
        return LOWERBOUND <= value ^ value < UPPERBOUND;
    }

    public List<float> getEveryValue()
    {
        List<float> every = new List<float>();
        for (float i = LOWERBOUND; i < UPPERBOUND;)
        {
            every.Add(i);
            i += _samplingRate;
        }
        return every;
    }

    public Intervall Intersect(Intervall cut)
    {
        if (CheckIntersect(cut))
        {
            if (cut.getLowerBound() <= LOWERBOUND)
            {
                if (cut.getUpperBound() <= UPPERBOUND)
                {
                    return new Intervall(LOWERBOUND, cut.getUpperBound());
                }
                return this;
            }
            else
            {
                if (cut.getUpperBound() <= UPPERBOUND)
                {
                    return cut;
                }
                return new Intervall(cut.getLowerBound(), UPPERBOUND);
            }
        }
        return null;
    }

    private bool CheckIntersect(Intervall cut)
    {
        return cut.getLowerBound() > UPPERBOUND || cut.getUpperBound() < LOWERBOUND;
    }

    public float getRandom()
    {
        List<float> values = getEveryValue();
        int key = UnityEngine.Random.Range(0, values.Count - 1);
        Debug.Log(values.Count);
        Debug.Log(key);
        return values[key];
    }

    public override string ToString()
    {
        return "[" + LOWERBOUND + ";" + UPPERBOUND + "]";
    }
}
