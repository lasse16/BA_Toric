using System;
using System.Collections.Generic;
using UnityEngine;
/**
 * a float intervall representation with a lower and upper bound and a optional smapling rate
 * The sampling rate means how many floats are going to be read out of this interval
 * (e.g sampling rate = 2 => add two on the lower bound , then 2 on that value,
 * then two again ... repeat , till you are at or above the upper bound)
 * 
 * @author Lasse Haffke 
 * @version 29.06.18
 */
public class Intervall
{
    private readonly float LOWERBOUND;
    private readonly float UPPERBOUND;
    private float _samplingRate;

    //static interval for the degrees in a circle
    public static Intervall DegreeInterval = new Intervall(0, 360, 1);

    //TODO better NaN handling
    public Intervall(float lBound, float uBound, float samplingRate = 0.05f)
    {
        //ensure lower bound is smaller then the upper bound
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

    /**
     * return the lower bound
     * @return the lower bound
     */
    public float getLowerBound()
    {
        return LOWERBOUND;
    }

    /**
    * return the upper bound
    * @return the upper bound
    */
    public float getUpperBound()
    {
        return UPPERBOUND;
    }

    /**
    * return the sampling rate 
    * @return the sampling rate
    */
    public float getSamplingRate()
    {
        return _samplingRate;
    }

    /**
     * returns true if the value is inside the interval(including lower bound, excluding upper bound)
     * @param float value the value to check
     * @return true if the value is inside
     */
    public Boolean IsInside(float value)
    {
        return LOWERBOUND <= value && value <= UPPERBOUND;
    }

    /**
     * Returns a list of every value inside
     * (e.g sampling rate = 2 => add two on the lower bound , then 2 on that value,
     * then two again ... repeat , till you are at or above the upper bound)
     * 
     * @return a list of float values calculated by the process described above
     */
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

    /**
     * returns the intersection of two intervals
     * @param cut the first interval
     * @return the resulting interval 
     */
    public Intervall Intersect(Intervall cut)
    {
        if (DoNotIntersect(cut))
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

    /**
     * returns true if the intervals have no intersection
     * 
     */
    private bool DoNotIntersect(Intervall cut)
    {
        return cut == null || cut.getLowerBound() > UPPERBOUND || cut.getUpperBound() < LOWERBOUND;
    }

    public float getRandom()
    {
        List<float> values = getEveryValue();
        int key = UnityEngine.Random.Range(0, values.Count);



        return values[key];
    }

    public override string ToString()
    {
        return "[" + LOWERBOUND + ";" + UPPERBOUND + "]";
    }

    public float Length()
    {
        return (UPPERBOUND - LOWERBOUND);
    }

    public float normedLength()
    {
        return (UPPERBOUND - LOWERBOUND) / _samplingRate;
    }

    public Vector2 toVector()
    {
        return new Vector2(LOWERBOUND, UPPERBOUND);
    }

    public static Intervall fromFloatArray(float[] res)
    {
        return new Intervall(Mathf.Max(res), Mathf.Min(res));
    }

    public static Intervall operator -(Intervall a, Intervall b)
    {
        return new Intervall(a.LOWERBOUND - b.UPPERBOUND, a.UPPERBOUND - b.LOWERBOUND);
    }

    public bool setSamplingRate(float newSamples)
    {
        if (newSamples < Length())
        {
            _samplingRate = newSamples;
            return true;
        }
        else return false;
    }
}
