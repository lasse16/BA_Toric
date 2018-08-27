using System;
using System.Collections.Generic;
using UnityEngine;

/**
 * 
 * class for a 2D ellipse centered at 0,0 with the major axis pointing up - y direction, can be recalculated to a 3D ellipse
 * 
 */
internal class Ellipse : conicSection
{
    //TODO Check 3D
    private Vector3 _orientation;
    private Vector2 _majorAxis;
    private Vector2 _minorAxis;
    private Vector3 _middlePointEllipse;
    private Vector2[] _vertices;
    private Vector2 _foci1;
    private Vector2 _foci2;
    private float e;

    public Ellipse(Vector3 majorAxis, Vector3 minorAxis, Vector3 middlePointEllipse)
    {
        if (Vector3.Dot(majorAxis, minorAxis) == 0)
        {
            _orientation = majorAxis;
            _majorAxis = majorAxis.magnitude * Vector2.up;
            _minorAxis = minorAxis.magnitude * Vector2.right;
            _middlePointEllipse = middlePointEllipse;

            calculateVertices();
            FindFoci();
        }
        throw new ArgumentException(); 
    }

    public Ellipse(float distanceMajor, float distanceMinor, Vector2 midpoint)
    {
        _majorAxis = 2 *  distanceMajor * Vector2.up;
        _minorAxis = 2 * distanceMinor * Vector2.right;
        _middlePointEllipse = new Vector3(midpoint.x, midpoint.y, 0);

        FindFoci();
    }

    private void calculateVertices()
    {
        Vector2[] res = new Vector2[4];

        
        Vector2 S1 = _majorAxis / 2;
        Vector2 S2 = -_majorAxis / 2;
        Vector2 S3 = _minorAxis / 2;
        Vector2 S4 = -_minorAxis / 2;

        res[0] = S1;
        res[1] = S2;
        res[2] = S3;
        res[3] = S4;

        _vertices = res;
    }

    private void FindFoci()
    {
        float a = _majorAxis.magnitude / 2;
        float b = _minorAxis.magnitude / 2;

        float fociDistance = Mathf.Sqrt(Mathf.Pow(a, 2) - Mathf.Pow(b, 2));

        _foci1 = _majorAxis * fociDistance;
        _foci2 = _majorAxis * -fociDistance;
        e = fociDistance;
    }

    public Boolean IsInside(Vector3 pointToTest)
    {
        Vector3 normalOfPlaneOfEllipse = Vector3.Cross(_majorAxis, _minorAxis);
        if (Vector3.Dot(normalOfPlaneOfEllipse, pointToTest) == 0)
        {
            Vector2 distanceToFoci = calculateDistanceToFoci(pointToTest);

            return distanceToFoci[0] + distanceToFoci[1] <= _majorAxis.magnitude;
        }
        else return false;


    }

    private Vector2 calculateDistanceToFoci(Vector3 pointToTest)
    {
        Vector3 normedPos = pointToTest - _middlePointEllipse;
        float distToFoci1 = Mathf.Sqrt(Mathf.Pow((normedPos.x - e), 2) + Mathf.Pow(normedPos.y, 2));
        float distToFoci2 = Mathf.Sqrt(Mathf.Pow((normedPos.x + e), 2) + Mathf.Pow(normedPos.y, 2));

        return new Vector2(distToFoci1, distToFoci2);
    }


    
    
}