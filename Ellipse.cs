using System;
using UnityEngine;

internal class Ellipse : conicSection
{
    //TODO Check 3D
    private Vector3 _majorAxis;
    private Vector3 _minorAxis;
    private Vector3 _middlePointEllipse;
    private Vector3[] _vertices;
    private Vector3 _foci1;
    private Vector3 _foci2;
    private float e;

    public Ellipse(Vector3 majorAxis, Vector3 minorAxis, Vector3 middlePointEllipse)
    {
        _majorAxis = majorAxis;
        _minorAxis = minorAxis;
       _middlePointEllipse = middlePointEllipse;

        calculateVertices();
        FindFoci();
    }

    private void calculateVertices()
    {
        Vector3[] res = new Vector3[4];

        //TODO check if correct in 3D
        Vector3 S1 = _middlePointEllipse + _majorAxis / 2;
        Vector3 S2 = _middlePointEllipse - _majorAxis / 2;
        Vector3 S3 = _middlePointEllipse + _minorAxis / 2;
        Vector3 S4 = _middlePointEllipse - _minorAxis / 2;

        res[0] = S1;
        res[1] = S2;
        res[2] = S3;
        res[3] = S4;

        _vertices =  res;
    }

    private void FindFoci()
    {
        float a = _majorAxis.magnitude/2;
        float b = _minorAxis.magnitude/2;

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
            Vector3 normedPos = pointToTest - _middlePointEllipse;
            float distToFoci1 = Mathf.Sqrt(Mathf.Pow((normedPos.x - e), 2) + Mathf.Pow(normedPos.y, 2));
            float distToFoci2 = Mathf.Sqrt(Mathf.Pow((normedPos.x + e), 2) + Mathf.Pow(normedPos.y, 2));

            return distToFoci1 + distToFoci2 <= _majorAxis.magnitude;
        }
        else return false;

    }

    public Vector3 IntersectEdge(conicSection otherForm)
    {
        throw new NotImplementedException();
    }
}