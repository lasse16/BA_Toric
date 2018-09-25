using System;
using UnityEngine;

/**
 * 
 * class for a 2D ellipse centered at 0,0 with the major axis pointing up - y direction, can be recalculated to a 3D ellipse
 * 
 */
internal class Ellipse : conicSection
{
    //TODO Check 3D
    private Vector3 _normal;
    private float _rotation;

    private Vector2 _majorAxis;
    private Vector2 _minorAxis;
    private Vector3 _middlePointEllipse;
    private Vector2[] _vertices;
    private Vector2 _foci1;
    private Vector2 _foci2;
    private float e;

    public Ellipse(Vector3 majorAxis, Vector3 minorAxis, Vector3 middlePointEllipse, float rotation = 0)
    {
        if (Vector3.Dot(majorAxis, minorAxis) == 0)
        {
            _normal =  Vector3.Cross(majorAxis,minorAxis).normalized;
            _majorAxis = majorAxis.magnitude * Vector2.up;
            _minorAxis = minorAxis.magnitude * Vector2.right;
            _middlePointEllipse = middlePointEllipse;
            _rotation = rotation;

            calculateVertices();
            FindFoci();
        }
        throw new ArgumentException(); 
    }

    public Ellipse(float distanceMajor, float distanceMinor, Vector3 midpoint, Vector3 orientation, float rotation = 0)
    {
        _majorAxis = 2 *  distanceMajor * Vector2.up;
        _minorAxis = 2 * distanceMinor * Vector2.right;
        _middlePointEllipse = midpoint;
        _normal = orientation.normalized;
        _rotation = rotation;

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
        
        if (Vector3.Dot(_normal, pointToTest) == 0)
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

    public void draw(Color color )
    {
        Vector3[] positions;
        int resolution = 16;
        float a = _majorAxis.magnitude / 2;
        float b = _minorAxis.magnitude / 2;

        

            positions = new Vector3[resolution + 1];
           
            Vector3 center = _middlePointEllipse;
        
            Quaternion q = Quaternion.AngleAxis(_rotation, _normal);
            q *= Quaternion.LookRotation(_normal);
            Debug.Log(_normal);

            for (int i = 0; i <= resolution; i++)
            {
                float angle = (float)i / (float)resolution * 2.0f * Mathf.PI;
                positions[i] = new Vector3(a * Mathf.Cos(angle), b * Mathf.Sin(angle), 0.0f);
                positions[i] = q * positions[i] + center;
            }

        Vector3 priorPoint = _middlePointEllipse;
            foreach (Vector3 point in positions)
            {
            Debug.DrawLine(priorPoint, point , color, Mathf.Infinity, false);
            priorPoint = point;
        }

        Debug.DrawLine(center + _normal, center, color, Mathf.Infinity, false);
    }

    public Vector3 getCenter()
    {
        return _middlePointEllipse;
    }
    
    
}