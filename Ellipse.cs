using System;
using UnityEngine;

/// <summary>
/// A representation of an 3D ellipse
/// </summary>
internal class Ellipse : conicSection
{
    /// <summary>
    /// The normal vector of the ellipse
    /// </summary>
    private Vector3 _normal;

    private float _rotation;

    private Vector2 _majorAxis;
    private Vector2 _minorAxis;

    public readonly Vector3 _middlePointEllipse;

    private Vector2 _foci1;
    private Vector2 _foci2;
    private float e;


    /// <summary>
    /// Sets up a 2D representation of the ellipse
    /// </summary>
    /// <param name="distanceMajor">the length along its major axis</param>
    /// <param name="distanceMinor">the length along its minor axis</param>
    /// <param name="midpoint">midpoint of the ellipse</param>
    /// <param name="normal">normal of the ellipse</param>
    /// <param name="rotation">rotation around the normal axis</param>
    public Ellipse(float distanceMajor, float distanceMinor, Vector3 midpoint, Vector3 normal, float rotation = 0)
    {
        _majorAxis = 2 * distanceMajor * Vector2.up;
        _minorAxis = 2 * distanceMinor * Vector2.right;
        _middlePointEllipse = midpoint;
        _normal = normal.normalized;
        _rotation = rotation;

        FindFoci();
    }


    /// <summary>
    /// Checks whether a 3D point is inside the 2D ellipse
    /// </summary>
    /// <param name="pointToTest">point to test</param>
    /// <returns></returns>
    public Boolean IsInside(Vector3 pointToTest)
    {
        if (Vector3.Dot(_normal, pointToTest) == 0)
        {
            Vector2 distanceToFoci = calculateDistanceToFoci(pointToTest);

            return distanceToFoci[0] + distanceToFoci[1] <= _majorAxis.magnitude;
        }
        else return false;
    }


    /// <summary>
    /// draws the ellipse in a desired color - using Debug.Draw
    /// </summary>
    /// <param name="color">desired color</param>
    public void Draw(Color color)
    {
        Vector3[] positions;
        int resolution = 32;
        float a = _majorAxis.magnitude / 2;
        float b = _minorAxis.magnitude / 2;

        positions = new Vector3[resolution + 1];

        Vector3 center = _middlePointEllipse;

        Quaternion q = Quaternion.AngleAxis(_rotation, _normal);
        q *= Quaternion.LookRotation(_normal);


        for (int i = 0; i <= resolution; i++)
        {
            float angle = (float)i / (float)resolution * 2.0f * Mathf.PI;
            positions[i] = new Vector3(a * Mathf.Cos(angle), b * Mathf.Sin(angle), 0.0f);
            positions[i] = q * positions[i] + center;
        }

        Vector3 priorPoint = positions[0];
        foreach (Vector3 point in positions)
        {
            Debug.DrawLine(priorPoint, point, color, Mathf.Infinity, false);
            priorPoint = point;
        }
    }


    //private methods


    private void FindFoci()
    {
        float a = _majorAxis.magnitude / 2;
        float b = _minorAxis.magnitude / 2;

        float fociDistance = Mathf.Sqrt(Mathf.Pow(a, 2) - Mathf.Pow(b, 2));

        _foci1 = _majorAxis * fociDistance;
        _foci2 = _majorAxis * -fociDistance;
        e = fociDistance;
    }

    private Vector2 calculateDistanceToFoci(Vector3 pointToTest)
    {
        Vector3 normedPos = pointToTest - _middlePointEllipse;
        float distToFoci1 = Mathf.Sqrt(Mathf.Pow((normedPos.x - e), 2) + Mathf.Pow(normedPos.y, 2));
        float distToFoci2 = Mathf.Sqrt(Mathf.Pow((normedPos.x + e), 2) + Mathf.Pow(normedPos.y, 2));

        return new Vector2(distToFoci1, distToFoci2);
    }
}