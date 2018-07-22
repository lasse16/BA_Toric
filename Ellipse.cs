using UnityEngine;

internal class Ellipse
{
    private Vector3 _majorAxis;
    private Vector3 _minorAxis;
    private Vector3 _middlePointEllipse;

    public Ellipse(Vector3 majorAxis, Vector3 minorAxis, Vector3 middlePointEllipse)
    {
        _majorAxis = majorAxis;
        _minorAxis = minorAxis;
       _middlePointEllipse = middlePointEllipse;
    }
}