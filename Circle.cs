using UnityEngine;
using System.Linq;

internal class Circle: conicSection
{
    private Vector3 _middlePointOfCirclePhi;
    private float _radius;
    private Vector3 _normal;

    public Circle(Vector3 middlePointOfCirclePhi, float v, Vector3 normal)
    {
        _middlePointOfCirclePhi = middlePointOfCirclePhi;
        _radius = v;
        _normal = normal;
    }

    public Vector3 IntersectEdge(conicSection otherForm)
    {
        throw new System.NotImplementedException();
    }

    public bool IsInside(Vector3 pointToTest)
    {

        if (Vector3.Dot(_normal, pointToTest) == 0)
        {
            Vector3 normedPos = pointToTest - _middlePointOfCirclePhi;
            return Vector3.Distance(normedPos, _middlePointOfCirclePhi) <= _radius;
        }
        else return false;
    }
}