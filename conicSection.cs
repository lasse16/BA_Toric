using UnityEngine;

internal interface conicSection
{
     bool IsInside(Vector3 pointToTest);
    Vector3 IntersectEdge(conicSection otherForm);
}