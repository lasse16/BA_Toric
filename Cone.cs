
using UnityEngine;

/// <summary>
/// A representation of a 3D cone defined by its apex, mid line and a deviation angle around the midline
/// </summary>
internal class Cone
{
    public readonly Vector3 _apex;
    public readonly Vector3 _midline;
    public readonly float _deviationAngle;


    public Cone(Vector3 midline, float deviationAngle, Vector3 apex)
    {
        _midline = midline;
        _deviationAngle = deviationAngle;
        _apex = apex;
    }

    public Cone(Vector3 prefferedVantageAngle, float deviationAngle) : this(prefferedVantageAngle, deviationAngle, Vector3.zero)
    { }

    /// <summary>
    /// Tests wether a 3D poin is inside the cone
    /// </summary>
    /// <param name="pointToTest">the point to test</param>
    /// <returns></returns>
    public bool IsInside(Vector3 pointToTest)
    {
        return Vector3.Angle(pointToTest - _apex, _midline) <= _deviationAngle;
    }


    public void Draw(Color color, float coneHeight = 2)
    {
        Vector3[] positions;
        int resolution = 16;
        float a = Mathf.Tan(_deviationAngle) * coneHeight;




        positions = new Vector3[resolution + 1];
        Quaternion q = Quaternion.LookRotation(_midline);
        Vector3 center = _apex + _midline * coneHeight;

        for (int i = 0; i <= resolution; i++)
        {
            float angle = (float)i / (float)resolution * 2.0f * Mathf.PI;
            positions[i] = new Vector3(a * Mathf.Cos(angle), a * Mathf.Sin(angle), 0.0f);
            positions[i] = q * positions[i] + center;
        }

        Vector3 priorPoint = center;
        foreach (Vector3 point in positions)
        {
            Debug.DrawLine(priorPoint, point, Color.green, Mathf.Infinity, false);
            Debug.DrawLine(_apex, point, Color.green, Mathf.Infinity, false);
            priorPoint = point;
        }
    }




}